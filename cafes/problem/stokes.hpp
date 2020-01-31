// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef PARTICLE_PROBLEM_STOKES_HPP_INCLUDED
#define PARTICLE_PROBLEM_STOKES_HPP_INCLUDED

#include <algorithm>
#include <fem/assembling.hpp>
#include <fem/matrixFree.hpp>
#include <fem/mesh.hpp>
#include <fem/rhs.hpp>
#include <iostream>
#include <petsc.h>
#include <problem/context.hpp>
#include <problem/options.hpp>
#include <problem/problem.hpp>

namespace cafes
{
    namespace problem
    {

#undef __FUNCT__
#define __FUNCT__ "setPMMSolver"
        PetscErrorCode setPMMSolver(KSP ksp, DM dm)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;
            PC pc, pc_i;
            KSP *sub_ksp;
            PetscInt MGlevels;
            DM dav, dap;
            DMDALocalInfo info;

            PetscFunctionBeginUser;

            ierr = DMCompositeGetEntries(dm, &dav, &dap);
            CHKERRQ(ierr);
            ierr = DMDAGetLocalInfo(dav, &info);
            CHKERRQ(ierr);

            ierr = KSPSetType(ksp, KSPGCR);
            CHKERRQ(ierr);
            ierr = KSPGetPC(ksp, &pc);
            CHKERRQ(ierr);
            ierr = PCSetType(pc, PCFIELDSPLIT);
            CHKERRQ(ierr);
            ierr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);
            CHKERRQ(ierr);
            ierr = PCFieldSplitSetOffDiagUseAmat(pc, PETSC_TRUE);
            CHKERRQ(ierr);
            ierr = PCFieldSplitSetSchurFactType(pc,
                                                PC_FIELDSPLIT_SCHUR_FACT_UPPER);
            CHKERRQ(ierr);
            ierr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER,
                                           PETSC_NULL);
            CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);
            CHKERRQ(ierr);

            ierr = PCFieldSplitGetSubKSP(pc, nullptr, &sub_ksp);
            CHKERRQ(ierr);

            /* Set MG solver on velocity field*/
            ierr = KSPGetPC(sub_ksp[0], &pc_i);
            CHKERRQ(ierr);
            ierr = KSPSetDM(sub_ksp[0], dav);
            CHKERRQ(ierr);
            ierr = KSPSetDMActive(sub_ksp[0], PETSC_FALSE);
            CHKERRQ(ierr);

            // ierr = KSPSetType(sub_ksp[0], KSPCG);CHKERRQ(ierr);
            ierr = KSPSetType(sub_ksp[0], KSPGCR);
            CHKERRQ(ierr);
            ierr = KSPSetTolerances(sub_ksp[0], 1e-2, PETSC_DEFAULT,
                                    PETSC_DEFAULT, PETSC_DEFAULT);
            CHKERRQ(ierr);
            ierr = PCSetType(pc_i, PCMG);
            CHKERRQ(ierr);

            auto i = (info.mx < info.my) ? info.mx : info.my;
            i = (info.mz == 1 || i < info.mz) ? i : info.mz;

            MGlevels = 1;
            while (i > 8)
            {
                i >>= 1;
                MGlevels++;
            }

            ierr = PCMGSetLevels(pc_i, MGlevels, PETSC_NULL);
            CHKERRQ(ierr);

            /* Set Jacobi preconditionner on pressure field*/
            ierr = KSPSetType(sub_ksp[1], KSPPREONLY);
            CHKERRQ(ierr);
            ierr = KSPSetTolerances(sub_ksp[1], 1e-1, PETSC_DEFAULT,
                                    PETSC_DEFAULT, PETSC_DEFAULT);
            CHKERRQ(ierr);
            ierr = KSPGetPC(sub_ksp[1], &pc_i);
            CHKERRQ(ierr);
            ierr = PCSetType(pc_i, PCJACOBI);
            CHKERRQ(ierr);

            // /* Set MG solver on pressure field*/
            // ierr = KSPGetPC(sub_ksp[1], &pc_i);CHKERRQ(ierr);
            // ierr = KSPSetDM(sub_ksp[0], dav);CHKERRQ(ierr);
            // ierr = KSPSetDMActive(sub_ksp[0], PETSC_FALSE);CHKERRQ(ierr);

            // ierr = KSPSetType(sub_ksp[1], KSPFGMRES);CHKERRQ(ierr);
            // ierr = KSPSetTolerances(sub_ksp[1], 1e-10, 1e-14, PETSC_DEFAULT,
            // PETSC_DEFAULT);CHKERRQ(ierr); ierr = PCSetType(pc_i,
            // PCMG);CHKERRQ(ierr);

            // ierr = DMDAGetLocalInfo(dap, &info);CHKERRQ(ierr);
            // i = (info.mx<info.my)? info.mx: info.my;
            // i = (info.mz == 1 || i<info.mz)? i: info.mz;

            // MGlevels = 1;
            // while(i > 8){
            //   i >>= 1;
            //   MGlevels++;
            // }

            // ierr = PCMGSetLevels(pc_i, MGlevels, PETSC_NULL);CHKERRQ(ierr);
            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "createLevelMatrices"
        template<typename CTX>
        PetscErrorCode createLevelMatrices(KSP ksp, Mat Alevel, Mat Plevel,
                                           void *ctx)
        {
            DM dm;
            PetscErrorCode ierr;
            CTX *ctx_ = (CTX *)ctx;
            int localsize, totalsize;
            PetscFunctionBeginUser;

            ierr = KSPGetDM(ksp, &dm);
            CHKERRQ(ierr);
            ctx_->dm = dm;

            ierr = fem::get_DM_sizes(dm, localsize, totalsize);
            CHKERRQ(ierr);

            ierr =
                MatSetSizes(Alevel, localsize, localsize, totalsize, totalsize);
            CHKERRQ(ierr);
            ierr = MatSetType(Alevel, MATSHELL);
            CHKERRQ(ierr);
            ierr = MatShellSetContext(Alevel, ctx);
            CHKERRQ(ierr);
            ierr = MatShellSetOperation(
                Alevel, MATOP_MULT, (void (*)())fem::diag_block_matrix<CTX>);
            CHKERRQ(ierr);
            ierr = MatShellSetOperation(
                Alevel, MATOP_GET_DIAGONAL,
                (void (*)())fem::diag_diag_block_matrix<CTX>);
            CHKERRQ(ierr);
            ierr = MatSetDM(Alevel, dm);
            CHKERRQ(ierr);
            ierr = MatSetFromOptions(Alevel);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

        //!
        //! Stokes class
        //!
        template<std::size_t Dimensions>
        struct stokes : public Problem<Dimensions>
        {
            // using parent = Problem<Dimensions>::Problem;
            using Ctx = context<Dimensions>;
            fem::rhs_conditions<Dimensions> rhsc_;
            options<Dimensions> opt{};

            Ctx *ctx; //!< Context for the PETSc solver
            Vec sol;  //!< The solution of Stokes problem
            Vec rhs;  //!< The RHS of Stokes problem
            Mat A;    //!< The matrix of Stokes problem
            Mat P;    //!< The preconditioner of Stokes problem
            KSP ksp;  //!< The solver of Stokes problem

            stokes(fem::dirichlet_conditions<Dimensions> bc,
                   fem::rhs_conditions<Dimensions> rhsc = {nullptr})
            {
                opt.process_options();

                DM mesh;
                fem::createMesh<Dimensions>(mesh, opt.mx, opt.xperiod);

                DMCreateGlobalVector(mesh, &sol);
                VecDuplicate(sol, &rhs);
                VecSet(rhs, 0.);

                std::array<double, Dimensions> hu;
                std::array<double, Dimensions> hp;
                for (std::size_t i = 0; i < Dimensions; ++i)
                {
                    hp[i] = opt.lx[i] / (opt.mx[i] - 1);
                    hu[i] = .5 * hp[i];
                }

                PetscErrorCode (*method)(
                    petsc::petsc_vec<Dimensions> &,
                    petsc::petsc_vec<Dimensions> &,
                    std::array<double, Dimensions> const &);

                if (opt.strain_tensor)
                    method = fem::strain_tensor_mult;
                else
                    method = fem::laplacian_mult;

                ctx = new Ctx{mesh, hu, method};
                ctx->set_dirichlet_bc(bc);

                DM dav, dap;
                DMCompositeGetEntries(mesh, &dav, &dap);

                rhsc_ = rhsc;

                // set Stokes matrix
                if (opt.assembling)
                {
                    PreallocateMat(ctx, opt, MATAIJ, &A, false);
                    PreallocateMat(ctx, opt, MATAIJ, &P, true);
                    MatSetDM(A, mesh);
                    MatSetDM(P, mesh);
                    MatSetFromOptions(A);
                    MatSetFromOptions(P);

                    if (opt.strain_tensor)
                    {
                        strain_tensor_assembling(dav, A, hu);
                        strain_tensor_assembling(dav, P, hu);
                    }
                    else
                    {
                        laplacian_assembling(dav, A, hu);
                        laplacian_assembling(dav, P, hu);
                    }

                    B_BT_assembling(mesh, A, hu);
                    diagonal_assembling<Dimensions>(mesh, A, 1.e-8); //0 à la place de 1e-8 (d'après H.)
                    cafes::mass_assembling(mesh, P, hp);

                    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
                    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
                    MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
                    MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);

                    cafes::fem::SetDirichletOnMat(A, bc);
                    cafes::fem::SetDirichletOnMat(P, bc);

                    // CreatePressureNullSpace(mesh); 
                }
                else
                {
                    // fix this to avoid raw pointer !!
                    A = fem::make_matrix<Ctx>(ctx, fem::stokes_matrix<Ctx>);
                    MatSetDM(A, mesh);
                    MatSetFromOptions(A);

                    // set preconditionner of Stokes matrix
                    DMSetMatType(dav, MATSHELL);
                    DMSetMatType(dap, MATSHELL);

                    Ctx *slap =
                        new Ctx{dav, hu, method, fem::diag_laplacian_mult};
                    slap->set_dirichlet_bc(bc);
                    auto A11 = fem::make_matrix<Ctx>(slap);

                    Ctx *smass =
                        new Ctx{dap, hp, fem::mass_mult, fem::diag_mass_mult};
                    auto A22 = fem::make_matrix<Ctx>(smass);

                    Mat bA[2][2];
                    bA[0][0] = A11;
                    bA[0][1] = PETSC_NULL;
                    bA[1][0] = PETSC_NULL;
                    bA[1][1] = A22;
                    MatCreateNest(PETSC_COMM_WORLD, 2, PETSC_NULL, 2,
                                  PETSC_NULL, &bA[0][0], &P);
                    MatSetDM(P, ctx->dm);
                    MatSetFromOptions(P);
                }
            }

#undef __FUNCT__
#define __FUNCT__ "setup_RHS"
            virtual PetscErrorCode setup_RHS() override
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                if (rhsc_.has_condition())
                {
                    ierr = fem::set_rhs<Dimensions, 1>(ctx->dm, rhs, rhsc_,
                                                       ctx->h);
                    CHKERRQ(ierr);
                }

                auto petsc_rhs = petsc::petsc_vec<Dimensions>(ctx->dm, rhs, 0);
                ierr = SetDirichletOnRHS(petsc_rhs, ctx->bc_, ctx->h);
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "setup_KSP"
            virtual PetscErrorCode setup_KSP() override
            {
                PetscErrorCode ierr;
                PC pc;
                PetscFunctionBeginUser;

                ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
                CHKERRQ(ierr);
                ierr = KSPSetDM(ksp, ctx->dm);
                CHKERRQ(ierr);
                ierr = KSPSetDMActive(ksp, PETSC_FALSE);
                CHKERRQ(ierr);
                ierr = KSPSetOptionsPrefix(ksp, "stokes_");
                CHKERRQ(ierr);
                ierr = KSPSetFromOptions(ksp);
                CHKERRQ(ierr);

                ierr = KSPSetOperators(ksp, A, P);
                CHKERRQ(ierr);

                if (opt.assembling)
                {
                    ierr = KSPSetUp(ksp);
                    CHKERRQ(ierr);
                    ierr = KSPGetPC(ksp, &pc);
                    CHKERRQ(ierr);

                    PetscBool ksp_preonly;
                    ierr = PetscObjectTypeCompare((PetscObject)ksp, KSPPREONLY,
                                                  &ksp_preonly);

                    PetscBool pc_lu;
                    ierr =
                        PetscObjectTypeCompare((PetscObject)pc, PCLU, &pc_lu);
                    if (ksp_preonly && pc_lu)
                    {
                        ierr = KSPSetOperators(ksp, A, A);
                        CHKERRQ(ierr);
                        ierr = PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
                        CHKERRQ(ierr);
                    }
                }
                else
                {
                    if (opt.pmm)
                    {
                        ierr = setPMMSolver(ksp, ctx->dm);
                        CHKERRQ(ierr);
                    }

                    ierr = KSPSetUp(ksp);
                    CHKERRQ(ierr);
                    ierr = KSPGetPC(ksp, &pc);
                    CHKERRQ(ierr);

                    PetscBool same;
                    ierr = PetscObjectTypeCompare((PetscObject)pc, PCFIELDSPLIT,
                                                  &same);

                    // If we use fieldsplit for the preconditionner
                    if (same)
                    {
                        KSP *subksp;
                        DM dav, dap;
                        ierr = PCFieldSplitGetSubKSP(pc, nullptr, &subksp);
                        CHKERRQ(ierr);
                        ierr = DMCompositeGetEntries(ctx->dm, &dav, &dap);
                        CHKERRQ(ierr);

                        ierr = KSPGetPC(subksp[0], &pc);
                        CHKERRQ(ierr);
                        ierr = KSPSetType(subksp[0], KSPGCR);
                        CHKERRQ(ierr);
                        ierr = KSPSetDM(subksp[0], dav);
                        CHKERRQ(ierr);
                        ierr = KSPSetDMActive(subksp[0], PETSC_FALSE);
                        CHKERRQ(ierr);
                        ierr = PetscObjectTypeCompare((PetscObject)pc, PCMG,
                                                      &same);
                        CHKERRQ(ierr);

                        // if MG is set for fieldsplit_0
                        if (same)
                        {
                            PetscErrorCode (*method)(
                                petsc::petsc_vec<Dimensions> &,
                                petsc::petsc_vec<Dimensions> &,
                                std::array<double, Dimensions> const &);

                            if (opt.strain_tensor)
                                method = fem::strain_tensor_mult;
                            else
                                method = fem::laplacian_mult;

                            PetscInt MGlevels;
                            KSP smoother;
                            ierr = PCMGGetLevels(pc, &MGlevels);
                            CHKERRQ(ierr);

                            for (std::size_t i = 0; i < MGlevels; ++i)
                            {
                                auto mg_h = ctx->h;
                                std::for_each(
                                    mg_h.begin(), mg_h.end(), [&](auto &x) {
                                        x *= (1 << (MGlevels - 1 - i));
                                    });

                                auto *mg_ctx =
                                    new Ctx{dav, mg_h, method,
                                            fem::diag_laplacian_mult};
                                mg_ctx->set_dirichlet_bc(ctx->bc_);

                                PC pcsmoother;
                                ierr = PCMGGetSmoother(pc, i, &smoother);
                                CHKERRQ(ierr);

                                ierr = KSPSetComputeOperators(
                                    smoother, createLevelMatrices<Ctx>,
                                    (void *)mg_ctx);
                                CHKERRQ(ierr);
                                // ierr = KSPSetType(smoother,
                                // KSPCG);CHKERRQ(ierr);
                                ierr = KSPSetType(smoother, KSPCG);
                                CHKERRQ(ierr);
                                ierr = KSPGetPC(smoother, &pcsmoother);
                                CHKERRQ(ierr);
                                ierr = PCSetType(pcsmoother, PCJACOBI);
                                CHKERRQ(ierr);
                            }

                            ierr = PCSetUp(pc);
                            CHKERRQ(ierr);
                        }

                        // ierr = KSPGetPC(subksp[1], &pc);CHKERRQ(ierr);
                        // ierr = KSPSetType(subksp[1],
                        // KSPFGMRES);CHKERRQ(ierr); ierr = KSPSetDM(subksp[1],
                        // dap);CHKERRQ(ierr); ierr = KSPSetDMActive(subksp[1],
                        // PETSC_FALSE);CHKERRQ(ierr); ierr =
                        // PetscObjectTypeCompare((PetscObject)pc, PCMG,
                        // &same);CHKERRQ(ierr);

                        // // if MG is set for fieldsplit_1
                        // if (same) {
                        //   PetscInt MGlevels;
                        //   KSP smoother;
                        //   ierr = PCMGGetLevels(pc, &MGlevels);CHKERRQ(ierr);

                        //   for(std::size_t i=0; i<MGlevels; ++i){
                        //     auto mg_h{ctx->h};
                        //     std::for_each(mg_h.begin(), mg_h.end(), [&](auto&
                        //     x){x*=2;});

                        //     std::for_each(mg_h.begin(), mg_h.end(), [&](auto&
                        //     x){x*=(1<<(MGlevels-1-i));});

                        //     auto *mg_ctx = new Ctx{dap, mg_h, fem::mass_mult,
                        //     fem::mass_mult_diag};
                        //     //mg_ctx->set_dirichlet_bc(ctx->bc_);

                        //     PC pcsmoother;
                        //     ierr = PCMGGetSmoother(pc, i,
                        //     &smoother);CHKERRQ(ierr);

                        //     ierr = KSPSetComputeOperators(smoother,
                        //     createLevelMatrices<Ctx>, (void *)
                        //     mg_ctx);CHKERRQ(ierr); ierr =
                        //     KSPSetType(smoother, KSPFGMRES);CHKERRQ(ierr);
                        //     ierr = KSPGetPC(smoother,
                        //     &pcsmoother);CHKERRQ(ierr); ierr =
                        //     PCSetType(pcsmoother, PCJACOBI);CHKERRQ(ierr);
                        //   }

                        //   ierr = PCSetUp(pc);CHKERRQ(ierr);
                        // }
                    }
                }
                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "CreatePressureNullSpace"
            PetscErrorCode CreatePressureNullSpace(DM mesh)
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                MatNullSpace nullSpacePres;

                Vec tmp;
                ierr = DMGetGlobalVector(mesh, &tmp);
                CHKERRQ(ierr);
                ierr = VecSet(tmp, 0.);
                CHKERRQ(ierr);

                auto lambda = [](const auto &dm, auto &v) {
                    auto xpetsc = petsc::petsc_vec<Dimensions>(dm, v, 1, false);

                    xpetsc.fill_global(1.);
                };

                lambda(mesh, tmp);
                ierr = VecNormalize(tmp, NULL);
                CHKERRQ(ierr);

                ierr = MatNullSpaceCreate(PetscObjectComm((PetscObject)mesh),
                                          PETSC_FALSE, 1, &tmp, &nullSpacePres);

                ierr = MatSetNullSpace(A, nullSpacePres);
                CHKERRQ(ierr);
                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "solve"
            virtual PetscErrorCode solve() override
            {
                PetscErrorCode ierr;
                PC pc;
                PetscFunctionBeginUser;

                ierr = KSPSolve(ksp, rhs, sol);
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }
        };
    } // namespace problem
    template<std::size_t Dimensions, std::size_t Ndof>
    problem::stokes<Dimensions>
    make_stokes(fem::dirichlet_conditions<Dimensions> const &dc,
                fem::rhs_conditions<Ndof> const &rhs)
    {
        return {dc, rhs};
    }

} // namespace cafes

#endif