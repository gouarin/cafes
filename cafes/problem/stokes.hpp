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
#include <petsc.h>
#include <problem/context.hpp>
#include <problem/options.hpp>
#include <problem/problem.hpp>

namespace cafes
{
    namespace problem
    {

        #undef __FUNCT__
        #define __FUNCT__ "createLevelMatrices"
        template<typename CTX, std::size_t Dim>
        PetscErrorCode createLevelMatrices(KSP ksp, Mat Alevel, Mat Plevel, void *ctx)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            CTX *ctx_ = (CTX *)ctx;

            DM dm;
            ierr = KSPGetDM(ksp, &dm);CHKERRQ(ierr);

            // DMView(dm, PETSC_VIEWER_STDOUT_WORLD);

            auto bounds = fem::get_global_bounds<Dim>(dm);
            auto *mg_ctx = new CTX(*ctx_);
            mg_ctx->dm = dm;
            for(std::size_t d = 0; d < Dim; ++d)
            {
                mg_ctx->h[d] = mg_ctx->domain_length[d]/(bounds[d] - 1);
            }

            PetscBool is_shell;
            MatType Atype;
            ierr = MatGetType(Alevel, &Atype);CHKERRQ(ierr);
            ierr = PetscStrcmp(Atype, MATSHELL, &is_shell);CHKERRQ(ierr);

            // std::cout << Atype << std::endl;
            if (is_shell)
            {
                ierr = MatShellSetContext(Alevel, mg_ctx);CHKERRQ(ierr);
                ierr = MatShellSetOperation(Alevel, MATOP_MULT, (void (*)())fem::diag_block_matrix<CTX>);CHKERRQ(ierr);
                ierr = MatShellSetOperation(Alevel, MATOP_GET_DIAGONAL,(void (*)())fem::diag_diag_block_matrix<CTX>);CHKERRQ(ierr);
            }
            else
            {
                ierr = laplacian_assembling(dm, Alevel, mg_ctx->h, mg_ctx->order);CHKERRQ(ierr);
                ierr = MatAssemblyBegin(Alevel, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
                ierr = MatAssemblyEnd(Alevel, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
                ierr = cafes::fem::SetDirichletOnMat(Alevel, mg_ctx->bc_);CHKERRQ(ierr);
            }

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "setPMMSolver"
        template<class Ctx, std::size_t Dim>
        PetscErrorCode setPMMSolver(KSP ksp, DM dm, Ctx *ctx)
        {
            PetscErrorCode ierr;
            PC pc, pc_i;
            KSP *sub_ksp;
            PetscInt MGlevels;
            DM dav, dap;
            DMDALocalInfo info;

            PetscFunctionBeginUser;


            ierr = KSPSetType(ksp, KSPGCR);CHKERRQ(ierr);
            ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
            ierr = PCSetType(pc, PCFIELDSPLIT);CHKERRQ(ierr);
            ierr = PCFieldSplitSetType(pc, PC_COMPOSITE_SCHUR);CHKERRQ(ierr);
            ierr = PCFieldSplitSetOffDiagUseAmat(pc, PETSC_TRUE);CHKERRQ(ierr);
            ierr = PCFieldSplitSetSchurFactType(pc, PC_FIELDSPLIT_SCHUR_FACT_UPPER);CHKERRQ(ierr);
            ierr = PCFieldSplitSetSchurPre(pc, PC_FIELDSPLIT_SCHUR_PRE_USER, PETSC_NULL);CHKERRQ(ierr);

            ierr = DMCompositeGetEntries(dm, &dav, &dap);CHKERRQ(ierr);
            // auto *mg_ctx = new Ctx{dav, ctx->h, ctx->order, ctx->apply, fem::diag_laplacian_mult};
            // mg_ctx->set_dirichlet_bc(ctx->bc_);
            auto *mg_ctx = new Ctx(*ctx);
            mg_ctx->dm = dav;
            mg_ctx->scale = 0.5;
            mg_ctx->apply_diag = fem::diag_laplacian_mult;

            ierr = DMKSPSetComputeOperators(dav, createLevelMatrices<Ctx, Dim>, (void *)mg_ctx);CHKERRQ(ierr);

            // ierr = PCSetUp(pc);CHKERRQ(ierr);

            // ierr = PCFieldSplitGetSubKSP(pc, nullptr, &sub_ksp);CHKERRQ(ierr);

            // // /* Set MG solver on velocity field*/
            // ierr = KSPGetPC(sub_ksp[0], &pc_i);CHKERRQ(ierr);
            // ierr = KSPSetDM(sub_ksp[0], dav);CHKERRQ(ierr);
            // ierr = KSPSetDMActive(sub_ksp[0], PETSC_FALSE);CHKERRQ(ierr);

            // ierr = KSPSetType(sub_ksp[0], KSPGCR);CHKERRQ(ierr);
            // ierr = KSPSetComputeOperators(sub_ksp[0], createLevelMatrices<Ctx, Dim>, (void *)mg_ctx);CHKERRQ(ierr);

            // ierr = KSPSetTolerances(sub_ksp[0], 1e-2, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
            // ierr = PCSetType(pc_i, PCMG);CHKERRQ(ierr);

            // ierr = DMDAGetLocalInfo(dav, &info);CHKERRQ(ierr);
            // auto i = (info.mx < info.my) ? info.mx : info.my;
            // i = (info.mz == 1 || i < info.mz) ? i : info.mz;

            // // MGlevels = 1;
            // // while (i > 8)
            // // {
            // //     i >>= 1;
            // //     MGlevels++;
            // // }
            // // loic: to remove
            // MGlevels = 3;

            // ierr = PCMGSetLevels(pc_i, MGlevels, PETSC_NULL);CHKERRQ(ierr);

            // for (std::size_t i = 0; i < MGlevels; ++i)
            // {
            //     KSP smoother;
            //     PC pcsmoother;
            //     ierr = PCMGGetSmoother(pc_i, i, &smoother);CHKERRQ(ierr);
            //     ierr = KSPSetType(smoother, KSPCG);CHKERRQ(ierr);
            //     ierr = KSPGetPC(smoother, &pcsmoother);CHKERRQ(ierr);
            //     ierr = PCSetType(pcsmoother, PCJACOBI);CHKERRQ(ierr);
            // }

            // PCSetUp(pc_i);
            // KSP coarse;
            // ierr = PCMGGetCoarseSolve(pc_i, &coarse);CHKERRQ(ierr);

            // PC pc_coarse;
            // ierr = KSPGetPC(coarse, &pc_coarse);CHKERRQ(ierr);
            // ierr = PCSetType(pc_coarse, PCTELESCOPE);CHKERRQ(ierr);
            // // ierr = PCTelescopeSetReductionFactor(pc_coarse, 4);CHKERRQ(ierr);

            // PCSetUp(pc_coarse);
            // KSP kspt;
            // PC pct;
            // ierr = PCTelescopeGetKSP(pc_coarse, &kspt);CHKERRQ(ierr);
            // // // ierr = KSPSetType(kspt, KSPCG);
            // ierr = KSPGetPC(kspt, &pct);
            // PCSetType(pct, PCLU);
            // // PCSetType(pct, PCMG);
            // // ierr = PCMGSetLevels(pct, 3, PETSC_NULL);CHKERRQ(ierr);
            // // ierr = KSPSetDM(kspt, dav);CHKERRQ(ierr);

            // // for (std::size_t i = 0; i < 2; ++i)
            // // {
            // //     KSP smoother;
            // //     PC pcsmoother;
            // //     ierr = PCMGGetSmoother(pct, i, &smoother);CHKERRQ(ierr);
            // //     ierr = KSPSetType(smoother, KSPCG);CHKERRQ(ierr);
            // //     ierr = KSPGetPC(smoother, &pcsmoother);CHKERRQ(ierr);
            // //     ierr = PCSetType(pcsmoother, PCJACOBI);CHKERRQ(ierr);
            // // }

            // // ierr = PCMGGetCoarseSolve(pct, &coarse);CHKERRQ(ierr);


            // // PCView(pct, PETSC_VIEWER_STDOUT_WORLD);
            // // Mat At, Pt, newAt;
            // // PCGetOperators(pct, &At, &Pt);

            // // DM dmc;
            // // MatGetDM(At, &dmc);
            // // auto bounds = fem::get_global_bounds<Dim>(dmc);
            // // auto *telescope_ctx = new Ctx(*ctx);
            // // telescope_ctx->dm = dmc;
            // // for(std::size_t d = 0; d < Dim; ++d)
            // // {
            // //     telescope_ctx->h[d] = ctx->domain_length[d]/(bounds[d] - 1);
            // // }

            // // {
            // //     int localsize, totalsize;

            // //     fem::get_DM_sizes(dmc, localsize, totalsize);
            // //     MatCreateShell(PETSC_COMM_WORLD, localsize, localsize, totalsize,
            // //                    totalsize, telescope_ctx, &newAt);
            // //     ierr = MatShellSetContext(newAt, telescope_ctx);CHKERRQ(ierr);
            // //     ierr = MatShellSetOperation(newAt, MATOP_MULT, (void (*)())fem::diag_block_matrix<Ctx>);CHKERRQ(ierr);
            // //     ierr = MatShellSetOperation(newAt, MATOP_GET_DIAGONAL,(void (*)())fem::diag_diag_block_matrix<Ctx>);CHKERRQ(ierr);
            // // }

            // // // PCSetOperators(pct, newAt, newAt);
            // // // // DMView(dmc, PETSC_VIEWER_STDOUT_WORLD);

            // // laplacian_assembling(dmc, At, telescope_ctx->h, telescope_ctx->order);
            // // MatAssemblyBegin(Pt, MAT_FINAL_ASSEMBLY);
            // // MatAssemblyEnd(Pt, MAT_FINAL_ASSEMBLY);
            // // MatCopy(At, Pt, SAME_NONZERO_PATTERN);
            // // MatView(Pt, PETSC_VIEWER_STDOUT_WORLD);
            // // cafes::fem::SetDirichletOnMat(At, ctx->bc_);


            // // DM dmc;
            // // KSPGetDM(coarse, &dmc);
            // // // DMView(dm, PETSC_VIEWER_STDOUT_WORLD);
            // // ierr = DMKSPSetComputeOperators(dmc, createLevelMatrices<Ctx, Dim>, (void *)mg_ctx);CHKERRQ(ierr);
            // // KSPSetUp(coarse);

            // // KSPView(coarse, PETSC_VIEWER_STDOUT_WORLD);
            // // DM dm;
            // // KSPGetDM(coarse, &dm);
            // // DMView(dm, PETSC_VIEWER_STDOUT_WORLD);
            // // Mat At, Pt;
            // // KSPGetOperators(coarse, &At, &Pt);
            // // MatSetType(At, MATSHELL);
            // // MatSetType(Pt, MATSHELL);
            //     // MatView(At, PETSC_VIEWER_STDERR_WORLD);

            // /* Set Jacobi preconditionner on pressure field*/
            // ierr = KSPSetType(sub_ksp[1], KSPPREONLY);CHKERRQ(ierr);
            // ierr = KSPSetTolerances(sub_ksp[1], 1e-1, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);CHKERRQ(ierr);
            // ierr = KSPGetPC(sub_ksp[1], &pc_i);CHKERRQ(ierr);
            // ierr = PCSetType(pc_i, PCJACOBI);CHKERRQ(ierr);

            // // /* Set MG solver on pressure field*/
            // // ierr = KSPGetPC(sub_ksp[1], &pc_i);CHKERRQ(ierr);
            // // ierr = KSPSetDM(sub_ksp[0], dav);CHKERRQ(ierr);
            // // ierr = KSPSetDMActive(sub_ksp[0], PETSC_FALSE);CHKERRQ(ierr);

            // // ierr = KSPSetType(sub_ksp[1], KSPFGMRES);CHKERRQ(ierr);
            // // ierr = KSPSetTolerances(sub_ksp[1], 1e-10, 1e-14, PETSC_DEFAULT,
            // // PETSC_DEFAULT);CHKERRQ(ierr); ierr = PCSetType(pc_i,
            // // PCMG);CHKERRQ(ierr);

            // // ierr = DMDAGetLocalInfo(dap, &info);CHKERRQ(ierr);
            // // i = (info.mx<info.my)? info.mx: info.my;
            // // i = (info.mz == 1 || i<info.mz)? i: info.mz;

            // // MGlevels = 1;
            // // while(i > 8){
            // //   i >>= 1;
            // //   MGlevels++;
            // // }

            // // ierr = PCMGSetLevels(pc_i, MGlevels, PETSC_NULL);CHKERRQ(ierr);
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
            static constexpr std::size_t Dim = Dimensions;
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
                    std::array<double, Dimensions> const &,
                    std::size_t);

                if (opt.strain_tensor)
                    method = fem::strain_tensor_mult;
                else
                    method = fem::laplacian_mult;

                ctx = new Ctx{mesh, opt.lx, 1., hu, opt.order, method};
                ctx->set_dirichlet_bc(bc);

                DM dav, dap;
                DMCompositeGetEntries(mesh, &dav, &dap);

                rhsc_ = rhsc;

                // set Stokes matrix
                if (opt.assembling)
                {
                    PreallocateMat<Dimensions>(ctx, opt, MATAIJ, &A, opt.order, false);
                    PreallocateMat<Dimensions>(ctx, opt, MATAIJ, &P, opt.order, true);

                    MatSetDM(A, mesh);
                    MatSetDM(P, mesh);
                    MatSetFromOptions(A);
                    MatSetFromOptions(P);

                    if (opt.strain_tensor)
                    {
                        strain_tensor_assembling(dav, A, hu, opt.order);
                        strain_tensor_assembling(dav, P, hu, opt.order);
                    }
                    else
                    {
                        laplacian_assembling(dav, A, hu, opt.order);
                        laplacian_assembling(dav, P, hu, opt.order);
                    }

                    B_BT_assembling(mesh, A, hu, opt.order);

                    diagonal_assembling<Dimensions>(mesh, A, 0); //0 à la place de 1e-8 (d'après H.)

                    DMDALocalInfo info;
                    DMDAGetLocalInfo(dav, &info);
                    int dec = info.dof * info.gxm * info.gym * info.gzm;
                    // cafes::mass_assembling(dap, A, hp, 1, dec, -1e-10);
                    cafes::mass_assembling(dap, P, hp, 1, dec);

                    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
                    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
                    MatAssemblyBegin(P, MAT_FINAL_ASSEMBLY);
                    MatAssemblyEnd(P, MAT_FINAL_ASSEMBLY);

                    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);

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

                    Ctx *slap = new Ctx{dav, opt.lx, 0.5, hu, opt.order, method, fem::diag_laplacian_mult};
                    slap->set_dirichlet_bc(bc);
                    auto A11 = fem::make_matrix<Ctx>(slap);

                    Ctx *smass = new Ctx{dap, opt.lx, 1, hp, 1, fem::mass_mult, fem::diag_mass_mult};
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
                    ierr = fem::set_rhs<Dimensions, 1>(ctx->dm, rhs, rhsc_, ctx->h, opt.order);CHKERRQ(ierr);
                }

                auto petsc_rhs = petsc::petsc_vec<Dimensions>(ctx->dm, rhs, 0);
                ierr = SetDirichletOnRHS(petsc_rhs, ctx->bc_, ctx->h);CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

            #undef __FUNCT__
            #define __FUNCT__ "setup_KSP"
            virtual PetscErrorCode setup_KSP() override
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;
                PC pc;

                ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
                ierr = KSPSetDM(ksp, ctx->dm);CHKERRQ(ierr);
                ierr = KSPSetDMActive(ksp, PETSC_FALSE);CHKERRQ(ierr);
                ierr = KSPSetOptionsPrefix(ksp, "stokes_");CHKERRQ(ierr);
                ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

                ierr = KSPSetOperators(ksp, A, P);CHKERRQ(ierr);

                if (opt.assembling)
                {
                    ierr = KSPSetUp(ksp);CHKERRQ(ierr);
                    ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);

                    PetscBool ksp_preonly;
                    ierr = PetscObjectTypeCompare((PetscObject)ksp, KSPPREONLY, &ksp_preonly);CHKERRQ(ierr);

                    PetscBool pc_lu;
                    ierr = PetscObjectTypeCompare((PetscObject)pc, PCLU, &pc_lu);CHKERRQ(ierr);

                    if (ksp_preonly && pc_lu)
                    {
                        ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
                        ierr = PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);CHKERRQ(ierr);
                    }
                    if (opt.pmm)
                    {
                        ierr = setPMMSolver<Ctx, Dim>(ksp, ctx->dm, ctx);CHKERRQ(ierr);
                    }
                }
                else
                {
                    if (opt.pmm)
                    {
                        ierr = setPMMSolver<Ctx, Dim>(ksp, ctx->dm, ctx);CHKERRQ(ierr);
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

                // Vec tmp;
                // ierr = DMGetGlobalVector(mesh, &tmp);
                // CHKERRQ(ierr);
                // ierr = VecSet(tmp, 0.);
                // CHKERRQ(ierr);

                // auto lambda = [](const auto &dm, auto &v) {
                //     auto xpetsc = petsc::petsc_vec<Dimensions>(dm, v, 1, false);

                //     xpetsc.fill_global(1.);
                // };

                // lambda(mesh, tmp);
                // ierr = VecNormalize(tmp, NULL);
                // CHKERRQ(ierr);

                // ierr = MatNullSpaceCreate(PetscObjectComm((PetscObject)mesh),
                //                           PETSC_FALSE, 1, &tmp, &nullSpacePres);

                ierr = MatNullSpaceCreate(PetscObjectComm((PetscObject)mesh),
                                          PETSC_TRUE, 0, PETSC_NULL, &nullSpacePres);
                CHKERRQ(ierr);
                ierr = MatSetNullSpace(A, nullSpacePres);
                CHKERRQ(ierr);
                PetscFunctionReturn(0);
            }

            #undef __FUNCT__
            #define __FUNCT__ "solve"
            virtual PetscErrorCode solve() override
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;
                PC pc;

                ierr = KSPSolve(ksp, rhs, sol);CHKERRQ(ierr);

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