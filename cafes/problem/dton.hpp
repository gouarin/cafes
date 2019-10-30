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

#ifndef CAFES_PROBLEM_DTON_HPP_INCLUDED
#define CAFES_PROBLEM_DTON_HPP_INCLUDED

#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>
#include <particle/geometry/vector.hpp>
#include <particle/particle.hpp>
#include <particle/singularity/add_singularity.hpp>
#include <problem/particle_operator.hpp>
#include <problem/problem.hpp>
#include <problem/stokes.hpp>

//#include <io/vtk.hpp>

#include <algorithm>
#include <iostream>
#include <memory>
#include <petsc.h>

namespace cafes
{
    namespace problem
    {

#undef __FUNCT__
#define __FUNCT__ "DtoN_matrix"
        template<std::size_t Dimensions, typename Ctx>
        PetscErrorCode DtoN_matrix(Mat A, Vec x, Vec y)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            Ctx *ctx;
            ierr = MatShellGetContext(A, (void **)&ctx);
            CHKERRQ(ierr);

            ierr = VecSet(ctx->problem.rhs, 0.);
            CHKERRQ(ierr);

            ierr = init_problem<Dimensions, Ctx>(*ctx, x);
            CHKERRQ(ierr);
        
            ierr = ctx->problem.solve();
            // if (!ctx->compute_rhs)
            // {
            //     ierr = cafes::io::save_VTK("Resultats", "two_part_ug_rhs",
            //                                 ctx->problem.rhs, ctx->problem.ctx->dm,
            //                                 ctx->problem.ctx->h);
            //     CHKERRQ(ierr);
            //     ierr = cafes::io::save_VTK("Resultats", "two_part_ug_sol",
            //                                 ctx->problem.sol, ctx->problem.ctx->dm,
            //                                 ctx->problem.ctx->h);
            //     CHKERRQ(ierr);
            //     std::exit(0);
            // }

            ierr = VecCopy(ctx->problem.sol, ctx->sol_tmp);
            CHKERRQ(ierr);

            std::vector<std::vector<geometry::vector<double, Dimensions>>> g;
            g.resize(ctx->particles.size());
            for (std::size_t ipart = 0; ipart < ctx->surf_points.size();
                 ++ipart)
                g[ipart].resize(ctx->surf_points[ipart].size());

            // interpolation
            ierr = interp_fluid_to_surf(*ctx, g, ctx->add_rigid_motion,
                                        ctx->compute_singularity);
            CHKERRQ(ierr);

            ierr = SL_to_Rhs(*ctx, g);
            CHKERRQ(ierr);

            ierr = ctx->problem.solve();
            CHKERRQ(ierr);

            ierr = compute_y<Dimensions, Ctx>(*ctx, y);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

        template<typename Shape, std::size_t Dimensions, typename Problem_type>
        struct DtoN : public Problem<Dimensions>
        {
            std::vector<particle<Shape>> parts_;
            Problem_type problem_;

            using position_type = geometry::position<double, Dimensions>;
            using position_type_i = geometry::position<int, Dimensions>;

            using Ctx = particle_context<Dimensions, Shape, Problem_type>;
            Ctx *ctx;

            std::vector<std::vector<std::pair<position_type_i, position_type>>>
                surf_points_;
            std::vector<std::vector<geometry::vector<double, Dimensions>>>
                radial_vec_;
            std::vector<int> nb_surf_points_;
            std::vector<int> num_;
            Vec sol;
            Vec rhs;
            Vec sol_rhs, sol_g, sol_tmp;
            Mat A;
            KSP ksp;
            std::size_t scale_ = 4;
            bool default_flags_ = true;
            bool use_sing = false;

            using dpart_type =
                typename std::conditional<Dimensions == 2, double,
                                          std::array<double, 2>>::type;
            dpart_type dpart_;

            DtoN(std::vector<particle<Shape>> &parts, Problem_type &p,
                 dpart_type dpart)
                : parts_{parts}, problem_{p}, dpart_{dpart}
            {
                problem_.setup_KSP();
                VecDuplicate(problem_.sol, &sol_tmp);
                VecDuplicate(problem_.sol, &sol_rhs);
                VecDuplicate(problem_.sol, &sol_g);
            }

#undef __FUNCT__
#define __FUNCT__ "create_Mat_and_Vec"
            PetscErrorCode create_Mat_and_Vec()
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                auto box = fem::get_DM_bounds<Dimensions>(problem_.ctx->dm, 0);
                auto &h = problem_.ctx->h;

                auto size = set_materials(parts_, surf_points_, radial_vec_,
                                          nb_surf_points_, num_, box, h, dpart_,
                                          scale_);

                ctx =
                    new Ctx{problem_,        parts_, surf_points_, radial_vec_,
                            nb_surf_points_, num_,   scale_,       false,
                            false,           false,  sol_tmp};

                ierr = MatCreateShell(PETSC_COMM_WORLD, size * Dimensions,
                                      size * Dimensions, PETSC_DECIDE,
                                      PETSC_DECIDE, ctx, &A);
                CHKERRQ(ierr);
                ierr = MatShellSetOperation(
                    A, MATOP_MULT,
                    (void (*)(void))DtoN_matrix<Dimensions, Ctx>);
                CHKERRQ(ierr);

                ierr = MatCreateVecs(A, &sol, &rhs);
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "setup_RHS"
            virtual PetscErrorCode setup_RHS() override
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                if (default_flags_)
                {
                    ctx->compute_rhs = true;
                    ctx->add_rigid_motion = true;
                    ctx->compute_singularity = false;
                    // ctx->compute_singularity = false;
                }

                ierr = MatMult(A, sol, rhs);
                CHKERRQ(ierr);
                ierr = VecScale(rhs, -1.);
                CHKERRQ(ierr);

                // ierr = VecCopy(sol_tmp, sol_rhs);
                // CHKERRQ(ierr);
                // ierr = cafes::io::save_VTK("Resultats", "two_part_u0",
                // sol_tmp,
                //                            problem_.ctx->dm,
                //                            problem_.ctx->h);
                // CHKERRQ(ierr);
                // ierr = cafes::io::save_VTK("Resultats", "two_part_w0",
                //                            problem_.sol, problem_.ctx->dm,
                //                            problem_.ctx->h);
                // CHKERRQ(ierr);
                // ierr = cafes::io::save_VTK("Resultats", "two_part_w0_rhs",
                //                            problem_.rhs, problem_.ctx->dm,
                //                            problem_.ctx->h);
                // CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "setup_KSP"
            virtual PetscErrorCode setup_KSP() override
            {
                PetscErrorCode ierr;
                PC pc;
                PetscFunctionBegin;

                ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
                CHKERRQ(ierr);
                ierr = KSPSetOptionsPrefix(ksp, "dton_");
                CHKERRQ(ierr);

                ierr = KSPSetOperators(ksp, A, A);
                CHKERRQ(ierr);
                ierr = KSPGetPC(ksp, &pc);
                CHKERRQ(ierr);
                ierr = PCSetType(pc, PCNONE);
                CHKERRQ(ierr);
                ierr = KSPSetFromOptions(ksp);
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "solve_last_problem"
            PetscErrorCode solve_last_problem()
            {
                PetscErrorCode ierr;
                PetscFunctionBegin;
                // solve the problem with the right control

                ierr = VecSet(ctx->problem.rhs, 0.);
                CHKERRQ(ierr);

                if (default_flags_)
                {
                    ctx->compute_rhs = false;
                    ctx->add_rigid_motion = false;
                    ctx->compute_singularity = false;
                }

                ierr = init_problem<Dimensions, Ctx>(*ctx, sol);
                CHKERRQ(ierr);
                // std::cout<<ctx->compute_rhs<<", "<<ctx->add_rigid_motion<<",
                // "<<ctx->compute_singularity<<"\n"; ierr =
                // cafes::io::save_VTK("Resultats", "two_part_reg", sol,
                // ctx->problem.ctx->dm, ctx->problem.ctx->h);CHKERRQ(ierr);

                ierr = ctx->problem.solve();
                CHKERRQ(ierr);
                // ierr = cafes::io::save_VTK("Resultats", "two_part_tilde_ug",
                // problem_.sol, problem_.ctx->dm,
                // problem_.ctx->h);CHKERRQ(ierr);

                VecAXPY(problem_.sol, 1, sol_rhs);
                // ierr = cafes::io::save_VTK("Resultats", "two_part_new_u",
                // problem_.sol, problem_.ctx->dm,
                // problem_.ctx->h);CHKERRQ(ierr);

                ierr = VecSet(ctx->problem.rhs, 0.);
                CHKERRQ(ierr);

                if (default_flags_)
                {
                    ctx->compute_rhs = true;
                    ctx->add_rigid_motion = true;
                    ctx->compute_singularity = true;
                }

                ierr = init_problem<Dimensions, Ctx>(*ctx, sol);
                CHKERRQ(ierr);
                // std::cout<<ctx->compute_rhs<<", "<<ctx->add_rigid_motion<<",
                // "<<ctx->compute_singularity<<"\n"; ierr =
                // cafes::io::save_VTK("Resultats", "two_part_reg", sol,
                // ctx->problem.ctx->dm, ctx->problem.ctx->h);CHKERRQ(ierr);

                ierr = ctx->problem.solve();
                CHKERRQ(ierr);
                // ierr = cafes::io::save_VTK("Resultats", "two_part_ureg",
                // problem_.sol, problem_.ctx->dm,
                // problem_.ctx->h);CHKERRQ(ierr);

                if (use_sing)
                {
                    ierr = singularity::add_singularity_to_last_sol<Dimensions,
                                                                    Ctx>(
                        *ctx, problem_.sol);
                    CHKERRQ(ierr);
                }
                // ierr = cafes::io::save_VTK("Resultats", "two_part_u",
                // problem_.sol, problem_.ctx->dm,
                // problem_.ctx->h);CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "solve"
            virtual PetscErrorCode solve() override
            {
                PetscErrorCode ierr;
                PetscFunctionBegin;

                if (default_flags_)
                {
                    ctx->compute_rhs = false;
                    ctx->add_rigid_motion = false;
                    // ctx->compute_singularity = false;
                    ctx->compute_singularity = true;
                }

                ierr = KSPSolve(ksp, rhs, sol);
                CHKERRQ(ierr);

                // if (default_flags_)
                // {
                //     std::cout << "Last Problem...\n";
                //     ierr = solve_last_problem();
                //     CHKERRQ(ierr);
                // }

                PetscFunctionReturn(0);
            }
        };
    } // namespace problem

    template<typename PL, typename Problem_type,
             typename Dimensions = typename PL::value_type::dimension_type>
    auto make_DtoN(
        PL &pt, Problem_type &p,
        typename std::conditional<Dimensions::value == 2, double,
                                  std::array<double, 2>>::type const &dpart)
    {
        using s_t = typename PL::value_type::shape_type;
        // using Dimensions = typename PL::value_type::dimension_type;
        return problem::DtoN<s_t, Dimensions::value, Problem_type>{pt, p,
                                                                   dpart};
    }

} // namespace cafes
#endif