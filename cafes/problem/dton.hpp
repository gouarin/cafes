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
#include <problem/options.hpp>

//#include <io/vtk.hpp>

#include <algorithm>
#include <iostream>
#include <memory>
#include <petsc.h>

namespace cafes
{
    namespace problem
    {

template <typename T, typename mesh_size> 
struct convergence_context
{
    Vec x0, x1;
    DM dm;
    mesh_size h;
    T particles;
};

template <typename T, typename mesh_size, std::size_t Dimensions, typename Ctx, typename bound_mesh_size> 
struct convergence_context2
{
    using position_type = geometry::position<double, Dimensions>;
    using position_type_i = geometry::position<int, Dimensions>;

    Vec sol, sol_rhs;
    bool compute_singularity;
    DM dm;
    mesh_size h;
    T particles;
    std::vector<std::vector<std::pair<position_type_i, position_type>>> surf_points;
    Ctx *ctx;
    bound_mesh_size dpart;
};

#undef __FUNCT__
#define __FUNCT__ "DtoN_matrix"
        template<std::size_t Dimensions, typename Ctx>
        PetscErrorCode DtoN_matrix(Mat A, Vec x, Vec y)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            std::cout << "dton x vecttor\n";
            VecView(x, PETSC_VIEWER_STDOUT_WORLD);

            Ctx *ctx;
            ierr = MatShellGetContext(A, (void **)&ctx);CHKERRQ(ierr);

            ierr = VecSet(ctx->problem.rhs, 0.);CHKERRQ(ierr);

            ierr = init_problem<Dimensions, Ctx>(*ctx, x);CHKERRQ(ierr);

            ierr = ctx->problem.solve();CHKERRQ(ierr);

            // std::cout << "solution\n";
            // VecView(ctx->problem.sol, PETSC_VIEWER_STDOUT_WORLD);
            if (!ctx->compute_rhs)
            {
                ierr = cafes::io::save_hdf5("Resultats", "two_part_ug_rhs",
                                            ctx->problem.rhs, ctx->problem.ctx->dm,
                                            ctx->problem.ctx->h);
                CHKERRQ(ierr);
                ierr = cafes::io::save_hdf5("Resultats", "two_part_ug_sol",
                                            ctx->problem.sol, ctx->problem.ctx->dm,
                                            ctx->problem.ctx->h);
                CHKERRQ(ierr);
                // std::exit(0);
            }

            ierr = VecCopy(ctx->problem.sol, ctx->sol_tmp);CHKERRQ(ierr);

            std::vector<std::vector<geometry::vector<double, Dimensions>>> g;
            g.resize(ctx->particles.size());
            for (std::size_t ipart = 0; ipart < ctx->surf_points.size(); ++ipart)
            {
                g[ipart].resize(ctx->surf_points[ipart].size());
            }

            // interpolation
            ierr = interp_fluid_to_surf(*ctx, g, ctx->add_rigid_motion,
                                        ctx->compute_singularity);
            CHKERRQ(ierr);

            // for (std::size_t ipart = 0; ipart < ctx->surf_points.size(); ++ipart)
            // {
            //     for(auto e: g[ipart])
            //     {
            //         std::cout << "g " << e << "\n";
            //     }
            // }

            ierr = SL_to_Rhs(*ctx, g);CHKERRQ(ierr);

            // std::cout << "rhs\n";
            // VecView(ctx->problem.rhs, PETSC_VIEWER_STDOUT_WORLD);

            ierr = ctx->problem.solve();CHKERRQ(ierr);

            ierr = compute_y<Dimensions, Ctx>(*ctx, y);CHKERRQ(ierr);
            std::cout << "dton y vecttor\n";
            ierr = VecView(y, PETSC_VIEWER_STDOUT_WORLD);
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
            Vec sol, sol_previous;
            //convergence_context<decltype(parts_), decltype(problem_.ctx->h)> *conv_ctx;// {NULL, NULL, problem_.ctx->dm, parts_};
            Vec rhs;
            Vec sol_rhs, sol_reg, sol_tmp, sol_last;
            Mat A;
            KSP ksp;
            PetscInt kspiter;
            PetscReal kspresnorm;
            std::size_t scale_ = 4;
            bool default_flags_ = true;
            bool use_sing = false;

            using dpart_type =
                typename std::conditional<Dimensions == 2, double,
                                          std::array<double, 2>>::type;
            dpart_type dpart_;
            convergence_context2<decltype(parts_), decltype(problem_.ctx->h), Dimensions, Ctx, decltype(dpart_)> *conv_ctx;

            DtoN(std::vector<particle<Shape>> &parts, Problem_type &p,
                 dpart_type dpart)
                : parts_{parts}, problem_{p}, dpart_{dpart}
            {
                problem_.setup_KSP();
                VecDuplicate(problem_.sol, &sol_tmp);
                VecDuplicate(problem_.sol, &sol_rhs);
                VecDuplicate(problem_.sol, &sol_reg);
            }

#undef __FUNCT__
#define __FUNCT__ "create_Mat_and_Vec"
            PetscErrorCode create_Mat_and_Vec()
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                auto box = fem::get_DM_bounds<Dimensions>(problem_.ctx->dm, 1);
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
                ierr = VecDuplicate(sol_tmp, &sol_previous);
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "setup_RHS"
            virtual PetscErrorCode setup_RHS() override
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;
                options<Dimensions> opt{};
                opt.process_options();
                std::cout << "COMPUTE_SINGULARITY = " << opt.compute_singularity << std::endl;

                if (default_flags_)
                {
                    ctx->compute_rhs = true;
                    ctx->add_rigid_motion = true;
                    ctx->compute_singularity = opt.compute_singularity;
                    // ctx->compute_singularity = false;
                }

                ierr = MatMult(A, sol, rhs);
                CHKERRQ(ierr);
                ierr = VecScale(rhs, -1.);
                CHKERRQ(ierr);

                ierr = VecCopy(sol_tmp, sol_rhs);
                CHKERRQ(ierr);

                // ierr = VecSet(rhs,0.);
                // CHKERRQ(ierr);
                // ierr = VecSet(sol_rhs,0.);
                // CHKERRQ(ierr);

                // ierr = cafes::io::save_VTK("Resultats", "two_part_u0", sol_tmp,
                //                            problem_.ctx->dm, problem_.ctx->h);
                // CHKERRQ(ierr);
                // ierr = cafes::io::save_hdf5("Resultats", "two_part_w0",
                //                            problem_.sol, problem_.ctx->dm,
                //                            problem_.ctx->h);
                // CHKERRQ(ierr);
                // ierr = cafes::io::save_hdf5("Resultats", "two_part_w0_rhs",
                //                            problem_.rhs, problem_.ctx->dm,
                //                            problem_.ctx->h);
                // CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

static PetscErrorCode convergeTest(KSP ksp, PetscInt it, PetscReal rnorm, KSPConvergedReason* reason, void* ctx)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    convergence_context<decltype(parts_), decltype(problem_.ctx->h)> *conv_ctx;
    KSPGetConvergenceContext(ksp, (void**) &conv_ctx);
    double norm, norm0;
    VecNorm(conv_ctx->x0, NORM_2, &norm0);
    DM dav, dap;
    DMDALocalInfo infop, infov;
    
    // Get dmda local information
    ierr = DMCompositeGetEntries(conv_ctx->dm, &dav, &dap);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dav, &infov);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);

    // Difference between to successive fluid velocities
    Vec diff;
    VecDuplicate(conv_ctx->x0, &diff);

    //VecCopy(conv_ctx->x0, diff);
    //VecAXPY(diff, -1, conv_ctx->x1);
    //VecNorm(diff, NORM_2, &norm);
    //VecView(conv_ctx->x0, PETSC_VIEWER_STDOUT_WORLD);
    //std::cout << norm << " " << it << "\n";

    // Extract fluid velocities and pressure from diff
    auto diffv = petsc::petsc_vec<Dimensions>(conv_ctx->dm, diff, 0, false);
    auto x0v = petsc::petsc_vec<Dimensions>(conv_ctx->dm, conv_ctx->x0, 0, true);
    auto x1v = petsc::petsc_vec<Dimensions>(conv_ctx->dm, conv_ctx->x1, 0, true);
    auto diffp = petsc::petsc_vec<Dimensions>(conv_ctx->dm, diff, 1, false);
    auto x0p = petsc::petsc_vec<Dimensions>(conv_ctx->dm, conv_ctx->x0, 1, true);
    auto x1p = petsc::petsc_vec<Dimensions>(conv_ctx->dm, conv_ctx->x1, 1, true);

    // Put zeros in diff velocity if point is inside particles
    for (std::size_t j =0; j< infov.my; j++ )
    {
        for (std::size_t i = 0; i <infov.mx; i++)
        {
            auto pos = cafes::geometry::position<int, Dimensions>{i,j};
            auto pts = cafes::geometry::position<double, Dimensions>{i*conv_ctx->h[0], j*conv_ctx->h[1]};
            auto u = diffv.at_g(pos);
            auto x0 = x0v.at_g(pos);
            auto x1 = x1v.at_g(pos);
            for(auto& p: conv_ctx->particles){
                if (p.contains(pts))
                {
                    for(std::size_t k=0; k<Dimensions; ++k)
                    {
                        u[k] = 0.;
                    }
                }
                else
                {
                    for(std::size_t k=0; k<Dimensions; ++k)
                    {
                        u[k] = x0[k] - x1[k];
                    }
                }
            }
        }
    }

    // Put zeros in diff pressure if point is inside particles
    for (std::size_t j =0; j< infop.my; j++ )
    {
        for (std::size_t i = 0; i <infop.mx; i++)
        {
            auto pos = cafes::geometry::position<int, Dimensions>{i,j};
            auto pts = cafes::geometry::position<double, Dimensions>{2*i*conv_ctx->h[0], 2*j*conv_ctx->h[1]};
            auto u = diffp.at_g(pos);
            u[0] = 0.;
            // auto x0 = x0p.at_g(pos);
            // auto x1 = x1p.at_g(pos);
            // for(auto& p: conv_ctx->particles){
            //     if (p.contains(pts))
            //     {
            //         u[0] = 0.;
            //     }
            //     else
            //     {
            //         u[0] = x0[0] - x1[0];
            //     }
            // }
        }
    }

    VecNorm(diff, NORM_2, &norm);
    if(norm0>1e-14)
    {
        norm /= norm0;
    }
    //VecView(conv_ctx->x0, PETSC_VIEWER_STDOUT_WORLD);
    std::cout << norm << " " << it << "\n";
    if(it == 0)
    {
        *reason = KSP_CONVERGED_ITERATING;
    }
    else
    {
        if(norm < 1e-10)
        {
            *reason = KSP_CONVERGED_RTOL_NORMAL;
        }
        else
        {
            *reason = KSP_CONVERGED_ITERATING;
        }
    }
    VecCopy(conv_ctx->x0, conv_ctx->x1);
    PetscFunctionReturn(0);
}

static PetscErrorCode convergeTest2(KSP ksp, PetscInt it, PetscReal rnorm, KSPConvergedReason* reason, void* ctx)
{
    PetscErrorCode ierr;
    PetscFunctionBegin;
    convergence_context2<decltype(parts_), decltype(problem_.ctx->h), Dimensions, Ctx, decltype(dpart_)> *conv_ctx;
    KSPGetConvergenceContext(ksp, (void**) &conv_ctx);
    double cost = 0, testperimetre = 0; 
    Vec sol, sol_reg;
    VecDuplicate(conv_ctx->sol, &sol);
    VecDuplicate(conv_ctx->sol_rhs, &sol_reg);
    ierr = KSPBuildSolution(ksp, NULL, &sol);
    ierr = init_problem<Dimensions, Ctx>(*conv_ctx->ctx, sol);
    CHKERRQ(ierr);
    ierr = conv_ctx->ctx->problem.solve();
    CHKERRQ(ierr);
    ierr = VecCopy(conv_ctx->ctx->problem.sol, sol_reg);
    CHKERRQ(ierr);
    // std::string stout1 = "two_part_test_interp_";
    // const char * stw1 = stout1.c_str();
    // ierr = cafes::io::save_VTK("Resultats", stw1, sol_reg, conv_ctx->dm, conv_ctx->h);CHKERRQ(ierr);
    VecAXPY(sol_reg, 1., conv_ctx->sol_rhs);
    auto p1 = conv_ctx->particles[0];
    using shape_type = typename decltype(p1)::shape_type;
    auto sing = cafes::singularity::singularity<shape_type, Dimensions>(conv_ctx->particles[0], conv_ctx->particles[1], conv_ctx->h[0]);
    // if (conv_ctx->compute_singularity)
    // {
    //     cafes::singularity::add_singularity_to_ureg(conv_ctx->dm, conv_ctx->h, sol_reg, conv_ctx->particles);
    // }

    // auto solreg = petsc::petsc_vec<Dimensions>(conv_ctx->dm, sol_reg, 0);

    // ierr = solreg.global_to_local(ADD_VALUES);CHKERRQ(ierr);
    // double U = 1.;
    // std::size_t ipart=0;
    // for(auto& spts: conv_ctx->surf_points){
    //     auto nb_surf_points = conv_ctx->surf_points[ipart].size();
    //     auto gammak = conv_ctx->particles[ipart].surface_area()/nb_surf_points;
    //     for(std::size_t i=0; i<spts.size(); ++i){
    //         auto bfunc = fem::P1_integration(get_position(spts[i]), conv_ctx->h);
    //         auto ielem = fem::get_element(get_index(spts[i]));
    //         double tmpx = 0., tmpy = 0.;
    //         for (std::size_t j=0; j<bfunc.size(); ++j){
    //             auto u = solreg.at(ielem[j]);
    //             tmpx += u[0]*bfunc[j];
    //             tmpy += u[1]*bfunc[j];
    //         }
    //         cost += ((tmpx-U)*(tmpx-U) + tmpy*tmpy)*gammak;
    //     }
    //     U = -1.;
    //     ipart++;
    // }

    auto solreg = petsc::petsc_vec<Dimensions>(conv_ctx->dm, sol_reg, 0);
    ierr = solreg.global_to_local(INSERT_VALUES);CHKERRQ(ierr);
    
    // PARTICLE 1
    for(std::size_t ipart=0; ipart<1; ++ipart)
    {
        //std::size_t ipart = 0; 
        double U = 1.;
        if (ipart == 1){U = -1.;}
        auto spts = conv_ctx->particles[ipart].surface(conv_ctx->dpart);
        auto nb_surf_points = conv_ctx->surf_points[ipart].size();
        auto gammak = conv_ctx->particles[ipart].surface_area()/nb_surf_points;
        for(std::size_t isurf=0; isurf<nb_surf_points; ++isurf){
            auto ind = get_index(conv_ctx->surf_points[ipart][isurf]);
            auto pos = spts[isurf];
            pos[0] -= ind[0]*conv_ctx->h[0];
            pos[1] -= ind[1]*conv_ctx->h[1];
            auto bfunc = fem::P1_integration(get_position(conv_ctx->surf_points[ipart][isurf]), conv_ctx->h); //conv_ctx->surf_points[ipart][isurf]), conv_ctx->h);
            auto ielem = fem::get_element(get_index(conv_ctx->surf_points[ipart][isurf]));
            auto Using = sing.get_u_sing(spts[isurf]);
            double tmpx = Using[0];
            double tmpy = Using[1];
            for (std::size_t j=0; j<bfunc.size(); ++j){
                auto u = solreg.at(ielem[j]);
                tmpx += u[0]*bfunc[j];
                tmpy += u[1]*bfunc[j];
            }
            //auto pos = get_position(conv_ctx->surf_points[ipart][isurf]);
            //cout << "ind : " << ind << " --> ux : " << tmpx  << " particle : " << spts[isurf] << " rel pos : " << pos << "\n";
            //cout << "ind : " << ind << " --> ux : " << tmpx << "   uy : " << tmpy << "\n";
            cost += ((tmpx-U)*(tmpx-U) + tmpy*tmpy)*gammak;
            //cost += std::abs(tmpx)*gammak;
            testperimetre += gammak;
        }
    }
    cost = std::sqrt(cost);
    //auto velocity = petsc::petsc_vec<Dimensions>(conv_ctx->dm, conv_ctx->sol, 0, false);
    // Vec sol, diff;
    // VecDuplicate(conv_ctx->sol_previous, &sol);
    // VecDuplicate(conv_ctx->sol_previous, &diff);
    // ierr = KSPBuildSolution(ksp, NULL, &sol);
    // VecCopy(sol, diff);
    // VecAXPY(diff, -1., conv_ctx->sol_previous);
    // VecNorm(diff, NORM_2, &norm);
    // VecCopy(sol, conv_ctx->sol_previous);

    const double pi = 3.14159265358979323846;
    double perimetre = 2*pi*conv_ctx->particles[0].shape_factors_[0];
    //std::cout << "--> integrale particle 1 : " << testperimetre << ", vrai perimetre : " <<  perimetre << "\n";
    std::cout << "--> COST : " << cost << "\n";
    if(it == 0)
    {
        *reason = KSP_CONVERGED_ITERATING;
    }
    else
    {
        if(cost < 1e-6)
        {
            *reason = KSP_CONVERGED_RTOL_NORMAL;
        }
        else
        {
            *reason = KSP_CONVERGED_ITERATING;
        }
    }
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
                //conv_ctx = new convergence_context<decltype(parts_), decltype(problem_.ctx->h)> {sol_tmp, sol_previous, problem_.ctx->dm, problem_.ctx->h, parts_};
                //ierr = KSPSetConvergenceTest(ksp, convergeTest, conv_ctx, PETSC_NULL);
                //conv_ctx = new convergence_context2<decltype(parts_), decltype(problem_.ctx->h), Dimensions, Ctx, decltype(dpart_)> {sol, sol_rhs, ctx->compute_singularity, problem_.ctx->dm, problem_.ctx->h, parts_, surf_points_, ctx, dpart_};
                //ierr = KSPSetConvergenceTest(ksp, convergeTest2, conv_ctx, PETSC_NULL);
                //CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "test"
            PetscErrorCode test()
            {
                PetscErrorCode ierr;
                PetscFunctionBegin;
                // solve the problem with the right control

                ierr = VecSet(ctx->problem.rhs, 0.);
                CHKERRQ(ierr);

                if (default_flags_)
                {
                    ctx->compute_rhs = true;
                    ctx->add_rigid_motion = false;
                    ctx->compute_singularity = false;
                }

                VecSet(sol, 0.);
                ierr = init_problem<Dimensions, Ctx>(*ctx, sol);
                CHKERRQ(ierr);

                ierr = ctx->problem.solve();
                CHKERRQ(ierr);
                VecView(ctx->problem.sol, PETSC_VIEWER_STDOUT_WORLD);
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

                ierr = ctx->problem.solve();
                CHKERRQ(ierr);
                ierr = VecDuplicate(ctx->problem.sol, &sol_last);
                CHKERRQ(ierr);
                ierr = VecCopy(ctx->problem.sol, sol_last);
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
                    ctx->compute_singularity = false;
                }
                //ierr = VecSet(sol, 0.1);
                //CHKERRQ(ierr);
                ierr = KSPSolve(ksp, rhs, sol);
                CHKERRQ(ierr);
                ierr = init_problem<Dimensions, Ctx>(*ctx, sol);
                CHKERRQ(ierr);
                ierr = ctx->problem.solve();
                CHKERRQ(ierr);
                ierr = VecCopy(ctx->problem.sol, sol_tmp);
                CHKERRQ(ierr);
                VecCopy(sol_tmp, sol_reg);
                VecAXPY(sol_reg, 1., sol_rhs);
                ierr = KSPGetIterationNumber(ksp, &kspiter);CHKERRQ(ierr);
                ierr = KSPGetResidualNorm(ksp, &kspresnorm);CHKERRQ(ierr);

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