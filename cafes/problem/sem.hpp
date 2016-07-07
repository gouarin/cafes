// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
// OF SUCH DAMAGE.

#ifndef PARTICLE_PROBLEM_SEM_HPP_INCLUDED
#define PARTICLE_PROBLEM_SEM_HPP_INCLUDED

#include <problem/particle_operator.hpp>
#include <problem/problem.hpp>
#include <problem/stokes.hpp>
#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <particle/particle.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>
#include <particle/geometry/vector.hpp>

#include <petsc.h>
#include <iostream>
#include <memory>
#include <algorithm>

namespace cafes
{
  namespace problem
  {

    #undef __FUNCT__
    #define __FUNCT__ "sem_matrix"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode sem_matrix(Mat A, Vec x, Vec y){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      Ctx *ctx;
      ierr = MatShellGetContext(A, (void**) &ctx);CHKERRQ(ierr);

      ierr = VecSet(ctx->problem.rhs, 0.);CHKERRQ(ierr);

      ierr = init_problem<Dimensions, Ctx>(*ctx, x, ctx->compute_rhs);CHKERRQ(ierr);
      //VecView(ctx->problem.rhs, 0);
      ierr = ctx->problem.solve();CHKERRQ(ierr);

      std::vector<std::vector<geometry::vector<double, Dimensions>>> g;
      g.resize(ctx->particles.size());
      for(std::size_t ipart=0; ipart<ctx->surf_points.size(); ++ipart)
        g[ipart].resize(ctx->surf_points[ipart].size());

      // interpolation
      ierr = interp_fluid_to_surf(*ctx, g);CHKERRQ(ierr);

      std::vector<geometry::vector<double, Dimensions>> mean(ctx->particles.size());
      using cross_type  = typename std::conditional<Dimensions==2,
                                                    double, 
                                                    geometry::vector<double, 3>>::type;
      std::vector<cross_type> cross_prod(ctx->particles.size());

      ierr = simple_layer(*ctx, g, mean, cross_prod);CHKERRQ(ierr);

      ierr = SL_to_Rhs(*ctx, g);CHKERRQ(ierr);
      
      ierr = ctx->problem.solve();CHKERRQ(ierr);

      ierr = projection<Dimensions, Ctx>(*ctx, mean, cross_prod);CHKERRQ(ierr);

      ierr = compute_y(*ctx, y, mean, cross_prod);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


    template<typename Shape, std::size_t Dimensions, typename Problem_type>
    struct SEM : public Problem<Dimensions>
    {
      std::vector<particle<Shape>> parts_;
      Problem_type problem_;

      using position_type   = geometry::position<double, Dimensions>;
      using position_type_i = geometry::position<int, Dimensions>;

      using Ctx = particle_context<Dimensions, Shape, Problem_type>;
      Ctx *ctx;

      std::vector<std::vector<std::pair<position_type_i, position_type>>> surf_points_;
      std::vector<std::vector<geometry::vector<double, Dimensions>>> radial_vec_;
      std::vector<int> nb_surf_points_;
      std::vector<int> num_;
      Vec sol;
      Vec rhs;
      Mat A;
      KSP ksp;
      std::size_t scale_ = 4;

      using dpart_type = typename std::conditional<Dimensions == 2, 
                                  double, 
                                  std::array<double, 2>>::type;
      dpart_type dpart_; 

      SEM(std::vector<particle<Shape>>const& parts, Problem_type& p, dpart_type dpart):
      parts_{parts}, problem_{p}, dpart_{dpart}
      {
        problem_.setup_KSP();
      }
      
      #undef __FUNCT__
      #define __FUNCT__ "create_Mat_and_Vec"
      PetscErrorCode create_Mat_and_Vec()
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;

        auto box = fem::get_DM_bounds<Dimensions>(problem_.ctx->dm, 0);
        auto& h = problem_.ctx->h;

        auto size = set_materials(parts_, surf_points_, radial_vec_,
                                  nb_surf_points_, num_, box,
                                  h, dpart_, scale_);

        ctx = new Ctx{problem_, parts_, surf_points_, radial_vec_, nb_surf_points_, num_, scale_, false, false};

        ierr = MatCreateShell(PETSC_COMM_WORLD, size*Dimensions, size*Dimensions, PETSC_DECIDE, PETSC_DECIDE, ctx, &A);CHKERRQ(ierr);
        ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))sem_matrix<Dimensions, Ctx>);CHKERRQ(ierr);

        ierr = MatCreateVecs(A, &sol, &rhs);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      
      #undef __FUNCT__
      #define __FUNCT__ "setup_RHS"
      virtual PetscErrorCode setup_RHS() override 
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;

        ctx->compute_rhs = true;

        ierr = MatMult(A, sol, rhs);CHKERRQ(ierr);
        ierr = VecScale(rhs, -1.);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "setup_KSP"
      virtual PetscErrorCode setup_KSP() override
      {
        PetscErrorCode ierr;
        PC             pc;
        PetscFunctionBegin;

        ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
        ierr = KSPSetOptionsPrefix(ksp, "sem_");CHKERRQ(ierr);

        ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
        ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
        ierr = PCSetType(pc, PCNONE);CHKERRQ(ierr);
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "solve"
      virtual PetscErrorCode solve() override
      {
        PetscErrorCode ierr;
        PetscFunctionBegin;

        ctx->compute_rhs = false;

        ierr = KSPSolve(ksp, rhs, sol);CHKERRQ(ierr);

        // solve the problem with the right control
        ctx->compute_rhs = true;
        ierr = init_problem<Dimensions, Ctx>(*ctx, sol, true);CHKERRQ(ierr);
        ierr = ctx->problem.solve();CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "get_new_velocities"
      PetscErrorCode get_new_velocities()
      {
        PetscErrorCode ierr;
        PetscFunctionBegin;

        std::vector<std::vector<geometry::vector<double, Dimensions>>> g;
        g.resize(ctx->particles.size());
        for(std::size_t ipart=0; ipart<ctx->surf_points.size(); ++ipart)
          g[ipart].resize(ctx->surf_points[ipart].size());

        // interpolation
        ierr = interp_fluid_to_surf(*ctx, g);CHKERRQ(ierr);

        std::vector<geometry::vector<double, Dimensions>> mean(ctx->particles.size());
        using cross_type  = typename std::conditional<Dimensions==2,
                                                      double, 
                                                      geometry::vector<double, 3>>::type;
        std::vector<cross_type> cross_prod(ctx->particles.size());

        // simplify this part: we only need to compute mean to set 
        // the new velocity
        ierr = simple_layer(*ctx, g, mean, cross_prod);CHKERRQ(ierr);

        for(std::size_t ipart=0; ipart<ctx->particles.size(); ++ipart)
          for(std::size_t d=0; d<Dimensions; ++d)
            parts_[ipart].velocity_[d] = mean[ipart][d];
          
        PetscFunctionReturn(0);
      }
    };
  }

  template<typename PL, typename Problem_type, typename Dimensions = typename PL::value_type::dimension_type> 
  auto make_SEM(PL const& pt, Problem_type& p, 
                typename std::conditional<Dimensions::value == 2, double, 
                                  std::array<double, 2>>::type const& dpart)
  {
    using s_t = typename PL::value_type::shape_type;
    //using Dimensions = typename PL::value_type::dimension_type;
    return problem::SEM<s_t, Dimensions::value, Problem_type>{pt, p, dpart};
  }

}
#endif