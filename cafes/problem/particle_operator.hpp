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

#ifndef PROBLEM_PARTICLE_OPERATOR_HPP_INCLUDED
#define PROBLEM_PARTICLE_OPERATOR_HPP_INCLUDED

#include <algorithm/iterate.hpp>
#include <problem/problem.hpp>
#include <problem/stokes.hpp>
#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <particle/particle.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/cross_product.hpp>
#include <particle/geometry/position.hpp>
#include <particle/singularity/add_singularity.hpp>
#include <particle/geometry/vector.hpp>
#include <petsc/vec.hpp>

#include <io/vtk.hpp>


#include <petsc.h>
#include <iostream>
#include <memory>
#include <algorithm>
#include <chrono>

namespace cafes
{
  namespace problem
  {
    template<typename T>
    decltype(auto) get_index(T&& c)
    {
      return std::forward<T>(c).first;
    }

    template<typename T>
    decltype(auto) get_position(T&& c)
    {
      return std::forward<T>(c).second;
    }

    #undef __FUNCT__
    #define __FUNCT__ "set_rhs_problem_impl"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode set_rhs_problem_impl(Ctx& ctx, Vec x, petsc::petsc_vec<Dimensions>& petsc_u, bool apply_forces)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      int num = 0;
      PetscScalar const *px;
      ierr = VecGetArrayRead(x, &px);CHKERRQ(ierr);

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      for(auto& p: ctx.particles){
        auto pbox = p.bounding_box(ctx.problem.ctx->h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::overlap_box(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, ctx.problem.ctx->h);
          for(auto& ind: pts)
          {
            auto u = petsc_u.at_g(ind);
            if (apply_forces)
              for(std::size_t i=0; i<Dimensions; ++i)
                u[i] = px[num++] + p.rho_*p.force_[i];
            else
              for(std::size_t i=0; i<Dimensions; ++i)
                u[i] = px[num++];
          }       
        }
      }

      ierr = VecRestoreArrayRead(x, &px);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "set_rhs_problem"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode set_rhs_problem(Ctx& ctx, Vec x, bool apply_forces=false)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      Vec forces, mass_mult;
      ierr = DMGetGlobalVector(ctx.problem.ctx->dm, &forces);CHKERRQ(ierr);
      ierr = DMGetGlobalVector(ctx.problem.ctx->dm, &mass_mult);CHKERRQ(ierr);
      ierr = VecSet(forces, 0);CHKERRQ(ierr);

      auto petsc_forces = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, forces, 0, false);
      
      ierr = set_rhs_problem_impl(ctx, x, petsc_forces, apply_forces);CHKERRQ(ierr);

      ierr = fem::apply_mass_matrix(ctx.problem.ctx->dm, forces, mass_mult, ctx.problem.ctx->h);CHKERRQ(ierr);
      ierr = VecAXPY(ctx.problem.rhs, 1., mass_mult);CHKERRQ(ierr);

      ierr = DMRestoreGlobalVector(ctx.problem.ctx->dm, &mass_mult);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(ctx.problem.ctx->dm, &forces);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "interp_rigid_motion_"
    template<typename Shape, std::size_t Dimensions>
    PetscErrorCode interp_rigid_motion_(std::vector<particle<Shape>> const& particles,
                                        std::vector<std::vector<geometry::vector<double, Dimensions>>> const& r, 
                                        std::vector<std::vector<geometry::vector<double, Dimensions>>>& g){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;
      
      std::cout<<"add rigid motion to surf...\n";

      std::size_t ipart=0;
      for(auto& rpts: r)
      {
        for(std::size_t i=0; i<rpts.size(); ++i){
          particles[ipart].velocity_ - geometry::cross_product(particles[ipart].angular_velocity_, r[ipart][i]);
          g[ipart][i] -= particles[ipart].velocity_ - geometry::cross_product(particles[ipart].angular_velocity_, r[ipart][i]);
        }
        ipart++;
      }
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "interp_fluid_to_surf"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode interp_fluid_to_surf(Ctx& ctx, std::vector<std::vector<geometry::vector<double, Dimensions>>>& g,
                                        bool rigid_motion = false, bool singularity = false)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.sol, 0);

      ierr = sol.global_to_local(INSERT_VALUES);CHKERRQ(ierr);

      std::size_t ipart=0;
      for(auto& spts: ctx.surf_points){
        for(std::size_t i=0; i<spts.size(); ++i){
          auto bfunc = fem::P1_integration(get_position(spts[i]), ctx.problem.ctx->h);
          auto ielem = fem::get_element(get_index(spts[i]));
          
          std::fill(g[ipart][i].begin(), g[ipart][i].end(), 0.);
          
          for (std::size_t j=0; j<bfunc.size(); ++j){
            auto u = sol.at(ielem[j]);
            for (std::size_t d=0; d<Dimensions; ++d)
              g[ipart][i][d] += u[d]*bfunc[j];
          }
        }
        ipart++;
      }

      if (rigid_motion)
      { 
        ierr = interp_rigid_motion_(ctx.particles, ctx.radial_vec, g);CHKERRQ(ierr);
      }
      // if (singularity)
      // { 
      //   ierr = singularity::add_singularity_to_surf<Dimensions, Ctx>(ctx, g);CHKERRQ(ierr);
      // }

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "simple_layer"
    template<std::size_t Dimensions, typename Ctx, typename cross_type>
    PetscErrorCode simple_layer(Ctx& ctx, std::vector<std::vector<geometry::vector<double, Dimensions>>>& g,
                                std::vector<geometry::vector<double, Dimensions>>& mean,
                                cross_type& cross_prod)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      for(std::size_t ipart=0; ipart<g.size(); ++ipart)
        for(std::size_t isurf=0; isurf<g[ipart].size(); ++isurf)
        {
          mean[ipart] += g[ipart][isurf];
          cross_prod[ipart] += geometry::cross_product(ctx.radial_vec[ipart][isurf], g[ipart][isurf]);
        }

      MPI_Allreduce(MPI_IN_PLACE, mean.data(), mean.size()*Dimensions, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, cross_prod.data(), cross_prod.size()*(Dimensions==2?1:3), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      
      for(std::size_t ipart=0; ipart<g.size(); ++ipart)
      {
        mean[ipart] /= ctx.nb_surf_points[ipart];
        cross_prod[ipart] *= ctx.particles[ipart].Cd_R()/ctx.nb_surf_points[ipart];
      }

      for(std::size_t ipart=0; ipart<g.size(); ++ipart)
      {
        for(std::size_t isurf=0; isurf<g[ipart].size(); ++isurf)
        {
          g[ipart][isurf] -= mean[ipart] + geometry::cross_product(cross_prod[ipart], ctx.radial_vec[ipart][isurf]);
          
        }
      }

      PetscFunctionReturn(0);
    }
    #undef __FUNCT__
    #define __FUNCT__ "SL_to_Rhs"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode SL_to_Rhs(Ctx& ctx, std::vector<std::vector<geometry::vector<double, Dimensions>>>const& g){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = VecSet(ctx.problem.rhs, 0.);CHKERRQ(ierr);

      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.rhs, 0, false);
      ierr = sol.fill(0.);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<g.size(); ++ipart){ 
        auto& radius = ctx.particles[ipart].shape_factors_[0];
        //auto gammak = ctx.particles[ipart].perimeter/ctx.nb_surf_points[ipart];  
        // remove this line !!
        auto gammak = ctx.particles[ipart].surface_area()/ctx.nb_surf_points[ipart];
        for(std::size_t isurf=0; isurf<ctx.surf_points[ipart].size(); ++isurf){
          auto bfunc = fem::P1_integration(get_position(ctx.surf_points[ipart][isurf]), ctx.problem.ctx->h);
          auto ielem = fem::get_element(get_index(ctx.surf_points[ipart][isurf]));
          for (std::size_t j=0; j<bfunc.size(); ++j){
            auto u = sol.at(ielem[j]);
            for (std::size_t d=0; d<Dimensions; ++d)
              u[d] += g[ipart][isurf][d]*bfunc[j]*gammak;
          }
        }
      }

      // if (ctx.compute_singularity)
      // { 
      //   ierr = singularity::add_singularity_to_surf<Dimensions, Ctx>(ctx, sol);CHKERRQ(ierr);
      // }

      ierr = sol.local_to_global(ADD_VALUES);CHKERRQ(ierr);

      auto petsc_rhs = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.rhs, 0);
      ierr = SetNullDirichletOnRHS(petsc_rhs, ctx.problem.ctx->bc_);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "init_problem"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode init_problem(Ctx& ctx, Vec x, bool apply_forces=false){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      if (ctx.compute_rhs)
      {
        ierr = ctx.problem.setup_RHS();CHKERRQ(ierr);
      }
      if (ctx.compute_singularity)
      {
        ierr = singularity::add_singularity_in_fluid<Dimensions, Ctx>(ctx);CHKERRQ(ierr);
      }

      ierr = set_rhs_problem<Dimensions, Ctx>(ctx, x, apply_forces);CHKERRQ(ierr);
      
      auto petsc_rhs = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.rhs, 0);

      if (ctx.compute_rhs)
      {
        ierr = SetDirichletOnRHS(petsc_rhs, ctx.problem.ctx->bc_, ctx.problem.ctx->h);CHKERRQ(ierr);
      }
      else
      {
        ierr = SetNullDirichletOnRHS(petsc_rhs, ctx.problem.ctx->bc_);CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }


    auto const kernel_projection = [](auto const& p, auto& sol, 
                            auto const& h, auto const& hs, 
                            auto& box_scale, auto& mean, auto& cross)
    {
      auto const kernel_pos = [&](auto const& pos){
        auto const kernel = [&](auto const& pos_scale)
        {
          using mean_type = typename std::remove_reference<decltype(mean)>::type;
          auto pts = pos*h + pos_scale*hs;
          if (p.contains(pts))
          {
            auto pts_corner = pos_scale*hs;
            auto bfunc = fem::P1_integration(pts_corner, h);
            auto ielem = fem::get_element(pos);
            mean_type tmp{};

            for(std::size_t ib=0; ib<bfunc.size(); ++ib)
            {
              auto u = sol.at(ielem[ib]);
              for (std::size_t d=0; d<pos.dimensions; ++d)
                tmp[d] += u[d]*bfunc[ib];
            }

            mean_type r;
            for (std::size_t d=0; d<pos.dimensions; ++d)
              r[d] = pts[d] - p.center_[d];

            mean += tmp;
            cross += geometry::cross_product(r, tmp);
          }
        };
        algorithm::iterate(box_scale, kernel);
      };
      return kernel_pos;
    };

    #undef __FUNCT__
    #define __FUNCT__ "projection_impl"
    template<typename Shape, typename cross_type, std::size_t Dimensions>
    PetscErrorCode projection_impl(std::vector<particle<Shape>> const& particles,
                              petsc::petsc_vec<Dimensions>& sol,
                              geometry::box<int, Dimensions> const& box,
                              std::vector<geometry::vector<double, Dimensions>>& mean,
                              cross_type& cross_prod,
                              std::vector<int> const& num,
                              std::size_t scale,
                              std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;


      std::array<double, Dimensions> hs;
      for (std::size_t d=0; d<Dimensions; ++d)
        hs[d] = h[d]/scale;

      geometry::position<std::size_t, Dimensions> p1, p2;
      p1.fill(0);
      p2.fill(scale);
      geometry::box<std::size_t, Dimensions> box_scale{ p1, p2};

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart){
        auto& p = particles[ipart];
        mean[ipart] = 0;
        cross_prod[ipart] = 0;
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          algorithm::iterate(new_box, kernel_projection(p, sol, h, hs, box_scale, mean[ipart], cross_prod[ipart]));
        }
        mean[ipart] /= num[ipart];
        cross_prod[ipart] /= num[ipart];
      }
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "projection"
    template<std::size_t Dimensions, typename Ctx, typename cross_type>
    PetscErrorCode projection(Ctx& ctx,
                              std::vector<geometry::vector<double, Dimensions>>& mean,
                              cross_type& cross_prod)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      using position_type = geometry::position<double, Dimensions>;
      using position_type_i = geometry::position<int, Dimensions>;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.sol, 0);
      sol.global_to_local(INSERT_VALUES);

      ierr = projection_impl(ctx.particles, sol, box, mean, cross_prod, ctx.num, ctx.scale, h);

      MPI_Allreduce(MPI_IN_PLACE, mean.data(), mean.size()*Dimensions, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, cross_prod.data(), cross_prod.size()*(Dimensions==2?1:3), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "compute_y"
    template<std::size_t Dimensions, typename Ctx, typename cross_type>
    PetscErrorCode compute_y(Ctx& ctx, Vec y, 
                             std::vector<geometry::vector<double, Dimensions>> const& mean, 
                             cross_type const& cross_prod)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;
      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.sol, 0);

      int num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<ctx.particles.size(); ++ipart){
        auto& p = ctx.particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts)
          {
            geometry::vector<double, Dimensions> r;
            for (std::size_t d=0; d<Dimensions; ++d)
              r[d] = ind[d]*h[d] - p.center_[d];
            auto tmp = mean[ipart] + p.Ci_R()*geometry::cross_product(cross_prod[ipart], r);
            auto usol = sol.at(ind);
            for (std::size_t d=0; d<Dimensions; ++d)
              py[num++] = usol[d] - tmp[d];
          }
        }
      }

      ierr = VecRestoreArray(y, &py);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "compute_y"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode compute_y(Ctx& ctx, Vec y)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;


      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.sol, 0);

      int num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<ctx.particles.size(); ++ipart){
        auto& p = ctx.particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts)
          {
            auto usol = sol.at(ind);
            for (std::size_t d=0; d<Dimensions; ++d)
              py[num++] = usol[d];
          }
        }
      }

      ierr = VecRestoreArray(y, &py);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }
}

#endif