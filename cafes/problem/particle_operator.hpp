#ifndef PROBLEM_PARTICLE_OPERATOR_HPP_INCLUDED
#define PROBLEM_PARTICLE_OPERATOR_HPP_INCLUDED

#include <problem/problem.hpp>
#include <problem/stokes.hpp>
#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <particle/particle.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/cross_product.hpp>
#include <particle/geometry/position.hpp>
#include <petsc/vec.hpp>

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
    #define __FUNCT__ "set_rhs_problem_"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode set_rhs_problem_(Ctx& ctx, Vec x, Vec out, bool apply_forces, std::integral_constant<int, 2>)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;      

      Vec outu;
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, out, &outu, NULL);CHKERRQ(ierr);

      int num = 0;
      PetscScalar const *px;
      ierr = VecGetArrayRead(x, &px);CHKERRQ(ierr);

      DM dau;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm , &dau, PETSC_NULL);CHKERRQ(ierr);

      PetscScalar ***pout;
      ierr = DMDAVecGetArrayDOF(dau, outu, &pout);CHKERRQ(ierr);
      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      for(auto& p: ctx.particles){
        auto pbox = p.bounding_box(ctx.problem.ctx->h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::overlap_box(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, ctx.problem.ctx->h);
          for(auto& ind: pts){
            if (apply_forces)
              for(std::size_t i=0; i<2; ++i)
                pout[ind[1]][ind[0]][i] = px[num++] + p.rho_*p.force_[i];
            else
              for(std::size_t i=0; i<2; ++i)
                pout[ind[1]][ind[0]][i] = px[num++];
          }       
        }
      }

      ierr = DMDAVecRestoreArrayDOF(dau, outu, &pout);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, out, &outu, NULL);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(x, &px);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "set_rhs_problem_"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode set_rhs_problem_(Ctx& ctx, Vec x, Vec out, bool apply_forces, std::integral_constant<int, 3>)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;      

      Vec outu;
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, out, &outu, NULL);CHKERRQ(ierr);

      int num = 0;
      PetscScalar const *px;
      ierr = VecGetArrayRead(x, &px);CHKERRQ(ierr);

      DM dau;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm , &dau, PETSC_NULL);CHKERRQ(ierr);

      PetscScalar ****pout;
      ierr = DMDAVecGetArrayDOF(dau, outu, &pout);CHKERRQ(ierr);
      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      for(auto& p: ctx.particles){
        auto pbox = p.bounding_box(ctx.problem.ctx->h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::overlap_box(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, ctx.problem.ctx->h);
          for(auto& ind: pts){
            if (apply_forces)
              for(std::size_t i=0; i<3; ++i)
                pout[ind[2]][ind[1]][ind[0]][i] = px[num++] + p.rho_*p.force_[i];
            else
              for(std::size_t i=0; i<3; ++i)
                pout[ind[2]][ind[1]][ind[0]][i] = px[num++];
          }       
        }
      }

      ierr = DMDAVecRestoreArrayDOF(dau, outu, &pout);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, out, &outu, NULL);CHKERRQ(ierr);
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
      ierr = set_rhs_problem_<Dimensions, Ctx>(ctx, x, forces, true, std::integral_constant<int, Dimensions>{});CHKERRQ(ierr);
      ierr = fem::apply_mass_matrix(ctx.problem.ctx->dm, forces, mass_mult, ctx.problem.ctx->h);CHKERRQ(ierr);
      ierr = VecAXPY(ctx.problem.rhs, 1., mass_mult);CHKERRQ(ierr);

      ierr = DMRestoreGlobalVector(ctx.problem.ctx->dm, &mass_mult);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(ctx.problem.ctx->dm, &forces);CHKERRQ(ierr);


      DM dav;
      Vec rhsv;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm, &dav, nullptr);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, ctx.problem.rhs, &rhsv, nullptr);CHKERRQ(ierr);
      ierr = SetDirichletOnRHS(dav, ctx.problem.ctx->bc_, rhsv, ctx.problem.ctx->h);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, ctx.problem.rhs, &rhsv, nullptr);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    using surf_2d = std::vector<std::vector<std::pair<geometry::position<2, int>, geometry::position<2, double>>>>;
    using surf_3d = std::vector<std::vector<std::pair<geometry::position<3, int>, geometry::position<3, double>>>>;
    using array_2d_g = std::vector<std::vector<std::array<double, 2>>>;
    using array_3d_g = std::vector<std::vector<std::array<double, 3>>>;
    #undef __FUNCT__
    #define __FUNCT__ "interp_rigid_motion_"
    template<typename Shape>
    PetscErrorCode interp_rigid_motion_(std::vector<particle<Shape>> const& particles,
                                        std::vector<std::vector<geometry::position<2, double>>> const& r, 
                                        array_2d_g& g){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      std::size_t ipart=0;
      for(auto& rpts: r){       
        for(std::size_t i=0; i<rpts.size(); ++i){
          //g[ipart][i] -= particles[ipart].velocity_ - geometry::cross_product(particles[ipart].angular_velocity_[2], r[ipart][i]);
          g[ipart][i][0] -= particles[ipart].velocity_[0] - particles[ipart].angular_velocity_[2]*r[ipart][i][1];
          g[ipart][i][1] -= particles[ipart].velocity_[1] + particles[ipart].angular_velocity_[2]*r[ipart][i][0];
        }
        ipart++;
      }
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "interp_rigid_motion_"
    template<typename Shape>
    PetscErrorCode interp_rigid_motion_(std::vector<particle<Shape>> const& particles,
                                        std::vector<std::vector<geometry::position<3, double>>> const& r, 
                                        array_3d_g& g){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      std::size_t ipart=0;
      for(auto& rpts: r){       
        for(std::size_t i=0; i<rpts.size(); ++i){
          g[ipart][i][0] -= particles[ipart].velocity_[0] 
                          - particles[ipart].angular_velocity_[2]*r[ipart][i][1]
                          + particles[ipart].angular_velocity_[1]*r[ipart][i][2];
          g[ipart][i][1] -= particles[ipart].velocity_[1] 
                          + particles[ipart].angular_velocity_[2]*r[ipart][i][0]
                          - particles[ipart].angular_velocity_[0]*r[ipart][i][2];
          g[ipart][i][2] -= particles[ipart].velocity_[2] 
                          - particles[ipart].angular_velocity_[1]*r[ipart][i][0]
                          + particles[ipart].angular_velocity_[0]*r[ipart][i][1];
          //g[ipart][i] -= particles[ipart].velocity_ + cross_product(r[ipart][i], particles[ipart].angular_velocity_);
        }
        ipart++;
      }
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "interp_fluid_to_surf"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode interp_fluid_to_surf(Ctx& ctx, std::vector<std::vector<std::array<double, Dimensions>>>& g,
                                        bool rigid_motion = false)
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

      if (rigid_motion){
        ierr = interp_rigid_motion_(ctx.particles, ctx.radial_vec, g);CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "simple_layer"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode simple_layer(Ctx& ctx, std::vector<std::vector<std::array<double, Dimensions>>>& g,
                                std::vector<double>& mean, std::vector<double>& cross_prod)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      std::vector<double> mean_local;
      std::vector<double> cross_prod_local;

      mean_local.resize(mean.size());
      cross_prod_local.resize(cross_prod.size());

      std::fill(mean_local.begin(), mean_local.end(), 0.);
      std::fill(cross_prod_local.begin(), cross_prod_local.end(), 0.);

      for(std::size_t ipart=0; ipart<g.size(); ++ipart)
        for(std::size_t isurf=0; isurf<g[ipart].size(); ++isurf){
          for (std::size_t d=0; d<Dimensions; ++d)
            mean_local[ipart*Dimensions + d] += g[ipart][isurf][d];

          if (Dimensions == 2)
            cross_prod_local[ipart] += geometry::cross_product(ctx.radial_vec[ipart][isurf], g[ipart][isurf]);

            //cross_prod_local[ipart] += ctx.radial_vec[ipart][isurf][0]*g[ipart][isurf][1]
            //                         - ctx.radial_vec[ipart][isurf][1]*g[ipart][isurf][0];
          else{
            cross_prod_local[ipart*Dimensions]     += ctx.radial_vec[ipart][isurf][1]*g[ipart][isurf][2]
                                                    - ctx.radial_vec[ipart][isurf][2]*g[ipart][isurf][1];
            cross_prod_local[ipart*Dimensions + 1] += ctx.radial_vec[ipart][isurf][1]*g[ipart][isurf][2]
                                                    - ctx.radial_vec[ipart][isurf][2]*g[ipart][isurf][1];
            cross_prod_local[ipart*Dimensions + 2] += ctx.radial_vec[ipart][isurf][0]*g[ipart][isurf][1]
                                                    - ctx.radial_vec[ipart][isurf][1]*g[ipart][isurf][0];
          }
        }

      MPI_Allreduce(mean_local.data(), mean.data(), mean.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(cross_prod_local.data(), cross_prod.data(), cross_prod.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

      for(std::size_t ipart=0; ipart<g.size(); ++ipart){
        auto& radius = ctx.particles[ipart].shape_factors_[0];
        for(std::size_t d=0; d<Dimensions; ++d)
          mean[ipart*Dimensions + d] /= ctx.nb_surf_points[ipart];
        // fix for non circle shape
        if (Dimensions == 2)
          cross_prod[ipart] /= radius*radius*ctx.nb_surf_points[ipart];
        else
          for(std::size_t d=0; d<Dimensions; ++d)
            cross_prod[ipart*Dimensions + d] *= 2./(5*radius*radius*ctx.nb_surf_points[ipart]);
      }

      if (Dimensions == 2)
        for(std::size_t ipart=0; ipart<g.size(); ++ipart)
          for(std::size_t isurf=0; isurf<g[ipart].size(); ++isurf){
              g[ipart][isurf][0] -= mean[ipart*Dimensions] - cross_prod[ipart]*ctx.radial_vec[ipart][isurf][1];
              g[ipart][isurf][1] -= mean[ipart*Dimensions + 1] + cross_prod[ipart]*ctx.radial_vec[ipart][isurf][0];
          }
      else
        for(std::size_t ipart=0; ipart<g.size(); ++ipart)
          for(std::size_t isurf=0; isurf<g[ipart].size(); ++isurf){
              g[ipart][isurf][0] -= mean[ipart*Dimensions] 
                                  + cross_prod[ipart*Dimensions + 1]*ctx.radial_vec[ipart][isurf][2]
                                  - cross_prod[ipart*Dimensions + 2]*ctx.radial_vec[ipart][isurf][1];
              g[ipart][isurf][1] -= mean[ipart*Dimensions + 1]
                                  + cross_prod[ipart*Dimensions + 2]*ctx.radial_vec[ipart][isurf][0]
                                  - cross_prod[ipart*Dimensions    ]*ctx.radial_vec[ipart][isurf][2];
              g[ipart][isurf][2] -= mean[ipart*Dimensions + 2] 
                                  + cross_prod[ipart*Dimensions    ]*ctx.radial_vec[ipart][isurf][1]
                                  - cross_prod[ipart*Dimensions + 1]*ctx.radial_vec[ipart][isurf][0];
          }

      PetscFunctionReturn(0);
    }
    #undef __FUNCT__
    #define __FUNCT__ "SL_to_Rhs"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode SL_to_Rhs(Ctx& ctx, std::vector<std::vector<std::array<double, Dimensions>>>const& g){
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
      ierr = sol.local_to_global(ADD_VALUES);CHKERRQ(ierr);

      DM dav;
      Vec rhs;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm, &dav, nullptr);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, ctx.problem.rhs, &rhs, nullptr);CHKERRQ(ierr);
      ierr = SetNullDirichletOnRHS(dav, ctx.problem.ctx->bc_, rhs, ctx.problem.ctx->h);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, ctx.problem.rhs, &rhs, nullptr);CHKERRQ(ierr);
      PetscFunctionReturn(0);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "init_problem_1"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode init_problem_1(Ctx& ctx, Vec x, bool apply_forces=false){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      if (ctx.compute_rhs){
        ierr = ctx.problem.setup_RHS();CHKERRQ(ierr);
      }
      
      ierr = set_rhs_problem<Dimensions, Ctx>(ctx, x, apply_forces);CHKERRQ(ierr);
      
      DM dav;
      Vec rhs;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm, &dav, nullptr);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, ctx.problem.rhs, &rhs, nullptr);CHKERRQ(ierr);

      if (ctx.compute_rhs){
        ierr = SetDirichletOnRHS(dav, ctx.problem.ctx->bc_, rhs, ctx.problem.ctx->h);CHKERRQ(ierr);
      }
      else{
        ierr = SetNullDirichletOnRHS(dav, ctx.problem.ctx->bc_, rhs, ctx.problem.ctx->h);CHKERRQ(ierr); 
      }
      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, ctx.problem.rhs, &rhs, nullptr);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "projection"
    template<typename Shape>
    PetscErrorCode projection(std::vector<particle<Shape>> const& particles,
                              DM dm, 
                              Vec sol,
                              geometry::box<2, int> const& box,
                              std::vector<double>& mean, 
                              std::vector<double>& cross_prod,
                              std::vector<int> const& num,
                              std::size_t scale,
                              std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 2;
      using position_type = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;

      std::array<double, Dimensions> hs = {{h[0]/scale, h[1]/scale}};

      PetscScalar ***psol;
      ierr = DMDAVecGetArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart){
        auto& p = particles[ipart];
        auto pbox = p.bounding_box(h);
        double meanu=0, meanv=0, cross=0;
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          for(std::size_t j=new_box.bottom_left[1]; j<new_box.upper_right[1]; ++j)
            for(std::size_t i=new_box.bottom_left[0]; i<new_box.upper_right[0]; ++i)
              for(std::size_t js=0; js<scale; ++js)
                for(std::size_t is=0; is<scale; ++is){
                  position_type pts = {i*h[0] + is*hs[0], j*h[1] + js*hs[1]};
                  position_type_i pts_i = {i, j};
                  if (p.contains(pts)){
                    position_type pts_corner = {is*hs[0], js*hs[1]};
                    auto bfunc = fem::P1_integration(pts_corner, h);
                    auto ielem = fem::get_element(pts_i);
                    auto solu = 0.;
                    auto solv = 0.;
                    for(std::size_t ib=0; ib<bfunc.size(); ++ib){
                      solu += psol[ielem[ib][1]][ielem[ib][0]][0]*bfunc[ib];
                      solv += psol[ielem[ib][1]][ielem[ib][0]][1]*bfunc[ib];
                    }
                    meanu += solu;
                    meanv += solv;
                    cross += (pts[0] - p.center_[0])*solv - (pts[1] - p.center_[1])*solu;
                  }
                }
          mean[ipart*Dimensions] = meanu/num[ipart];
          mean[ipart*Dimensions + 1] = meanv/num[ipart];
          cross_prod[ipart] = cross/num[ipart];
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "projection"
    template<typename Shape>
    PetscErrorCode projection(std::vector<particle<Shape>> const& particles,
                              DM dm, 
                              Vec sol,
                              geometry::box<3, int> const& box,
                              std::vector<double>& mean, 
                              std::vector<double>& cross_prod,
                              std::vector<int> const& num,
                              std::size_t scale,
                              std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 3;
      using position_type = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;

      std::array<double, Dimensions> hs = {{h[0]/scale, h[1]/scale, h[2]/scale}};

      PetscScalar ****psol;
      ierr = DMDAVecGetArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart){
        auto& p = particles[ipart];
        auto pbox = p.bounding_box(h);
        std::array<double, Dimensions> mean_tmp{0};
        std::array<double, Dimensions> cross_tmp{0};
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);

          for(std::size_t k=new_box.bottom_left[2]; k<new_box.upper_right[2]; ++k)
            for(std::size_t j=new_box.bottom_left[1]; j<new_box.upper_right[1]; ++j)
              for(std::size_t i=new_box.bottom_left[0]; i<new_box.upper_right[0]; ++i)
                for(std::size_t ks=0; ks<scale; ++ks)
                  for(std::size_t js=0; js<scale; ++js)
                    for(std::size_t is=0; is<scale; ++is){
                      position_type pts = {i*h[0] + is*hs[0], j*h[1] + js*hs[1], k*h[2] + ks*hs[2]};
                      position_type_i pts_i = {i, j, k};
                      if (p.contains(pts)){
                        position_type pts_corner = {is*hs[0], js*hs[1], ks*hs[2]};
                        auto bfunc = fem::P1_integration(pts_corner, h);
                        auto ielem = fem::get_element(pts_i);
                        std::array<double, Dimensions> tmp{0.};
                        for(std::size_t ib=0; ib<bfunc.size(); ++ib)
                          for(std::size_t d=0; d<Dimensions; ++d)
                            tmp[d] += psol[ielem[ib][2]][ielem[ib][1]][ielem[ib][0]][d]*bfunc[ib];
                        for(std::size_t d=0; d<Dimensions; ++d)
                          mean_tmp[d] += tmp[d];
                        cross_tmp[0] += (pts[1] - p.center_[1])*mean_tmp[2] - (pts[2] - p.center_[2])*mean_tmp[1];
                        cross_tmp[1] += (pts[2] - p.center_[2])*mean_tmp[0] - (pts[0] - p.center_[0])*mean_tmp[2];
                        cross_tmp[2] += (pts[0] - p.center_[0])*mean_tmp[1] - (pts[1] - p.center_[1])*mean_tmp[0];
                      }
                    }
          for(std::size_t d=0; d<Dimensions; ++d){
            mean[ipart*Dimensions + d ] = mean_tmp[d]/num[ipart];
            cross_prod[ipart*Dimensions + d ] = cross_tmp[d]/num[ipart];
          }
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "projection"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode projection(Ctx& ctx, std::vector<double>& mean, std::vector<double>& cross_prod)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      using position_type = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      std::vector<double> mean_local;
      std::vector<double> cross_prod_local;

      mean_local.resize(mean.size());
      cross_prod_local.resize(cross_prod.size());

      std::fill(mean_local.begin(), mean_local.end(), 0.);
      std::fill(cross_prod_local.begin(), cross_prod_local.end(), 0.);

      DM dav;
      Vec sol;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm, &dav, nullptr);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, ctx.problem.sol, &sol, nullptr);CHKERRQ(ierr);
      
      Vec sol_local;
      ierr = DMGetLocalVector(dav, &sol_local);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(dav, sol, INSERT_VALUES, sol_local);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dav, sol, INSERT_VALUES, sol_local);CHKERRQ(ierr);

      ierr = projection(ctx.particles, dav, sol_local, box, mean_local, cross_prod_local, ctx.num, ctx.scale, h);

      MPI_Allreduce(mean_local.data(), mean.data(), mean.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(cross_prod_local.data(), cross_prod.data(), cross_prod.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

      ierr = DMRestoreLocalVector(dav, &sol_local);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, ctx.problem.sol, &sol, nullptr);CHKERRQ(ierr);      
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "compute_y"
    template<typename Shape>
    PetscErrorCode compute_y(std::vector<particle<Shape>> const& particles,
                             DM dm, 
                             Vec y, 
                             Vec sol,
                             geometry::box<2, int> const& box,
                             std::vector<double> const& mean, 
                             std::vector<double> const& cross_prod,
                             std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 2;

      PetscScalar ***psol;
      ierr = DMDAVecGetArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);

      int num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart){
        auto& p = particles[ipart];
        auto pbox = p.bounding_box(h);
        auto radius = p.shape_factors_[0];
        auto Cr = 2./(radius*radius);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts){
            py[num++] = psol[ind[1]][ind[0]][0] - mean[ipart*Dimensions  ] + Cr*(ind[1]*h[1] - p.center_[1])*cross_prod[ipart];
            py[num++] = psol[ind[1]][ind[0]][1] - mean[ipart*Dimensions+1] - Cr*(ind[0]*h[0] - p.center_[0])*cross_prod[ipart];
          }
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);
      ierr = VecRestoreArray(y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }
    
    #undef __FUNCT__
    #define __FUNCT__ "compute_y"
    template<typename Shape>
    PetscErrorCode compute_y(std::vector<particle<Shape>> const& particles,
                             DM dm, 
                             Vec y, 
                             Vec sol,
                             geometry::box<3, int> const& box,
                             std::vector<double> const& mean, 
                             std::vector<double> const& cross_prod,
                             std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 3;
      using position_type = geometry::position<Dimensions, double>;

      PetscScalar ***psol;
      ierr = DMDAVecGetArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);

      int num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart){
        auto& p = particles[ipart];
        auto pbox = p.bounding_box(h);
        auto radius = p.shape_factors_[0];
        auto Cr = 2./(5*radius*radius);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts){
            position_type pos { ind[0]*h[0] - p.center_[0], ind[1]*h[1] - p.center_[1], ind[2]*h[2] - p.center_[2] };
            py[num++] = psol[ind[2]][ind[1]][ind[0]][0] - mean[ipart*Dimensions] 
                      - Cr*( cross_prod[ipart*Dimensions + 1]*pos[2] 
                           - cross_prod[ipart*Dimensions + 2]*pos[1]);
            py[num++] = psol[ind[2]][ind[1]][ind[0]][1] - mean[ipart*Dimensions+1] 
                      - Cr*( cross_prod[ipart*Dimensions    ]*pos[2] 
                           - cross_prod[ipart*Dimensions + 2]*pos[0]);
            py[num++] = psol[ind[2]][ind[1]][ind[0]][2] - mean[ipart*Dimensions+2] 
                      - Cr*( cross_prod[ipart*Dimensions + 1]*pos[2] 
                           - cross_prod[ipart*Dimensions + 2]*pos[1]);
          }
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);
      ierr = VecRestoreArray(y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "compute_y"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode compute_y(Ctx& ctx, Vec y, std::vector<double> const& mean, std::vector<double> const& cross_prod)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      Vec sol;
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, ctx.problem.sol, &sol, NULL);CHKERRQ(ierr);
      DM dav;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm , &dav, PETSC_NULL);CHKERRQ(ierr);

      compute_y(ctx.particles, dav, y, sol, box, mean, cross_prod, h);

      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, ctx.problem.sol, &sol, NULL);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "compute_y"
    template<typename Shape>
    PetscErrorCode compute_y(std::vector<particle<Shape>> const& particles,
                             DM dm, 
                             Vec y, 
                             Vec sol,
                             geometry::box<2, int> const& box,
                             std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 2;

      PetscScalar ***psol;
      ierr = DMDAVecGetArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);

      int num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart){
        auto& p = particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts){
            py[num++] = psol[ind[1]][ind[0]][0];
            py[num++] = psol[ind[1]][ind[0]][1];
          }
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);
      ierr = VecRestoreArray(y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }
    
    #undef __FUNCT__
    #define __FUNCT__ "compute_y"
    template<typename Shape>
    PetscErrorCode compute_y(std::vector<particle<Shape>> const& particles,
                             DM dm, 
                             Vec y, 
                             Vec sol,
                             geometry::box<3, int> const& box,
                             std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 3;
      using position_type = geometry::position<Dimensions, double>;

      PetscScalar ****psol;
      ierr = DMDAVecGetArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);

      int num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart){
        auto& p = particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::box_inside(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts){
            py[num++] = psol[ind[2]][ind[1]][ind[0]][0];
            py[num++] = psol[ind[2]][ind[1]][ind[0]][1];
            py[num++] = psol[ind[2]][ind[1]][ind[0]][2];
          }
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dm, sol, &psol);CHKERRQ(ierr);
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

      Vec sol;
      ierr = DMCompositeGetAccess(ctx.problem.ctx->dm, ctx.problem.sol, &sol, NULL);CHKERRQ(ierr);
      DM dav;
      ierr = DMCompositeGetEntries(ctx.problem.ctx->dm , &dav, PETSC_NULL);CHKERRQ(ierr);

      compute_y(ctx.particles, dav, y, sol, box, h);

      ierr = DMCompositeRestoreAccess(ctx.problem.ctx->dm, ctx.problem.sol, &sol, NULL);CHKERRQ(ierr);
      PetscFunctionReturn(0);
    }
  }
}

#endif