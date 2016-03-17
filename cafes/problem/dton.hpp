#ifndef CAFES_PROBLEM_DTON_HPP_INCLUDED
#define CAFES_PROBLEM_DTON_HPP_INCLUDED

#include <problem/particle_operator.hpp>
#include <problem/problem.hpp>
#include <problem/stokes.hpp>
#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <particle/particle.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>

#include <petsc.h>
#include <iostream>
#include <memory>
#include <algorithm>

namespace cafes
{
  namespace problem
  {

    #undef __FUNCT__
    #define __FUNCT__ "DtoN_matrix"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode DtoN_matrix(Mat A, Vec x, Vec y){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      Ctx *ctx;
      ierr = MatShellGetContext(A, (void**) &ctx);CHKERRQ(ierr);

      ierr = VecSet(ctx->problem.rhs, 0.);CHKERRQ(ierr);

      ierr = init_problem_1<Dimensions, Ctx>(*ctx, x);CHKERRQ(ierr);
      ierr = ctx->problem.solve();CHKERRQ(ierr);

      std::vector<std::vector<std::array<double, Dimensions>>> g;
      g.resize(ctx->particles.size());
      for(std::size_t ipart=0; ipart<ctx->surf_points.size(); ++ipart)       
        g[ipart].resize(ctx->surf_points[ipart].size());

      // interpolation
      ierr = interp_fluid_to_surf(*ctx, g, ctx->compute_rhs);CHKERRQ(ierr);

      ierr = SL_to_Rhs(*ctx, g);CHKERRQ(ierr);
      
      ierr = ctx->problem.solve();CHKERRQ(ierr);

      ierr = compute_y<Dimensions, Ctx>(*ctx, y);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }



    template<typename Shape, std::size_t Dimensions, typename Problem_type>
    struct DtoN : public Problem<Dimensions>
    {
      std::vector<particle<Shape>> parts_;
      Problem_type problem_;

      using position_type   = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;

      using Ctx = particle_context<Dimensions, Shape, Problem_type>;
      Ctx *ctx;

      std::vector<std::vector<std::pair<position_type_i, position_type>>> surf_points_;
      std::vector<std::vector<position_type>> radial_vec_;
      std::vector<int> nb_surf_points_;
      std::vector<int> num_;
      Vec sol;
      Vec rhs;
      Mat A;
      KSP ksp;
      std::size_t scale_=4;

      using dpart_type = typename std::conditional<Dimensions == 2, 
                                  double, 
                                  std::array<double, 2>>::type;
      dpart_type dpart_; 

      DtoN(std::vector<particle<Shape>>const& parts, Problem_type& p, dpart_type dpart):
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


        surf_points_.resize(parts_.size());
        radial_vec_.resize(parts_.size());
        nb_surf_points_.resize(parts_.size());

        std::vector<int> nb_surf_points_local;
        nb_surf_points_local.resize(parts_.size());

        auto box = fem::get_DM_bounds<Dimensions>(problem_.ctx->dm, 0);
        std::size_t size = 0;
        std::size_t ipart = 0;

        std::vector<int> num_local;
        num_local.resize(parts_.size());
        num_.resize(parts_.size());

        auto& h = problem_.ctx->h;
        std::array<double, Dimensions> hs;
        for(std::size_t d=0; d<Dimensions; ++d)
          hs[d] = h[d]/scale_;

        for(auto& p: parts_){
          auto pbox = p.bounding_box(h);
          if (geometry::intersect(box, pbox)){
            auto new_box = geometry::box_inside(box, pbox);
            auto pts = find_fluid_points_insides(p, new_box, h);
            size += pts.size();

            auto spts = p.surface(dpart_);   
            auto spts_valid = find_surf_points_insides(spts, new_box, h);
            surf_points_[ipart].assign(spts_valid.begin(), spts_valid.end());
            nb_surf_points_local[ipart] = surf_points_[ipart].size();
            
            auto radial_valid = find_radial_surf_points_insides(spts, new_box, h, p.center_);
            radial_vec_[ipart].assign(radial_valid.begin(), radial_valid.end());
            
            if (Dimensions == 2)
            {
              for(std::size_t j=new_box.bottom_left[1]; j<new_box.upper_right[1]; ++j)
                for(std::size_t i=new_box.bottom_left[0]; i<new_box.upper_right[0]; ++i)
                  for(std::size_t js=0; js<scale_; ++js)
                    for(std::size_t is=0; is<scale_; ++is){
                      position_type pts_2 = {i*h[0] + is*hs[0], j*h[1] + js*hs[1]};
                        if (p.contains(pts_2)) 
                          num_local[ipart]++;
                      }
            }
            // else
            // {
            //   for(std::size_t k=new_box.bottom_left[2]; k<new_box.upper_right[2]; ++k)
            //     for(std::size_t j=new_box.bottom_left[1]; j<new_box.upper_right[1]; ++j)
            //       for(std::size_t i=new_box.bottom_left[0]; i<new_box.upper_right[0]; ++i)
            //         for(std::size_t ks=0; ks<scale_; ++ks)
            //           for(std::size_t js=0; js<scale_; ++js)
            //             for(std::size_t is=0; is<scale_; ++is){
            //               position_type pts_2 = {i*h[0] + is*hs[0], j*h[1] + js*hs[1], k*h[2] + ks*hs[2]};
            //                 if (p.contains(pts_2)) 
            //                   num_local[ipart]++;
            //               }
            // }
          }
          ipart++;
        }

        MPI_Allreduce(nb_surf_points_local.data(), nb_surf_points_.data(), parts_.size(), MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce(num_local.data(), num_.data(), parts_.size(), MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

        ctx = new Ctx{problem_, parts_, surf_points_, radial_vec_, nb_surf_points_, num_, scale_, false, false};

        ierr = MatCreateShell(PETSC_COMM_WORLD, size*Dimensions, size*Dimensions, PETSC_DECIDE, PETSC_DECIDE, ctx, &A);CHKERRQ(ierr);
        ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))DtoN_matrix<Dimensions, Ctx>);CHKERRQ(ierr);

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
        ctx->compute_sol = false;

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
        ierr = KSPSetOptionsPrefix(ksp, "DtoN_");CHKERRQ(ierr);

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
        ctx->compute_sol = true;

        ierr = KSPSolve(ksp, rhs, sol);CHKERRQ(ierr);

        // solve the problem with the right control
        ctx->compute_rhs = true;
        ctx->compute_sol = true;
        ierr = init_problem_1<Dimensions, Ctx>(*ctx, sol);CHKERRQ(ierr);
        ierr = ctx->problem.solve();CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      };
  }

  template<typename PL, typename Problem_type, typename Dimensions = typename PL::value_type::dimension_type> 
  auto make_DtoN(PL const& pt, Problem_type& p, 
                typename std::conditional<Dimensions::value == 2, double, 
                                  std::array<double, 2>>::type const& dpart)
  {
    using s_t = typename PL::value_type::shape_type;
    //using Dimensions = typename PL::value_type::dimension_type;
    return problem::DtoN<s_t, Dimensions::value, Problem_type>{pt, p, dpart};
  }

}
#endif