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

#ifndef CAFES_FEM_BC_HPP_INCLUDED
#define CAFES_FEM_BC_HPP_INCLUDED

#include <algorithm/iterate.hpp>
#include <fem/mesh.hpp>
#include <petsc/vec.hpp>

#include <array>

namespace cafes
{
  using condition_fn = void(*)(const PetscReal*, PetscScalar*);
  namespace fem
  {
    template<std::size_t Dimensions>
    struct dirichlet_conditions
    {
      dirichlet_conditions( std::initializer_list<std::array<condition_fn, Dimensions>> il )
      {
        std::copy(il.begin(), il.end(), conditions_.begin());
      }

      std::array<std::array<condition_fn, Dimensions>, 2*Dimensions> conditions_;

      dirichlet_conditions(dirichlet_conditions const&) = default;
      dirichlet_conditions(dirichlet_conditions&&)      = default;

      dirichlet_conditions& operator=(dirichlet_conditions const&) = default;
      dirichlet_conditions& operator=(dirichlet_conditions&&)      = default;

    };

    auto const kernel_vec_on_bc = [](auto const& x, auto& y, auto const& bc_conditions, auto const& h){
      auto const kernel_pos = [&](auto const& pos){
        auto ux = x.at_g(pos);
        auto uy = y.at_g(pos);
        for(std::size_t dof=0; dof<y.dof_; ++dof)
          if (bc_conditions[dof])
            uy[dof] = ux[dof];
      };
      return kernel_pos;
    };

    auto const kernel_const_on_bc = [](auto const& x, auto& y, auto const& bc_conditions, auto const& h){
      auto const kernel_pos = [&](auto const& pos){
        auto uy = y.at_g(pos);
        for(std::size_t dof=0; dof<y.dof_; ++dof)
          if (bc_conditions[dof])
            uy[dof] = x;
      };
      return kernel_pos;
    };

    auto const kernel_func_on_bc = [](auto const& x, auto& y, auto const& bc_conditions, auto const& h){
      auto const kernel_pos = [&](auto const& pos){
        auto uy = y.at_g(pos);
        double coord[pos.dimensions];

        for(std::size_t i=0; i<pos.dimensions; ++i)
          coord[i] = pos[i]*h[i];

        for(std::size_t dof=0; dof<y.dof_; ++dof)
          if (bc_conditions[dof])
            bc_conditions[dof](coord, &uy[dof]);
      };
      return kernel_pos;
    };

    #undef __FUNCT__
    #define __FUNCT__ "set_dirichlet_impl"
    template<std::size_t Dimensions, typename Function, typename x_type, typename h_type>
    PetscErrorCode set_dirichlet_impl(x_type&& x, 
                                      petsc::petsc_vec<Dimensions>& y,
                                      dirichlet_conditions<Dimensions> bc,
                                      h_type&& h,
                                      Function&& kernel)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto bd_type = get_boundary_type<Dimensions>(y.dm_);
      auto gbounds = get_global_bounds<Dimensions>(y.dm_);
      auto box = get_DM_bounds<Dimensions>(y.dm_, false);

      for(std::size_t d=0; d<Dimensions; ++d)
      {
        if (bd_type[d] != DM_BOUNDARY_PERIODIC)
        {
          if (box.bottom_left[d] == 0)
          {
            auto new_box{box};
            new_box.upper_right[d] = 1;
            algorithm::iterate(new_box, kernel(x, y, bc.conditions_[2*d], h));
          }

          if (box.upper_right[d] == gbounds[d])
          {
            auto new_box{box};
            new_box.bottom_left[d] = gbounds[d]-1;
            algorithm::iterate(new_box, kernel(x, y, bc.conditions_[2*d+1], h));
          }

        }
      }
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetDirichletOnVec"
    template<std::size_t Dimensions>
    PetscErrorCode SetDirichletOnVec(petsc::petsc_vec<Dimensions>& x, 
                                     petsc::petsc_vec<Dimensions>& y,
                                     dirichlet_conditions<Dimensions> bc)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      set_dirichlet_impl(x, y, bc, nullptr, kernel_vec_on_bc);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetDirichletOnRHS"
    template<std::size_t Dimensions>
    PetscErrorCode SetDirichletOnRHS(petsc::petsc_vec<Dimensions>& x,
                                     dirichlet_conditions<Dimensions> bc,
                                     std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      set_dirichlet_impl(nullptr, x, bc, h, kernel_func_on_bc);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetNullDirichletOnRHS"
    template<std::size_t Dimensions>
    PetscErrorCode SetNullDirichletOnRHS(petsc::petsc_vec<Dimensions>& x,
                                         dirichlet_conditions<Dimensions> bc)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      set_dirichlet_impl(0., x, bc, nullptr, kernel_const_on_bc);

      PetscFunctionReturn(0);
    }
  }
  template<std::size_t Dimensions>
  fem::dirichlet_conditions<Dimensions> make_bc(std::initializer_list<std::array<condition_fn, Dimensions>> il){
    return {il};
  }
}
#endif