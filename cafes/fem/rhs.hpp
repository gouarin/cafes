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

#ifndef CAFES_FEM_RHS_HPP_INCLUDED
#define CAFES_FEM_RHS_HPP_INCLUDED

#include <algorithm/iterate.hpp>
#include <fem/matrixFree.hpp>
#include <particle/geometry/position.hpp>
#include <petsc.h>
#include <array>
#include <type_traits>

namespace cafes
{
  using condition_fn = void(*)(const PetscReal*, PetscScalar*);
  namespace fem
  {
    template<std::size_t Ndof>
    struct rhs_conditions
    {
      std::array<condition_fn, Ndof> conditions_;

      auto has_condition(){
        return !std::all_of(conditions_.cbegin(), conditions_.cend(), [](auto i){ return i == nullptr; });
      }
    };

    auto const kernel_func_on_rhs = [](auto& x, auto const& rhs_conditions, auto const& h){
      auto const kernel_pos = [&](auto const& pos){
        auto u = x.at_g(pos);
        double coord[pos.dimensions];

        for(std::size_t i=0; i<pos.dimensions; ++i)
          coord[i] = pos[i]*h[i];

        for(std::size_t dof=0; dof<x.dof_; ++dof)
          if (rhs_conditions.conditions_[dof])
            rhs_conditions.conditions_[dof](coord, &u[dof]);
      };
      return kernel_pos;
    };

    #undef __FUNCT__
    #define __FUNCT__ "set_rhs_impl"
    template<std::size_t Dimensions, typename Function>
    PetscErrorCode set_rhs_impl(petsc::petsc_vec<Dimensions>& x,
                                rhs_conditions<Dimensions> rhs,
                                std::array<double, Dimensions> const& h,
                                Function&& kernel)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      // warning: we assume that both u and v on a border have a Dirichlet condition
      auto box = get_DM_bounds<Dimensions>(x.dm_, false);
      algorithm::iterate(box, kernel(x, rhs, h));

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "set_rhs"
    template<std::size_t Ndof, std::size_t Ndm, std::size_t Dimensions>
    typename std::enable_if<Ndm==1, PetscErrorCode>::type
    set_rhs(DM dm, Vec rhs, rhs_conditions<Ndof> const& rhs_cond, std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;
    
      Vec tmp;
      ierr = DMGetGlobalVector(dm, &tmp);CHKERRQ(ierr);
      ierr = VecSet(tmp, 0.);CHKERRQ(ierr);

      auto xpetsc = petsc::petsc_vec<Dimensions>(dm, tmp, 0, false);
      ierr = set_rhs_impl(xpetsc, rhs_cond, h, kernel_func_on_rhs);CHKERRQ(ierr);

      auto ypetsc = petsc::petsc_vec<Dimensions>(dm, rhs, 0, false);
      ierr = xpetsc.global_to_local(INSERT_VALUES);CHKERRQ(ierr);
      ierr = ypetsc.fill(0.);CHKERRQ(ierr);
      ierr = mass_mult(xpetsc, ypetsc, h);CHKERRQ(ierr);
      ierr = ypetsc.local_to_global(ADD_VALUES);CHKERRQ(ierr);

      ierr = DMRestoreGlobalVector(dm, &tmp);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }
  }

  template<std::size_t Ndof>
  fem::rhs_conditions<Ndof> make_rhs(std::array<condition_fn, Ndof> conditions)
  {
    return {conditions};
  }
}

#endif