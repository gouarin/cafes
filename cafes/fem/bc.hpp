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
        PetscFunctionReturn(0);
      };
      return kernel_pos;
    };

    auto const kernel_const_on_bc = [](auto const& x, auto& y, auto const& bc_conditions, auto const& h){
      auto const kernel_pos = [&](auto const& pos){
        auto uy = y.at_g(pos);
        for(std::size_t dof=0; dof<y.dof_; ++dof)
          if (bc_conditions[dof])
            uy[dof] = x;
        PetscFunctionReturn(0);
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
        PetscFunctionReturn(0);
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
            auto new_box = box;
            new_box.upper_right[d] = 1;
            algorithm::iterate(new_box, kernel(x, y, bc.conditions_[2*d], h));
          }

          if (box.upper_right[d] == gbounds[d])
          {
            auto new_box = box;
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

    #undef __FUNCT__
    #define __FUNCT__ "set_dirichlet_on_mat_impl"
    PetscErrorCode set_dirichlet_on_mat_impl(Mat& A,
                                     dirichlet_conditions<2> bc)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      PetscInt       i, j, d, ind, nx, ny, size;
      PetscInt       local_id, nbcs;
      PetscInt       *bc_id1, *bc_id2;
      DMDALocalInfo  info;
      IS is;

      // warning: we assume that both u and v on a border have a Dirichlet condition
      DM dm;
      PetscBool is_composite;
      ierr = MatGetDM(A, &dm);CHKERRQ(ierr);
      ierr = PetscObjectTypeCompare((PetscObject)dm, DMCOMPOSITE, &is_composite);CHKERRQ(ierr);
      if (is_composite)
      {
        ierr = DMCompositeGetEntries(dm, &dm, PETSC_NULL);CHKERRQ(ierr);
      }
      auto bd_type = get_boundary_type<2>(dm);

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);

      ierr = PetscMalloc(sizeof(PetscInt)*2*info.gym, &bc_id1);CHKERRQ(ierr);
      ierr = PetscMalloc(sizeof(PetscInt)*2*info.gxm, &bc_id2);CHKERRQ(ierr);

      MPI_Comm_size(PETSC_COMM_WORLD, &size);
      nx = (size == 1)? info.mx: info.gxm;
      ny = (size == 1)? info.my: info.gym;

      //left border
      nbcs = 0;
      i = 0;
      ind = 0;
      for(j=0; j<ny; j++){
        local_id = info.dof*(i + j*nx);
        bc_id1[ind] = local_id;
        bc_id1[ind + 1] = local_id + 1;
        ind += 2;
      }

      if (bd_type[0] != DM_BOUNDARY_PERIODIC)
        if(info.xs == 0 && info.bx != DM_BOUNDARY_PERIODIC){
          nbcs = 2*ny;
        }

      ierr = MatZeroRowsLocal(A, nbcs, bc_id1, 1.0, NULL, NULL);CHKERRQ(ierr);

      // right border
      nbcs = 0;
      i = (size == 1) ? info.mx-1: info.gxm-1;
      ind = 0;
      for(j=0; j<ny; j++){
        local_id = info.dof*(i + j*nx);
        bc_id1[ind] = local_id;
        bc_id1[ind+1] = local_id + 1;
        ind += 2;
      }

      if (bd_type[0] != DM_BOUNDARY_PERIODIC)
        if(info.xs + info.xm == info.mx && info.bx != DM_BOUNDARY_PERIODIC){
          nbcs = 2*ny;
        }

      ierr = MatZeroRowsLocal(A, nbcs, bc_id1, 1.0, NULL, NULL);CHKERRQ(ierr);

      nx = (size == 1) ? info.mx: info.gxm;

      // bottom border
      nbcs = 0;
      j = 0;
      ind = 0;
      for(i=0; i<nx; i++){
        local_id = info.dof*(i + j*nx);
        bc_id2[ind] = local_id;
        bc_id2[ind+1] = local_id + 1;
        ind += 2;
      }

      if (bd_type[1] != DM_BOUNDARY_PERIODIC)
        if(info.ys == 0 && info.by != DM_BOUNDARY_PERIODIC){
          nbcs =  2*nx;
        }

      ierr = MatZeroRowsLocal(A, nbcs, bc_id2, 1.0, NULL, NULL);CHKERRQ(ierr);

      // top border
      nbcs = 0;
      j = (size == 1) ? info.my-1:  info.gym-1;
      ind = 0;
      for(i=0; i<nx; i++){
        local_id = info.dof*(i + j*nx);
        bc_id2[ind] = local_id;
        bc_id2[ind+1] = local_id + 1;
        ind += 2;
      }

      if (bd_type[1] != DM_BOUNDARY_PERIODIC)
        if(info.ys + info.ym == info.my && info.by != DM_BOUNDARY_PERIODIC){
          nbcs = 2*nx;
        }

      ierr = MatZeroRowsLocal(A, nbcs, bc_id2, 1.0, NULL, NULL);CHKERRQ(ierr);

      ierr = PetscFree(bc_id1);CHKERRQ(ierr);
      ierr = PetscFree(bc_id2);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "SetDirichletOnMat"
    template<std::size_t Dimensions>
    PetscErrorCode SetDirichletOnMat(Mat& A,
                                     dirichlet_conditions<Dimensions> bc)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      set_dirichlet_on_mat_impl(A, bc);

      PetscFunctionReturn(0);
    }

  }

  template<std::size_t Dimensions>
  fem::dirichlet_conditions<Dimensions> make_bc(std::initializer_list<std::array<condition_fn, Dimensions>> il){
    return {il};
  }
}
#endif