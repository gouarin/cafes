#ifndef CAFES_FEM_RHS_HPP_INCLUDED
#define CAFES_FEM_RHS_HPP_INCLUDED

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
                                auto const& h,
                                Function&& kernel)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      // warning: we assume that both u and v on a border have a Dirichlet condition
      auto box = get_DM_bounds<Dimensions>(x.dm_, false);
      iterate(box, kernel(x, rhs, h));
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