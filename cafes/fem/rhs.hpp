#ifndef CAFES_FEM_RHS_HPP_INCLUDED
#define CAFES_FEM_RHS_HPP_INCLUDED

#include <fem/matrixFree.hpp>
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

    #undef __FUNCT__
    #define __FUNCT__ "set_rhs"
    template<std::size_t Ndof, std::size_t Ndm>
    typename std::enable_if<Ndm==1, PetscErrorCode>::type
    set_rhs(DM dm, Vec rhs, rhs_conditions<Ndof> const& rhs_cond, std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      DM             dav, dap;
      DMDALocalInfo  infov, infop;
      Vec            bv, bp;
      PetscScalar    ***pbv, **pbp;
      PetscReal      x[2];
      PetscFunctionBeginUser;
    
      Vec tmp;
      ierr = DMGetGlobalVector(dm, &tmp);CHKERRQ(ierr);
      ierr = VecSet(tmp, 0.);CHKERRQ(ierr);
      ierr = DMCompositeGetEntries(dm, &dav, &dap);CHKERRQ(ierr);
      ierr = DMDAGetLocalInfo(dav, &infov);CHKERRQ(ierr);
      ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(dm, tmp, &bv, &bp);CHKERRQ(ierr);

      // set RHS for the velocity field
      ierr = DMDAVecGetArrayDOF(dav, bv, &pbv);CHKERRQ(ierr);
      for (std::size_t j=infov.ys; j<infov.ys+infov.ym; ++j){
        x[1] = j*h[1];
        for (std::size_t i=infov.xs; i<infov.xs+infov.xm; ++i){
          x[0] = i*h[0];
          for (std::size_t d=0; d<infov.dof;++d)
            if (rhs_cond.conditions_[d])
              rhs_cond.conditions_[d](x, &pbv[j][i][d]);
        }
      }
      ierr = DMDAVecRestoreArrayDOF(dav, bv, &pbv);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(dm, tmp, &bv, &bp);CHKERRQ(ierr);

      auto hp{h};
      std::for_each(hp.begin(), hp.end(), [](auto x){x*=2;});

      using ctx = problem::context<2, 2>;
      ctx s{dm, {{h, hp }}, mass_mult};
      Mat A = make_matrix<ctx>(&s, diag_block_matrix<ctx>);

      ierr = MatMult(A, tmp, rhs);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(dm, &tmp);CHKERRQ(ierr);
      ierr = MatDestroy(&A);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


    #undef __FUNCT__
    #define __FUNCT__ "set_rhs"
    template<std::size_t Ndof, std::size_t Ndm>
    typename std::enable_if<Ndm==1, PetscErrorCode>::type
    set_rhs(DM dm, Vec rhs, rhs_conditions<Ndof> const& rhs_cond, std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      DM             dav, dap;
      DMDALocalInfo  infov, infop;
      Vec            bv, bp;
      PetscScalar    ****pbv, ***pbp;
      PetscReal      x[3];
      PetscFunctionBeginUser;
    
      Vec tmp;
      ierr = DMGetGlobalVector(dm, &tmp);CHKERRQ(ierr);
      ierr = VecSet(tmp, 0.);CHKERRQ(ierr);
      ierr = DMCompositeGetEntries(dm, &dav, &dap);CHKERRQ(ierr);
      ierr = DMDAGetLocalInfo(dav, &infov);CHKERRQ(ierr);
      ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);
      ierr = DMCompositeGetAccess(dm, tmp, &bv, &bp);CHKERRQ(ierr);

      // set RHS for the velocity field
      ierr = DMDAVecGetArrayDOF(dav, bv, &pbv);CHKERRQ(ierr);
      for (std::size_t k=infov.zs; k<infov.zs+infov.zm; ++k){
        x[2] = k*h[2];
        for (std::size_t j=infov.ys; j<infov.ys+infov.ym; ++j){
          x[1] = j*h[1];
          for (std::size_t i=infov.xs; i<infov.xs+infov.xm; ++i){
            x[0] = i*h[0];
            for (std::size_t d=0; d<infov.dof;++d)
              if (rhs_cond.conditions_[d])
                rhs_cond.conditions_[d](x, &pbv[k][j][i][d]);
          }
        }
      }
      ierr = DMDAVecRestoreArrayDOF(dav, bv, &pbv);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(dm, tmp, &bv, &bp);CHKERRQ(ierr);

      auto hp{h};
      std::for_each(hp.begin(), hp.end(), [](auto x){x*=2;});

      using ctx = problem::context<3, 2>;
      ctx s{dm, {{h, hp }}, mass_mult};
      Mat A = make_matrix<ctx>(&s, diag_block_matrix<ctx>);

      ierr = MatMult(A, tmp, rhs);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(dm, &tmp);CHKERRQ(ierr);
      ierr = MatDestroy(&A);CHKERRQ(ierr);

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