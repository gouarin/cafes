#ifndef CAFES_PETSC_VEC_HPP_INCLUDED
#define CAFES_PETSC_VEC_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <petsc.h>

#include <iostream>
#include <type_traits>

namespace cafes
{
  namespace petsc
  {
    DM get_DM(DM dm, int entry)
    {
      int nDM;
      DMCompositeGetNumberDM(dm, &nDM);

      DM dm_array[nDM];
      DMCompositeGetEntriesArray(dm, dm_array);
      return dm_array[entry];
    }

    Vec get_Vec(DM dm, Vec v, int entry)
    {
      Vec vtmp[1];
      int entries[1] = {entry};
      DMCompositeGetAccessArray(dm, v, 1, entries, vtmp);
      return vtmp[0];
    }

    void retore_Vec(DM dm, Vec v, Vec v_entry, int entry)
    {
      Vec vtmp[1] = {v_entry};
      int entries[1] = {entry};
      DMCompositeRestoreAccessArray(dm, v, 1, entries, vtmp);
    }

    template<std::size_t Dimensions>
    struct petsc_vec{
      DM dm_;
      DM dm_global_;
      Vec v_;
      Vec v_global_;
      Vec v_entry_;
      std::size_t entry_;
      int dof_;

      bool readonly_;
      PetscBool iscomposite;

      typename std::conditional<Dimensions==2, PetscScalar ***, PetscScalar ****>::type pv_;
      typename std::conditional<Dimensions==2, PetscScalar ***, PetscScalar ****>::type pvg_;

      petsc_vec(petsc_vec const &) = default;
      petsc_vec(petsc_vec &&) = default;

      petsc_vec(DM dm, Vec v, std::size_t entry, bool readonly=true):
      dm_global_{dm}, v_global_{v}, entry_{entry}, readonly_{readonly}
      {
        PetscObjectTypeCompare((PetscObject)dm,DMCOMPOSITE,&iscomposite);
        if (iscomposite)
        {
          v_entry_ = get_Vec(dm, v, entry_);
          dm_ = get_DM(dm, entry);
        }
        else
        {
          v_entry_ = v;
          dm_ = dm;
        }

        DMGetLocalVector(dm_, &v_);

        DMDAGetDof(dm_, &dof_);

        if (readonly_)
        {
          DMDAVecGetArrayDOFRead(dm_, v_, &pv_);
          DMDAVecGetArrayDOFRead(dm_, v_entry_, &pvg_);
        }
        else
        {
          DMDAVecGetArrayDOF(dm_, v_, &pv_);
          DMDAVecGetArrayDOF(dm_, v_entry_, &pvg_);
        }
      }

      double* at(std::array<int, 2> indices){
        return pv_[indices[1]][indices[0]];
      }

      double* at(std::array<int, 3> indices){
        return pv_[indices[2]][indices[1]][indices[0]];
      }

      double const* at(std::array<int, 2> indices) const{
        return pv_[indices[1]][indices[0]];
      }

      double const* at(std::array<int, 3> indices) const{
        return pv_[indices[2]][indices[1]][indices[0]];
      }

      double* at(geometry::position<2, int> indices){
        return pv_[indices[1]][indices[0]];
      }

      double* at(geometry::position<3, int> indices){
        return pv_[indices[2]][indices[1]][indices[0]];
      }

      double const* at(geometry::position<2, int> indices) const{
        return pv_[indices[1]][indices[0]];
      }

      double const* at(geometry::position<3, int> indices) const{
        return pv_[indices[2]][indices[1]][indices[0]];
      }

      double* at_g(geometry::position<2, int> indices){
        return pvg_[indices[1]][indices[0]];
      }

      double* at_g(geometry::position<3, int> indices){
        return pvg_[indices[2]][indices[1]][indices[0]];
      }

      double const* at_g(geometry::position<2, int> indices) const{
        return pvg_[indices[1]][indices[0]];
      }

      double const* at_g(geometry::position<3, int> indices) const{
        return pvg_[indices[2]][indices[1]][indices[0]];
      }
      
      #undef __FUNCT__
      #define __FUNCT__ "local_to_global"
      PetscErrorCode local_to_global(InsertMode iora)
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
        ierr=DMLocalToGlobalBegin(dm_, v_, iora, v_entry_);CHKERRQ(ierr);
        ierr=DMLocalToGlobalEnd(dm_, v_, iora, v_entry_);CHKERRQ(ierr);
        PetscFunctionReturn(0);

      }

      #undef __FUNCT__
      #define __FUNCT__ "global_to_local"
      PetscErrorCode global_to_local(InsertMode iora)
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
        ierr=DMGlobalToLocalBegin(dm_, v_entry_, iora, v_);CHKERRQ(ierr);
        ierr=DMGlobalToLocalEnd(dm_, v_entry_, iora, v_);CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "fill"
      PetscErrorCode fill(double value)
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
        ierr=VecSet(v_, value);CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "fill_global"
      PetscErrorCode fill_global(double value)
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
        ierr=VecSet(v_entry_, value);CHKERRQ(ierr);
        PetscFunctionReturn(0);
      }

      ~petsc_vec()
      {
        if (readonly_)
        {
          DMDAVecRestoreArrayDOFRead(dm_, v_, &pv_);
          DMDAVecRestoreArrayDOFRead(dm_, v_entry_, &pvg_);
        }
        else
        {
          DMDAVecRestoreArrayDOF(dm_, v_, &pv_);
          DMDAVecRestoreArrayDOF(dm_, v_entry_, &pvg_);
        }

        DMRestoreLocalVector(dm_, &v_);
        if (iscomposite)
        {
          retore_Vec(dm_global_, v_global_, v_entry_, entry_);
        } 
      }
    };
  }
}
#endif