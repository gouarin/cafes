#ifndef PARTICLE_PROBLEM_OPTIONS_HPP_INCLUDED
#define PARTICLE_PROBLEM_OPTIONS_HPP_INCLUDED

#include <array>
#include <petsc.h>

namespace cafes
{
  namespace problem
  {
    template<std::size_t Dimensions>
    struct options{
      std::array<PetscInt, Dimensions> mx;
      std::array<double, Dimensions> lx;
      std::array<PetscBool, Dimensions> xperiod;
      PetscBool strain_tensor = PETSC_FALSE;
      PetscBool pmm = PETSC_FALSE;

      options(){
        mx.fill(17);
        lx.fill(1.);
        xperiod.fill(PETSC_FALSE);
      }

      #undef __FUNCT__
      #define __FUNCT__ "process_options_"
      PetscErrorCode process_options_(std::integral_constant<int,2> const&){
        PetscErrorCode ierr;
        PetscFunctionBeginUser;

        ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Problem Options", "");CHKERRQ(ierr);
        ierr = PetscOptionsInt("-mx", "The number of points in x direction", "options.hpp", mx[0], &mx[0], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsInt("-my", "The number of points in y direction", "options.hpp", mx[1], &mx[1], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsReal("-Lx", "The size of x", "options.hpp", lx[0], &lx[0], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsReal("-Ly", "The size of y", "options.hpp", lx[1], &lx[1], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-xperiod", "set periodic condition in x direction", "options.hpp", xperiod[0], &xperiod[0], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-yperiod", "set periodic condition in y direction", "options.hpp", xperiod[1], &xperiod[1], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-strain_tensor", "Use strain tensor term instead of the Laplacian", "options.hpp", strain_tensor, &strain_tensor, nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-PMM", "MG with Pressure Mass Matrix preconditionner", "options.hpp", pmm, &pmm, nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsEnd();CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "process_options_"
      PetscErrorCode process_options_(std::integral_constant<int,3> const&){
        PetscErrorCode ierr;
        PetscFunctionBeginUser;

        ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "", "Problem Options", "");CHKERRQ(ierr);
        ierr = PetscOptionsInt("-mx", "The number of points in x direction", "options.hpp", mx[0], &mx[0], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsInt("-my", "The number of points in y direction", "options.hpp", mx[1], &mx[1], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsInt("-mz", "The number of points in z direction", "options.hpp", mx[2], &mx[2], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsReal("-Lx", "The size of x", "options.hpp", lx[0], &lx[0], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsReal("-Ly", "The size of y", "options.hpp", lx[1], &lx[1], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsReal("-Lz", "The size of z", "options.hpp", lx[2], &lx[2], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-xperiod", "set periodic condition in x direction", "options.hpp", xperiod[0], &xperiod[0], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-yperiod", "set periodic condition in y direction", "options.hpp", xperiod[1], &xperiod[1], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-zperiod", "set periodic condition in z direction", "options.hpp", xperiod[2], &xperiod[2], nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-strain_tensor", "Use strain tensor term instead of the Laplacian", "options.hpp", strain_tensor, &strain_tensor, nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsBool("-PMM", "MG with Pressure Mass Matrix preconditionner", "options.hpp", pmm, &pmm, nullptr);CHKERRQ(ierr);
        ierr = PetscOptionsEnd();CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "process_options"
      PetscErrorCode process_options(){
        PetscErrorCode ierr;
        PetscFunctionBeginUser;
      
        ierr = process_options_(std::integral_constant<int, Dimensions>{});CHKERRQ(ierr);
      
        PetscFunctionReturn(0);
      }

    };
  }
}
#endif
