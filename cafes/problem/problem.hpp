#ifndef PARTICLE_PROBLEM_PROBLEM_HPP_INCLUDED
#define PARTICLE_PROBLEM_PROBLEM_HPP_INCLUDED

#include <problem/context.hpp>
#include <petsc.h>
#include <memory>

namespace cafes
{
  namespace problem
  {
    template<std::size_t Dimensions>
    struct Problem
    {
      virtual PetscErrorCode setup_RHS() {return 0;};
      virtual PetscErrorCode setup_KSP() {return 0;};
      virtual PetscErrorCode solve() {return 0;};
    };
  }
}
#endif