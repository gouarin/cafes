#ifndef PARTICLE_PROBLEM_CONTEXT_HPP_INCLUDED
#define PARTICLE_PROBLEM_CONTEXT_HPP_INCLUDED

#include<fem/bc.hpp>
#include <array>
#include <iostream>
#include <type_traits>
#include <petsc.h>

namespace cafes
{
  namespace problem
  {
    template<std::size_t Dimensions, std::size_t Ndm=1>
    struct context{
      template<std::size_t N> using int_ = std::integral_constant<std::size_t, N>;
      using ndm_type = int_<Ndm>;
      DM dm;
      using array1d = std::array<double, Dimensions>;
      using array2d = std::array<std::array<double, Dimensions>, Ndm>;
      typename std::conditional<Ndm == 1, array1d, array2d>::type h;
      PetscErrorCode(*apply)(DM, Vec, Vec, std::array<double, Dimensions> const&);
      PetscErrorCode(*apply_diag)(DM, Vec, std::array<double, Dimensions> const&) = nullptr;
      fem::dirichlet_conditions<Dimensions> bc_{};
      bool set_bc_ = false;

      context(context const&) = default;
      context(context&&)      = default;

      context& operator=(context const&) = default;
      context& operator=(context&&)      = default;

      void set_dirichlet_bc(fem::dirichlet_conditions<Dimensions> bc)
      {
        set_bc_ = true;
        bc_ = bc;
      }
    };

  }
}

#endif
