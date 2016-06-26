#ifndef PARTICLE_PROBLEM_CONTEXT_HPP_INCLUDED
#define PARTICLE_PROBLEM_CONTEXT_HPP_INCLUDED

#include <fem/bc.hpp>
#include <particle/particle.hpp>
#include <problem/problem.hpp>
#include <particle/geometry/position.hpp>
#include <petsc/vec.hpp>
#include <array>
#include <iostream>
#include <vector>
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
      using dimension_type = int_<Dimensions>;

      DM dm;
      using array1d = std::array<double, Dimensions>;
      using array2d = std::array<std::array<double, Dimensions>, Ndm>;
      typename std::conditional<Ndm == 1, array1d, array2d>::type h;
      PetscErrorCode(*apply)(petsc::petsc_vec<Dimensions>&,
                             petsc::petsc_vec<Dimensions>&,
                             std::array<double, Dimensions> const&);
      PetscErrorCode(*apply_diag)(petsc::petsc_vec<Dimensions>&,
                                  std::array<double, Dimensions> const&) = nullptr;
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

    template<std::size_t Dimensions, typename Shape, typename Problem_type>
    struct particle_context{
      template<std::size_t N> using int_ = std::integral_constant<std::size_t, N>;
      using dimension_type = int_<Dimensions>;
 
      using position_type   = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;
      using force_type   = physics::force<Dimensions>;

      Problem_type& problem;
      std::vector<particle<Shape>>& particles;
      std::vector<std::vector<std::pair<position_type_i, position_type>>> const& surf_points;
      std::vector<std::vector<position_type>> const& radial_vec;
      std::vector<int> const& nb_surf_points;
      std::vector<int> const& num;
      std::size_t scale;
      bool compute_rhs;
      bool add_rigid_motion;
      bool compute_singularity;
      Vec sol_tmp;
    };
  }
}

#endif
