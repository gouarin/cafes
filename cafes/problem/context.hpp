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
      static constexpr std::size_t dim = Dimensions;
      template<std::size_t N> using int_ = std::integral_constant<std::size_t, N>;
      using ndm_type = int_<Ndm>;
      using dimension_type = int_<Dimensions>;

      DM dm;
      using array1d = std::array<double, Dimensions>;
      using array2d = std::array<std::array<double, Dimensions>, Ndm>;
      typename std::conditional<Ndm == 1, array1d, array2d>::type h;
      int order;
      PetscErrorCode(*apply)(petsc::petsc_vec<Dimensions>&,
                             petsc::petsc_vec<Dimensions>&,
                             std::array<double, Dimensions> const&,
                             std::size_t);
      PetscErrorCode(*apply_diag)(petsc::petsc_vec<Dimensions>&,
                                  std::array<double, Dimensions> const&, std::size_t) = nullptr;
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
 
      using position_type   = geometry::position<double, Dimensions>;
      using position_type_i = geometry::position<int, Dimensions>;
      using force_type   = physics::force<Dimensions>;

      Problem_type& problem;
      std::vector<particle<Shape>>& particles;
      std::vector<std::vector<std::pair<position_type_i, position_type>>> const& surf_points;
      std::vector<std::vector<geometry::vector<double, Dimensions>>> const& radial_vec;
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
