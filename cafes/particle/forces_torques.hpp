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

#ifndef CAFES_PARTICLE_FORCES_TORQUES_HPP_INCLUDED
#define CAFES_PARTICLE_FORCES_TORQUES_HPP_INCLUDED

#include <particle/singularity/add_singularity.hpp>

namespace cafes
{
    void set_torques(std::vector<double>& torques, std::size_t ipart,
                     cafes::geometry::position<double, 3> const& center,
                     cafes::geometry::position<int, 3> const& ind,
                     std::array<double, 3> const& g,
                     std::array<double, 3> const& h
                     )
    {
      cafes::geometry::position<double, 3> pos { ind[0]*h[0] - center[0], ind[1]*h[1] - center[1], ind[2]*h[2] - center[2] };
      torques[ipart*3    ] += pos[1]*g[2] - pos[2]*g[1];
      torques[ipart*3 + 1] += pos[2]*g[0] - pos[0]*g[2];
      torques[ipart*3 + 2] += pos[0]*g[1] - pos[1]*g[0];
    }

    void set_torques(std::vector<double>& torques, std::size_t ipart,
                     cafes::geometry::position<double, 2> const& center,
                     cafes::geometry::position<int, 2> const& ind,
                     std::array<double, 2> const& g,
                     std::array<double, 2> const& h
                     )
    {
      cafes::geometry::position<double, 2> pos { ind[0]*h[0] - center[0], ind[1]*h[1] - center[1]};
      torques[ipart] += pos[0]*g[1] - pos[1]*g[0];
    }

    #undef __FUNCT__
    #define __FUNCT__ "forces_torques_with_control"
    template<std::size_t Dimensions, typename Shape>
    PetscErrorCode forces_torques_with_control(std::vector<particle<Shape>> const& particles,
                                               Vec control,
                                               geometry::box<int, Dimensions> const& box,
                                               std::vector<double>& forces,
                                               std::vector<double>& torques,
                                               std::vector<int> const& integration_points_size,
                                               std::array<double, Dimensions> const& h,
                                               bool compute_singularity)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      std::vector<double> local_forces;
      std::vector<double> local_torques;

      local_forces.resize(forces.size());
      local_torques.resize(torques.size());

      std::fill(local_forces.begin(), local_forces.end(), 0.);
      std::fill(local_torques.begin(), local_torques.end(), 0.);

      std::size_t icontrol = 0;
      PetscScalar const *pcontrol;
      ierr = VecGetArrayRead(control, &pcontrol);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<particles.size(); ++ipart)
      {
        auto p = particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::overlap_box(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts){
            std::array<double, Dimensions> g;
            for(std::size_t d=0; d<Dimensions; ++d){
              local_forces[ipart*Dimensions + d] += pcontrol[icontrol];
              g[d] = pcontrol[icontrol++];
            }
            set_torques(local_torques, ipart, p.center_, ind, g, h);
          }       
        }
        for(std::size_t d=0; d<Dimensions; ++d){
          local_forces[ipart*Dimensions + d] *= p.volume()/integration_points_size[ipart];
          local_torques[ipart*Dimensions + d] *= p.volume()/integration_points_size[ipart];
        }

      }
      ierr = VecRestoreArrayRead(control, &pcontrol);CHKERRQ(ierr);

      MPI_Allreduce(local_forces.data(), forces.data(), forces.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(local_torques.data(), torques.data(), torques.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

      if (compute_singularity)
      {
        std::vector<physics::force<Dimensions>> sing_forces;
        sing_forces.resize(particles.size());

        ierr = singularity::compute_singular_forces(particles, sing_forces, h, 100);CHKERRQ(ierr);

        for(std::size_t ipart=0; ipart<particles.size(); ++ipart)
        {
          for(std::size_t d=0; d<Dimensions; ++d)
            forces[ipart*Dimensions + d] += sing_forces[ipart][d];
          // TODO: add singular torques if needed
        }
      }
      PetscFunctionReturn(0);
    }
}

#endif
