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

#ifndef POSTTRAITEMENT_INTERPOLATION_HPP_INCLUDED
#define POSTTRAITEMENT_INTERPOLATION_HPP_INCLUDED

//#include <algorithm/iterate.hpp>
//#include <problem/problem.hpp>
//#include <problem/stokes.hpp>
#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <particle/particle.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/cross_product.hpp>
#include <particle/geometry/position.hpp>
#include <particle/singularity/add_singularity.hpp>
#include <particle/geometry/vector.hpp>
#include <petsc/vec.hpp>

//#include <io/vtk.hpp>


#include <petsc.h>
#include <iostream>
#include <memory>
#include <algorithm>
#include <chrono>

namespace cafes
{
	namespace posttraitement
  {
  	#undef __FUNCT__
    #define __FUNCT__ "linear_interpolation"
    // --------------------------------------------------------------------------
    // Linear interpolation of the solution of a Stokes problem on a rafined mesh                                                
    // --------------------------------------------------------------------------
    //
    //                                    DM dm :   mesh of the solution
    //                                  Vec sol :   solution of the Stokes problem
    //                             DM dm_refine :   refined mesh on wich sol will be interpolated (construct within the function)
    //                           Vec sol_refine :   refined solution (construct within the function)
    // const std::array<int, Dimensions> refine :   refinement factor of the mesh in each direction
    //  const std::array<double, Dimensions>& h :   size of the mesh dm
    template<std::size_t Dimensions>
    PetscErrorCode linear_interpolation(DM& dm, Vec& sol, DM& dm_refine, Vec& sol_refine, const std::array<int, Dimensions>& refine, const std::array<double, Dimensions>& h)
    {
    	PetscErrorCode ierr;
    	PetscFunctionBeginUser;
    	DM dav, dap, dav_refine, dap_refine;
    	Vec solv_refine, solp_refine, coord_dav;
      PetscInt xstart, xend, ystart, yend, zstart, zend;
      PetscInt localsize;

    	DMCompositeGetEntries(dm, &dav, &dap);
      DMClone(dav, &dav_refine);
      DMClone(dap, &dap_refine);
    	ierr = DMDASetRefinementFactor(dav, refine[0], refine[1], PETSC_NULL);CHKERRQ(ierr);
      ierr = DMRefine(dav, PETSC_COMM_WORLD, &dav_refine);
      ierr = DMDASetRefinementFactor(dap, refine[0], refine[1], PETSC_NULL);CHKERRQ(ierr);
      ierr = DMRefine(dap, PETSC_COMM_WORLD, &dap_refine);
      DMCompositeCreate(PETSC_COMM_WORLD, &dm_refine);
      DMCompositeAddDM(dm_refine, dav_refine);
      DMCompositeAddDM(dm_refine, dap_refine);
    	ierr = DMCreateGlobalVector(dm_refine, &sol_refine);CHKERRQ(ierr);
    	auto solv = petsc::petsc_vec<Dimensions>(dm, sol, 0);
      auto solv_ref = petsc::petsc_vec<Dimensions>(dm_refine, sol_refine, 0, false);
      auto solp = petsc::petsc_vec<Dimensions>(dm, sol, 1);
      auto solp_ref = petsc::petsc_vec<Dimensions>(dm_refine, sol_refine, 1, false);

      // RECOVER VECTORS FROM PETSC_VECTORS
    	ierr = solv.global_to_local(INSERT_VALUES);CHKERRQ(ierr);
      ierr = solp.global_to_local(INSERT_VALUES);CHKERRQ(ierr);

      // INTERPOLATION KERNEL
      auto const kernel = [](const auto& solv, auto& solv_ref, const auto& refine, const auto& h)
      {
        auto const kernel_pos = [&](auto const &pos)
        {
          auto ielem = fem::get_element(pos);
          for(std::size_t j=0; j<refine[0]+1; j++)
          {
            for(std::size_t i =0; i<refine[1]+1; i++)
            {
              auto refcoord = pos*h + geometry::position<double, Dimensions>{i*h[0]/refine[0], j*h[1]/refine[1]};
              auto bfunc = fem::P1_integration(geometry::position<double,Dimensions>{refcoord[0] - h[0]*pos[0], refcoord[1] - h[1]*pos[1]}, h);
              auto refpos = geometry::position<int, Dimensions>{i+pos[0]*refine[0],j+pos[1]*refine[1]};
              auto uref = solv_ref.at(refpos);
              for (std::size_t d=0; d<Dimensions; d++)
              {
                uref[d] = 0;
              }
              for (std::size_t k=0; k<bfunc.size(); k++)
              {
                auto u = solv.at(ielem[k]);
                for (std::size_t d=0; d<Dimensions; d++)
                {
                  uref[d] += u[d]*bfunc[k];
                }
              }
            }
          }
        };
        return kernel_pos;
      };

      // VELOCITY INTERPOLATION
      auto box = fem::get_DM_bounds<Dimensions>(dav);
      algorithm::iterate(box, kernel(solv, solv_ref, refine, h));

      // PRESSURE INTERPOLATION
      auto boxp = fem::get_DM_bounds<Dimensions>(dap);
      algorithm::iterate(boxp, kernel(solp, solp_ref, refine, h));

      // LOCAL TO GLOBAL
      ierr = solv_ref.local_to_global(INSERT_VALUES);CHKERRQ(ierr);
      ierr = solp_ref.local_to_global(INSERT_VALUES);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  }
 
}
#endif