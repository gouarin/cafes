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
    // Linear interpolation of the solution of a Stokes problem on a refined mesh
    // --------------------------------------------------------------------------
    //
    //                                    DM dm :   mesh of the solution
    //                                  Vec sol :   solution of the Stokes problem
    //                             DM dm_refine :   refined mesh on which sol will be interpolated (construct within the function)
    //                           Vec sol_refine :   refined solution (filled within the function)
    // const std::array<int, Dimensions> refine :   refinement factor of the mesh in each direction
    //  const std::array<double, Dimensions>& h :   size of the mesh dm
    template<std::size_t Dimensions>
    PetscErrorCode linear_interpolation(DM& dm, Vec& sol, DM& dm_refine, Vec& sol_refine, const std::array<int, Dimensions>& refine, const std::array<double, Dimensions>& h)
    {
    	PetscErrorCode ierr;
    	DM dav, dap, dav_refine, dap_refine;
    	Vec solv_refine, solp_refine, coord_dav;
      PetscInt xstart, xend, ystart, yend, zstart, zend;
      PetscInt localsize;
      PetscInt         mpx, mpy, npx, npy, mvx, mvy;
      const PetscInt   *lx, *ly;
      PetscInt         *lxu, *lyu;

      PetscFunctionBeginUser;

      int rank;
      ierr=MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

      /*

      DMCompositeGetEntries(dm, &dav, &dap);
      std::cout << "View dap" << std::endl;
      DMView(dap, PETSC_VIEWER_STDOUT_WORLD);
      DMClone(dap, &dap_refine);
      ierr = DMDASetRefinementFactor(dap, refine[0], refine[1], PETSC_NULL);CHKERRQ(ierr);
      ierr = DMRefine(dap, PETSC_COMM_WORLD, &dap_refine);
      ierr = DMDASetFieldName(dap_refine, 0, "p");CHKERRQ(ierr);
      std::cout << "View dap_refine" << std::endl;
      DMView(dap_refine, PETSC_VIEWER_STDOUT_WORLD);



      ierr = DMDAGetInfo(dap_refine, NULL, &mpx, &mpy, NULL, &npx, &npy, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);
      ierr = DMDAGetOwnershipRanges(dap_refine, &lx, &ly, 0);CHKERRQ(ierr);

      lxu = fem::set_lxu(PETSC_FALSE, npx, lx);
      lyu = fem::set_lxu(PETSC_FALSE, npy, ly);

      mvx = 2*mpx - 1;
      mvy = 2*mpy - 1;

      ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,
                          mvx, mvy, npx, npy,
                          2, 1, lxu, lyu, &dav_refine);CHKERRQ(ierr);
      ierr = DMSetUp(dav_refine);CHKERRQ(ierr);
      std::cout << "View dav" << std::endl;
      DMView(dav, PETSC_VIEWER_STDOUT_WORLD);
      std::cout << "View dav_refine" << std::endl;
      DMView(dav_refine, PETSC_VIEWER_STDOUT_WORLD);
      ierr = DMDASetFieldName(dav_refine, 0, "u");CHKERRQ(ierr);
      ierr = DMDASetFieldName(dav_refine, 1, "v");CHKERRQ(ierr);

      ierr = DMCompositeCreate(PETSC_COMM_WORLD, &dm_refine);CHKERRQ(ierr);
      ierr = DMCompositeAddDM(dm_refine, dav_refine);CHKERRQ(ierr);
      ierr = DMCompositeAddDM(dm_refine, dap_refine);CHKERRQ(ierr);
      ierr = DMSetFromOptions(dm_refine);CHKERRQ(ierr);

      // ierr = DMDestroy(&dap);CHKERRQ(ierr);
      // ierr = DMDestroy(&dav);CHKERRQ(ierr);
      // ierr = PetscFree2(lxu, lyu);CHKERRQ(ierr);

      */



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
          using position_type = geometry::position<double, Dimensions>;
          using position_type_i = geometry::position<int, Dimensions>;
          std::array<double, Dimensions> hr = {{h[0]/refine[0], h[1]/refine[1]}};

          auto ielem = fem::get_element(pos);
          for(std::size_t j=0; j<refine[1]; ++j)
          {
            for(std::size_t i =0; i<refine[0]; ++i)
            {
              position_type_i pts_i = {pos[0]*refine[0] + i, pos[1]*refine[1] + j};
              position_type pts_loc = {i*hr[0], j*hr[1]};
              auto bfunc = fem::P1_integration(pts_loc, h);
              auto uref = solv_ref.at_g(pts_i);

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


      auto const kernelbis = [](const auto& solv, auto& solv_ref, const auto& refine, const auto& h)
      {
        auto const kernel_pos = [&](auto const &pos)
        {
          using position_type = geometry::position<double, Dimensions>;
          using position_type_i = geometry::position<int, Dimensions>;
          std::array<double, Dimensions> hr = {{h[0]/refine[0], h[1]/refine[1]}};
          position_type_i pos_coarse = {std::floor(pos[0]/refine[0]), std::floor(pos[1]/refine[1])};
          auto ielem_coarse = fem::get_element(pos_coarse, 1);
          // std::cout << ielem_coarse[0][0] << " " << ielem_coarse[0][1] << " " << ielem_coarse[1][0] << " " << ielem_coarse[1][1] << " " << ielem_coarse[2][0] << " " << ielem_coarse[2][1] << " "<< ielem_coarse[3][0] << " " << ielem_coarse[3][1]  << std::endl;
          auto ielem = fem::get_element(pos, 1);

          for (std::size_t i=0; i<ielem.size(); ++i)
          {
            auto uref = solv_ref.at(ielem[i]);
            position_type pts_loc = {ielem[i][0]*hr[0]-pos_coarse[0]*h[0], ielem[i][1]*hr[1]-pos_coarse[1]*h[1]};
            auto bfunc = fem::P1_integration(pts_loc, h);
            for (std::size_t d=0; d<Dimensions; d++)
            {
              uref[d] = 0;
            }
            for (std::size_t j=0; j<ielem_coarse.size(); j++)
            {
              auto u = solv.at(ielem_coarse[j]);
              for (std::size_t d=0; d<Dimensions; d++)
              {
                uref[d] += u[d]*bfunc[j];
              }
            }
          }
        };
        return kernel_pos;
      };


      // VELOCITY INTERPOLATION
      auto box = fem::get_DM_bounds<Dimensions>(dav_refine);
      algorithm::iterate(box, kernelbis(solv, solv_ref, refine, h));
      // std::cout << "ici" << " " << rank << std::endl;

      // PRESSURE INTERPOLATION
      std::array<double, Dimensions> hp {2*h[0], 2*h[1]};
      auto boxp = fem::get_DM_bounds<Dimensions>(dap_refine);
      algorithm::iterate(boxp, kernelbis(solp, solp_ref, refine, hp));

      // LOCAL TO GLOBAL
      ierr = solv_ref.local_to_global(INSERT_VALUES);CHKERRQ(ierr);
      ierr = solp_ref.local_to_global(INSERT_VALUES);CHKERRQ(ierr);



      PetscFunctionReturn(0);
    }

  }

}
#endif