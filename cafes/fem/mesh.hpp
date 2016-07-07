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

#ifndef CAFES_FEM_MESH_HPP_INCLUDED
#define CAFES_FEM_MESH_HPP_INCLUDED

#include <particle/geometry/box.hpp>
#include <petsc.h>
#include <iostream>
#include <array>
#include <algorithm>

namespace cafes
{
  namespace fem
  {
    auto get_dof(DM const& dm)
    {
      DMDALocalInfo info;
      DMDAGetLocalInfo(dm, &info);
      return info.dof;
    }

    auto get_DM_bounds_(DM const& dm, std::integral_constant<int, 2>, bool remove_final_points)
    {
      DMDALocalInfo info;
      DMDAGetLocalInfo(dm, &info);
      
      geometry::box<int, 2> box{ { info.xs          , info.ys           },
                                 { info.xs + info.xm, info.ys + info.ym }
                                };

      // For the not periodic boundary condition.
      if (remove_final_points)
      {
        if(box.upper_right[1] == info.my && info.by != DM_BOUNDARY_PERIODIC)
          --box.upper_right[1];
        if(box.upper_right[0] == info.mx && info.bx != DM_BOUNDARY_PERIODIC)
          --box.upper_right[0];
      }
      return box;
    }

    auto get_DM_bounds_(DM const& dm, std::integral_constant<int, 3>, bool remove_final_points)
    {
      DMDALocalInfo info;
      DMDAGetLocalInfo(dm, &info);

      geometry::box<int, 3> box{ { info.xs          , info.ys          , info.zs           },
                                 { info.xs + info.xm, info.ys + info.ym, info.zs + info.zm }
                                };
      
      // For the not periodic boundary condition.
      if (remove_final_points)
      {
        if(box.upper_right[2] == info.mz && info.bz != DM_BOUNDARY_PERIODIC)
          --box.upper_right[2];
        if(box.upper_right[1] == info.my && info.by != DM_BOUNDARY_PERIODIC)
          --box.upper_right[1];
        if(box.upper_right[0] == info.mx && info.bx != DM_BOUNDARY_PERIODIC)
          --box.upper_right[0];
      }
      return box;
    }

    template<int Dimensions>
    auto get_DM_bounds(DM const& dm, bool remove_final_points=true)
    {
      return get_DM_bounds_(dm, std::integral_constant<int, Dimensions>{}, remove_final_points);
    }

    template<int Dimensions>
    auto get_DM_bounds(DM const& dm, int i, bool remove_final_points=true)
    {
      int ndm;
      DMCompositeGetNumberDM(dm, &ndm);
      DM dmc[ndm];
      DMCompositeGetEntriesArray(dm, dmc);
      return get_DM_bounds_(dmc[i], std::integral_constant<int, Dimensions>{}, remove_final_points);
    }

    std::array<int, 2> get_global_bounds_(DM const& dm, std::integral_constant<int, 2>)
    {
      DMDALocalInfo info;
      DMDAGetLocalInfo(dm, &info);
      return {{info.mx, info.my}};
    }

    std::array<int, 3> get_global_bounds_(DM const& dm, std::integral_constant<int, 3>)
    {
      DMDALocalInfo info;
      DMDAGetLocalInfo(dm, &info);
      return {{info.mx, info.my, info.mz}};
    }

    template<int Dimensions>
    auto get_global_bounds(DM const& dm)
    {
      return get_global_bounds_(dm, std::integral_constant<int, Dimensions>{});
    }

    std::array<DMBoundaryType, 2> get_boundary_type_(DM const& dm, std::integral_constant<int, 2>)
    {
      DMDALocalInfo info;
      DMDAGetLocalInfo(dm, &info);
      return {{info.bx, info.by}};
    }

    std::array<DMBoundaryType, 3> get_boundary_type_(DM const& dm, std::integral_constant<int, 3>)
    {
      DMDALocalInfo info;
      DMDAGetLocalInfo(dm, &info);
      return {{info.bx, info.by, info.bz}};
    }

    template<int Dimensions>
    auto get_boundary_type(DM const& dm)
    {
      return get_boundary_type_(dm, std::integral_constant<int, Dimensions>{});
    }

    #undef __FUNCT__
    #define __FUNCT__ "get_DM_sizes"
    PetscErrorCode get_DM_sizes(DM const& dm, int & localsize, int & totalsize)
    {
        PetscErrorCode ierr;
        DMDALocalInfo  info;
        PetscBool      dm_composite;
        PetscFunctionBeginUser;

        localsize = 0;
        totalsize = 0;

        ierr = PetscObjectTypeCompare((PetscObject)dm, DMCOMPOSITE, &dm_composite);CHKERRQ(ierr);
        if (dm_composite){
            int ndm;
            ierr = DMCompositeGetNumberDM(dm, &ndm);CHKERRQ(ierr);

            DM dmc[ndm];
            ierr = DMCompositeGetEntriesArray(dm, dmc);CHKERRQ(ierr);

            for(std::size_t i=0; i<ndm; ++i){
                ierr = DMDAGetLocalInfo(dmc[i], &info);CHKERRQ(ierr);
                localsize += info.dof*info.xm*info.ym*info.zm;
                totalsize += info.dof*info.mx*info.my*info.mz;
            }

        }
        else{
            ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
            localsize = info.dof*info.xm*info.ym*info.zm;
            totalsize = info.dof*info.mx*info.my*info.mz;
        }

        PetscFunctionReturn(0);
    }

    PetscInt* set_lxu(bool period, int nx, const PetscInt* lx)
    {
        PetscInt *lxu;
        PetscMalloc1(nx, &lxu);
        for(std::size_t i=0; i<nx; ++i) lxu[i] = 2*lx[i];
        if (!period)
          lxu[nx-1]--;
        return lxu;
    }

    #undef __FUNCT__
    #define __FUNCT__ "createDMDA"
    PetscErrorCode createDMDA(DM& dm, 
                              std::array<int, 2> const& mpres, 
                              std::array<int, 2> const& mvel, 
                              std::array<DMBoundaryType, 2> const& b_type, 
                              std::array<PetscBool, 2> const& period)
    {
        PetscErrorCode   ierr;
        DM               DAPressure, DAVelocity;
        PetscInt         npx, npy;
        const PetscInt   *lx, *ly;
        PetscInt         *lxu, *lyu;

        PetscFunctionBeginUser;

        ierr = DMDACreate2d(PETSC_COMM_WORLD, b_type[0], b_type[1], DMDA_STENCIL_BOX,
                            mpres[0], mpres[1], PETSC_DECIDE, PETSC_DECIDE,
                            1, 1, 0, 0, &DAPressure);CHKERRQ(ierr);
        ierr = DMDASetFieldName(DAPressure, 0, "p");CHKERRQ(ierr);

        ierr = DMDAGetInfo(DAPressure, NULL, NULL, NULL, NULL, &npx, &npy, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);
        ierr = DMDAGetOwnershipRanges(DAPressure, &lx, &ly, 0);CHKERRQ(ierr);

        /* Compute the number of points for the velacity field using the decomposition
           used by the pressure and set by PETSc */
        lxu = set_lxu(period[0], npx, lx);
        lyu = set_lxu(period[1], npy, ly);

        /* set velocity grid*/
        ierr = DMDACreate2d(PETSC_COMM_WORLD, b_type[0], b_type[1], DMDA_STENCIL_BOX,
                            mvel[0], mvel[1], npx, npy,
                            2, 1, lxu, lyu, &DAVelocity);CHKERRQ(ierr);
        ierr = DMDASetFieldName(DAVelocity, 0, "u");CHKERRQ(ierr);
        ierr = DMDASetFieldName(DAVelocity, 1, "v");CHKERRQ(ierr);
        
        ierr = DMCompositeAddDM(dm, DAVelocity);CHKERRQ(ierr);
        ierr = DMCompositeAddDM(dm, DAPressure);CHKERRQ(ierr);
        ierr = DMSetFromOptions(dm);CHKERRQ(ierr);

        ierr = DMDestroy(&DAPressure);CHKERRQ(ierr);
        ierr = DMDestroy(&DAVelocity);CHKERRQ(ierr);
        ierr = PetscFree2(lxu, lyu);CHKERRQ(ierr);

        PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "createDMDA"
    PetscErrorCode createDMDA(DM& dm, 
                              std::array<int, 2> const& mvel, 
                              std::array<DMBoundaryType, 2> const& b_type, 
                              std::array<PetscBool, 2> const& period)
    {
        PetscErrorCode   ierr;
        PetscFunctionBeginUser;

        ierr = DMDACreate2d(PETSC_COMM_WORLD, b_type[0], b_type[1], DMDA_STENCIL_BOX,
                            mvel[0], mvel[1], PETSC_DECIDE, PETSC_DECIDE,
                            1, 1, 0, 0, &dm);CHKERRQ(ierr);
        ierr = DMDASetFieldName(dm, 0, "u");CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "createDMDA"
    PetscErrorCode createDMDA(DM& dm, 
                              std::array<int, 3> const& mpres, 
                              std::array<int, 3> const& mvel, 
                              std::array<DMBoundaryType, 3> const& b_type, 
                              std::array<PetscBool, 3> const& period)
    {
        PetscErrorCode   ierr;
        DM               DAPressure, DAVelocity;
        PetscInt         npx, npy, npz;
        const PetscInt   *lx, *ly, *lz;
        PetscInt         *lxu, *lyu, *lzu;

        PetscFunctionBeginUser;

        ierr = DMDACreate3d(PETSC_COMM_WORLD, b_type[0], b_type[1], b_type[2], DMDA_STENCIL_BOX,
                            mpres[0], mpres[1], mpres[2], PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                            1, 1, 0, 0, 0, &DAPressure);CHKERRQ(ierr);
        ierr = DMDASetFieldName(DAPressure, 0, "p");CHKERRQ(ierr);

        ierr = DMDAGetInfo(DAPressure, NULL, NULL, NULL, NULL, &npx, &npy, &npz, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);
        ierr = DMDAGetOwnershipRanges(DAPressure, &lx, &ly, &lz);CHKERRQ(ierr);

        /* Compute the number of points for the velacity field using the decomposition
           used by the pressure and set by PETSc */
        lxu = set_lxu(period[0], npx, lx);
        lyu = set_lxu(period[1], npy, ly);
        lzu = set_lxu(period[2], npz, lz);

        /* set velocity grid*/
        ierr = DMDACreate3d(PETSC_COMM_WORLD, b_type[0], b_type[1], b_type[2], DMDA_STENCIL_BOX,
                            mvel[0], mvel[1], mvel[2], npx, npy, npz,
                            3, 1, lxu, lyu, lzu, &DAVelocity);CHKERRQ(ierr);
        ierr = DMDASetFieldName(DAVelocity, 0, "u");CHKERRQ(ierr);
        ierr = DMDASetFieldName(DAVelocity, 1, "v");CHKERRQ(ierr);
        ierr = DMDASetFieldName(DAVelocity, 2, "w");CHKERRQ(ierr);

        
        ierr = DMCompositeAddDM(dm, DAVelocity);CHKERRQ(ierr);
        ierr = DMCompositeAddDM(dm, DAPressure);CHKERRQ(ierr);
        ierr = DMSetFromOptions(dm);CHKERRQ(ierr);

        ierr = DMDestroy(&DAPressure);CHKERRQ(ierr);
        ierr = DMDestroy(&DAVelocity);CHKERRQ(ierr);
        ierr = PetscFree3(lxu, lyu, lzu);CHKERRQ(ierr);

        PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "createDMDA"
    PetscErrorCode createDMDA(DM& dm, 
                              std::array<int, 3> const& mvel, 
                              std::array<DMBoundaryType, 3> const& b_type, 
                              std::array<PetscBool, 3> const& period)
    {
        PetscErrorCode   ierr;
        PetscFunctionBeginUser;

        ierr = DMDACreate3d(PETSC_COMM_WORLD, b_type[0], b_type[1], b_type[2], DMDA_STENCIL_BOX,
                            mvel[0], mvel[1], mvel[2], PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                            1, 1, 0, 0, 0, &dm);CHKERRQ(ierr);
        ierr = DMDASetFieldName(dm, 0, "u");CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "createMesh"
    template<int Dimensions>
    PetscErrorCode createMesh(DM& dm, std::array<int, Dimensions> const& npoints, std::array<PetscBool, Dimensions> const& period)
    {
        PetscErrorCode   ierr;
        auto mpres{npoints};
        auto mvel{npoints};
        std::array<DMBoundaryType, Dimensions> b_type;
        PetscFunctionBeginUser;

        b_type.fill(DM_BOUNDARY_NONE);

        std::for_each(mvel.begin(), mvel.end(), [](auto & i){i=2*i-1;});

        for(std::size_t i=0; i<Dimensions; ++i)
            if(period[i]){
                b_type[i] = DM_BOUNDARY_PERIODIC;
                --mpres[i];
                --mvel[i];
            }

        ierr = DMCompositeCreate(PETSC_COMM_WORLD, &dm);CHKERRQ(ierr);
        ierr = createDMDA(dm, mpres, mvel, b_type, period);CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "createLaplacianMesh"
    template<int Dimensions>
    PetscErrorCode createLaplacianMesh(DM& dm, std::array<int, Dimensions> const& npoints, std::array<PetscBool, Dimensions> const& period)
    {
        PetscErrorCode   ierr;
        auto mvel{npoints};
        std::array<DMBoundaryType, Dimensions> b_type;
        PetscFunctionBeginUser;

        b_type.fill(DM_BOUNDARY_NONE);

        for(std::size_t i=0; i<Dimensions; ++i)
            if(period[i]){
                b_type[i] = DM_BOUNDARY_PERIODIC;
                --mvel[i];
            }

        ierr = createDMDA(dm, mvel, b_type, period);CHKERRQ(ierr);
        PetscFunctionReturn(0);
    }
  }
}

#endif