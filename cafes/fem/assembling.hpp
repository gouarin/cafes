// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef CAFES_FEM_ASSEMBLING_HPP_INCLUDED
#define CAFES_FEM_ASSEMBLING_HPP_INCLUDED

#include <petsc.h>

#include "../algorithm/iterate.hpp"
#include "matElem.hpp"
#include "mesh.hpp"
#include "quadrature.hpp"

inline PetscInt DMDALocalIndex2D(DMDALocalInfo *info, PetscInt i, PetscInt j)
{
    return info->dof * ((j - info->gys) * info->gxm + i - info->gxs);
}

inline PetscInt DMDALocalIndex3D(DMDALocalInfo *info, PetscInt i, PetscInt j,
                                 PetscInt k)
{
    return info->dof *
           (((k - info->gzs) * info->gym + j - info->gys) * info->gxm + i -
            info->gxs);
}

#undef __FUNCT__
#define __FUNCT__ "PreallocateMat"
template<std::size_t Dimensions, typename Context, typename Options>
PetscErrorCode PreallocateMat(Context *user, Options opt, const MatType mtype,
                              Mat *J, std::size_t order, bool precond)
{
    PetscFunctionBegin;
    PetscErrorCode ierr;
    PetscInt *dnz, *onz;
    MPI_Comm comm;
    Mat A;
    ISLocalToGlobalMapping ltog;

    DM pack = user->dm, dav, dap;
    DMDALocalInfo infov, infop;

    ierr = DMCompositeGetEntries(pack, &dav, &dap);CHKERRQ(ierr);

    ierr = DMDAGetLocalInfo(dav, &infov);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);

    PetscInt localsize = infov.dof * infov.xm * infov.ym * infov.zm + infop.xm * infop.ym * infop.zm;
    PetscInt totalsize = infov.dof * infov.mx * infov.my * infov.mz + infop.mx * infop.my * infop.mz;
    PetscInt dec = infov.dof * infov.gxm * infov.gym * infov.gzm;

    ierr = PetscObjectGetComm((PetscObject)pack, &comm);CHKERRQ(ierr);
    ierr = MatCreate(comm, &A);CHKERRQ(ierr);
    ierr = MatSetSizes(A, localsize, localsize, totalsize, totalsize);CHKERRQ(ierr);
    ierr = MatSetType(A, (const MatType)mtype);CHKERRQ(ierr);
    ierr = MatSetDM(A, pack);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);

    // Get the local to global mapping
    ierr = DMGetLocalToGlobalMapping(pack, &ltog);CHKERRQ(ierr);

    // determine the matrix preallocation information
    ierr = MatPreallocateInitialize(comm, localsize, localsize, dnz, onz);CHKERRQ(ierr);

    // preallocation  of the velocity -> block (1, 1)
    if (opt.strain_tensor)
    {
        auto box = cafes::fem::get_DM_bounds<Dimensions>(dav);

        auto const kernel = [&]()
        {
            auto const kernel_pos = [&](auto const &pos)
            {
                auto ielem = cafes::fem::get_indices_tensor(infov, pos, order);
                MatPreallocateSetLocal(ltog, ielem.size(), ielem.data(), ltog,
                                       ielem.size(), ielem.data(), dnz, onz);
            };
            return kernel_pos;
        };

        cafes::algorithm::iterate(box, kernel(), order);

    }
    else
    {
        auto box = cafes::fem::get_DM_bounds<Dimensions>(dav);

        auto kernel = [&]()
        {
            auto const kernel_pos = [&](auto const &pos)
            {
                for (std::size_t d = 0; d < infov.dof; ++d)
                {
                    auto ielem = cafes::fem::get_indices(infov, pos, d, order);
                    MatPreallocateSetLocal(ltog, ielem.size(), ielem.data(), ltog,
                                           ielem.size(), ielem.data(), dnz, onz);
                }
            };
            return kernel_pos;
        };

        cafes::algorithm::iterate(box, kernel(), order);
    }

    {
        // preallocation of the pressure -> block (2, 2)
        auto box = cafes::fem::get_DM_bounds<Dimensions>(dap);

        auto kernel = [&]()
        {
            auto const kernel_pos = [&](auto const &pos)
            {
                auto ielem = cafes::fem::get_indices(infop, pos, dec, 1);
                MatPreallocateSetLocal(ltog, ielem.size(), ielem.data(), ltog,
                                       ielem.size(), ielem.data(), dnz, onz);
            };
            return kernel_pos;
        };
        cafes::algorithm::iterate(box, kernel());
    }

    if (!precond)
    {

        // preallocation of the pressure-velocity interaction -> block (2, 1)
        auto box = cafes::fem::get_DM_bounds<Dimensions>(dap);

        auto kernel = [&]()
        {
            auto const kernel_pos = [&](auto const &pos)
            {
                auto ielem_p = cafes::fem::get_indices(infop, pos, dec, 1);
                auto ielem_v = cafes::fem::get_indices_tensor(infov, 2*pos, 2);

                MatPreallocateSetLocal(ltog, ielem_p.size(), ielem_p.data(), ltog,
                                       ielem_v.size(), ielem_v.data(), dnz, onz);
                MatPreallocateSetLocal(ltog, ielem_v.size(), ielem_v.data(), ltog,
                                       ielem_p.size(), ielem_p.data(), dnz, onz);
            };
            return kernel_pos;
        };
        cafes::algorithm::iterate(box, kernel());
    }

    // Preallocation
    ierr = MatSeqAIJSetPreallocation(A, 0, dnz);CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A, 0, dnz, 0, onz);CHKERRQ(ierr);
    ierr = MatPreallocateFinalize(dnz, onz);CHKERRQ(ierr);

    ierr = MatSetLocalToGlobalMapping(A, ltog, ltog);CHKERRQ(ierr);

    *J = A;
    PetscFunctionReturn(0);
}

namespace cafes
{
    auto const kernel_diag_block_A = [](Mat &A, auto &info, auto const &matelem,
                                      const int order,
                                      const int dec = 0) {
        auto const kernel_pos = [&](auto const &pos) {
            for (std::size_t d = 0; d < info.dof; ++d)
            {
                auto ielem = fem::get_indices(info, pos, d + dec, order);
                MatSetValuesLocal(A, ielem.size(), ielem.data(), ielem.size(),
                                  ielem.data(), matelem.data(), ADD_VALUES);
            }
        };
        return kernel_pos;
    };

    template<std::size_t Dimensions>
    PetscErrorCode laplacian_assembling(DM &dv, Mat &A, const std::array<double, Dimensions> &h, std::size_t order)
    {
        PetscErrorCode ierr;

        DMDALocalInfo info;
        DMDAGetLocalInfo(dv, &info);

        auto box = fem::get_DM_bounds<Dimensions>(dv);

        auto matelem = getMatElemLaplacian(h, order);

        algorithm::iterate(box, kernel_diag_block_A(A, info, matelem, order), order);

        return ierr;
    }

    template<std::size_t Dimensions>
    PetscErrorCode mass_assembling(DM &dm, Mat &A,
                                   const std::array<double, Dimensions> &h,
                                   std::size_t order, int dec, double coeff=1)
    {
        PetscErrorCode ierr;
        DMDALocalInfo info;
        DMDAGetLocalInfo(dm, &info);

        auto box = fem::get_DM_bounds<Dimensions>(dm);

        auto matelem = getMatElemMass(h, order);

        for(std::size_t i=0; i<matelem.size(); ++i)
            matelem[i] *= coeff;

        algorithm::iterate(box, kernel_diag_block_A(A, info, matelem, order, dec), order);

        return ierr;
    }

    auto const kernel_tensor_block_A = [](Mat &A, auto &info,
                                          auto const &matelem,
                                          std::size_t order) {
        auto const kernel_pos = [&](auto const &pos) {
            auto ielem = fem::get_indices_tensor(info, pos, order);
            MatSetValuesLocal(A, ielem.size(), ielem.data(), ielem.size(),
                              ielem.data(), matelem.data(), ADD_VALUES);
        };
        return kernel_pos;
    };

    template<std::size_t Dimensions>
    void strain_tensor_assembling(DM &dv, Mat &A,
                                  const std::array<double, Dimensions> &h, std::size_t order)
    {
        DMDALocalInfo info;
        DMDAGetLocalInfo(dv, &info);

        auto box = fem::get_DM_bounds<Dimensions>(dv);

        auto matelem = getMatElemStrainTensor(h, order);

        algorithm::iterate(box, kernel_tensor_block_A(A, info, matelem, order), order);
    }

    template<std::size_t Dimensions>
    void diagonal_assembling(DM &dm, Mat &A, double value)
    {
        DMDALocalInfo infov, infop;
        DM dav, dap;

        DMCompositeGetEntries(dm, &dav, &dap);

        DMDAGetLocalInfo(dav, &infov);
        DMDAGetLocalInfo(dap, &infop);

        auto box = fem::get_DM_bounds<Dimensions>(dap);

        std::array<double, 16> matelem;
        matelem.fill(value);

        int dec = infov.dof * infov.gxm * infov.gym * infov.gzm;

        algorithm::iterate(box, kernel_diag_block_A(A, infop, matelem, dec));
    }

    auto const kernel_off_diag_block_A = [](auto &A, auto &dec, auto &infov,
                                            auto &infop, const auto &matelemB,
                                            const auto &matelemBT) {
        auto const kernel_pos = [&](auto const &pos) {
            auto ielem_p = fem::get_indices(infop, pos, dec, 1);
            auto ielem_v = fem::get_indices_tensor(infov, 2*pos, 2);

            MatSetValuesLocal(A, ielem_p.size(), ielem_p.data(), ielem_v.size(),
                              ielem_v.data(), matelemB.data(), ADD_VALUES);
            MatSetValuesLocal(A, ielem_v.size(), ielem_v.data(), ielem_p.size(),
                              ielem_p.data(), matelemBT.data(), ADD_VALUES);
        };
        return kernel_pos;
    };

    template<std::size_t Dimensions>
    PetscErrorCode B_BT_assembling(DM &dm, Mat &A,
                                   const std::array<double, Dimensions> &h,
                                   std::size_t order)
    {
        PetscErrorCode ierr;
        DMDALocalInfo infov, infop;
        DM dav, dap;

        ierr = DMCompositeGetEntries(dm, &dav, &dap);CHKERRQ(ierr);

        DMDAGetLocalInfo(dav, &infov);
        DMDAGetLocalInfo(dap, &infop);

        auto box = fem::get_DM_bounds<Dimensions>(dap);

        auto matelemB = getMatElemPressure(h, order);
        auto matelemBT = getMatElemPressureT(h, order);

        int dec = infov.dof * infov.gxm * infov.gym * infov.gzm;

        algorithm::iterate(box, kernel_off_diag_block_A(A, dec, infov, infop,
                                                        matelemB, matelemBT));

        return ierr;
    }
} // namespace cafes
#endif