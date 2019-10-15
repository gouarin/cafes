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

#define SetInRange(i, m) ((i < 0) ? m + i : ((i >= m) ? i - m : i))

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
template<typename Context, typename Options>
PetscErrorCode PreallocateMat(Context *user, Options opt, const MatType mtype,
                              Mat *J, bool precond)
{
    PetscFunctionBegin;
    PetscErrorCode ierr;
    PetscInt *dnz, *onz;
    MPI_Comm comm;
    Mat A;
    ISLocalToGlobalMapping ltog;
    PetscScalar *values;

    DM pack = user->dm, dav, dap;
    DMDALocalInfo infov, infop;
    PetscInt localsize, totalsize;

    PetscInt slotu, slotp, colu, colp, cnt;
    PetscInt stencilu = 2, stencilp = 1, s = 1;
    PetscInt lstart, lend, pstart, pend, rstart, rend;
    PetscInt i, j, k, kd, l, p, r;
    PetscInt *rowsu = PETSC_NULL, *colsu = PETSC_NULL, rowsuu;
    PetscInt *rowsp = PETSC_NULL, *colsp = PETSC_NULL;
    PetscInt dec;
    PetscInt indxEnd, indyEnd, indzEnd;

    ierr = DMCompositeGetEntries(pack, &dav, &dap);
    CHKERRQ(ierr);

    ierr = DMDAGetLocalInfo(dav, &infov);
    CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dap, &infop);
    CHKERRQ(ierr);

    localsize = infov.dof * infov.xm * infov.ym * infov.zm +
                infop.xm * infop.ym * infop.zm;
    totalsize = infov.dof * infov.mx * infov.my * infov.mz +
                infop.mx * infop.my * infop.mz;
    dec = infov.dof * infov.gxm * infov.gym * infov.gzm;

    // allocation de la matrice du probleme
    ierr = PetscObjectGetComm((PetscObject)pack, &comm);
    CHKERRQ(ierr);
    ierr = MatCreate(comm, &A);
    CHKERRQ(ierr);
    ierr = MatSetSizes(A, localsize, localsize, totalsize, totalsize);
    CHKERRQ(ierr);
    ierr = MatSetType(A, (const MatType)mtype);
    CHKERRQ(ierr);
    ierr = MatSetDM(A, pack);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);

    // Get the local to global mapping
    ierr = DMGetLocalToGlobalMapping(pack, &ltog);
    CHKERRQ(ierr);
    // MGK ierr = DMGetLocalToGlobalMappingBlock(pack, &ltogb);CHKERRQ(ierr);

    // determine the matrix preallocation information
    ierr = MatPreallocateInitialize(comm, localsize, localsize, dnz, onz);
    CHKERRQ(ierr);

    /*
      A  tB
      B  0
    */

    colu = 2 * stencilu + 1;
    colp = 2 * stencilp + 1;
    for (i = 1; i < user->dim; i++)
    {
        colu *= 2 * stencilu + 1;
        colp *= 2 * stencilp + 1;
    }

    ierr = PetscMalloc2(infov.dof, &rowsu, infov.dof * colu, &colsu);
    CHKERRQ(ierr);

    indxEnd = infov.xs + infov.xm;
    indyEnd = infov.ys + infov.ym;
    indzEnd = infov.zs + infov.zm;

    // preallocation  of the velocity -> block (1, 1)
    if (opt.strain_tensor)
    {
        for (i = infov.xs; i < indxEnd; i++)
        {
            pstart =
                (infov.bx == DM_BOUNDARY_PERIODIC) ? -s : (PetscMax(-s, -i));
            pend = (infov.bx == DM_BOUNDARY_PERIODIC)
                       ? s
                       : (PetscMin(s, infov.mx - i - 1));

            for (j = infov.ys; j < indyEnd; j++)
            {
                lstart = (infov.by == DM_BOUNDARY_PERIODIC)
                             ? -s
                             : (PetscMax(-s, -j));
                lend = (infov.by == DM_BOUNDARY_PERIODIC)
                           ? s
                           : (PetscMin(s, infov.my - j - 1));

                for (k = infov.zs; k < indzEnd; k++)
                {
                    rstart = (infov.bz == DM_BOUNDARY_PERIODIC)
                                 ? -s
                                 : (PetscMax(-s, -k));
                    rend = (infov.bz == DM_BOUNDARY_PERIODIC)
                               ? s
                               : (PetscMin(s, infov.mz - k - 1));
                    cnt = 0;
                    for (kd = 0; kd < infov.dof; kd++)
                    {
                        for (r = rstart; r < rend + 1; r++)
                        {
                            for (l = lstart; l < lend + 1; l++)
                            {
                                for (p = pstart; p < pend + 1; p++)
                                {
                                    colsu[cnt++] =
                                        kd + DMDALocalIndex3D(&infov, i + p,
                                                              j + l, k + r);
                                }
                            }
                        }
                        rowsu[kd] = kd + DMDALocalIndex3D(&infov, i, j, k);
                    }
                    ierr = MatPreallocateSetLocal(ltog, infov.dof, rowsu, ltog,
                                                  cnt, colsu, dnz, onz);
                    CHKERRQ(ierr);
                }
            }
        }
    }
    else
    {
        for (i = infov.xs; i < indxEnd; i++)
        {
            pstart =
                (infov.bx == DM_BOUNDARY_PERIODIC) ? -s : (PetscMax(-s, -i));
            pend = (infov.bx == DM_BOUNDARY_PERIODIC)
                       ? s
                       : (PetscMin(s, infov.mx - i - 1));

            for (j = infov.ys; j < indyEnd; j++)
            {
                lstart = (infov.by == DM_BOUNDARY_PERIODIC)
                             ? -s
                             : (PetscMax(-s, -j));
                lend = (infov.by == DM_BOUNDARY_PERIODIC)
                           ? s
                           : (PetscMin(s, infov.my - j - 1));

                for (k = infov.zs; k < indzEnd; k++)
                {
                    rstart = (infov.bz == DM_BOUNDARY_PERIODIC)
                                 ? -s
                                 : (PetscMax(-s, -k));
                    rend = (infov.bz == DM_BOUNDARY_PERIODIC)
                               ? s
                               : (PetscMin(s, infov.mz - k - 1));

                    for (kd = 0; kd < infov.dof; kd++)
                    {
                        cnt = 0;
                        for (r = rstart; r < rend + 1; r++)
                            for (l = lstart; l < lend + 1; l++)
                                for (p = pstart; p < pend + 1; p++)
                                    colsu[cnt++] =
                                        kd + DMDALocalIndex3D(&infov, i + p,
                                                              j + l, k + r);
                        rowsuu = kd + DMDALocalIndex3D(&infov, i, j, k);
                        ierr = MatPreallocateSetLocal(ltog, 1, &rowsuu, ltog,
                                                      cnt, colsu, dnz, onz);
                        CHKERRQ(ierr);
                    }
                }
            }
        }
    }
    ierr = PetscFree2(rowsu, colsu);
    CHKERRQ(ierr);

    ierr = PetscMalloc2(infov.dof * colu, &rowsu, infov.dof * colu, &colsu);
    CHKERRQ(ierr);
    ierr = PetscMalloc2(colp, &rowsp, colp, &colsp);
    CHKERRQ(ierr);

    // if (precond)
    {
        // preallocation of the pressure -> block (2, 2)
        indxEnd = infop.xs + infop.xm;
        indyEnd = infop.ys + infop.ym;
        indzEnd = infop.zs + infop.zm;

        for (i = infop.xs; i < indxEnd; i++)
        {
            pstart =
                (infop.bx == DM_BOUNDARY_PERIODIC) ? -s : (PetscMax(-s, -i));
            pend = (infop.bx == DM_BOUNDARY_PERIODIC)
                       ? s
                       : (PetscMin(s, infop.mx - i - 1));

            for (j = infop.ys; j < indyEnd; j++)
            {
                lstart = (infop.by == DM_BOUNDARY_PERIODIC)
                             ? -s
                             : (PetscMax(-s, -j));
                lend = (infop.by == DM_BOUNDARY_PERIODIC)
                           ? s
                           : (PetscMin(s, infop.my - j - 1));

                for (k = infop.zs; k < indzEnd; k++)
                {
                    rstart = (infop.bz == DM_BOUNDARY_PERIODIC)
                                 ? -s
                                 : (PetscMax(-s, -k));
                    rend = (infop.bz == DM_BOUNDARY_PERIODIC)
                               ? s
                               : (PetscMin(s, infop.mz - k - 1));

                    cnt = 0;
                    for (r = rstart; r < rend + 1; r++)
                        for (l = lstart; l < lend + 1; l++)
                            for (p = pstart; p < pend + 1; p++)
                                colsp[cnt++] =
                                    dec + DMDALocalIndex3D(&infop, i + p, j + l,
                                                           k + r);
                    rowsp[0] = dec + DMDALocalIndex3D(&infop, i, j, k);
                    ierr = MatPreallocateSetLocal(ltog, 1, rowsp, ltog, cnt,
                                                  colsp, dnz, onz);
                    CHKERRQ(ierr);
                }
            }
        }
    }
    // else
    if (!precond)
    {
        // preallocation of the pressure-velocity interaction -> block (2, 1)
        PetscInt pst, lst, rst, icnt, ii, kk, jj, rank;
        PetscInt px, py, pz;
        PetscInt npx, npy, npz;

        ierr = DMDAGetInfo(dav, NULL, NULL, NULL, NULL, &npx, &npy, &npz, NULL,
                           NULL, NULL, NULL, NULL, NULL);
        CHKERRQ(ierr);

        px = (npx == 1) ? 0 : infov.gxs;
        py = (npy == 1) ? 0 : infov.gys;
        pz = (npz == 1) ? 0 : infov.gzs;

        for (k = infop.zs; k < infop.zs + infop.zm; k++)
        {
            rstart = (infop.bz == DM_BOUNDARY_PERIODIC)
                         ? -2 * s
                         : PetscMax(-2 * s, -2 * k);
            rend = (infop.bz == DM_BOUNDARY_PERIODIC)
                       ? 2 * s
                       : PetscMin(2 * s, infov.mz - 2 * k - 1);

            for (j = infop.ys; j < infop.ys + infop.ym; j++)
            {
                lstart = (infop.by == DM_BOUNDARY_PERIODIC)
                             ? -2 * s
                             : PetscMax(-2 * s, -2 * j);
                lend = (infop.by == DM_BOUNDARY_PERIODIC)
                           ? 2 * s
                           : PetscMin(2 * s, infov.my - 2 * j - 1);

                for (i = infop.xs; i < infop.xs + infop.xm; i++)
                {
                    pstart = (infop.bx == DM_BOUNDARY_PERIODIC)
                                 ? -2 * s
                                 : PetscMax(-2 * s, -2 * i);
                    pend = (infop.bx == DM_BOUNDARY_PERIODIC)
                               ? 2 * s
                               : PetscMin(2 * s, infov.mx - 2 * i - 1);

                    cnt = 0;
                    for (r = rstart; r < rend + 1; r++)
                    {
                        kk = (infop.bz == DM_BOUNDARY_PERIODIC)
                                 ? SetInRange(2 * k + r - infov.gzs, infov.mz)
                                 : 2 * k + r - infov.gzs;
                        for (l = lstart; l < lend + 1; l++)
                        {
                            jj = (infop.by == DM_BOUNDARY_PERIODIC)
                                     ? SetInRange(2 * j + l - infov.gys,
                                                  infov.my)
                                     : 2 * j + l - infov.gys;
                            for (p = pstart; p < pend + 1; p++)
                            {
                                ii = (infop.bx == DM_BOUNDARY_PERIODIC)
                                         ? SetInRange(2 * i + p - infov.gxs,
                                                      infov.mx)
                                         : 2 * i + p - infov.gxs;

                                if (ii + px < infov.xs ||
                                    ii + px > infov.xs + infov.xm ||
                                    jj + py < infov.ys ||
                                    jj + py > infov.ys + infov.ym ||
                                    kk + pz < infov.zs ||
                                    kk + pz > infov.zs + infov.zm)
                                {
                                    for (kd = 0; kd < infov.dof; kd++)
                                        colsu[cnt++] = -1;
                                }
                                else
                                {
                                    slotu =
                                        (kk * infov.gym + jj) * infov.gxm + ii;
                                    for (kd = 0; kd < infov.dof; kd++)
                                        colsu[cnt++] = infov.dof * slotu + kd;
                                }
                            }
                        }
                    }
                    rowsp[0] = dec + DMDALocalIndex3D(&infop, i, j, k);

                    ierr = MatPreallocateSetLocal(ltog, 1, rowsp, ltog, cnt,
                                                  colsu, dnz, onz);
                    CHKERRQ(ierr);
                }
            }
        }

        indxEnd = infov.xs + infov.xm;
        indyEnd = infov.ys + infov.ym;
        indzEnd = infov.zs + infov.zm;

        // preallocation of the velocity-pressure interaction -> block (1, 2)
        for (i = infov.xs; i < indxEnd; i++)
        {
            pstart = (infop.bx == DM_BOUNDARY_PERIODIC)
                         ? -(1 - (i & 1))
                         : (PetscMax(-(1 - (i & 1)), -i / 2));
            pend = (infop.bx == DM_BOUNDARY_PERIODIC)
                       ? s
                       : (PetscMin(s, infop.mx - i / 2 - 1));
            for (j = infov.ys; j < indyEnd; j++)
            {
                lstart = (infop.by == DM_BOUNDARY_PERIODIC)
                             ? -(1 - (j & 1))
                             : (PetscMax(-(1 - (j & 1)), -j / 2));
                lend = (infop.by == DM_BOUNDARY_PERIODIC)
                           ? s
                           : (PetscMin(s, infop.my - j / 2 - 1));

                for (k = infov.zs; k < indzEnd; k++)
                {
                    rstart = (infop.bz == DM_BOUNDARY_PERIODIC)
                                 ? -(1 - (k & 1))
                                 : (PetscMax(-(1 - (k & 1)), -k / 2));
                    rend = (infop.bz == DM_BOUNDARY_PERIODIC)
                               ? s
                               : (PetscMin(s, infop.mz - k / 2 - 1));

                    slotu = ((k - infov.gzs) * infov.gym + j - infov.gys) *
                                infov.gxm +
                            (i - infov.gxs);
                    slotp =
                        ((k / 2 - infop.gzs) * infop.gym + j / 2 - infop.gys) *
                            infop.gxm +
                        (i / 2 - infop.gxs);

                    cnt = 0;
                    for (r = rstart; r < rend + 1; r++)
                        for (l = lstart; l < lend + 1; l++)
                            for (p = pstart; p < pend + 1; p++)
                                colsp[cnt++] = dec + slotp +
                                               infop.gxm * (l + infop.gym * r) +
                                               p;

                    for (kd = 0; kd < infov.dof; kd++)
                        rowsu[kd] = kd + infov.dof * slotu;

                    ierr = MatPreallocateSetLocal(ltog, infov.dof, rowsu, ltog,
                                                  cnt, colsp, dnz, onz);
                    CHKERRQ(ierr);
                }
            }
        }
    }
    ierr = PetscFree2(rowsu, colsu);
    CHKERRQ(ierr);
    ierr = PetscFree2(rowsp, colsp);
    CHKERRQ(ierr);

    // Preallocation
    ierr = MatSeqAIJSetPreallocation(A, 0, dnz);
    CHKERRQ(ierr);
    ierr = MatMPIAIJSetPreallocation(A, 0, dnz, 0, onz);
    CHKERRQ(ierr);
    ierr = MatPreallocateFinalize(dnz, onz);
    CHKERRQ(ierr);

    ierr = MatSetLocalToGlobalMapping(A, ltog, ltog);
    CHKERRQ(ierr);
    *J = A;
    PetscFunctionReturn(0);
}

namespace cafes
{
    auto const kernel_diag_block_A = [](Mat &A, auto &info, auto const &matelem,
                                        int const &dec = 0) {
        auto const kernel_pos = [&](auto const &pos) {
            for (std::size_t d = 0; d < info.dof; ++d)
            {
                auto ielem = fem::get_indices(info, pos, d + dec);
                MatSetValuesLocal(A, ielem.size(), ielem.data(), ielem.size(),
                                  ielem.data(), matelem.data(), ADD_VALUES);
            }
        };
        return kernel_pos;
    };

    template<std::size_t Dimensions>
    PetscErrorCode laplacian_assembling(DM &dv, Mat &A,
                                        const std::array<double, Dimensions> &h)
    {
        PetscErrorCode ierr;

        DMDALocalInfo info;
        DMDAGetLocalInfo(dv, &info);

        auto box = fem::get_DM_bounds<Dimensions>(dv);

        auto matelem = getMatElemLaplacianA(h);

        algorithm::iterate(box, kernel_diag_block_A(A, info, matelem));

        return ierr;
    }

    template<std::size_t Dimensions>
    PetscErrorCode mass_assembling(DM &dm, Mat &A,
                                   const std::array<double, Dimensions> &h)
    {
        PetscErrorCode ierr;
        DMDALocalInfo infov, infop;
        DM dav, dap;

        ierr = DMCompositeGetEntries(dm, &dav, &dap);
        CHKERRQ(ierr);

        DMDAGetLocalInfo(dav, &infov);
        DMDAGetLocalInfo(dap, &infop);

        auto box = fem::get_DM_bounds<Dimensions>(dap);

        auto matelem = getMatElemMassA(h);

        int dec = infov.dof * infov.gxm * infov.gym * infov.gzm;

        algorithm::iterate(box, kernel_diag_block_A(A, infop, matelem, dec));

        return ierr;
    }

    auto const kernel_tensor_block_A = [](Mat &A, auto &info,
                                          auto const &matelem) {
        auto const kernel_pos = [&](auto const &pos) {
            auto ielem = fem::get_indices_tensor(info, pos);
            MatSetValuesLocal(A, ielem.size(), ielem.data(), ielem.size(),
                              ielem.data(), matelem.data(), ADD_VALUES);
        };
        return kernel_pos;
    };

    template<std::size_t Dimensions>
    void strain_tensor_assembling(DM &dv, Mat &A,
                                  const std::array<double, Dimensions> &h)
    {
        DMDALocalInfo info;
        DMDAGetLocalInfo(dv, &info);

        auto box = fem::get_DM_bounds<Dimensions>(dv);

        auto matelem = getMatElemStrainTensorA(h);

        algorithm::iterate(box, kernel_tensor_block_A(A, info, matelem));
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
        matelem.fill(0);

        int dec = infov.dof * infov.gxm * infov.gym * infov.gzm;

        algorithm::iterate(box, kernel_diag_block_A(A, infop, matelem, dec));
    }

    auto const kernel_off_diag_block_A = [](auto &A, auto &dec, auto &infov,
                                            auto &infop, const auto &matelemB,
                                            const auto &matelemBT) {
        auto const kernel_pos = [&](auto const &pos) {
            auto ielem_p = fem::get_indices(infop, pos, dec);
            auto ielem_v = fem::get_indices_4Q1(infov, pos);

            MatSetValuesLocal(A, ielem_p.size(), ielem_p.data(), ielem_v.size(),
                              ielem_v.data(), matelemB.data(), ADD_VALUES);
            MatSetValuesLocal(A, ielem_v.size(), ielem_v.data(), ielem_p.size(),
                              ielem_p.data(), matelemBT.data(), ADD_VALUES);
        };
        return kernel_pos;
    };

    template<std::size_t Dimensions>
    PetscErrorCode B_BT_assembling(DM &dm, Mat &A,
                                   const std::array<double, Dimensions> &h)
    {
        PetscErrorCode ierr;
        DMDALocalInfo infov, infop;
        DM dav, dap;

        ierr = DMCompositeGetEntries(dm, &dav, &dap);
        CHKERRQ(ierr);

        DMDAGetLocalInfo(dav, &infov);
        DMDAGetLocalInfo(dap, &infop);

        auto box = fem::get_DM_bounds<Dimensions>(dap);

        auto matelemB = getMatElemPressureA(h);
        auto matelemBT = getMatElemPressureAT(h);

        int dec = infov.dof * infov.gxm * infov.gym * infov.gzm;

        algorithm::iterate(box, kernel_off_diag_block_A(A, dec, infov, infop,
                                                        matelemB, matelemBT));

        return ierr;
    }
} // namespace cafes
#endif