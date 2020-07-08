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

#ifndef CAFES_FEM_OPERATOR_HPP_INCLUDED
#define CAFES_FEM_OPERATOR_HPP_INCLUDED

#include <algorithm/iterate.hpp>
#include <array>
#include <cmath>
#include <fem/matElem.hpp>
#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <iostream>
#include <particle/geometry/position.hpp>
#include <petsc.h>
#include <petsc/vec.hpp>

namespace cafes
{
    namespace fem
    {
        auto const kernel_diag_block = [](auto const &x, auto &y,
                                          auto const &matelem,
                                          std::size_t order) {
            auto const kernel_pos = [&](auto const &pos) {
                auto const ielem = get_element(pos, order);
                auto const nbasis = (order + 1)*(order + 1);

                for (std::size_t k1 = 0; k1 < nbasis; ++k1)
                {
                    auto uy = y.at(ielem[k1]);
                    for (std::size_t k2 = 0; k2 < nbasis; ++k2)
                    {
                        auto ux = x.at(ielem[k2]);
                        for (std::size_t d = 0; d < x.dof_; ++d)
                        {
                            uy[d] += ux[d] * matelem[k1*nbasis + k2];
                        }
                    }
                }
            };
            return kernel_pos;
        };

        auto const kernel_tensor_block = [](auto const &x, auto &y,
                                            auto const &matelem,
                                            std::size_t order) {
            auto const kernel_pos = [&](auto const &pos) {
                auto const ielem = get_element(pos, order);
                auto const nbasis = (order + 1)*(order + 1);

                for (std::size_t k1 = 0; k1 < nbasis; ++k1)
                {
                    auto uy = y.at(ielem[k1]);
                    for (std::size_t k2 = 0; k2 < nbasis; ++k2)
                    {
                        auto ux = x.at(ielem[k2]);
                        for (std::size_t d1 = 0; d1 < x.dof_; ++d1)
                        {
                            for (std::size_t d2 = 0; d2 < x.dof_; ++d2)
                            {
                                uy[d1] += ux[d2] * matelem[(k1 + d1*nbasis)*nbasis + k2 + d2*nbasis];
                            }
                        }
                    }
                }
            };
            return kernel_pos;
        };

        auto const kernel_diag_diag_block = [](auto &x, auto const &matelem, std::size_t order) {
            auto const kernel_pos = [&](auto const &pos) {
                auto const ielem = get_element(pos, order);
                auto const nbasis = (order + 1)*(order + 1);

                for (std::size_t ie = 0; ie < nbasis; ++ie)
                {
                    auto u = x.at(ielem[ie]);
                    for (std::size_t d = 0; d < x.dof_; ++d)
                    {
                        u[d] += matelem[ie*order + ie];
                    }
                }
            };
            return kernel_pos;
        };

        auto const kernel_off_diag_block = [](auto const &x1, auto const &x2,
                                              auto &y1, auto &y2,
                                              auto const &matelem) {
            auto const kernel_pos = [&](auto const &pos) {
                auto const ielem_p = get_element(pos, 1);
                auto const ielem_v = get_element(2*pos, 2);
                auto const nbasis_p = ielem_p.size();
                auto const nbasis_v = ielem_v.size();

                for (std::size_t ie_v = 0; ie_v < nbasis_v; ++ie_v)
                {
                    auto ux = x1.at(ielem_v[ie_v]);
                    auto uy = y1.at(ielem_v[ie_v]);
                    for (std::size_t ie_p = 0; ie_p < nbasis_p; ++ie_p)
                    {
                        auto uxp = x2.at(ielem_p[ie_p]);
                        auto uyp = y2.at(ielem_p[ie_p]);
                        for (std::size_t d = 0; d < x1.dof_; ++d)
                        {
                            uyp[0] += ux[d] * matelem[ie_p*nbasis_v + ie_v + d*nbasis_p*nbasis_v];
                            uy[d] -= uxp[0] * matelem[ie_p*nbasis_v + ie_v + d*nbasis_p*nbasis_v];
                        }
                    }
                }
            };
            return kernel_pos;
        };

        #undef __FUNCT__
        #define __FUNCT__ "diag_block_mult"
        template<typename MatElem, typename Function, std::size_t Dimensions>
        PetscErrorCode diag_block_mult(petsc::petsc_vec<Dimensions> &x,
                                       petsc::petsc_vec<Dimensions> &y,
                                       MatElem const &matelem,
                                       std::size_t order,
                                       Function &&kernel)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            auto box = get_DM_bounds<Dimensions>(x.dm_);

            algorithm::iterate(box, kernel(x, y, matelem, order), order);

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "diag_diag_block_mult"
        template<typename MatElem, typename Function, std::size_t Dimensions>
        PetscErrorCode diag_diag_block_mult(petsc::petsc_vec<Dimensions> &x,
                                            MatElem const &matelem,
                                            std::size_t order,
                                            Function &&kernel)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            auto box = get_DM_bounds<Dimensions>(x.dm_);

            algorithm::iterate(box, kernel(x, matelem, order), order);

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "off_diag_block_mult"
        template<typename MatElem, typename Function, std::size_t Dimensions>
        PetscErrorCode
        off_diag_block_mult(petsc::petsc_vec<Dimensions> const &x1,
                            petsc::petsc_vec<Dimensions> const &x2,
                            petsc::petsc_vec<Dimensions> &y1,
                            petsc::petsc_vec<Dimensions> &y2,
                            MatElem const &matelem, Function &&kernel)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            auto box = get_DM_bounds<Dimensions>(x2.dm_);

            algorithm::iterate(box, kernel(x1, x2, y1, y2, matelem));

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "laplacian_mult"
        template<std::size_t Dimensions>
        PetscErrorCode laplacian_mult(petsc::petsc_vec<Dimensions> &x,
                                      petsc::petsc_vec<Dimensions> &y,
                                      std::array<double, Dimensions> const &h,
                                      std::size_t order)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            ierr = diag_block_mult(x, y, getMatElemLaplacian(h, order), order,  
                                   kernel_diag_block);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "diag_laplacian_mult"
        template<std::size_t Dimensions>
        PetscErrorCode
        diag_laplacian_mult(petsc::petsc_vec<Dimensions> &x,
                            std::array<double, Dimensions> const &h,
                            std::size_t order)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            ierr = diag_diag_block_mult(x, getMatElemLaplacian(h, order), order,
                                        kernel_diag_diag_block);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "mass_mult"
        template<std::size_t Dimensions>
        PetscErrorCode mass_mult(petsc::petsc_vec<Dimensions> &x,
                                 petsc::petsc_vec<Dimensions> &y,
                                 std::array<double, Dimensions> const &h,
                                 std::size_t order)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            ierr = diag_block_mult(x, y, getMatElemMass(h, order), order, kernel_diag_block);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "diag_mass_mult"
        template<std::size_t Dimensions>
        PetscErrorCode diag_mass_mult(petsc::petsc_vec<Dimensions> &x,
                                      std::array<double, Dimensions> const &h,
                                      std::size_t order)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            ierr = diag_diag_block_mult(x, getMatElemMass(h, order), order,
                                        kernel_diag_diag_block);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "strain_tensor_mult"
        template<std::size_t Dimensions>
        PetscErrorCode
        strain_tensor_mult(petsc::petsc_vec<Dimensions> &x,
                           petsc::petsc_vec<Dimensions> &y,
                           std::array<double, Dimensions> const &h,
                           std::size_t order)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            ierr = diag_block_mult(x, y, getMatElemStrainTensor(h, order), order,
                                   kernel_tensor_block);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

        #undef __FUNCT__
        #define __FUNCT__ "B_and_BT_mult"
        template<std::size_t Dimensions>
        PetscErrorCode B_and_BT_mult(petsc::petsc_vec<Dimensions> const &x1,
                                     petsc::petsc_vec<Dimensions> const &x2,
                                     petsc::petsc_vec<Dimensions> &y1,
                                     petsc::petsc_vec<Dimensions> &y2,
                                     std::array<double, Dimensions> const &h,
                                     std::size_t order)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            ierr = off_diag_block_mult(x1, x2, y1, y2, getMatElemPressure(h, order),
                                       kernel_off_diag_block);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }
    } // namespace fem
} // namespace cafes
#endif