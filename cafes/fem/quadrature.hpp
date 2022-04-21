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

#ifndef PARTICLE_FEM_POSITION_HPP_INCLUDED
#define PARTICLE_FEM_POSITION_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <tuple>

namespace cafes
{
    namespace fem
    {
        auto P1_integration(geometry::position<double, 2> const &x,
                            std::array<double, 2> const &h)
        {
            std::array<double, 4> that;
            double c = 1. / (h[0] * h[1]);
            that[0] = c * (x[0] - h[0]) * (x[1] - h[1]);
            that[1] = c * x[0] * (h[1] - x[1]);
            that[2] = c * (h[0] - x[0]) * x[1];
            that[3] = c * x[0] * x[1];
            return that;
        }

        auto P2_integration(geometry::position<double, 2> const &x,
                            std::array<double, 2> const &h)
        {
            std::array<double, 9> that;
            double c = 1. / (h[0]*h[0]*h[1]*h[1]);
            that[0] = .25*c*(-2*h[0] + x[0])*(-h[0] + x[0])*(-2*h[1] + x[1])*(-h[1] + x[1]);
            that[1] = -.5*c*x[0]*(-2*h[0] + x[0])*(-2*h[1] + x[1])*(-h[1] + x[1]);
            that[2] = .25*c*x[0]*(-h[0] + x[0])*(-2*h[1] + x[1])*(-h[1] + x[1]);
            that[3] = -.5*c*x[1]*(-2*h[0] + x[0])*(-h[0] + x[0])*(-2*h[1] + x[1]);
            that[4] =  c*x[0]*x[1]*(-2*h[0] + x[0])*(-2*h[1] + x[1]);
            that[5] = -.5*c*x[0]*x[1]*(-h[0] + x[0])*(-2*h[1] + x[1]);
            that[6] = .25*c*x[1]*(-2*h[0] + x[0])*(-h[0] + x[0])*(-h[1] + x[1]);
            that[7] = -.5*c*x[0]*x[1]*(-2*h[0] + x[0])*(-h[1] + x[1]);
            that[8] = .25*c*x[0]*x[1]*(-h[0] + x[0])*(-h[1] + x[1]);
            return that;
        }

         auto P1_integration_sing(geometry::position<double, 2> const &x,
                            std::array<double, 2> const &h)
        {
            std::array<double, 4> that;
            that[0] = (x[0] - h[0]) * (x[1] - h[1]);
            that[1] = x[0] * (h[1] - x[1]);
            that[2] = (h[0] - x[0]) * x[1];
            that[3] = x[0] * x[1];
            return that;
        }

        auto P1_integration_grad(geometry::position<double, 2> const &x,
                                 std::array<double, 2> const &h)
        {
            std::array<std::array<double, 2>, 4> that;
            // double c = 1./(h[0]*h[1]);
            that[0] = {{x[1] - h[1], x[0] - h[0]}};
            that[1] = {{h[1] - x[1], -x[0]}};
            that[2] = {{-x[1], h[0] - x[0]}};
            that[3] = {{x[1], x[0]}};
            return that;
        }

        auto P1_integration(geometry::position<double, 3> const &x,
                            std::array<double, 3> const &h)
        {
            std::array<double, 8> that;
            double c = 1. / (h[0] * h[1] * h[2]);

            that[0] = c * (h[0] - x[0]) * (h[1] - x[1]) * (h[2] - x[2]);
            that[1] = c * x[0] * (h[1] - x[1]) * (h[2] - x[2]);
            that[2] = c * (h[0] - x[0]) * x[1] * (h[2] - x[2]);
            that[3] = c * x[0] * x[1] * (h[2] - x[2]);
            that[4] = c * (h[0] - x[0]) * (h[1] - x[1]) * x[2];
            that[5] = c * x[0] * (h[1] - x[1]) * x[2];
            that[6] = c * (h[0] - x[0]) * x[1] * x[2];
            that[7] = c * x[0] * x[1] * x[2];
            return that;
        }

        auto P1_integration_grad(geometry::position<double, 3> const &x,
                                 std::array<double, 3> const &h)
        {
            std::array<std::array<double, 3>, 8> that;

            that[0] = {{-(h[1] - x[1]) * (h[2] - x[2]),
                        -(h[0] - x[0]) * (h[2] - x[2]),
                        -(h[0] - x[0]) * (h[1] - x[1])}};
            that[1] = {{(h[1] - x[1]) * (h[2] - x[2]), -x[0] * (h[2] - x[2]),
                        -x[0] * (h[1] - x[1])}};
            that[2] = {{-x[1] * (h[2] - x[2]), (h[0] - x[0]) * (h[2] - x[2]),
                        -(h[0] - x[0]) * x[1]}};
            that[3] = {
                {x[1] * (h[2] - x[2]), x[0] * (h[2] - x[2]), -x[0] * x[1]}};
            that[4] = {{-(h[1] - x[1]) * x[2], -(h[0] - x[0]) * x[2],
                        (h[0] - x[0]) * (h[1] - x[1])}};
            that[5] = {
                {(h[1] - x[1]) * x[2], -x[0] * x[2], x[0] * (h[1] - x[1])}};
            that[6] = {
                {-x[1] * x[2], (h[0] - x[0]) * x[2], (h[0] - x[0]) * x[1]}};
            that[7] = {{x[1] * x[2], x[0] * x[2], x[0] * x[1]}};
            return that;
        }

        auto get_element(geometry::position<PetscInt, 2> const &ix, std::size_t order, std::size_t step = 1)
        {
            std::size_t size = order + 1;
            std::vector<std::array<PetscInt, 2>> that(size*size);

            std::size_t ind = 0;

            for(int j=0; j<size; ++j)
            {
                for(int i=0; i<size; ++i)
                {
                    that[ind][0] = ix[0] + i*step;
                    that[ind][1] = ix[1] + j*step;
                    ind++;
                }
            }
            return that;
        }

        inline PetscInt DMDALocalIndex2D(DMDALocalInfo &info, int i, int j)
        {
            return info.dof * ((j - info.gys) * info.gxm + i - info.gxs);
        }

        auto get_indices(DMDALocalInfo &info,
                         geometry::position<PetscInt, 2> const &ix,
                         int dof,
                         std::size_t order)
        {
            std::size_t size = order + 1;
            std::vector<PetscInt> that(size*size);

            int ind = 0;
            for(int j=0; j<size; ++j)
            {
                for(int i=0; i<size; ++i)
                {
                    that[ind++] = dof + DMDALocalIndex2D(info, ix[0] + i, ix[1] + j);
                }
            }
            return that;
        }

        auto get_range(const DMBoundaryType& bt, int i, int min, int max, int stencil_size)
        {
            return std::make_pair<PetscInt, int>(
                (bt == DM_BOUNDARY_PERIODIC)? -stencil_size : PetscMax(-stencil_size, -i + min),
                (bt == DM_BOUNDARY_PERIODIC)? stencil_size: PetscMin(stencil_size, max - i - 1)
            );
        }

        auto get_col_indices(DMDALocalInfo &info,
                             geometry::position<PetscInt, 2> const &ix,
                             int dof,
                             std::size_t order)
        {
            int istart, iend;
            int jstart, jend;

            std::tie(istart, iend) = get_range(info.bx, ix[0], 0, info.mx, order);
            std::tie(jstart, jend) = get_range(info.by, ix[1], 0, info.my, order);

            // std::cout << ix[0] << " " << ix[1] << " " << istart << " " << iend << " " << jstart << " " << jend << std::endl;
            std::vector<PetscInt> that((iend - istart + 1)*(jend - jstart + 1));

            int ind = 0;
            for(int j = jstart; j < jend + 1; ++j)
            {
                for(int i = istart; i < iend + 1; ++i)
                {
                    that[ind++] = dof + DMDALocalIndex2D(info, ix[0] + i, ix[1] + j);
                }
            }
            return that;
        }

        auto get_row_indices(DMDALocalInfo &info,
                             geometry::position<PetscInt, 2> const &ix,
                             int dof,
                             std::size_t order)
        {
            return dof + DMDALocalIndex2D(info, ix[0], ix[1]);
        }

        auto get_row_indices_tensor(DMDALocalInfo &info,
                                    geometry::position<PetscInt, 2> const &ix,
                                    std::size_t order)
        {
            int dof = 2;
            std::vector<PetscInt> that(dof);

            std::size_t ind = 0;
            for(std::size_t d = 0; d < dof; ++d)
            {
                that[ind++] = d + DMDALocalIndex2D(info, ix[0], ix[1]);
            }
            return that;
        }

        auto get_indices_tensor(DMDALocalInfo &info,
                                geometry::position<PetscInt, 2> const &ix,
                                std::size_t order)
        {
            std::size_t size = order + 1;
            std::size_t dof = 2;
            std::vector<PetscInt> that(dof*size*size);

            std::size_t ind = 0;
            for(std::size_t d = 0; d < dof; ++d)
            {
                for(int j=0; j<size; ++j)
                {
                    for(int i=0; i<size; ++i)
                    {
                        that[ind++] = d + DMDALocalIndex2D(info, ix[0] + i, ix[1] + j);
                    }
                }
            }
            return that;
        }

        auto get_col_indices_tensor(DMDALocalInfo &info,
                             geometry::position<PetscInt, 2> const &ix,
                             std::size_t order)
        {
            std::size_t dof = 2;

            int istart, iend;
            int jstart, jend;

            std::tie(istart, iend) = get_range(info.bx, ix[0], 0, info.mx, order);
            std::tie(jstart, jend) = get_range(info.by, ix[1], 0, info.my, order);

            std::vector<PetscInt> that(dof*(iend - istart + 1)*(jend - jstart + 1));

            int ind = 0;

            for(int j = jstart; j < jend + 1; ++j)
            {
                for(int i = istart; i < iend + 1; ++i)
                {
                    auto index = DMDALocalIndex2D(info, ix[0] + i, ix[1] + j);
                    for(std::size_t d=0; d<dof; ++d)
                    {
                        that[ind++] = d + index;
                    }
                }
            }
            return that;
        }

        auto get_element(geometry::position<PetscInt, 3> const &ix)
        {
            std::array<std::array<PetscInt, 3>, 8> that;

            that[0] = {ix[0], ix[1], ix[2]};
            that[1] = {ix[0] + 1, ix[1], ix[2]};
            that[2] = {ix[0], ix[1] + 1, ix[2]};
            that[3] = {ix[0] + 1, ix[1] + 1, ix[2]};
            that[4] = {ix[0], ix[1], ix[2] + 1};
            that[5] = {ix[0] + 1, ix[1], ix[2] + 1};
            that[6] = {ix[0], ix[1] + 1, ix[2] + 1};
            that[7] = {ix[0] + 1, ix[1] + 1, ix[2] + 1};
            return that;
        }

        auto get_element_4Q1(geometry::position<PetscInt, 2> const &ix)
        {
            std::array<std::array<PetscInt, 2>, 9> that;

            that[0] = {2 * ix[0], 2 * ix[1]};
            that[1] = {2 * ix[0] + 1, 2 * ix[1]};
            that[2] = {2 * ix[0] + 2, 2 * ix[1]};
            that[3] = {2 * ix[0], 2 * ix[1] + 1};
            that[4] = {2 * ix[0] + 1, 2 * ix[1] + 1};
            that[5] = {2 * ix[0] + 2, 2 * ix[1] + 1};
            that[6] = {2 * ix[0], 2 * ix[1] + 2};
            that[7] = {2 * ix[0] + 1, 2 * ix[1] + 2};
            that[8] = {2 * ix[0] + 2, 2 * ix[1] + 2};
            return that;
        }

        auto get_indices_4Q1(DMDALocalInfo &info,
                             geometry::position<PetscInt, 2> const &ix)
        {
            std::array<PetscInt, 18> that;

            for (std::size_t d = 0; d < 2; ++d)
            {
                that[0 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0]    ,     2 * ix[1]);
                that[1 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0] + 1,     2 * ix[1]);
                that[2 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0] + 2,     2 * ix[1]);
                that[3 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0]    , 2 * ix[1] + 1);
                that[4 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0] + 1, 2 * ix[1] + 1);
                that[5 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0] + 2, 2 * ix[1] + 1);
                that[6 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0]    , 2 * ix[1] + 2);
                that[7 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0] + 1, 2 * ix[1] + 2);
                that[8 + d * 9] = d + DMDALocalIndex2D(info, 2 * ix[0] + 2, 2 * ix[1] + 2);
            }
            return that;
        }

        auto get_element_4Q1(geometry::position<PetscInt, 3> const &ix)
        {
            std::array<std::array<PetscInt, 3>, 27> that;

            that[0] = {2 * ix[0], 2 * ix[1], 2 * ix[2]};
            that[1] = {2 * ix[0] + 1, 2 * ix[1], 2 * ix[2]};
            that[2] = {2 * ix[0] + 2, 2 * ix[1], 2 * ix[2]};
            that[3] = {2 * ix[0], 2 * ix[1] + 1, 2 * ix[2]};
            that[4] = {2 * ix[0] + 1, 2 * ix[1] + 1, 2 * ix[2]};
            that[5] = {2 * ix[0] + 2, 2 * ix[1] + 1, 2 * ix[2]};
            that[6] = {2 * ix[0], 2 * ix[1] + 2, 2 * ix[2]};
            that[7] = {2 * ix[0] + 1, 2 * ix[1] + 2, 2 * ix[2]};
            that[8] = {2 * ix[0] + 2, 2 * ix[1] + 2, 2 * ix[2]};
            that[9] = {2 * ix[0], 2 * ix[1], 2 * ix[2] + 1};
            that[10] = {2 * ix[0] + 1, 2 * ix[1], 2 * ix[2] + 1};
            that[11] = {2 * ix[0] + 2, 2 * ix[1], 2 * ix[2] + 1};
            that[12] = {2 * ix[0], 2 * ix[1] + 1, 2 * ix[2] + 1};
            that[13] = {2 * ix[0] + 1, 2 * ix[1] + 1, 2 * ix[2] + 1};
            that[14] = {2 * ix[0] + 2, 2 * ix[1] + 1, 2 * ix[2] + 1};
            that[15] = {2 * ix[0], 2 * ix[1] + 2, 2 * ix[2] + 1};
            that[16] = {2 * ix[0] + 1, 2 * ix[1] + 2, 2 * ix[2] + 1};
            that[17] = {2 * ix[0] + 2, 2 * ix[1] + 2, 2 * ix[2] + 1};
            that[18] = {2 * ix[0], 2 * ix[1], 2 * ix[2] + 2};
            that[19] = {2 * ix[0] + 1, 2 * ix[1], 2 * ix[2] + 2};
            that[20] = {2 * ix[0] + 2, 2 * ix[1], 2 * ix[2] + 2};
            that[21] = {2 * ix[0], 2 * ix[1] + 1, 2 * ix[2] + 2};
            that[22] = {2 * ix[0] + 1, 2 * ix[1] + 1, 2 * ix[2] + 2};
            that[23] = {2 * ix[0] + 2, 2 * ix[1] + 1, 2 * ix[2] + 2};
            that[24] = {2 * ix[0], 2 * ix[1] + 2, 2 * ix[2] + 2};
            that[25] = {2 * ix[0] + 1, 2 * ix[1] + 2, 2 * ix[2] + 2};
            that[26] = {2 * ix[0] + 2, 2 * ix[1] + 2, 2 * ix[2] + 2};

            return that;
        }
    } // namespace fem
} // namespace cafes

#endif