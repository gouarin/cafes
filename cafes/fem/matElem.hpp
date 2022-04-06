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

#ifndef CAFES_FEM_MATELEM_HPP_INCLUDED
#define CAFES_FEM_MATELEM_HPP_INCLUDED

#include <array>
#include <petsc.h>

/* Elementary matrices for the 2D problem */

/* \int dx(\phi_i) dx(\phi_j) */
std::vector<double> matElem2d_dxudxu(std::size_t order)
{
    std::vector<double> out;
    if (order == 1)
    {
        out = {1. / 3, -1. / 3, 1. / 6, -1. / 6,
              -1. / 3, 1. / 3, -1. / 6, 1. / 6,
               1. / 6, -1. / 6, 1. / 3, -1. / 3,
              -1. / 6, 1. / 6, -1. / 3, 1. / 3};
    }
    else if (order == 2)
    {
        out = {14./45, -16./45, 2./45, 7./45, -8./45, 1./45, -7./90, 4./45, -1./90,
                -16./45, 32./45, -16./45, -8./45, 16./45, -8./45, 4./45, -8./45, 4./45,
                2./45, -16./45, 14./45, 1./45, -8./45, 7./45, -1./90, 4./45, -7./90,
                7./45, -8./45, 1./45, 56./45, -64./45, 8./45, 7./45, -8./45, 1./45,
                -8./45, 16./45, -8./45, -64./45, 128./45, -64./45, -8./45, 16./45, -8./45,
                1./45, -8./45, 7./45, 8./45, -64./45, 56./45, 1./45, -8./45, 7./45,
                -7./90, 4./45, -1./90, 7./45, -8./45, 1./45, 14./45, -16./45, 2./45,
                4./45, -8./45, 4./45, -8./45, 16./45, -8./45, -16./45, 32./45, -16./45,
                -1./90, 4./45, -7./90, 1./45, -8./45, 7./45, 2./45, -16./45, 14./45};
    }
    return out;
}

std::vector<double> matElem2d_dyudyu(std::size_t order)
{
    std::vector<double> out;

    if (order == 1)
    {
        out = {1. / 3, 1. / 6, -1. / 3, -1. / 6,
                1. / 6, 1. / 3, -1. / 6, -1. / 3,
                -1. / 3, -1. / 6, 1. / 3, 1. / 6,
                -1. / 6, -1. / 3, 1. / 6, 1. / 3};
    }
    else if (order == 2)
    {
        out = {14./45, 7./45, -7./90, -16./45, -8./45, 4./45, 2./45, 1./45, -1./90,
                7./45, 56./45, 7./45, -8./45, -64./45, -8./45, 1./45, 8./45, 1./45,
                -7./90, 7./45, 14./45, 4./45, -8./45, -16./45, -1./90, 1./45, 2./45,
                -16./45, -8./45, 4./45, 32./45, 16./45, -8./45, -16./45, -8./45, 4./45,
                -8./45, -64./45, -8./45, 16./45, 128./45, 16./45, -8./45, -64./45, -8./45,
                4./45, -8./45, -16./45, -8./45, 16./45, 32./45, 4./45, -8./45, -16./45,
                2./45, 1./45, -1./90, -16./45, -8./45, 4./45, 14./45, 7./45, -7./90,
                1./45, 8./45, 1./45, -8./45, -64./45, -8./45, 7./45, 56./45, 7./45,
                -1./90, 1./45, 2./45, 4./45, -8./45, -16./45, -7./90, 7./45, 14./45};
    }
    return out;
}

std::vector<double> matElem2d_dxudyu(std::size_t order)
{
     std::vector<double> out;

    if (order == 1)
    {
        out = {0.25, -0.25, 0.25, -0.25,
               0.25, -0.25, 0.25, -0.25,
               -0.25, 0.25, -0.25, 0.25,
               -0.25, 0.25, -0.25, 0.25};
    }
    else if (order == 2)
    {
        out = {1./4, 1./3, -1./12, -1./3, -4./9, 1./9, 1./12, 1./9, -1./36,
                -1./3, 0, 1./3, 4./9, 0, -4./9, -1./9, 0, 1./9,
                1./12, -1./3, -1./4, -1./9, 4./9, 1./3, 1./36, -1./9, -1./12,
                1./3, 4./9, -1./9, 0, 0, 0, -1./3, -4./9, 1./9,
                -4./9, 0, 4./9, 0, 0, 0, 4./9, 0, -4./9,
                1./9, -4./9, -1./3, 0, 0, 0, -1./9, 4./9, 1./3,
                -1./12, -1./9, 1./36, 1./3, 4./9, -1./9, -1./4, -1./3, 1./12,
                1./9, 0, -1./9, -4./9, 0, 4./9, 1./3, 0, -1./3,
                -1./36, 1./9, 1./12, 1./9, -4./9, -1./3, -1./12, 1./3, 1./4};
    }
    return out;
}

std::vector<double> matElemMass2d(std::size_t order)
{
    std::vector<double> out;

    if (order == 1)
    {
        out = {1. / 9, 1. / 18, 1. / 18, 1. / 36,
                1. / 18, 1. / 9, 1. / 36, 1. / 18,
                1. / 18, 1. / 36, 1. / 9, 1. / 18,
                1. / 36, 1. / 18, 1. / 18, 1. / 9};
    }
    else if (order == 2)
    {
        out = {16./225, 8./225, -4./225, 8./225, 4./225, -2./225, -4./225, -2./225, 1./225,
                8./225, 64./225, 8./225, 4./225, 32./225, 4./225, -2./225, -16./225, -2./225,
                -4./225, 8./225, 16./225, -2./225, 4./225, 8./225, 1./225, -2./225, -4./225,
                8./225, 4./225, -2./225, 64./225, 32./225, -16./225, 8./225, 4./225, -2./225,
                4./225, 32./225, 4./225, 32./225, 256./225, 32./225, 4./225, 32./225, 4./225,
                -2./225, 4./225, 8./225, -16./225, 32./225, 64./225, -2./225, 4./225, 8./225,
                -4./225, -2./225, 1./225, 8./225, 4./225, -2./225, 16./225, 8./225, -4./225,
                -2./225, -16./225, -2./225, 4./225, 32./225, 4./225, 8./225, 64./225, 8./225,
                1./225, -2./225, -4./225, -2./225, 4./225, 8./225, -4./225, 8./225, 16./225};
    }
    return out;
}

std::vector<double> matElemMass_pu_2d(std::size_t order)
{
    std::vector<double> out;

    if (order == 2)
    {
        out = {1./9, 0, 0, 0,
                2./9, 2./9, 0, 0,
                0, 1./9, 0, 0,
                2./9, 0, 2./9, 0,
                4./9, 4./9, 4./9, 4./9,
                0, 2./9, 0, 2./9,
                0, 0, 1./9, 0,
                0, 0, 2./9, 2./9,
                0, 0, 0, 1./9};
    }
    return out;
}

std::vector<double> matElem2d_pdxu(std::size_t order)
{
    std::vector<double> out;

    if (order == 1)
    {
        out = {15. / 48, -10. / 48, -5. / 48, 18. / 48, -12. / 48, -6. / 48, 3. / 48, -2. / 48, -1. / 48,
                5. / 48, 10. / 48, -15. / 48, 6. / 48, 12. / 48, -18. / 48, 1. / 48, 2. / 48, -3. / 48,
                3. / 48, -2. / 48, -1. / 48, 18. / 48, -12. / 48, -6. / 48, 15. / 48, -10. / 48, -5. / 48,
                1. / 48, 2. / 48, -3. / 48, 6. / 48, 12. / 48, -18. / 48, 5. / 48, 10. / 48, -15. / 48};
    }
    else if (order == 2)
    {
        out = { -5./18, 2./9, 1./18, -5./9, 4./9, 1./9, 0, 0, 0,
                 -1./18, -2./9, 5./18, -1./9, -4./9, 5./9, 0, 0, 0,
                 0, 0, 0, -5./9, 4./9, 1./9, -5./18, 2./9, 1./18,
                 0, 0, 0, -1./9, -4./9, 5./9, -1./18, -2./9, 5./18 };
    }
    return out;
}

std::vector<double> matElem2d_pdyu(std::size_t order)
{
    std::vector<double> out;

    if (order == 1)
    {
        out = {15. / 48, 18. / 48, 3. / 48, -10. / 48, -12. / 48, -2. / 48, -5. / 48, -6. / 48, -1. / 48,
                3. / 48, 18. / 48, 15. / 48, -2. / 48, -12. / 48, -10. / 48, -1. / 48, -6. / 48, -5. / 48,
                5. / 48, 6. / 48, 1. / 48, 10. / 48, 12. / 48, 2. / 48, -15. / 48, -18. / 48, -3. / 48,
                1. / 48, 6. / 48, 5. / 48, 2. / 48, 12. / 48, 10. / 48, -3. / 48, -18. / 48, -15. / 48};
    }
    else if (order == 2)
    {
        out = {-5./18, -5./9, 0, 2./9, 4./9, 0, 1./18, 1./9, 0,
                0, -5./9, -5./18, 0, 4./9, 2./9, 0, 1./9, 1./18,
                -1./18, -1./9, 0, -2./9, -4./9, 0, 5./18, 5./9, 0,
                0, -1./9, -1./18, 0, -4./9, -2./9, 0, 5./9, 5./18};
    }
    return out;
}

auto getMatElemLaplacian(std::array<double, 2> const &h, std::size_t order)
{
    std::size_t size = (order + 1)*(order + 1);
    std::vector<double> MatElem(size*size);
    double hxy = h[0] / h[1], hyx = h[1] / h[0];

    auto dxudxu = matElem2d_dxudxu(order);
    auto dyudyu = matElem2d_dyudyu(order);
    for (std::size_t i = 0; i < size; ++i)
    {
        for (std::size_t j = 0; j < size; ++j)
        {
            MatElem[i*size + j] = hyx * dxudxu[i*size + j] + hxy * dyudyu[i*size + j];
        }
    }
    return MatElem;
}

auto getMatElemStrainTensor(std::array<double, 2> const &h, std::size_t order)
{
    std::size_t size = (order + 1)*(order + 1);
    std::vector<double> MatElem(4*size*size);
    double hxy = h[0] / h[1], hyx = h[1] / h[0];

    auto dxudxu = matElem2d_dxudxu(order);
    auto dyudyu = matElem2d_dyudyu(order);
    auto dxudyu = matElem2d_dxudyu(order);

    for (std::size_t i = 0; i < size; ++i)
    {
        for (std::size_t j = 0; j < size; ++j)
        {
            MatElem[i*2*size + j] = 2 * hyx * dxudxu[i*size + j] + hxy * dyudyu[i*size + j];
            MatElem[i*2*size + j + size] = dxudyu[i*size + j];
            MatElem[(i + size)*2*size + j] = dxudyu[j*size + i];
            MatElem[(i + size)*2*size + j + size] = hyx * dxudxu[i*size + j] + 2 * hxy * dyudyu[i*size + j];
        }
    }

    return MatElem;
}

auto getMatElemPressure(std::array<double, 2> const &h, std::size_t order)
{
    std::vector<double> MatElem(2*9*4);

    auto pdxu = matElem2d_pdxu(order);
    auto pdyu = matElem2d_pdyu(order);

    for (std::size_t i = 0; i < 4; ++i)
    {
        for (std::size_t j = 0; j < 9; ++j)
        {
            MatElem[i*2*9 + j] = h[1] * pdxu[i*9 + j];
            MatElem[i*2*9 + j + 9] = h[0] * pdyu[i*9 + j];
        }
    }

    return MatElem;
}

auto getMatElemPressureT(std::array<double, 2> const &h, std::size_t order)
{
    std::vector<double> MatElem(2*9*4);

    auto pdxu = matElem2d_pdxu(order);
    auto pdyu = matElem2d_pdyu(order);

    for (std::size_t i = 0; i < 9; ++i)
    {
        for (std::size_t j = 0; j < 4; ++j)
        {
            MatElem[i*4 + j] = -h[1] * pdxu[j*9 + i];
            MatElem[(i + 9)*4 + j] = -h[0] * pdyu[j*9 + i];
        }
    }

    return MatElem;
}

auto getMatElemMass(std::array<double, 2> const &h, std::size_t order)
{
    std::size_t size = (order + 1)*(order + 1);
    std::vector<double> MatElem(size*size);

    auto mass = matElemMass2d(order);

    for (std::size_t i = 0; i < size; ++i)
    {
        for (std::size_t j = 0; j < size; ++j)
        {
            MatElem[i*size + j] = h[0] * h[1] * mass[i*size + j];
        }
    }

    return MatElem;
}

auto getMatElemMassP2U(std::array<double, 2> const &h, std::size_t order)
{
    std::vector<double> MatElem(4*9);

    auto mass = matElemMass_pu_2d(order);

    for (std::size_t i = 0; i < 9; ++i)
    {
        for (std::size_t j = 0; j < 4; ++j)
        {
            MatElem[i*4 + j] = h[0] * h[1] * mass[i*4 + j];
        }
    }

    return MatElem;
}
#endif