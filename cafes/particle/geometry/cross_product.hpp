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

#ifndef PARTICLE_GEOMETRY_CROSS_PRODUCT_HPP_INCLUDED
#define PARTICLE_GEOMETRY_CROSS_PRODUCT_HPP_INCLUDED

#include <particle/physics/velocity.hpp>

namespace cafes
{
  namespace geometry
  {

    template<template<typename, std::size_t> class T, typename U>
    U cross_product(T<U, 2> const& x, T<U, 2> const& y)
    {
       return x[0]*y[1] - x[1]*y[0];
    }

    template<template<typename, std::size_t> class T1, template<typename, std::size_t> class T2, typename U>
    U cross_product(T1<U, 2> const& x, T2<U, 2> const& y)
    {
       return x[0]*y[1] - x[1]*y[0];
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(double rz, T<U, 2> const& x){
      return { -rz*x[1],
                rz*x[0]
             };
    }


    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(T<U, 1> const& rz, T<U, 2> const& x){
      return { -rz[0]*x[1],
                rz[0]*x[0]
             };
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(T<U, 2> const& x, double rz){
      return {  rz*x[1],
               -rz*x[0]
             };
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(T<U, 2> const& x, T<U, 1> const & rz){
      return {  rz[0]*x[1],
               -rz[0]*x[0]
             };
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 3> cross_product(T<U, 3> const& x, T<U, 3> const& y){
      return { x[1]*y[2] - x[2]*y[1],
               x[2]*y[0] - x[0]*y[2],
               x[0]*y[1] - x[1]*y[0]
             };
    }

  }
}
#endif