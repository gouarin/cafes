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

#ifndef PARTICLE_GEOMETRY_VECTOR_HPP_INCLUDED
#define PARTICLE_GEOMETRY_VECTOR_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/physics/velocity.hpp>
#include <array>
#include <iostream>
#include <type_traits>

namespace cafes
{
  namespace geometry
  {
    template<typename T, std::size_t Dimensions>
    struct vector : private std::array<T, Dimensions>
    {
      using parent = std::array<T, Dimensions>;
      using parent::operator[];
      using parent::fill;
      using parent::begin;
      using parent::end;

      vector() = default;
      vector(vector const&) = default;

      template<template<typename, std::size_t> class T1, typename U>
      vector(T1<U, Dimensions> const& p)
      {
        for(std::size_t i=0; i<Dimensions; ++i)
          (*this)[i] = static_cast<T>(p[i]);
      } 

      vector(physics::velocity<Dimensions> const& v)
      {
        for(std::size_t i=0; i<Dimensions; ++i)
          (*this)[i] = static_cast<T>(v[i]);
      } 

      template<typename... Ts>
      vector( Ts... ts ) : parent{{ static_cast<T>(ts)...}}
      {
         static_assert( sizeof...(Ts) <= Dimensions
                      , "Size mismatch in vector variadic constructor"
                      );
      }

      friend std::ostream &operator<<( std::ostream &output, 
                                       vector const& vec )
      { 
        output << "vector( ";

        for(std::size_t i=0; i<Dimensions; ++i)
          output << vec[i] << " ";

        output << ")";

        return output;            
      }

      // vector& operator*=(quaternion const& q)
      // {
      //     //todo: MATH!
      //     return *this;
      // }

      vector& operator*=(vector<T, Dimensions> const& v)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= v[i];
          return *this;
      }

      vector& operator*=(T const& c)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= c;
          return *this;
      }

      vector& operator+=(vector<T, Dimensions> const& v)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
          return *this;
      }

      vector& operator-=(vector<T, Dimensions> const& v)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= v[i];
          return *this;
      }

      vector& operator+=(physics::velocity<Dimensions> const& v)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= v[i];
          return *this;
      }

      vector& operator/=(T const& c)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] /= c;
          return *this;
      }

      // vector& operator=(vector<Dimensions, T> const& other)
      // {
      //   for(std::size_t i=0;i<Dimensions;++i) (*this)[i] = other[i];
      //   return *this;
      // }

    };
    
    template<typename T, std::size_t Dimensions>
    auto operator*(vector<T, Dimensions> const& v, T const& c)
    {
      vector<T, Dimensions> that{v};
      return that*=c;
    }

    template<typename T, std::size_t Dimensions>
    auto operator*(T const& c, vector<T, Dimensions> const& v)
    {
      return v*c;
    }

    template<typename T, std::size_t Dimensions>
    auto operator+(vector<T, Dimensions> const& v1, vector<T, Dimensions> const& v2)
    {
      vector<T, Dimensions> that{v1};
      return that+=v2;
    }

    template<typename T, std::size_t Dimensions>
    auto operator-(vector<T, Dimensions> const& v)
    {
      vector<T, Dimensions> that{v};
      return that *= -1;
    }

    template<typename T, std::size_t Dimensions>
    vector<T, Dimensions> operator-( vector<T, Dimensions> const& p
                                   , physics::velocity<Dimensions> const& v
                                  )
    {
        vector<T, Dimensions> that{p};
        return that -= v;
    }

    template<typename T, std::size_t Dimensions>
    vector<T, Dimensions> operator-( physics::velocity<Dimensions> const& v
                                   , vector<T, Dimensions> const& p
                                  )
    {
        vector<T, Dimensions> that{p};
        return that -= v;
    }

    template<typename T, std::size_t Dimensions>
    auto operator/(vector<T, Dimensions> const& v, T const& c)
    {
      vector<T, Dimensions> that{v};
      return that/=c;
    }

    template<typename T, std::size_t Dimensions>
    auto operator/(vector<T, Dimensions> const& v, int const& c)
    {
      vector<T, Dimensions> that{v};
      return that/=c;
    }

  }
}
#endif