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

#ifndef PARTICLE_GEOMETRY_POSITION_HPP_INCLUDED
#define PARTICLE_GEOMETRY_POSITION_HPP_INCLUDED

#include <particle/physics/velocity.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>

namespace cafes
{
  namespace geometry
  {
      template<typename T, std::size_t Dimensions>
      struct position : private std::array<T, Dimensions>
      {
        using parent = std::array<T, Dimensions>;
        using parent::operator[];
        using parent::begin;
        using parent::end;
        using parent::fill;
        using parent::data;

        static constexpr std::size_t dimensions = Dimensions;

        position() = default;
        position(position const&) = default;

        template<typename U> position(position<U, Dimensions> const& p)
        {
          for(std::size_t i=0; i<Dimensions; ++i)
            (*this)[i] = static_cast<T>(p[i]);
        } 

        template<typename... Ts>
        position( Ts... ts ) : parent{{ static_cast<T>(ts)...}}
        {
           static_assert( sizeof...(Ts) <= Dimensions
                        , "Size mismatch in position variadic constructor"
                        );
        }

        friend std::ostream &operator<<( std::ostream &output, 
                                         position const& pos )
        { 
          output << "position( ";

          for(std::size_t i=0; i<Dimensions; ++i)
            output << pos[i] << " ";

          output << ")";

          return output;            
        }

        // position& operator*=(quaternion const& q)
        // {
        //     //todo: MATH!
        //     return *this;
        // }

        position& operator*=(position<T, Dimensions> const& p)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= p[i];
            return *this;
        }

        position& operator*=(std::array<T, Dimensions> const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= v[i];
            return *this;
        }

        position& operator-=(position<T, Dimensions> const& p)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= p[i];
            return *this;
        }

        position& operator+=(position<T, Dimensions> const& p)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += p[i];
            return *this;
        }

        position& operator+=(physics::velocity<Dimensions> const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
            return *this;
        }

        position& operator+=(T const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v;
            return *this;
        }

        position& operator-=(T const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= v;
            return *this;
        }

        position& operator/=(std::array<T, Dimensions> const& v)
        {
            static_assert(std::is_floating_point<T>::value, "invalid operation in division");
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] /= v[i];
            return *this;
        }

        position& operator/=(T const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] /= v;
            return *this;
        }

      };

      template<typename T, std::size_t Dimensions>
      position<T, Dimensions> operator+( position<T, Dimensions> const& p
                                       , physics::velocity<Dimensions> const& v
                                    )
      {
          position<T, Dimensions> that{p};
          return that += v;
      }

      template<typename T, std::size_t Dimensions>
      position<T, Dimensions> operator+( position<T, Dimensions> const& p
                                       , T const& v
                                    )
      {
          position<T, Dimensions> that{p};
          return that += v;
      }

      template<typename T, std::size_t Dimensions>
      position<T, Dimensions> operator+( position<T, Dimensions> const& p1
                                       , position<T, Dimensions> const& p2
                                       )
      {
          position<T, Dimensions> that{p1};
          return that += p2;
      }

      template<typename T, std::size_t Dimensions>
      position<T, Dimensions> operator-( position<T, Dimensions> const& p
                                    , T const& v
                                    )
      {
          position<T, Dimensions> that{p};
          return that -= v;
      }

      // template<std::size_t Dimensions, typename T>
      // position<Dimensions, T> operator*( position<Dimensions, T> const& p
      //                               , quaternion const& q
      //                               )
      // {
      //     position<Dimensions, T> that{p};
      //     return that *= q;
      // }

      template<typename T, typename U, std::size_t Dimensions>
      position<double, Dimensions> operator*( position<T, Dimensions> const& p1
                                         ,    position<U, Dimensions> const& p2
                                    )
      {
          position<double, Dimensions> that{p1};
          position<double, Dimensions> tmp{p2};
          return that *= p2;
      }

      template<typename T, typename U, std::size_t Dimensions>
      position<double, Dimensions> operator*( position<T, Dimensions> const& p
                                            , std::array<U, Dimensions> const& v
                                            )
      {
          position<double, Dimensions> that{p};
          return that *= v;
      }

      template<typename T, std::size_t Dimensions>
      position<double, Dimensions> operator-( position<T, Dimensions> const& p1
                                            , position<T, Dimensions> const& p2
                                            )
      {
          position<T, Dimensions> that{p1};
          return that -= p2;
      }

      template<typename T, std::size_t Dimensions>
      position<T, Dimensions> operator*( position<T, Dimensions> const& p1
                                    ,    position<T, Dimensions> const& p2
                                    )
      {
          position<T, Dimensions> that{p1};
          return that *= p2;
      }

      template<typename T, std::size_t Dimensions>
      position<T, Dimensions> operator/( position<T, Dimensions> const& p
                                       , std::array<T, Dimensions> const& v
                                       )
      {
          position<T, Dimensions> that{p};
          return that /= v;
      }

      template<typename T, std::size_t Dimensions>
      position<T, Dimensions> operator/( position<T, Dimensions> const& p
                                       , T const& v
                                       )
      {
          position<T, Dimensions> that{p};
          return that /= v;
      }

      template<std::size_t Dimensions>
      position<int, Dimensions> ceil(position<double, Dimensions> const& pos)
      {
          position<double, Dimensions> that;
          std::transform(pos.begin(), pos.end(), that.begin(), [](auto& p){return std::ceil(p);});
          return that;
      }

      template<std::size_t Dimensions>
      position<int, Dimensions> floor(position<double, Dimensions> const& pos)
      {
          position<double, Dimensions> that;
          std::transform(pos.begin(), pos.end(), that.begin(), [](auto& p){return std::floor(p);});
          return that;
      }

  }
}

#endif
