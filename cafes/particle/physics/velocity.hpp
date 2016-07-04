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

#ifndef PARTICLE_PHYSICS_VELOCITY_HPP_INCLUDED
#define PARTICLE_PHYSICS_VELOCITY_HPP_INCLUDED

//#include <particle/geometry/vector.hpp>
#include <array>

namespace cafes
{
  namespace physics
  {
      template<std::size_t Dimensions>
      struct velocity : private std::array<double,Dimensions>
      {
        using parent = std::array<double,Dimensions>;
        using parent::operator[];
        using parent::begin;
        using parent::end;

        velocity() = default;
        velocity(velocity const&) = default;

        template<typename... Ts>
        velocity( Ts... ts ) : parent{{ static_cast<double>(ts)...}}
        {
           static_assert( sizeof...(Ts) <= Dimensions
                        , "Size mismatch in velocity variadic constructor"
                        );
        }

        velocity& operator+=(velocity const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
            return *this;
        }

        // velocity& operator-=(geometry::vector<double, Dimensions> const& v)
        // {
        //     for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
        //     return *this;
        // }
      };

      template<std::size_t Dimensions>
      velocity<Dimensions> operator+( velocity<Dimensions> const& v1
                                    , velocity<Dimensions> const& v2
                                    )
      {
          velocity<Dimensions> that{v1};
          return that += v2;
      }

      // template<std::size_t Dimensions>
      // velocity<Dimensions> operator-( velocity<Dimensions> const& v1
      //                               , geometry::vector<double, Dimensions> const& v2
      //                               )
      // {
      //     velocity<Dimensions> that{v1};
      //     return that -= v2;
      // }
  }
}

#endif
