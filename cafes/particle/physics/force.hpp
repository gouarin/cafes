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

#ifndef PARTICLE_PHYSICS_FORCE_HPP_INCLUDED
#define PARTICLE_PHYSICS_FORCE_HPP_INCLUDED

#include <array>
#include <algorithm>
#include <initializer_list>

namespace cafes
{
  namespace physics
  {
      template<std::size_t Dimensions>
      struct force : private std::array<double, Dimensions>
      {
        using parent = std::array<double, Dimensions>;
        using parent::operator[];

        force() = default;
        force(force const&) = default;

        template<typename... Ts>
        force( Ts... ts ) : parent{{ static_cast<double>(ts)...}}
        {
           static_assert( sizeof...(Ts) <= Dimensions
                        , "Size mismatch in force variadic constructor"
                        );
        }

        friend std::ostream &operator<<( std::ostream &output, 
                                         force const& f )
        { 
          output << "force( ";

          for(std::size_t i=0; i<Dimensions; ++i)
            output << f[i] << " ";

          output << ")";

          return output;            
        }

        force& operator*=(double scale)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= scale;
            return *this;
        }

        force& operator+=(double scale)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += scale;
            return *this;
        }

        force& operator+=(force<Dimensions> const& f)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += f[i];
            return *this;
        }

        force& operator-=(force<Dimensions> const& f)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= f[i];
            return *this;
        }

        void fill(double value)
        {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] = value;
        }
      };

      template<std::size_t Dimensions>
      force<Dimensions> operator*(force<Dimensions> const& f, double scale)
      {
          force<Dimensions> that{f};
          return that *= scale;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator*(double scale, force<Dimensions> const& f)
      {
          force<Dimensions> that{f};
          return that *= scale;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator+(force<Dimensions> const& f, double scale)
      {
          force<Dimensions> that{f};
          return that += scale;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator+(double scale, force<Dimensions> const& f)
      {
          force<Dimensions> that{f};
          return that += scale;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator+(force<Dimensions> const& f1, force<Dimensions> const& f2)
      {
          force<Dimensions> that{f1};
          return that += f2;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator-(force<Dimensions> const& f1, force<Dimensions> const& f2)
      {
          force<Dimensions> that{f1};
          return that -= f2;
      }
  }
}

#endif
