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

#ifndef PARTICLE_GEOMETRY_CIRCLE_HPP_INCLUDED
#define PARTICLE_GEOMETRY_CIRCLE_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/super_ellipsoid.hpp>

namespace cafes
{
  namespace geometry
  {

    template<typename T=double>
    struct circle : public super_ellipsoid<T, 2>
    {
      using position_type = position<T, 2>;

      using parent = super_ellipsoid<T, 2>;
      using parent::perimeter;
      using parent::shape_factors_;

      position_type center;

      circle(position_type const& c, double const& radius, quaternion q={}):
             super_ellipsoid<T, 2>(c, {{radius, radius}}, 2, q), center{c}
      {
      }
      
      double surface_area() const
      {
        return 2*M_PI*shape_factors_[0];
      }

      double volume() const
      {
        return M_PI*shape_factors_[0]*shape_factors_[0];
      }

      double Cd_R() const
      {
        return 1./(shape_factors_[0]*shape_factors_[0]);
      }

      double Ci_R() const
      {
        return 2./(shape_factors_[0]*shape_factors_[0]);
      }

      std::vector<position_type>
      surface(double const& k, double tol=1e-2) const
      {
        std::size_t npoints = surface_area()/k;
        double step = 2*M_PI/npoints;
        std::vector<position<double, 2>> that;
        for(double theta=0; theta<2*M_PI; theta+=step)
        {
            that.push_back({shape_factors_[0]*std::cos(theta),
                            shape_factors_[0]*std::sin(theta)});
        }
        std::for_each(that.begin(), that.end(),[&](auto& p){p += center;});
        return that;
      }
      
    };
  }

  geometry::circle<> make_circle(geometry::position<double, 2> const& center, double const& radius, geometry::quaternion q={})
  {
    return {center, radius, q};
  }
}

#endif