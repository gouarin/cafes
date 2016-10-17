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

#ifndef PARTICLE_GEOMETRY_QUATERNION_HPP_INCLUDED
#define PARTICLE_GEOMETRY_QUATERNION_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/vector.hpp>
#include <particle/geometry/cross_product.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numeric>

namespace cafes
{
  namespace geometry
  {
    // template<std::size_t Dimensions>
    // struct quaternion;

    // template<std::size_t Dimensions>
    // auto operator*( quaternion<Dimensions> const& q1
    //               , quaternion<Dimensions> const& q2
    //               );

    struct quaternion
    {
      using vector_type   = vector<double, 3>;

      vector<double, 4> components_;

      quaternion() = default;

      quaternion(double angle, vector_type axes)
      {
        components_ = { std::sin(angle/2)*axes[0]
                      , std::sin(angle/2)*axes[1]
                      , std::sin(angle/2)*axes[2]
                      , std::cos(angle/2)};
        normalize();
      }

      quaternion(double angle)
      {
        *this = {angle, {0, 0, 1}};
      }

      friend std::ostream &operator<<( std::ostream &output, 
                                       quaternion const& q )
      { 
        output << "quaternion( ";

        for(std::size_t i=0; i<4; ++i)
          output << q.components_[i] << " ";

        output << ")";

        return output;            
      }


      void normalize()
      {
        auto norm = std::inner_product(components_.begin(), components_.end(), components_.begin(), 0.);
        components_ /= norm;
      }

      quaternion conj() const
      {
       quaternion that{*this}; 
       for(std::size_t i=0; i<3; ++i) that.components_[i] = -components_[i];
       return that; 
      }


      // quaternion& operator*=(quaternion<Dimensions> const& q)
      // {
      //   auto tmp = coef_*q.coef_ - std::inner_product(components_.begin(), components_.end(), q.components_.begin(), 0);
      //   components_ = coef_*q.components_ + q.coef_*components_ + cross_product(components_, q.components_);
      //   coef_ = tmp;
      //   //normalize();
      // }

      //todo : proper interface for quaternion compositon and application
      position<double, 2> rotate(position<double, 2> pos) const
      {
        auto& x = components_[0];
        auto& y = components_[1];
        auto& z = components_[2];
        auto& w = components_[3];
        return { (1-2*y*y-2*z*z)*pos[0] +   (2*x*y-2*z*w)*pos[1],
                   (2*x*y+2*z*w)*pos[0] + (1-2*x*x-2*z*z)*pos[1]
               };
      }

      position<double, 3> rotate(position<double, 3> pos) const
      {
        auto& x = components_[0];
        auto& y = components_[1];
        auto& z = components_[2];
        auto& w = components_[3];
        return { (1-2*y*y-2*z*z)*pos[0] +   (2*x*y-2*z*w)*pos[1] +   (2*x*z+2*y*w)*pos[2],
                   (2*x*y+2*z*w)*pos[0] + (1-2*x*x-2*z*z)*pos[1] +   (2*y*z-2*x*w)*pos[2],
                   (2*x*z-2*y*w)*pos[0] +   (2*y*z+2*x*w)*pos[1] + (1-2*x*x-2*y*y)*pos[2]
               };
      }

      bool is_rotate() const
      {
        return std::any_of(components_.begin(), components_.end(), [](double c){ return c!=0; });
      }
    };

    // template<std::size_t Dimensions>
    // auto operator*( quaternion<Dimensions> const& q1
    //               , quaternion<Dimensions> const& q2
    //               )
    // {
    //     quaternion<Dimensions> that{q1};
    //     // fix this !!
    //     // return that *= q2;
    //     that.coef_ = q1.coef_*q2.coef_ - std::inner_product(q1.components_.begin(), q1.components_.end(), q2.components_.begin(), 0);
    //     that.components_ = q1.coef_*q2.components_ + q2.coef_*q1.components_ + cross_product(q1.components_, q2.components_);
    //     that.normalize();

    //     return that;
    // }

  }
}

#endif
