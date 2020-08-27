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

#ifndef PARTICLE_PARTICLE_HPP_INCLUDED
#define PARTICLE_PARTICLE_HPP_INCLUDED

#include <array>
#include <algorithm>
#include <initializer_list>
#include <numeric>
#include <particle/physics/force.hpp>
#include <particle/physics/velocity.hpp>
#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>
#include <particle/geometry/vector.hpp>

namespace cafes
{
  template<typename Shape>
  struct particle : private Shape
  {
      using shape_type             = Shape;
      using dimension_type         = typename Shape::dimension_type;
      using velocity_type          = physics::velocity<dimension_type::value>;
      using angular_velocity_type  = typename std::conditional<dimension_type::value==2,
                                                               double, 
                                                               geometry::vector<double, 3>>::type;
      using force_type             = physics::force<dimension_type::value>;

      force_type force_{};
      double rho_;
      particle( Shape const& s, force_type const& f, double r )
      : Shape(s), force_(f), rho_(r)
      {}

      particle( Shape const& s, velocity_type const& v, angular_velocity_type const& ang_v)
      : Shape(s), velocity_(v), angular_velocity_(ang_v)
      {}

      using Shape::surface;
      using Shape::contains;
      using Shape::bounding_box;
      using Shape::center_;
      using Shape::shape_factors_;
      using Shape::surface_area;
      using Shape::volume;
      using Shape::Cd_R;
      using Shape::Ci_R;

      velocity_type velocity_;
      angular_velocity_type angular_velocity_;

      static constexpr std::size_t dimensions = dimension_type::value;
      // std::conditional<dimension_type::value==2,
      //                  double, 
      //                  geometry::vector<double, 3>> angular_velocity_;
  };

  template<typename Shape>
  auto find_fluid_points_insides( particle<Shape> const& p, cafes::geometry::box<int, 2> const& b, std::array<double, 2> const& h)
  {
      std::vector<geometry::position<int, 2>> that;
      that.reserve(b.length());

      double x, y;
      std::size_t ix, iy;
      for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
          for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
          {
              cafes::geometry::position<double, 2> pt {x,y};
              if(p.contains(pt)) that.push_back({ix, iy});
          }

      return that;
  }

  template<typename Shape>
  auto find_FE_insides( particle<Shape> const& p, cafes::geometry::box<int, 2> const& b, std::array<double, 2> const& h)
  {
      std::vector<geometry::position<int, 2>> that;
      that.reserve(b.length());

      double x, y;
      std::size_t ix, iy;
      for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
      {
          for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
          {
              // bool found = true;
              // for(std::size_t i=0; i<2; ++i)
              // {
              //   for(std::size_t j=0; j<2; ++j)
              //   {
              //     cafes::geometry::position<double, 2> pt {x + i*h[0], y + j*h[1]};
              //     if(!p.contains(pt))
              //     {
              //         found = false;
              //         break;
              //     }
              //   }
              // }
              // if (found)
              // {
                that.push_back({ix, iy});
              // }
          }
      }
      return that;
  }

  template<typename Shape>
  auto find_fluid_points_insides( particle<Shape> const& p, cafes::geometry::box<int, 3> const& b, std::array<double, 3> const& h)
  {
    std::vector<geometry::position<int, 3>> that;
    that.reserve(b.length());

    double x, y, z;
    std::size_t ix, iy, iz;
    for(z=b.bottom_left[2]*h[2], iz=b.bottom_left[2]; iz<b.upper_right[2]; z+=h[2], ++iz)
        for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
            for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
            {
               cafes::geometry::position<double, 3> pt {x, y, z};
               if(p.contains(pt)) that.push_back({ix, iy, iz});
            }   

    return that;
  }

  template<typename Shape>
  auto find_FE_insides( particle<Shape> const& p, cafes::geometry::box<int, 3> const& b, std::array<double, 3> const& h)
  {
      std::vector<geometry::position<int, 2>> that;
      that.reserve(b.length());

      double x, y, z;
      std::size_t ix, iy, iz;
      for(z=b.bottom_left[2]*h[2], iz=b.bottom_left[2]; iz<b.upper_right[2]; z+=h[2], ++iz)
      {
        for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
        {
            for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
            {
                bool found = true;
                for(std::size_t i=0; i<2; ++i)
                {
                    for(std::size_t j=0; j<2; ++j)
                    {
                        for(std::size_t k=0; k<2; ++k)
                        {
                            cafes::geometry::position<double, 3> pt {x + i*h[0], y + j*h[1], z + k*h[2]};
                            if(!p.contains(pt))
                            {
                                found = false;
                                break;
                            }
                        }
                    }
                }
                if (found)
                {
                    that.push_back({ix, iy});
                }
            }
        }
      }
      return that;
  }

  template<std::size_t Dimensions>
  auto find_surf_points_insides( std::vector<geometry::position<double, Dimensions>> const& surf_p, cafes::geometry::box<int, Dimensions> const& b, std::array<double, Dimensions> const& h)
  {
    std::vector<std::pair<geometry::position<int, Dimensions>, geometry::position<double, Dimensions>>> that;
    that.reserve(surf_p.size());
    
    for(std::size_t i=0; i<surf_p.size(); ++i){
        auto surf_pi = static_cast<geometry::position<int, Dimensions>>(surf_p[i]/h);
        if (cafes::geometry::point_inside(b, surf_pi)){
            that.push_back(std::make_pair(surf_pi, surf_p[i] - surf_pi*h));
        }
    }
    return that;
  }

  template<std::size_t Dimensions>
  auto find_radial_surf_points_insides( std::vector<geometry::position<double, Dimensions>> const& surf_p, 
                                       cafes::geometry::box<int, Dimensions> const& b, 
                                       std::array<double,Dimensions> const& h,
                                       geometry::position<double, Dimensions> center)
  {
    std::vector<geometry::vector<double, Dimensions>> that;
    that.reserve(surf_p.size());

    for(std::size_t i=0; i<surf_p.size(); ++i){
        auto surf_pi = static_cast<geometry::position<int, Dimensions>>(surf_p[i]/h);
        if (cafes::geometry::point_inside(b, surf_pi)){
            that.push_back(surf_p[i] - center);
        }
    }
    return that;
  }

  template<typename Shape, std::size_t Dimensions>
  auto position_diff(particle<Shape> const& p1, particle<Shape> const& p2)
  {
    geometry::position<double, Dimensions> pos_diff;
    auto diff = [](double x, double y){return x-y;};
    std::transform(p1.center_.begin(), p1.center_.end(), p2.center_.begin(), pos_diff.begin(), diff);
    return pos_diff;
  }


  template<typename Shape, std::size_t Dimensions>
  auto velocity_diff(particle<Shape> const& p1, particle<Shape> const& p2)
  {
    physics::velocity<Dimensions> vel_diff;
    auto diff = [](double x, double y){return x-y;};
    std::transform(p1.velocity_.begin(), p1.velocity_.end(), p2.velocity_.begin(), vel_diff.begin(), diff);
    return vel_diff;
  }

  template<typename Shape, std::size_t Dimensions>
  auto distance(particle<Shape> const& p1, particle<Shape> const& p2)
  {
    double dist = 0.;
    auto pos_diff = position_diff<Shape, Dimensions>(p1, p2);
    
    dist = std::inner_product(pos_diff.begin(), pos_diff.end(), pos_diff.begin(), 0.);
    return std::sqrt(dist);
  }

  // auto kernel_num_count = [](auto const& p, auto const& h, auto const& hs, auto scale, auto& num)
  // {
  //   auto const kernel_pos = [&](auto const& pos){
  //     auto const kernel = [&](auto const& pos_scale)
  //     {
  //       using position_type = geometry::position<double, pos_scale.dimensions>;
  //       position_type pts = pos*h + pos_scale*hs;
  //       if (p.contains(pts))
  //         num++;
  //     };
  //     using position_type = geometry::position<int, pos.dimensions>;
  //     position_type p1, p2;
  //     p1.fill(0);
  //     p2.fill(scale);
  //     geometry::box<int, pos.dimensions> box{ p1, p2};
  //     algorithm::iterate(box, kernel);
  //   };
  //   return kernel_pos;
  // };

  auto kernel_num_count = [](auto const& p, auto const& h, auto const& hs, auto& box_scale, auto& num)
  {
    auto const kernel_pos = [&](auto const& pos){
      auto const kernel = [&](auto const& pos_scale)
      {
        auto pts = pos*h + pos_scale*hs;
        if (p.contains(pts))
          num++;
      };
      algorithm::iterate(box_scale, kernel);
    };
    return kernel_pos;
  };

  template<std::size_t Dimensions,
           typename part_type,
           typename surf_type,
           typename radial_type,
           typename nb_type,
           typename num_type,
           typename box_type,
           typename dpart_type>
  auto set_materials(part_type& parts, surf_type& surf_points, radial_type& radial_vec,
                     nb_type& nb_surf_points, num_type& num, box_type const& box,
                     std::array<double, Dimensions> const &h, dpart_type const& dpart, std::size_t const scale)
  {
    surf_points.resize(parts.size());
    radial_vec.resize(parts.size());
    nb_surf_points.resize(parts.size());
    num.resize(parts.size());

    std::size_t size = 0;
    std::size_t ipart = 0;

    geometry::position<std::size_t, Dimensions> p1, p2;
    p1.fill(0);
    p2.fill(scale);
    geometry::box<std::size_t, Dimensions> box_scale{ p1, p2};

    std::array<double, Dimensions> hp;
    for(std::size_t i = 0; i<h.size(); ++i)
    {
      hp[i] = 2*h[i];
    }

    std::array<double, Dimensions> hs;
    for(std::size_t d=0; d<Dimensions; ++d)
      hs[d] = hp[d]/scale;

    for(auto& p: parts){
      auto pbox = p.bounding_box(hp);
      if (geometry::intersect(box, pbox)){
        auto new_box = geometry::box_inside(box, pbox);
        auto pts = find_fluid_points_insides(p, new_box, hp);
        size += pts.size();

        // auto spts = p.surface(dpart);
        auto spts = p.surface(hs[0]);
        auto spts_valid = find_surf_points_insides(spts, new_box, hp);
        surf_points[ipart].assign(spts_valid.begin(), spts_valid.end());
        nb_surf_points[ipart] = surf_points[ipart].size();

        auto radial_valid = find_radial_surf_points_insides(spts, new_box, hp, p.center_);
        radial_vec[ipart].assign(radial_valid.begin(), radial_valid.end());
        
        algorithm::iterate(new_box, kernel_num_count(p, hp, hs, box_scale, num[ipart]));

      }
      ipart++;
    }

    MPI_Allreduce(MPI_IN_PLACE, nb_surf_points.data(), parts.size(), MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, num.data(), parts.size(), MPI_INT, MPI_SUM, PETSC_COMM_WORLD);

    return size;
  }

  template<typename ST>
  particle<ST> make_particle_with_force(ST const& se, physics::force<ST::dimension_type::value> const& a, double d)
  {
    return {se,a,d};
  }

  template<typename ST>
  particle<ST> make_particle_with_velocity(ST const& se,
                                           physics::velocity<ST::dimension_type::value> const& v,
                                           typename std::conditional<ST::dimension_type::value==2,
                                                                     double, 
                                                                     physics::velocity<3>>::type const& ang_v)
  {
    return {se, v, ang_v};
  }
}
#endif
