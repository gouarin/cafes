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

namespace cafes
{
  template<typename Shape>
  struct particle : private Shape
  {
      using shape_type     = Shape;
      using dimension_type = typename Shape::dimension_type;
      using velocity_type  = physics::velocity<dimension_type::value>;
      using force_type     = physics::force<dimension_type::value>;

      force_type force_{};
      double rho_;
      particle( Shape const& s, force_type const& f, double r )
      : Shape(s), force_(f), rho_(r)
      {}

      particle( Shape const& s, velocity_type const& v, physics::velocity<3> const& ang_v)
      : Shape(s), velocity_(v), angular_velocity_(ang_v)
      {}

      using Shape::surface;
      using Shape::contains;
      using Shape::bounding_box;
      using Shape::center_;
      using Shape::shape_factors_;
      using Shape::surface_area;

      velocity_type velocity_;
      physics::velocity<3> angular_velocity_{};
      
  };

  template<typename Shape>
  auto find_fluid_points_insides( particle<Shape> const& p, cafes::geometry::box<2, int> const& b, std::array<double,2> const& h)
  {
      std::vector<geometry::position<2, int>> that;
      that.reserve(b.length());

      double x, y;
      std::size_t ix, iy;
      for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
          for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
          {
              cafes::geometry::position<2, double> pt {x,y};
              if(p.contains(pt)) that.push_back({ix, iy});
          }

      return that;
  }

  template<typename Shape>
  auto find_fluid_points_insides( particle<Shape> const& p, cafes::geometry::box<3, int> const& b, std::array<double,3> const& h)
  {
    std::vector<geometry::position<3, int>> that;
    that.reserve(b.length());

    double x, y, z;
    std::size_t ix, iy, iz;
    for(z=b.bottom_left[2]*h[2], iz=b.bottom_left[2]; iz<b.upper_right[2]; z+=h[2], ++iz)
        for(y=b.bottom_left[1]*h[1], iy=b.bottom_left[1]; iy<b.upper_right[1]; y+=h[1], ++iy)
            for(x=b.bottom_left[0]*h[0], ix=b.bottom_left[0]; ix<b.upper_right[0]; x+=h[0], ++ix)
            {
               cafes::geometry::position<3, double> pt {x, y, z};
               if(p.contains(pt)) that.push_back({ix, iy, iz});
            }   

    return that;
  }

  template<std::size_t Dimensions>
  auto find_surf_points_insides( std::vector<geometry::position<Dimensions, double>> const& surf_p, cafes::geometry::box<Dimensions, int> const& b, std::array<double,Dimensions> const& h)
  {
    std::vector<std::pair<geometry::position<Dimensions, int>, geometry::position<Dimensions, double>>> that;
    that.reserve(surf_p.size());

    for(std::size_t i=0; i<surf_p.size(); ++i){
        auto surf_pi = static_cast<geometry::position<Dimensions, int>>(surf_p[i]/h);
        if (cafes::geometry::point_inside(b, surf_pi)){
            that.push_back(std::make_pair(surf_pi, surf_p[i] - surf_pi*h));
        }
    }
    return that;
  }

  template<std::size_t Dimensions>
  auto find_radial_surf_points_insides( std::vector<geometry::position<Dimensions, double>> const& surf_p, 
                                       cafes::geometry::box<Dimensions, int> const& b, 
                                       std::array<double,Dimensions> const& h,
                                       geometry::position<Dimensions, double> center)
  {
    std::vector<geometry::position<Dimensions, double>> that;
    that.reserve(surf_p.size());

    for(std::size_t i=0; i<surf_p.size(); ++i){
        auto surf_pi = static_cast<geometry::position<Dimensions, int>>(surf_p[i]/h);
        if (cafes::geometry::point_inside(b, surf_pi)){
            that.push_back(surf_p[i] - center);
        }
    }
    return that;
  }

  template<typename Shape, std::size_t Dimensions>
  auto position_diff(particle<Shape> const& p1, particle<Shape> const& p2)
  {
    geometry::position<Dimensions, double> pos_diff;
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

  template<typename ST>
  particle<ST> make_particle(ST const& se, physics::force<ST::dimension_type::value> const& a, double d)
  {
    return {se,a,d};
  }

  template<typename ST>
  particle<ST> make_particle(ST const& se, physics::velocity<ST::dimension_type::value> const& v, physics::velocity<3> const& ang_v)
  {
    return {se, v, ang_v};
  }
}



#endif
