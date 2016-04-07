#ifndef PARTICLE_GEOMETRY_CIRCLE_HPP_INCLUDED
#define PARTICLE_GEOMETRY_CIRCLE_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/super_ellipsoid.hpp>

namespace cafes
{
  namespace geometry
  {

    template<typename T=double>
    struct circle : public super_ellipsoid<2, T>
    {
      using position_type = position<2, T>;

      using parent = super_ellipsoid<2, T>;
      using parent::perimeter;
      using parent::shape_factors_;

      circle(position_type const& c, double const& radius, quaternion q={}):
             super_ellipsoid<2, T>(c, {radius, radius}, 2, q)
      {
      }
      
      double surface_area() const
      {
        return 2*M_PI*shape_factors_[0];
      }

      double voulume() const
      {
        return M_PI*shape_factors_[0]*shape_factors_[0];
      }
    };
  }

  geometry::circle<> make_circle(geometry::position<2,double> const& center, double const& radius, geometry::quaternion q)
  {
    return {center, radius, q};
  }
}

#endif