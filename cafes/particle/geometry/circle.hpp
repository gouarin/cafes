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

      circle(position_type const& c, double const& radius, quaternion q={}):
             super_ellipsoid<T, 2>(c, {radius, radius}, 2, q)
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
      
    };
  }

  geometry::circle<> make_circle(geometry::position<double, 2> const& center, double const& radius, geometry::quaternion q={})
  {
    return {center, radius, q};
  }
}

#endif