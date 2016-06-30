#ifndef PARTICLE_GEOMETRY_SPHERE_HPP_INCLUDED
#define PARTICLE_GEOMETRY_SPHERE_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/super_ellipsoid.hpp>

namespace cafes
{
  namespace geometry
  {

    template<typename T=double>
    struct sphere : public super_ellipsoid<T, 3>
    {
      using position_type = position<T, 3>;

      using parent = super_ellipsoid<T, 3>;
      using parent::perimeter;
      using parent::shape_factors_;

      sphere(position_type const& c, double const& radius, quaternion q={}):
             super_ellipsoid<T, 3>(c, {radius, radius, radius}, 2, 2, q)
      {
      }
      
      double surface_area() const
      {
        return 4*M_PI*shape_factors_[0]*shape_factors_[0];
      }

      double volume() const
      {
        return 4./3*M_PI*shape_factors_[0]*shape_factors_[0]*shape_factors_[0];
      }

      double Cd_R() const
      {
        return 1.5/(shape_factors_[0]*shape_factors_[0]);
      }

      double Ci_R() const
      {
        return 2.5/(shape_factors_[0]*shape_factors_[0]);
      }

    };
  }

  geometry::sphere<> make_sphere(geometry::position<double, 3> const& center, double const& radius, geometry::quaternion q={})
  {
    return {center, radius, q};
  }
}

#endif