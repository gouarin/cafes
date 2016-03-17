#ifndef PARTICLE_GEOMETRY_SPHERE_HPP_INCLUDED
#define PARTICLE_GEOMETRY_SPHERE_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/super_ellipsoid.hpp>

namespace cafes
{
  namespace geometry
  {

    template<typename T=double>
    struct sphere : public super_ellipsoid<3, T>
    {
      using position_type = position<3, T>;

      using parent = super_ellipsoid<3, T>;
      using parent::perimeter;
      using parent::shape_factors_;

      sphere(position_type const& c, double const& radius, quaternion q={}):
             super_ellipsoid<3, T>(c, {radius, radius, radius}, 2, 2, q)
      {
      }
      
      void surface_perimeter()
      {
        perimeter = 2*M_PI*shape_factors_[0];
      }
    };
  }

  geometry::sphere<> make_sphere(geometry::position<3,double> const& center, double const& radius, geometry::quaternion q)
  {
    return {center, radius, q};
  }
}

#endif