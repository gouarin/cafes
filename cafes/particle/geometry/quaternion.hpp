#ifndef PARTICLE_GEOMETRY_QUATERNION_HPP_INCLUDED
#define PARTICLE_GEOMETRY_QUATERNION_HPP_INCLUDED

#include <array>

namespace cafes
{
  namespace geometry
  {
      struct quaternion
      {
        std::array<double,4> components_;

        //todo : proper interface for quaternion compositon and application
      };
  }
}

#endif
