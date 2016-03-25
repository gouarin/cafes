#ifndef PARTICLE_PHYSICS_VELOCITY_HPP_INCLUDED
#define PARTICLE_PHYSICS_VELOCITY_HPP_INCLUDED

#include <particle/geometry/vector.hpp>
#include <array>

namespace cafes
{
  namespace physics
  {
      template<std::size_t Dimensions>
      struct velocity : private std::array<double,Dimensions>
      {
        using parent = std::array<double,Dimensions>;
        using parent::operator[];
        using parent::begin;
        using parent::end;

        velocity() = default;
        velocity(velocity const&) = default;

        template<typename... Ts>
        velocity( Ts... ts ) : parent{{ static_cast<double>(ts)...}}
        {
           static_assert( sizeof...(Ts) <= Dimensions
                        , "Size mismatch in velocity variadic constructor"
                        );
        }

        velocity& operator+=(velocity const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
            return *this;
        }

        velocity& operator-=(geometry::vector<Dimensions, double> const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
            return *this;
        }
      };

      template<std::size_t Dimensions>
      velocity<Dimensions> operator+( velocity<Dimensions> const& v1
                                    , velocity<Dimensions> const& v2
                                    )
      {
          velocity<Dimensions> that{v1};
          return that += v2;
      }

      template<std::size_t Dimensions>
      velocity<Dimensions> operator-( velocity<Dimensions> const& v1
                                    , geometry::vector<Dimensions, double> const& v2
                                    )
      {
          velocity<Dimensions> that{v1};
          return that -= v2;
      }
  }
}

#endif
