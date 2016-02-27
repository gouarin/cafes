#ifndef PARTICLE_GEOMETRY_POSITION_HPP_INCLUDED
#define PARTICLE_GEOMETRY_POSITION_HPP_INCLUDED

#include <particle/physics/velocity.hpp>
#include <particle/geometry/quaternion.hpp>
#include <array>
#include <type_traits>

namespace cafes
{
  namespace geometry
  {
      template<std::size_t Dimensions, typename T>
      struct position : private std::array<T, Dimensions>
      {
        using parent = std::array<T, Dimensions>;
        using parent::operator[];

        position() = default;
        position(position const&) = default;

        template<typename U> position(position<Dimensions, U> const& p)
        {
          for(std::size_t i=0; i<Dimensions; ++i)
            (*this)[i] = static_cast<T>(p[i]);
        } 

        template<typename... Ts>
        position( Ts... ts ) : parent{{ static_cast<T>(ts)...}}
        {
           static_assert( sizeof...(Ts) <= Dimensions
                        , "Size mismatch in position variadic constructor"
                        );
        }

        position& operator*=(quaternion const& q)
        {
            //todo: MATH!
            return *this;
        }

        position& operator+=(physics::velocity<Dimensions> const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
            return *this;
        }

        position& operator/=(std::array<T, Dimensions> const& v)
        {
            static_assert(std::is_floating_point<T>::value, "invalid operation in division");
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] /= v[i];
            return *this;
        }

      };

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator+( position<Dimensions, T> const& p
                                    , physics::velocity<Dimensions> const& v
                                    )
      {
          position<Dimensions, T> that{p};
          return that += v;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator*( position<Dimensions, T> const& p
                                    , quaternion const& q
                                    )
      {
          position<Dimensions, T> that{p};
          return that *= q;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator/( position<Dimensions, T> const& p
                                    , std::array<T, Dimensions> const& v
                                    )
      {
          position<Dimensions, T> that{p};
          return that /= v;
      }

  }
}

#endif
