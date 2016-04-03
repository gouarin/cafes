#ifndef PARTICLE_PHYSICS_FORCE_HPP_INCLUDED
#define PARTICLE_PHYSICS_FORCE_HPP_INCLUDED

#include <array>
#include <algorithm>
#include <initializer_list>

namespace cafes
{
  namespace physics
  {
      template<std::size_t Dimensions>
      struct force : private std::array<double,Dimensions>
      {
        using parent = std::array<double,Dimensions>;
        using parent::operator[];

        force() = default;
        force(force const&) = default;

        template<typename... Ts>
        force( Ts... ts ) : parent{{ static_cast<double>(ts)...}}
        {
           static_assert( sizeof...(Ts) <= Dimensions
                        , "Size mismatch in force variadic constructor"
                        );
        }

        force& operator*=(double scale)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= scale;
            return *this;
        }

        force& operator+=(double scale)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += scale;
            return +this;
        }
      };

      template<std::size_t Dimensions>
      force<Dimensions> operator*(force<Dimensions> const& f, double scale)
      {
          force<Dimensions> that{f};
          return that *= scale;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator*(double scale, force<Dimensions> const& f)
      {
          force<Dimensions> that{f};
          return that *= scale;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator+(force<Dimensions> const& f, double scale)
      {
          force<Dimensions> that{f};
          return that += scale;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator+(double scale, force<Dimensions> const& f)
      {
          force<Dimensions> that{f};
          return that += scale;
      }
  }
}

#endif
