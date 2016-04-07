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

        friend std::ostream &operator<<( std::ostream &output, 
                                         force const& f )
        { 
          output << "force( ";

          for(std::size_t i=0; i<Dimensions; ++i)
            output << f[i] << " ";

          output << ")";

          return output;            
        }

        force& operator*=(double scale)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= scale;
            return *this;
        }

        force& operator+=(double scale)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += scale;
            return *this;
        }

        force& operator+=(force<Dimensions> const& f)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += f[i];
            return *this;
        }

        force& operator-=(force<Dimensions> const& f)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= f[i];
            return *this;
        }

        void fill(double value)
        {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] = value;
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

      template<std::size_t Dimensions>
      force<Dimensions> operator+(force<Dimensions> const& f1, force<Dimensions> const& f2)
      {
          force<Dimensions> that{f1};
          return that += f2;
      }

      template<std::size_t Dimensions>
      force<Dimensions> operator-(force<Dimensions> const& f1, force<Dimensions> const& f2)
      {
          force<Dimensions> that{f1};
          return that -= f2;
      }
  }
}

#endif
