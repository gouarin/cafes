#ifndef PARTICLE_GEOMETRY_POSITION_HPP_INCLUDED
#define PARTICLE_GEOMETRY_POSITION_HPP_INCLUDED

#include <particle/physics/velocity.hpp>
#include <array>
#include <iostream>
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
        using parent::begin;
        using parent::end;

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

        friend std::ostream &operator<<( std::ostream &output, 
                                         position const& pos )
        { 
          output << "position( ";

          for(std::size_t i=0; i<Dimensions; ++i)
            output << pos[i] << " ";

          output << ")";

          return output;            
        }

        // position& operator*=(quaternion const& q)
        // {
        //     //todo: MATH!
        //     return *this;
        // }

        position& operator*=(position<Dimensions, T> const& p)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= p[i];
            return *this;
        }

        position& operator*=(std::array<T, Dimensions> const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= v[i];
            return *this;
        }

        position& operator-=(position<Dimensions, T> const& p)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= p[i];
            return *this;
        }

        position& operator+=(position<Dimensions, T> const& p)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += p[i];
            return *this;
        }

        position& operator+=(physics::velocity<Dimensions> const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
            return *this;
        }

        position& operator+=(T const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v;
            return *this;
        }

        position& operator-=(T const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] -= v;
            return *this;
        }

        position& operator/=(std::array<T, Dimensions> const& v)
        {
            static_assert(std::is_floating_point<T>::value, "invalid operation in division");
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] /= v[i];
            return *this;
        }

        position& operator/=(T const& v)
        {
            for(std::size_t i=0;i<Dimensions;++i) (*this)[i] /= v;
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
      position<Dimensions, T> operator+( position<Dimensions, T> const& p
                                    , T const& v
                                    )
      {
          position<Dimensions, T> that{p};
          return that += v;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator+( position<Dimensions, T> const& p1
                                       , position<Dimensions, T> const& p2
                                       )
      {
          position<Dimensions, T> that{p1};
          return that += p2;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator-( position<Dimensions, T> const& p
                                    , T const& v
                                    )
      {
          position<Dimensions, T> that{p};
          return that -= v;
      }

      // template<std::size_t Dimensions, typename T>
      // position<Dimensions, T> operator*( position<Dimensions, T> const& p
      //                               , quaternion const& q
      //                               )
      // {
      //     position<Dimensions, T> that{p};
      //     return that *= q;
      // }

      template<std::size_t Dimensions, typename T, typename U>
      position<Dimensions, double> operator*( position<Dimensions, T> const& p1
                                         ,    position<Dimensions, U> const& p2
                                    )
      {
          position<Dimensions, double> that{p1};
          position<Dimensions, double> tmp{p2};
          return that *= p2;
      }

      template<std::size_t Dimensions, typename T, typename U>
      position<Dimensions, double> operator*( position<Dimensions, T> const& p
                                    ,    std::array<U, Dimensions> const& v
                                    )
      {
          position<Dimensions, double> that{p};
          return that *= v;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, double> operator-( position<Dimensions, T> const& p1
                                         ,    position<Dimensions, T> const& p2
                                    )
      {
          position<Dimensions, T> that{p1};
          return that -= p2;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator*( position<Dimensions, T> const& p1
                                    ,    position<Dimensions, T> const& p2
                                    )
      {
          position<Dimensions, T> that{p1};
          return that *= p2;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator/( position<Dimensions, T> const& p
                                    , std::array<T, Dimensions> const& v
                                    )
      {
          position<Dimensions, T> that{p};
          return that /= v;
      }

      template<std::size_t Dimensions, typename T>
      position<Dimensions, T> operator/( position<Dimensions, T> const& p
                                    , T const& v
                                    )
      {
          position<Dimensions, T> that{p};
          return that /= v;
      }

  }
}

#endif
