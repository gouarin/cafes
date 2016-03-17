#ifndef PARTICLE_GEOMETRY_VECTOR_HPP_INCLUDED
#define PARTICLE_GEOMETRY_VECTOR_HPP_INCLUDED

#include <particle/physics/velocity.hpp>
#include <array>
#include <iostream>
#include <type_traits>

namespace cafes
{
  namespace geometry
  {
    template<std::size_t Dimensions, typename T>
    struct vector : private std::array<T, Dimensions>
    {
      using parent = std::array<T, Dimensions>;
      using parent::operator[];
      using parent::begin;
      using parent::end;

      vector() = default;
      vector(vector const&) = default;


      template<typename U> vector(vector<Dimensions, U> const& v)
      {
        for(std::size_t i=0; i<Dimensions; ++i)
          (*this)[i] = static_cast<T>(v[i]);
      } 

      template<typename... Ts>
      vector( Ts... ts ) : parent{{ static_cast<T>(ts)...}}
      {
         static_assert( sizeof...(Ts) <= Dimensions
                      , "Size mismatch in vector variadic constructor"
                      );
      }

      friend std::ostream &operator<<( std::ostream &output, 
                                       vector const& vec )
      { 
        output << "vector( ";

        for(std::size_t i=0; i<Dimensions; ++i)
          output << vec[i] << " ";

        output << ")";

        return output;            
      }

      // vector& operator*=(quaternion const& q)
      // {
      //     //todo: MATH!
      //     return *this;
      // }

      vector& operator*=(vector<Dimensions, T> const& v)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= v[i];
          return *this;
      }

      vector& operator*=(T const& c)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] *= c;
          return *this;
      }

      vector& operator+=(vector<Dimensions, T> const& v)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] += v[i];
          return *this;
      }

      vector& operator/=(T const& c)
      {
          for(std::size_t i=0;i<Dimensions;++i) (*this)[i] /= c;
          return *this;
      }

      // vector& operator=(vector<Dimensions, T> const& other)
      // {
      //   for(std::size_t i=0;i<Dimensions;++i) (*this)[i] = other[i];
      //   return *this;
      // }

    };
    
    template<std::size_t Dimensions, typename T>
    auto operator*(vector<Dimensions, T> const& v, T const& c)
    {
      vector<Dimensions, T> that{v};
      return that*=c;
    }

    template<std::size_t Dimensions, typename T>
    auto operator*(T const& c, vector<Dimensions, T> const& v)
    {
      return v*c;
    }

    template<std::size_t Dimensions, typename T>
    auto operator+(vector<Dimensions, T> const& v1, vector<Dimensions, T> const& v2)
    {
      vector<Dimensions, T> that{v1};
      return that+=v2;
    }

    template<std::size_t Dimensions, typename T>
    auto operator-(vector<Dimensions, T> const& v)
    {
      vector<Dimensions, T> that{v};
      return that *= -1;
    }

  }
}
#endif