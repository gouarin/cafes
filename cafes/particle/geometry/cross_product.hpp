#ifndef PARTICLE_GEOMETRY_CROSS_PRODUCT_HPP_INCLUDED
#define PARTICLE_GEOMETRY_CROSS_PRODUCT_HPP_INCLUDED

#include <particle/physics/velocity.hpp>

namespace cafes
{
  namespace geometry
  {

    template<template<typename, std::size_t> class T, typename U>
    U cross_product(T<U, 2> const& x, T<U, 2> const& y)
    {
       return x[0]*y[1] - x[1]*y[0];
    }

    template<template<typename, std::size_t> class T1, template<typename, std::size_t> class T2, typename U>
    U cross_product(T1<U, 2> const& x, T2<U, 2> const& y)
    {
       return x[0]*y[1] - x[1]*y[0];
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(double rz, T<U, 2> const& x){
      return { -rz*x[1],
                rz*x[0]
             };
    }


    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(T<U, 1> const& rz, T<U, 2> const& x){
      return { -rz[0]*x[1],
                rz[0]*x[0]
             };
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(T<U, 2> const& x, double rz){
      return {  rz*x[1],
               -rz*x[0]
             };
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 2> cross_product(T<U, 2> const& x, T<U, 1> const & rz){
      return {  rz[0]*x[1],
               -rz[0]*x[0]
             };
    }

    template<template<typename, std::size_t> class T, typename U>
    T<U, 3> cross_product(T<U, 3> const& x, T<U, 3> const& y){
      return { x[1]*y[2] - x[2]*y[1],
               x[2]*y[0] - x[0]*y[2],
               x[0]*y[1] - x[1]*y[0]
             };
    }

  }
}
#endif