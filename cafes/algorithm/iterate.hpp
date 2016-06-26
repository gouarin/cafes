#ifndef CAFES_ALGORITHM_ITERATE_HPP_INCLUDED
#define CAFES_ALGORITHM_ITERATE_HPP_INCLUDED

#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>

namespace cafes
{
  namespace algorithm
  {

    template<typename Box, typename Function, typename Position>
    void iterate_impl(Box const& b, Function&& f, Position& p, std::integral_constant<std::size_t,0> const&)
    {
      std::forward<Function>(f)(p);
    }

    template<typename Box, typename Function, typename Position, typename Index>
    void iterate_impl(Box const& b, Function&& f, Position& p, Index const&)
    {
      static constexpr std::size_t n = Index::value-1;

      for( p[n] = b.bottom_left[n]; p[n] < b.upper_right[n]; ++p[n] )
      {
        iterate_impl(b, std::forward<Function>(f), p, std::integral_constant<std::size_t,n>{});
      }
    }

    template<typename Box, typename Function>
    void iterate(Box const& b, Function&& f)
    {
      typename Box::position_type pos;
      iterate_impl(b, f, pos, std::integral_constant<std::size_t,Box::dimensions>{});
    }
  }
}

#endif