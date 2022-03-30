// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef CAFES_ALGORITHM_ITERATE_HPP_INCLUDED
#define CAFES_ALGORITHM_ITERATE_HPP_INCLUDED

#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>

namespace cafes
{
    namespace algorithm
    {

        template<typename Box, typename Function, typename Position>
        void iterate_impl(Box const &b, Function &&f, Position &p, std::size_t step,
                          std::integral_constant<std::size_t, 0> const &)
        {
            std::forward<Function>(f)(p);
        }

        template<typename Box, typename Function, typename Position, typename Index>
        void iterate_impl(Box const &b, Function &&f, Position &p, std::size_t step, Index const &)
        {
            static constexpr std::size_t n = Index::value - 1;

            for (p[n] = b.bottom_left[n]; p[n] < b.upper_right[n]; p[n] += step)
            {
                iterate_impl(b, std::forward<Function>(f), p, step, std::integral_constant<std::size_t, n>{});
            }
        }

        template<typename Box, typename Function>
        void iterate(Box const &b, Function &&f, std::size_t step=1)
        {
            typename Box::position_type pos;
            iterate_impl(b, f, pos, step, std::integral_constant<std::size_t, Box::dimensions>{});
        }
    } // namespace algorithm
} // namespace cafes

#endif