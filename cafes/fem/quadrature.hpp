// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
// OF SUCH DAMAGE.

#ifndef PARTICLE_FEM_POSITION_HPP_INCLUDED
#define PARTICLE_FEM_POSITION_HPP_INCLUDED

#include <particle/geometry/position.hpp>

namespace cafes
{
    namespace fem
    {
        auto P1_integration(geometry::position<double, 2> const& x, std::array<double, 2> const& h)
        {
            std::array<double, 4> that;    
            double c = 1./(h[0]*h[1]);
            that[0] = c*(x[0]-h[0])*(x[1]-h[1]);
            that[1] = c*x[0]*(h[1]-x[1]);
            that[2] = c*(h[0]-x[0])*x[1];
            that[3] = c*x[0]*x[1];
            return that;
        }

        auto P1_integration_grad(geometry::position<double, 2> const& x, std::array<double, 2> const& h)
        {
            std::array< std::array<double, 2>, 4> that;
            //double c = 1./(h[0]*h[1]);
            that[0] = {{ x[1]-h[1], x[0]-h[0] }};
            that[1] = {{ h[1]-x[1],     -x[0] }};
            that[2] = {{     -x[1], h[0]-x[0] }};
            that[3] = {{      x[1],      x[0] }};
            return that;
        }

        auto P1_integration(geometry::position<double, 3> const& x, std::array<double, 3> const& h)
        {
            std::array<double, 8> that; 
            double c = 1./(h[0]*h[1]*h[2]);

            that[0] = c*(h[0]-x[0])*(h[1]-x[1])*(h[2]-x[2]);
            that[1] = c*x[0]*(h[1]-x[1])*(h[2]-x[2]);
            that[2] = c*(h[0]-x[0])*x[1]*(h[2]-x[2]);
            that[3] = c*x[0]*x[1]*(h[2]-x[2]);
            that[4] = c*(h[0]-x[0])*(h[1]-x[1])*x[2];
            that[5] = c*x[0]*(h[1]-x[1])*x[2];
            that[6] = c*(h[0]-x[0])*x[1]*x[2];
            that[7] = c*x[0]*x[1]*x[2];
            return that;
        }

        auto P1_integration_grad(geometry::position<double, 3> const& x, std::array<double, 3> const& h)
        {
            std::array<std::array<double, 3>, 8> that; 

            that[0] = {{ -(h[1]-x[1])*(h[2]-x[2]), -(h[0]-x[0])*(h[2]-x[2]), -(h[0]-x[0])*(h[1]-x[1]) }};
            that[1] = {{  (h[1]-x[1])*(h[2]-x[2]),        -x[0]*(h[2]-x[2]),        -x[0]*(h[1]-x[1]) }};
            that[2] = {{        -x[1]*(h[2]-x[2]),  (h[0]-x[0])*(h[2]-x[2]),        -(h[0]-x[0])*x[1] }};
            that[3] = {{         x[1]*(h[2]-x[2]),         x[0]*(h[2]-x[2]),               -x[0]*x[1] }};
            that[4] = {{        -(h[1]-x[1])*x[2],        -(h[0]-x[0])*x[2],  (h[0]-x[0])*(h[1]-x[1]) }};
            that[5] = {{         (h[1]-x[1])*x[2],               -x[0]*x[2],         x[0]*(h[1]-x[1]) }};
            that[6] = {{               -x[1]*x[2],         (h[0]-x[0])*x[2],         (h[0]-x[0])*x[1] }};
            that[7] = {{                x[1]*x[2],                x[0]*x[2],                x[0]*x[1] }};
            return that;
        }

        auto get_element(geometry::position<int, 2> const& ix){
            std::array<std::array<int, 2>, 4> that;

            that[0] = {ix[0]  , ix[1]  };
            that[1] = {ix[0]+1, ix[1]  };
            that[2] = {ix[0]  , ix[1]+1};
            that[3] = {ix[0]+1, ix[1]+1};
            return that;
        }

        inline PetscInt DMDALocalIndex2D(DMDALocalInfo& info, int i, int j)
        {
            return info.dof*((j - info.gys)*info.gxm + i - info.gxs);
        }

        auto get_indices(DMDALocalInfo& info, geometry::position<int, 2> const& ix, int dof){
            std::array<int, 4> that;

            that[0] = dof + DMDALocalIndex2D(info, ix[0], ix[1]);
            that[1] = dof + DMDALocalIndex2D(info, ix[0]+1, ix[1]);
            that[2] = dof + DMDALocalIndex2D(info, ix[0], ix[1]+1);
            that[3] = dof + DMDALocalIndex2D(info, ix[0]+1, ix[1]+1);
            return that;
        }

        auto get_element(geometry::position<int, 3> const& ix){
            std::array<std::array<int, 3>, 8> that;

            that[0] = {ix[0]  , ix[1]  , ix[2]  };
            that[1] = {ix[0]+1, ix[1]  , ix[2]  };
            that[2] = {ix[0]  , ix[1]+1, ix[2]  };
            that[3] = {ix[0]+1, ix[1]+1, ix[2]  };
            that[4] = {ix[0]  , ix[1]  , ix[2]+1};
            that[5] = {ix[0]+1, ix[1]  , ix[2]+1};
            that[6] = {ix[0]  , ix[1]+1, ix[2]+1};
            that[7] = {ix[0]+1, ix[1]+1, ix[2]+1};
            return that;
        }

        auto get_element_4Q1(geometry::position<int, 2> const& ix){
            std::array<std::array<int, 2>, 9> that;

            that[0] = {2*ix[0]  , 2*ix[1]  };
            that[1] = {2*ix[0]+1, 2*ix[1]  };
            that[2] = {2*ix[0]+2, 2*ix[1]  };
            that[3] = {2*ix[0]  , 2*ix[1]+1};
            that[4] = {2*ix[0]+1, 2*ix[1]+1};
            that[5] = {2*ix[0]+2, 2*ix[1]+1};
            that[6] = {2*ix[0]  , 2*ix[1]+2};
            that[7] = {2*ix[0]+1, 2*ix[1]+2};
            that[8] = {2*ix[0]+2, 2*ix[1]+2};
            return that;
        }

        auto get_indices_4Q1(DMDALocalInfo& info, geometry::position<int, 2> const& ix){
            std::array<int, 18> that;

            for(std::size_t d=0; d<2; ++d)
            {
                that[0 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]  , 2*ix[1]  );
                that[1 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]+1, 2*ix[1]  );
                that[2 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]+2, 2*ix[1]  );
                that[3 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]  , 2*ix[1]+1);
                that[4 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]+1, 2*ix[1]+1);
                that[5 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]+2, 2*ix[1]+1);
                that[6 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]  , 2*ix[1]+2);
                that[7 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]+1, 2*ix[1]+2);
                that[8 + d*9] = d + DMDALocalIndex2D(info, 2*ix[0]+2, 2*ix[1]+2);
            }
            return that;
        }

        auto get_element_4Q1(geometry::position<int, 3> const& ix){
            std::array<std::array<int, 3>, 27> that;

            that[0]  = {2*ix[0]  , 2*ix[1]  , 2*ix[2]  };
            that[1]  = {2*ix[0]+1, 2*ix[1]  , 2*ix[2]  };
            that[2]  = {2*ix[0]+2, 2*ix[1]  , 2*ix[2]  };
            that[3]  = {2*ix[0]  , 2*ix[1]+1, 2*ix[2]  };
            that[4]  = {2*ix[0]+1, 2*ix[1]+1, 2*ix[2]  };
            that[5]  = {2*ix[0]+2, 2*ix[1]+1, 2*ix[2]  };
            that[6]  = {2*ix[0]  , 2*ix[1]+2, 2*ix[2]  };
            that[7]  = {2*ix[0]+1, 2*ix[1]+2, 2*ix[2]  };
            that[8]  = {2*ix[0]+2, 2*ix[1]+2, 2*ix[2]  };
            that[9]  = {2*ix[0]  , 2*ix[1]  , 2*ix[2]+1};
            that[10] = {2*ix[0]+1, 2*ix[1]  , 2*ix[2]+1};
            that[11] = {2*ix[0]+2, 2*ix[1]  , 2*ix[2]+1};
            that[12] = {2*ix[0]  , 2*ix[1]+1, 2*ix[2]+1};
            that[13] = {2*ix[0]+1, 2*ix[1]+1, 2*ix[2]+1};
            that[14] = {2*ix[0]+2, 2*ix[1]+1, 2*ix[2]+1};
            that[15] = {2*ix[0]  , 2*ix[1]+2, 2*ix[2]+1};
            that[16] = {2*ix[0]+1, 2*ix[1]+2, 2*ix[2]+1};
            that[17] = {2*ix[0]+2, 2*ix[1]+2, 2*ix[2]+1};
            that[18] = {2*ix[0]  , 2*ix[1]  , 2*ix[2]+2};
            that[19] = {2*ix[0]+1, 2*ix[1]  , 2*ix[2]+2};
            that[20] = {2*ix[0]+2, 2*ix[1]  , 2*ix[2]+2};
            that[21] = {2*ix[0]  , 2*ix[1]+1, 2*ix[2]+2};
            that[22] = {2*ix[0]+1, 2*ix[1]+1, 2*ix[2]+2};
            that[23] = {2*ix[0]+2, 2*ix[1]+1, 2*ix[2]+2};
            that[24] = {2*ix[0]  , 2*ix[1]+2, 2*ix[2]+2};
            that[25] = {2*ix[0]+1, 2*ix[1]+2, 2*ix[2]+2};
            that[26] = {2*ix[0]+2, 2*ix[1]+2, 2*ix[2]+2};

            return that;
        }
    }
}

#endif