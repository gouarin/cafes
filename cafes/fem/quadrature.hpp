#ifndef PARTICLE_FEM_POSITION_HPP_INCLUDED
#define PARTICLE_FEM_POSITION_HPP_INCLUDED

#include <particle/geometry/position.hpp>

namespace cafes
{
    namespace fem
    {
        auto P1_integration(geometry::position<2, double> const& x, std::array<double, 2> const& h)
        {
            std::array<double, 4> that;    
            double c = 1./(h[0]*h[1]);
            that[0] = c*(x[0]-h[0])*(x[1]-h[1]);
            that[1] = c*x[0]*(h[1]-x[1]);
            that[2] = c*(h[0]-x[0])*x[1];
            that[3] = c*x[0]*x[1];
            return that;
        }

        auto P1_integration(geometry::position<3, double> const& x, std::array<double, 3> const& h)
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

        auto get_element(geometry::position<2, int> const& ix){
            std::array<std::array<int, 2>, 4> that;

            that[0] = {ix[0]  , ix[1]  };
            that[1] = {ix[0]+1, ix[1]  };
            that[2] = {ix[0]  , ix[1]+1};
            that[3] = {ix[0]+1, ix[1]+1};
            return that;
        }

        auto get_element(geometry::position<3, int> const& ix){
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
    }
}

#endif