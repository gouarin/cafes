#ifndef PARTICLE_GEOMETRY_CROSS_PRODUCT_HPP_INCLUDED
#define PARTICLE_GEOMETRY_CROSS_PRODUCT_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <particle/geometry/vector.hpp>

#include <array>

namespace cafes
{
  namespace geometry
  {

    double cross_product(position<2, double> x, std::array<double, 2> y){
      return x[0]*y[1] - x[1]*y[0];
    }

    double cross_product(vector<2, double> x, vector<2, double> y){
      return x[0]*y[1] - x[1]*y[0];
    }

    double cross_product(std::array<double, 2> x, std::array<double, 2> y){
      return x[0]*y[1] - x[1]*y[0];
    }

    vector<2, double> cross_product(double rz, vector<2, double> x){
      return { -rz*x[1],
                rz*x[0]
             };
    }

    vector<2, double> cross_product(vector<2, double> x, double rz){
      return {  rz*x[1],
               -rz*x[0]
             };
    }

    vector<2, double> cross_product(double rz, position<2, double> x){
      return { -rz*x[1],
                rz*x[0]
             };
    }

    position<2, double> cross_product(position<2, double> x, double rz){
      return {  rz*x[1],
               -rz*x[0]
             };
    }

    vector<3, double> cross_product(vector<3, double> x, vector<3, double> y){
      return { x[1]*y[2] - x[2]*y[1],
               x[2]*y[0] - x[0]*y[2],
               x[0]*y[1] - x[1]*y[0]
             };
    }

  }
}
#endif