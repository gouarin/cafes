#ifndef PARTICLE_GEOMETRY_UNIFORM_SAMPLING_HPP_INCLUDED
#define PARTICLE_GEOMETRY_UNIFORM_SAMPLING_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

namespace cafes
{
  namespace geometry
  {

    double theta(double const& r1, double const& r2, double const& theta, double const& eps, 
                 double const & K, double const& tol)
    {
      if (theta >= tol){
        double denom = r1*r1*std::pow(std::cos(theta), 2*eps)*std::pow(std::sin(theta), 4)
                     + r2*r2*std::pow(std::sin(theta), 2*eps)*std::pow(std::cos(theta), 4);
        return K/eps*std::cos(theta)*std::sin(theta)/std::sqrt(denom);
      }
      else if (theta < tol && theta > 0)
        return std::pow(std::abs(K/r2  + std::pow(theta, eps)), 1./eps) - theta;
      else
        return std::pow(std::abs(K/r2), 1./eps);
    }


    template<std::size_t Dimensions>
    std::vector<position<Dimensions, double>> uniform_sampling(std::array<double, Dimensions> r, 
                                                          double const& eps, double const& k, 
                                                          double const& tol)
    {
      std::vector<position<Dimensions, double>> that;
      std::vector<double> save_eta;
      double eta = 0.;

      while (eta < M_PI/4){
        save_eta.push_back(eta);
        auto x = r[0]*std::pow(std::cos(eta), eps);
        auto y = r[1]*std::pow(std::sin(eta), eps);    

        auto update_eta = theta(1., 1., eta, eps, k, tol);
        eta += update_eta;

        that.push_back({x, y});
      }
      auto func = [&](auto eta){
        auto x = r[0]*std::pow(std::sin(eta), eps);
        auto y = r[1]*std::pow(std::cos(eta), eps);

        that.push_back({x, y});
      };

      std::for_each(save_eta.crbegin(), save_eta.crend(), func);
      auto size = that.size();
      for(std::size_t i = 1; i<size; ++i)
        that.push_back({-that[size-i-1][0], that[size-i-1][1]});

      for(std::size_t i = 1; i<size; ++i)
        that.push_back({-that[i][0], -that[i][1]});
      
      for(std::size_t i = 1; i<size; ++i)
        that.push_back({that[size-i-1][0], -that[size-i-1][1]});
      
      that.erase(that.begin());
      return that;
    }

    template<std::size_t Dimensions>
    std::vector<position<Dimensions, double>> uniform_sampling(std::array<double, Dimensions> r, 
                                                          double const& eps, double const& k,
                                                          double const& z, double const& phi, 
                                                          double const& tol)
    {
      std::vector<position<Dimensions, double>> that;
      std::vector<double> save_eta;
      double eta = 0.;

      while (eta < M_PI/4){
        save_eta.push_back(eta);
        auto x = r[0]*std::pow(std::cos(eta), eps)*phi;
        auto y = r[1]*std::pow(std::sin(eta), eps)*phi;    

        if (phi == 0)
          eta = M_PI/4;
        else{
          auto update_eta = theta(phi, phi, eta, eps, k, tol);
          eta += update_eta;
        }

        // FIX THIS !!!
        //that.push_back({x, y}, {-x, y}, {x, -y}, {-x, -y});
        that.push_back({x, y, z});
      }
      auto func = [&](auto eta){
        auto x = r[0]*std::pow(std::sin(eta), eps)*phi;
        auto y = r[1]*std::pow(std::cos(eta), eps)*phi;

        // FIX THIS !!!
        //that.push_back({x, y}, {-x, y}, {x, -y}, {-x, -y});
        that.push_back({x, y, z});
      };

      std::for_each(save_eta.crbegin(), save_eta.crend(), func);
      auto size = that.size();
      for(std::size_t i = 1; i<size; ++i)
        that.push_back({-that[size-i-1][0], that[size-i-1][1], that[size-i-1][2]});

      for(std::size_t i = 1; i<size; ++i)
        that.push_back({-that[i][0], -that[i][1], that[i][2]});
      
      for(std::size_t i = 1; i<size; ++i)
        that.push_back({that[size-i-1][0], -that[size-i-1][1], that[size-i-1][2]});
      
      that.erase(that.begin());
      return that;
    }
  }
}
#endif