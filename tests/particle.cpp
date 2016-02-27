#include <particle/particle.hpp>
#include <fem/quadrature.hpp>
#include <particle/geometry/super_ellipsoid.hpp>
#include <particle/geometry/linspace.hpp>
#include <cassert>
#include <iostream>
#include <algorithm>

template<typename T>
decltype(auto) index(T&& c)
{
  return std::forward<T>(c).first;
}

template<typename T>
decltype(auto) position(T&& c)
{
  return std::forward<T>(c).second;
}

int main()
{
  {
    cafes::geometry::super_ellipsoid<2> se{ {1.,1.}, {.5,.5}, 1. };
    cafes::particle<cafes::geometry::super_ellipsoid<2>> pt(se, {0.,150.}, 0.25);

    double p = M_PI/5.;
    std::vector<double> sf{0.,p,2*p,3*p,4*p,5*p,6*p,7*p,8*p,9*p};

    //auto s = pt.surface( sf );
    //std::for_each(s.begin(),s.end(), [](auto p) { std::cout << p[0] << ", " << p[1] << "\n";} );

    cafes::geometry::box<2, double> b{{0., 0.}, {2., 2.}};
    std::array<double, 2> h {.1, .1};
    cafes::geometry::box<2, int> bi{b.bottom_left/h, b.upper_right/h};

    // auto indices = cafes::find_insides(pt, b, h);
    // std::for_each(indices.begin(),indices.end(), [](auto p) { std::cout << p[0] << ", " << p[1] << "\n";} );

    auto surf_points = pt.surface(cafes::geometry::linspace(0., 2.*M_PI, 1000, false));
    std::cout << "size " << surf_points.size() << "\n";
    //std::for_each(surf_points.begin(),surf_points.end(), [](auto p) { std::cout << p[0] << ", " << p[1] << "\n";} );
    //assert( se.is_outside({1. ,1.}) );

    auto pts = cafes::find_surf_points_insides(surf_points, bi, h);
    std::cout << "size " << pts.size() << "\n";
    //std::for_each(pts.cbegin(), pts.cend(), [h](auto p) { std::cout << p.first[0] << ", " << p.first[1] << "\n";} );

    std::vector<std::vector<double>> sol(2./h[0]+1, std::vector<double>(2./h[1]+1));

    for(std::size_t k=0; k<10000; ++k){
      std::vector<double> g;
      g.resize(pts.size());

      for(std::size_t i=0; i<pts.size(); ++i){
        auto bfunc = cafes::fem::P1_integration(position(pts[i]), h);
        auto ielem = cafes::fem::get_element(index(pts[i]));
        for (std::size_t j=0; j<bfunc.size(); ++j){
          g[i] += sol[ielem[j][0]][ielem[j][1]]*bfunc[j];
        }
        auto sum = std::accumulate(g.cbegin(), g.cend(), 0);
      }

    }
    

  }

  return 0;
}
