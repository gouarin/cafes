#include <particle/geometry/super_ellipsoid.hpp>
#include <cassert>
#include <iostream>

int main()
{
  {
    using pt = cafes::geometry::super_ellipsoid<2>::position_type;
    using st = cafes::geometry::super_ellipsoid<2>::shapes_type;

    cafes::geometry::super_ellipsoid<2> se{ {0.,0.}, {1.,1.}, 1. };

    //for(std::size_t y = -1;y<1)
    std::cout << se.implicit({0.5 ,0.5 }) << "\n";
    std::cout << se.implicit({0. ,1. }) << "\n";
    std::cout << se.implicit({1.1, 0. }) << "\n";
    //assert( se.is_inside({0.99 ,0.99 }) );

    double p = M_PI/5.;
    std::vector<double> sf{0.,p,2*p,3*p,4*p,5*p,6*p,7*p,8*p,9*p};

    auto s = se.surface( sf );
    std::for_each(s.begin(),s.end(), [](auto p) { std::cout << p[0] << ", " << p[1] << "\n";} );
    
    std::array<double, 2> h {.1, .1};
    auto bi = se.bounding_box(h);
    std::cout << bi.bottom_left[0] << ", " << bi.bottom_left[1] << "\n";
    std::cout << bi.upper_right[0] << ", " << bi.upper_right[1] << "\n";

    auto b = se.bounding_box();
    std::cout << b.bottom_left[0] << ", " << b.bottom_left[1] << "\n";
    std::cout << b.upper_right[0] << ", " << b.upper_right[1] << "\n";

    //assert( se.is_outside({1. ,1.}) );
  }

  return 0;
}
