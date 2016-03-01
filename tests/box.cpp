#include <particle/geometry/box.hpp>
#include <particle/geometry/super_ellipsoid.hpp>
#include <cassert>
#include <iostream>

int main()
{
  {
    using pt = cafes::geometry::super_ellipsoid<2>::position_type;
    using st = cafes::geometry::super_ellipsoid<2>::shapes_type;
    cafes::geometry::super_ellipsoid<2> se{ {0.,0.}, {1.,1.}, 1. };

    cafes::geometry::box<2, double> b1=se.bounding_box(), b2{{-2.,-2.}, {5., .5}}, b3;


    std::cout << b1.bottom_left[0] << ", " << b1.bottom_left[1] << "\n";
    std::cout << b1.upper_right[0] << ", " << b1.upper_right[1] << "\n";

    std::cout << b2.bottom_left[0] << ", " << b2.bottom_left[1] << "\n";
    std::cout << b2.upper_right[0] << ", " << b2.upper_right[1] << "\n";

    b3 = overlap_box(b1, b2);
    std::cout << b3.bottom_left[0] << ", " << b3.bottom_left[1] << "\n";
    std::cout << b3.upper_right[0] << ", " << b3.upper_right[1] << "\n";

    cafes::geometry::box<2, int> bi1{{0, 0}, {2, 2}}, bi2{{1, 1}, {3, 3}}, bi3{{-2, -2}, {-1, -1}};
    assert(cafes::geometry::intersect(bi1, bi2));
    assert(!cafes::geometry::intersect(bi1, bi3));

    //assert( se.is_outside({1. ,1.}) );
  }

  return 0;
}
