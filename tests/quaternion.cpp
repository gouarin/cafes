#include <particle/geometry/position.hpp>
#include <particle/geometry/quaternion.hpp>
#include <iostream>

int main()
{
  cafes::geometry::quaternion q1{M_PI/4}, q3{};
  cafes::geometry::quaternion q2{M_PI/2, {0, 0, 1}};
  cafes::geometry::position<2, double> pos{1, 0};

  std::cout << q3 << "\n";
  std::cout << q3.is_rotate() << " " << q1.is_rotate();
  //std::cout << q1.apply(pos) << "\n";
  auto pos1 = q1.rotate(pos);
  std::cout << pos1 << "\n";
  std::cout << q1.conj().rotate(pos1) << "\n";
  
  // auto q3 = q1*q1;
  // std::cout << q1 << "\n";
  // std::cout << q2 << "\n";

  // std::cout << q3 << "\n";
  // std::cout << q3.apply(pos) << "\n";
  // q1 *= q2;
  // std::cout << q1.apply(pos) << "\n";

  return 0;
}