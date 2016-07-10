#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include <cafes.hpp>

TEST_CASE( "Box limits", "[box]" ) {
  cafes::geometry::box<double, 2> b1{ {0, 0}, {1, 1}};

  REQUIRE( (b1.bottom_left[0] == 0 && b1.bottom_left[1] == 0) );
  REQUIRE( (b1.upper_right[0] == 1 && b1.upper_right[1] == 1) );
}