#include <particle/physics/force.hpp>
#include <cassert>

int main()
{
  {
    cafes::physics::force<2> f{};
    assert( f[0] == 0. );
    assert( f[1] == 0. );
  }

  {
    cafes::physics::force<2> f{1., 2.};
    assert( f[0] == 1. );
    assert( f[1] == 2. );
  }

  {
    cafes::physics::force<2> f{1., 2.};
    f *= 3.1;
    assert( f[0] == 3.1 );
    assert( f[1] == 6.2 );
  }

  {
    cafes::physics::force<2> f{1., 2.}, g;
    g = f * 3.1;
    assert( g[0] == 3.1 );
    assert( g[1] == 6.2 );
  }

  {
    cafes::physics::force<3> f{};
    assert( f[0] == 0. );
    assert( f[1] == 0. );
    assert( f[2] == 0. );
  }

  {
    cafes::physics::force<3> f{1., 2.,3.};
    assert( f[0] == 1. );
    assert( f[1] == 2. );
    assert( f[2] == 3. );
  }

    {
      cafes::physics::force<3> f{1., 2.,3.};
      f *= 3.1;
      assert( f[0] == 3.1 );
      assert( f[1] == 6.2 );
      assert( f[2] == 9.3 );
    }

    {
      cafes::physics::force<3> f{1., 2.,3.}, g;
      g = f * 3.1;
      assert( g[0] == 3.1 );
      assert( g[1] == 6.2 );
      assert( g[2] == 9.3 );
    }

    {
      cafes::physics::force<3> f{1., 2.,3.}, g;
      g = 3.1*f;
      assert( g[0] == 3.1 );
      assert( g[1] == 6.2 );
      assert( g[2] == 9.3 );
    }

    return 0;
}
