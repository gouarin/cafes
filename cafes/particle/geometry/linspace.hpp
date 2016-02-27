#ifndef PARTICLE_GEOMETRY_LINSPACE_HPP_INCLUDED
#define PARTICLE_GEOMETRY_LINSPACE_HPP_INCLUDED

#include <cassert>
#include <vector>

namespace cafes
{
  namespace geometry
  {
    std::vector<double> linspace(double beg, double end, int num, bool include_end=true)
    {
        assert( beg <= end
                  && "beg must be lower than end"
                  );
        std::vector<double> that;
        that.reserve(num);
        int ns = (include_end)? num - 1: num;
        double step = (end - beg)/ns;

        for(size_t i=0; i<num; ++i)
            that.push_back(beg + i*step);

        return that;
    }
  }
}
#endif