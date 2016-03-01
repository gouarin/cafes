#ifndef PARTICLE_GEOMETRY_BOX_HPP_INCLUDED
#define PARTICLE_GEOMETRY_BOX_HPP_INCLUDED

#include <particle/geometry/position.hpp>
#include <iostream>

namespace cafes
{
  namespace geometry
  {
      template<std::size_t Dimensions, typename T> struct box
      {
        using position_type = position<Dimensions, T>;

        position_type bottom_left, upper_right;

        double length(std::size_t i) const { return std::abs(upper_right[i] - bottom_left[i]); }
        
        double length() const
        {
            auto result=1.;
            for(std::size_t i=0; i<Dimensions; ++i)
                result *= length(i);
            return result;
        }

        friend std::ostream &operator<<( std::ostream &output, 
                                         box const& b )
        { 
          output << "box : ( ";

          for(std::size_t i=0; i<Dimensions; ++i)
            output << b.bottom_left[i] << " ";

          output << "), ( ";

          for(std::size_t i=0; i<Dimensions; ++i)
            output << b.upper_right[i] << " ";

          output << ")";

          return output;            
        }

      };


      template<std::size_t Dimensions, typename T>
      bool intersect(box<Dimensions, T> const& b1, box<Dimensions, T> const& b2)
      {
        for(std::size_t i=0; i<Dimensions; ++i){
          if (b1.bottom_left[i] > b2.upper_right[i]) return false;
          if (b1.upper_right[i] < b2.bottom_left[i]) return false;
        }
        return true;
      }

      template<std::size_t Dimensions, typename T>
      box<Dimensions, T> overlap_box(box<Dimensions, T> const& b1, box<Dimensions, T> const& b2)
      {

          box<Dimensions, T> bout{};
          for(std::size_t i=0; i<Dimensions; ++i){
              bout.bottom_left[i] = std::max( b1.bottom_left[i], b2.bottom_left[i]);
              bout.upper_right[i] = std::min( b1.upper_right[i], b2.upper_right[i]);
          }
          return bout;
      }

      template<std::size_t Dimensions, typename T>
      box<Dimensions, T> box_inside(box<Dimensions, T> const& b1, box<Dimensions, T> const& b2)
      {

          box<Dimensions, T> bout{b1};
          for(std::size_t i=0; i<Dimensions; ++i){
            if (b2.bottom_left[i] > b1.bottom_left[i])
              bout.bottom_left[i] = b2.bottom_left[i];
            if (b2.upper_right[i] < b1.upper_right[i])
              bout.upper_right[i] = b2.upper_right[i];
          }
          return bout;
      }

      template<std::size_t Dimensions, typename T>
      bool point_inside(box<Dimensions, T> const& b, position<Dimensions, T> p){
        for(std::size_t i=0; i<Dimensions; ++i)
          if (p[i] < b.bottom_left[i] || p[i] >= b.upper_right[i]) // don't take the last point (fix for the border domain)
            return false;
        return true;
      }

      template<std::size_t Dimensions, typename T>
      bool check_box_inside(box<Dimensions, T> const& b1, box<Dimensions, T> const& b2)
      // check if a part of box b2 is inside the box b1
      {
          int count = 0;
          for(std::size_t i=0; i<Dimensions; ++i){
            if (b2.bottom_left[i] >= b1.bottom_left[i] && b2.bottom_left[i] <= b1.upper_right[i])
              count++;
            if (b2.upper_right[i] >= b1.bottom_left[i] && b2.upper_right[i] <= b1.upper_right[i])
              count++;
          }
          return (count>0)?true:false;
      }

  }
}

#endif
