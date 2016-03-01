#ifndef PARTICLE_GEOMETRY_TRIANGULATION_HPP_INCLUDED
#define PARTICLE_GEOMETRY_TRIANGULATION_HPP_INCLUDED

#include <particle/geometry/position.hpp>

#include <vector>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/QhullPoint.h>


namespace cafes
{
  namespace geometry
  {

    std::vector<std::array<int, 3>> triangulation_(std::vector<position<2, double>> points){
      orgQhull::Qhull qhull;
      std::vector<std::array<int, 3>> triangles;

      std::vector<double> coords;
      coords.resize(points.size()*2);
      
      for(std::size_t i=0; i<points.size(); ++i){
          coords[2*i]     = points[i][0];
          coords[2*i + 1] = points[i][1];
      }

      qhull.runQhull("", 2, points.size(), coords.data(), "d Qt Qbb Qz Qc");
      triangles.resize(qhull.facetCount());

      std::size_t ind = 0;
      for(auto &f: qhull.facetList()){
        if (f.isTriCoplanar()){
          orgQhull::QhullVertexSet ps= f.vertices();

          orgQhull::QhullVertexSetIterator i = ps;
          std::size_t ns = 0;
          while(i.hasNext()){
              const orgQhull::QhullVertex v= i.next();
              const orgQhull::QhullPoint p= v.point();
              triangles[ind][ns++] = p.id();
          }
          ind++;
        }
      }
      return triangles;
    }

    std::vector<std::array<int, 3>> triangulation_(std::vector<position<3, double>> points){
      orgQhull::Qhull qhull;
      std::vector<std::array<int, 3>> triangles;

      std::vector<double> coords;
      coords.resize(points.size()*3);
      
      for(std::size_t i=0; i<points.size(); ++i){
          coords[3*i]     = points[i][0];
          coords[3*i + 1] = points[i][1];
          coords[3*i + 2] = points[i][2];
      }

      qhull.runQhull("", 3, points.size(), coords.data(), "i Qt");
      triangles.resize(qhull.facetCount());

      std::size_t ind = 0;
      for(auto &f: qhull.facetList()){
          orgQhull::QhullVertexSet ps= f.vertices();

          orgQhull::QhullVertexSetIterator i = ps;
          std::size_t ns = 0;
          while(i.hasNext()){
              const orgQhull::QhullVertex v= i.next();
              const orgQhull::QhullPoint p= v.point();
              triangles[ind][ns++] = p.id();
          }
          ind++;
      }
      return triangles;
    }

    template<std::size_t Dimensions>
    std::vector<std::array<int, 3>> triangulation(std::vector<position<Dimensions, double>> points){
      return triangulation_(points);
    }
  }
}
#endif