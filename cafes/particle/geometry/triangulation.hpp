// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
// OF SUCH DAMAGE.

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

    std::vector<std::array<int, 3>> triangulation_(std::vector<position<double, 2>> points){
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

    std::vector<std::array<int, 3>> triangulation_(std::vector<position<double, 3>> points){
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
    std::vector<std::array<int, 3>> triangulation(std::vector<position<double, Dimensions>> points){
      return triangulation_(points);
    }
  }
}
#endif