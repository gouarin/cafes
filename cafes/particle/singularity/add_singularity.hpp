// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef CAFES_PARTICLE_SINGULARITY_ADD_SINGULARITY_HPP_INCLUDED
#define CAFES_PARTICLE_SINGULARITY_ADD_SINGULARITY_HPP_INCLUDED

#include <particle/geometry/box.hpp>
#include <particle/geometry/position.hpp>
#include <particle/particle.hpp>
#include <particle/singularity/UandPNormal.hpp>
#include <particle/singularity/singularity.hpp>
#include <petsc/vec.hpp>
#include <particle/singularity/truncation.hpp>

#include <cmath>
#include <iostream>
#include <petsc.h>
#include <sstream>

// #include "vtkDoubleArray.h"
// #include "vtkPoints.h"
// #include "vtkPointData.h"
// #include "vtkStructuredGrid.h"
// #include "vtkXMLStructuredGridWriter.h"

namespace cafes
{
    namespace singularity
    {
#undef __FUNCT__
#define __FUNCT__ "computesingularST"
        template<typename Shape>
        PetscErrorCode
        computesingularST(singularity<Shape, 2> sing, particle<Shape> const &p1,
                          particle<Shape> const &p2, petsc::petsc_vec<2> &sol,
                          geometry::box<int, 2> box,
                          std::array<double, 2> const &h)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            const int Dimensions = 2;
            using position_type = geometry::position<double, Dimensions>;
            using position_type_i = geometry::position<int, Dimensions>;

              if (!p1.contains(pts) && !p2.contains(pts))
              {
                auto pos_ref_part = sing.get_pos_in_part_ref(pts);

            for (std::size_t j = box.bottom_left[1]; j < box.upper_right[1];
                 ++j)
            {
                for (std::size_t i = box.bottom_left[0]; i < box.upper_right[0];
                     ++i)
                {
                    position_type_i pts_i = {i, j};
                    auto ielem = fem::get_element(pts_i);

                  auto gradUsing = sing.get_grad_u_sing(pts);
                  auto psing = sing.get_p_sing(pts);
                  //auto chix = sing.get_x_truncation(pts);
                  
                  for (std::size_t je=0; je<bfunc.size(); ++je)
                  {
                    auto u = sol.at(ielem[je]);

                    for (std::size_t d1=0; d1<Dimensions; ++d1)
                    {
                      for (std::size_t d2=0; d2<Dimensions; ++d2)
                        u[d1] -= coef*gradUsing[d1][d2]*bfunc[je][d2];//*chix;
                      u[d1] += coef*psing*bfunc[je][d1];//*chix;
                    }
                  }
                }
              }
              
              // Test with cut function
              // Inside P1
              else if (p1.contains(pts))
              {
                auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                {
                  position_type pts_polar = {sqrt( (pts[0]- p1.center_[0])*(pts[0]- p1.center_[0]) + (pts[1]- p1.center_[1])*(pts[1]- p1.center_[1]) ), atan2(pts[1]- p1.center_[1], pts[0]- p1.center_[0])};
                  position_type pts_border = {p1.center_[0] + p1.shape_factors_[0]*cos(pts_polar[1]), p1.center_[1] + p1.shape_factors_[0]*sin(pts_polar[1])};
                  position_type pts_loc = {is*hs[0], js*hs[1]};
                  auto bfunc = fem::P1_integration_grad(pts_loc, h);

                  auto gradUsing = sing.get_grad_u_sing(pts_border);
                  auto psing = sing.get_p_sing(pts_border);
                  auto chi = cafes::singularity::singTrunc(pts_polar[0]/p1.shape_factors_[0]);
                  //auto chix = sing.get_x_truncation(pts_border);
                  
                  for (std::size_t je=0; je<bfunc.size(); ++je)
                  {
                    auto u = sol.at(ielem[je]);

                    for (std::size_t d1=0; d1<Dimensions; ++d1)
                    {
                      for (std::size_t d2=0; d2<Dimensions; ++d2)
                        u[d1] -= coef*gradUsing[d1][d2]*bfunc[je][d2]*chi;//*chix;
                      u[d1] += coef*psing*bfunc[je][d1]*chi;//*chix;
                    }
                  }
                }
              }
              // Inside P2
              else if (p2.contains(pts))
              {
                auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                {
                  position_type pts_polar = {sqrt( (pts[0]- p2.center_[0])*(pts[0]- p2.center_[0]) + (pts[1]- p2.center_[1])*(pts[1]- p2.center_[1]) ), atan2(pts[1]- p2.center_[1], pts[0]- p2.center_[0])};
                  position_type pts_border = {p2.center_[0] + p2.shape_factors_[0]*cos(pts_polar[1]), p2.center_[1] + p2.shape_factors_[0]*sin(pts_polar[1])};
                  position_type pts_loc = {is*hs[0], js*hs[1]};
                  auto bfunc = fem::P1_integration_grad(pts_loc, h);

                  auto gradUsing = sing.get_grad_u_sing(pts_border);
                  auto psing = sing.get_p_sing(pts_border);
                  auto chi = cafes::singularity::singTrunc(pts_polar[0]/p2.shape_factors_[0]);
                  //auto chix = sing.get_x_truncation(pts_border);
                  
                  for (std::size_t je=0; je<bfunc.size(); ++je)
                  {
                    auto u = sol.at(ielem[je]);

                    for (std::size_t d1=0; d1<Dimensions; ++d1)
                    {
                      for (std::size_t d2=0; d2<Dimensions; ++d2)
                        u[d1] -= coef*gradUsing[d1][d2]*bfunc[je][d2]*chi;//*chix;
                      u[d1] += coef*psing*bfunc[je][d1]*chi;//*chix;
                    }
                }
            }
            PetscFunctionReturn(0);
        }

        // #undef __FUNCT__
        // #define __FUNCT__ "computesingularST_pressure"
        // template<typename Shape>
        // PetscErrorCode computesingularST_pressure(singularity<2, Shape> sing,
        //                                  particle<Shape> const& p1,
        //                                  particle<Shape> const& p2,
        //                                  petsc::petsc_vec<2>& sol,
        //                                  geometry::box<2, int> box,
        //                                  std::array<double, 2> const& h)
        // {
        //   PetscErrorCode ierr;
        //   PetscFunctionBeginUser;

        //   const int Dimensions = 2;
        //   using position_type = geometry::position<Dimensions, double>;
        //   using position_type_i = geometry::position<Dimensions, int>;

        //   std::array<double, Dimensions> hs = {{h[0]/(2*sing.scale),
        //   h[1]/(2*sing.scale)}}; double coef = 1./(4*sing.scale*sing.scale);

        //   for(std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j)
        //   {
        //     for(std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i)
        //     {
        //       position_type_i pts_i = {i, j};
        //       auto ielem = fem::get_element(pts_i);

        //       for(std::size_t js=0; js<2*sing.scale; ++js)
        //       {
        //         for(std::size_t is=0; is<2*sing.scale; ++is)
        //         {
        //           position_type pts = {i*h[0] + is*hs[0], j*h[1] + js*hs[1]};

        //           if (!p1.contains(pts) && !p2.contains(pts))
        //           {
        //             auto pos_ref_part = sing.get_pos_in_part_ref(pts);

        //             if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
        //             {
        //               position_type pts_loc = {is*hs[0], js*hs[1]};

        //               auto bfunc = fem::P1_integration(pts_loc, h);

        //               auto divUsing = sing.get_divu_sing(pts);

        //               for (std::size_t je=0; je<bfunc.size(); ++je)
        //               {
        //                 auto u = sol.at(ielem[je]);

        //                 u[0] -= coef*divUsing*bfunc[je];
        //               }
        //             }
        //           }
        //         }
        //       }
        //     }
        //   }
        //   PetscFunctionReturn(0);
        // }

#undef __FUNCT__
#define __FUNCT__ "computesingularST"
        template<typename Shape>
        PetscErrorCode
        computesingularST(singularity<Shape, 3> sing, particle<Shape> const &p1,
                          particle<Shape> const &p2, petsc::petsc_vec<3> &sol,
                          geometry::box<int, 3> box,
                          std::array<double, 3> const &h)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            const int Dimensions = 3;
            using position_type = geometry::position<double, Dimensions>;
            using position_type_i = geometry::position<int, Dimensions>;

                if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                {
                  position_type pts_loc = {is*hs[0], js*hs[1]};

                  auto bfunc = fem::P1_integration_sing(pts_loc, h);
                  //auto chix = sing.get_x_truncation(pts);

                  auto divUsing = sing.get_divu_sing(pts);
                  //auto chix = sing.get_x_truncation(pts);
                  //auto gradUsing = sing.get_grad_u_sing(pts);
                  
                  for (std::size_t je=0; je<bfunc.size(); ++je)
                  {
                    auto u = sol.at(ielem[je]);
                    u[0] -= coef*divUsing*bfunc[je];//*chix;
                    //u[0] -= coef*(gradUsing[0][0] + gradUsing[Dimensions][Dimensions])*bfunc[je];
                  }
                }
              }
              
              // Test with cut function
              // \bar{u}^{sing}(x,y) = u^{sing}(X,Y) \chi(r)
              // with :
              // X(\theta) = R\cos(\theta)
              // Y(\theta) = R \sin(\theta)
              // \theta(x,y) = atan2(y - y_part, x - x_part)
              // r(x,y) = \sqrt{ (x - x_part)^2 + (y - y_part)^2 }
              
              // Inside P1
              else if (p1.contains(pts))
              {
                auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                {
                  position_type pts_loc = {is*hs[0], js*hs[1]};

                  // (r, \theta)
                  position_type pts_polar = {sqrt( (pts[0]- p1.center_[0])*(pts[0]- p1.center_[0]) + (pts[1]- p1.center_[1])*(pts[1]- p1.center_[1]) ), atan2(pts[1]- p1.center_[1], pts[0]- p1.center_[0])};

                  // (X, Y)
                  position_type pts_border = {p1.center_[0] + p1.shape_factors_[0]*cos(pts_polar[1]), p1.center_[1] + p1.shape_factors_[0]*sin(pts_polar[1])};

                  // Scale for \partial (r, \theta)
                  double scale_theta = ( (pts[0] - p1.center_[0])*(pts[0] - p1.center_[0]) + (pts[1] - p1.center_[1])*(pts[1] - p1.center_[1]) );
                  if (scale_theta < 1e-14) { scale_theta = 1.;}
                  double scale_r = sqrt(scale_theta);

                  // \partial_x (r, \theta)
                  position_type dx_pts_polar = { (pts[0] - p1.center_[0])/scale_r, -1.*(pts[1] - p1.center_[1])/scale_theta};

                  // \partial_y (r, theta)
                  position_type dy_pts_polar = { (pts[1] - p1.center_[1])/scale_r, (pts[0] - p1.center_[0])/scale_theta};

                  // \partial_\theta (X, Y)
                  position_type dtheta_pts_border = { -1.*p1.shape_factors_[0]*sin(pts_polar[1]), p1.shape_factors_[0]*cos(pts_polar[1]) };
            
                  // Finite element function
                  auto bfunc = fem::P1_integration_sing(pts_loc, h);

                  // Get u^{sing}, its partial derivatives and \chi
                  auto gradUsing = sing.get_grad_u_sing(pts_border);
                  auto Using = sing.get_u_sing(pts_border);
                  auto chi = cafes::singularity::singTrunc(pts_polar[0]/p1.shape_factors_[0]);
                  auto drChi = cafes::singularity::dsingTrunc(pts_polar[0]/p1.shape_factors_[0]);
                  
                  for (std::size_t je=0; je<bfunc.size(); ++je)
                  {
                    auto u = sol.at(ielem[je]);
                    u[0] -= coef*(
                      // \partial_x u_x^{sing}(X,Y) * \chi(r)
                      (dx_pts_polar[1]*dtheta_pts_border[0]*gradUsing[0][0] + dx_pts_polar[1]*dtheta_pts_border[1]*gradUsing[1][0])*chi
                      // \partial_X u_y^{sing}(X,Y) * \chi(r)
                      + (dy_pts_polar[1]*dtheta_pts_border[0]*gradUsing[0][1] + dy_pts_polar[1]*dtheta_pts_border[1]*gradUsing[1][1])*chi
                      // \u_x^{sing}(X,Y) * \partial_x \chi(r)
                      + Using[0]*dx_pts_polar[0]*drChi
                      // \u_y^{sing}(X,Y) * \partial_y \chi(r)
                      + Using[1]*dy_pts_polar[0]*drChi
                    )*bfunc[je];
                  }
                }
              } // End Inside P1
              
              // Inside P2
              else if (p2.contains(pts))
              {
                auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                {
                  position_type pts_loc = {is*hs[0], js*hs[1]};

                  // (r, \theta)
                  position_type pts_polar = {sqrt( (pts[0]- p2.center_[0])*(pts[0]- p2.center_[0]) + (pts[1]- p2.center_[1])*(pts[1]- p2.center_[1]) ), atan2(pts[1]- p2.center_[1], pts[0]- p2.center_[0])};

                  // (X, Y)
                  position_type pts_border = {p2.center_[0] + p2.shape_factors_[0]*cos(pts_polar[1]), p2.center_[1] + p2.shape_factors_[0]*sin(pts_polar[1])};

                  // Scale for \partial (r, \theta)
                  double scale_theta = ( (pts[0] - p2.center_[0])*(pts[0] - p2.center_[0]) + (pts[1] - p2.center_[1])*(pts[1] - p2.center_[1]) );
                  if (scale_theta < 1e-14) { scale_theta = 1.;}
                  double scale_r = sqrt(scale_theta);

                  // \partial_x (r, \theta)
                  position_type dx_pts_polar = { (pts[0] - p2.center_[0])/scale_r, -1.*(pts[1] - p2.center_[1])/scale_theta};

                  // \partial_y (r, theta)
                  position_type dy_pts_polar = { (pts[1] - p2.center_[1])/scale_r, (pts[0] - p2.center_[0])/scale_theta};

                  // \partial_\theta (X, Y)
                  position_type dtheta_pts_border = { -1.*p2.shape_factors_[0]*sin(pts_polar[1]), p2.shape_factors_[0]*cos(pts_polar[1]) };
            
                  // Finite element function
                  auto bfunc = fem::P1_integration_sing(pts_loc, h);

                  // Get u^{sing}, its partial derivatives and \chi
                  auto gradUsing = sing.get_grad_u_sing(pts_border);
                  auto Using = sing.get_u_sing(pts_border);
                  auto chi = cafes::singularity::singTrunc(pts_polar[0]/p2.shape_factors_[0]);
                  auto drChi = cafes::singularity::dsingTrunc(pts_polar[0]/p2.shape_factors_[0]);
                  
                  for (std::size_t je=0; je<bfunc.size(); ++je)
                  {
                    auto u = sol.at(ielem[je]);
                    u[0] -= coef*(
                      // \partial_x u_x^{sing}(X,Y) * \chi(r)
                      (dx_pts_polar[1]*dtheta_pts_border[0]*gradUsing[0][0] + dx_pts_polar[1]*dtheta_pts_border[1]*gradUsing[1][0])*chi
                      // \partial_X u_y^{sing}(X,Y) * \chi(r)
                      + (dy_pts_polar[1]*dtheta_pts_border[0]*gradUsing[0][1] + dy_pts_polar[1]*dtheta_pts_border[1]*gradUsing[1][1])*chi
                      // \u_x^{sing}(X,Y) * \partial_x \chi(r)
                      + Using[0]*dx_pts_polar[0]*drChi
                      // \u_y^{sing}(X,Y) * \partial_y \chi(r)
                      + Using[1]*dy_pts_polar[0]*drChi
                    )*bfunc[je];
                  }
                }
              } // End Inside P2
            }
          }
        }
      }
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "computesingularST"
    template<typename Shape>
    PetscErrorCode computesingularST(singularity<Shape, 3> sing, 
                                     particle<Shape> const& p1, particle<Shape> const& p2,
                                     petsc::petsc_vec<3>& sol,
                                     geometry::box<int, 3> box,
                                     std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 3;
      using position_type = geometry::position<double, Dimensions>;
      using position_type_i = geometry::position<int, Dimensions>;

      std::array<double, Dimensions> hs = {{ h[0]/sing.scale, 
                                             h[1]/sing.scale,
                                             h[2]/sing.scale }};

      double coef = 1./(sing.scale*sing.scale*sing.scale);

      for(std::size_t k=box.bottom_left[2]; k<box.upper_right[2]; ++k)
      {
        for(std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j)
        {
          for(std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i)
          {
            position_type_i pts_i = {i, j, k};
            auto ielem = fem::get_element(pts_i);

            for(std::size_t ks=0; ks<sing.scale; ++ks)
            {
                for (std::size_t j = box.bottom_left[1]; j < box.upper_right[1];
                     ++j)
                {
                    for (std::size_t i = box.bottom_left[0];
                         i < box.upper_right[0]; ++i)
                    {
                        position_type_i pts_i = {i, j, k};
                        auto ielem = fem::get_element(pts_i);

                        for (std::size_t ks = 0; ks < sing.scale; ++ks)
                        {
                            for (std::size_t js = 0; js < sing.scale; ++js)
                            {
                                for (std::size_t is = 0; is < sing.scale; ++is)
                                {
                                    position_type pts = {i * h[0] + is * hs[0],
                                                         j * h[1] + js * hs[1],
                                                         k * h[2] + ks * hs[2]};

                                    if (!p1.contains(pts) && !p2.contains(pts))
                                    {
                                        position_type pts_loc = {
                                            is * hs[0], js * hs[1], ks * hs[2]};

                                        auto bfunc = fem::P1_integration_grad(
                                            pts_loc, h);

                                        auto gradUsing =
                                            sing.get_grad_u_sing(pts);
                                        auto psing = sing.get_p_sing(pts);

                                        for (std::size_t je = 0;
                                             je < bfunc.size(); ++je)
                                        {
                                            auto u = sol.at(ielem[je]);

                                            for (std::size_t d1 = 0;
                                                 d1 < Dimensions; ++d1)
                                            {
                                                for (std::size_t d2 = 0;
                                                     d2 < Dimensions; ++d2)
                                                    u[d1] -= coef *
                                                             gradUsing[d1][d2] *
                                                             bfunc[je][d2];
                                                u[d1] += coef * psing *
                                                         bfunc[je][d1];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "computesingularBC"
        template<typename Shape>
        PetscErrorCode computesingularBC(singularity<Shape, 2> &sing,
                                         petsc::petsc_vec<2> &sol,
                                         std::array<double, 2> const &h,
                                         geometry::box<int, 2> &box)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            using position_type = geometry::position<double, 2>;

            double theta = std::asin(sing.cutoff_dist_ * sing.H1_);
            int N = sing.cutoff_dist_ / h[0] * sing.scale;

            // first particle
            for (int i = -N; i < N; ++i)
            {
                double cos_t = std::cos(i * theta / N);
                double sin_t = std::sin(i * theta / N);

                position_type pos_ref{(cos_t - 1.) / sing.H1_,
                                      sin_t / sing.H1_};
                auto pos = sing.get_pos_from_part_ref(pos_ref);
                auto pos_i = static_cast<geometry::position<int, 2>>(pos / h);
                // std::cout<< "pos_ref"<<pos_ref<<"pos"<<pos<<"\n";
                // std::cout<< "pos_i"<<pos_i<<"\n";
                if (geometry::point_inside(box, pos_i))
                {
                    auto Using_ref = sing.get_u_sing_ref(pos_ref);
                    std::array<double, 2> Using{};

                    for (std::size_t d1 = 0; d1 < 2; ++d1)
                    {
                        Using[d1] = 0.;
                        for (std::size_t d2 = 0; d2 < 2; ++d2)
                            Using[d1] += Using_ref[d2] * sing.base_[d2][d1];
                    }

                    position_type spts = pos - pos_i * h;
                    // std::cout<< "spts"<<spts<<"\n";

                    auto bfunc = fem::P1_integration(spts, h);
                    auto ielem = fem::get_element(pos_i);

                    for (std::size_t je = 0; je < bfunc.size(); ++je)
                    {
                        auto u = sol.at(ielem[je]);
                        for (std::size_t d = 0; d < 2; ++d)
                            u[d] += Using[d] * bfunc[je] * theta / N / sing.H1_;
                    }
                }
            }

            theta = std::asin(sing.cutoff_dist_ * sing.H2_);

            // second particle
            for (int i = -N; i < N; ++i)
            {
                double cos_t = std::cos(i * theta / N);
                double sin_t = std::sin(i * theta / N);

                position_type pos_ref{sing.contact_length_ +
                                          (1. - cos_t) / sing.H2_,
                                      sin_t / sing.H2_};
                auto pos = sing.get_pos_from_part_ref(pos_ref);
                auto pos_i = static_cast<geometry::position<int, 2>>(pos / h);

                if (geometry::point_inside(box, pos_i))
                {
                    auto Using_ref = sing.get_u_sing_ref(pos_ref);
                    std::array<double, 2> Using{};

                    for (std::size_t d1 = 0; d1 < 2; ++d1)
                    {
                        Using[d1] = 0.;
                        for (std::size_t d2 = 0; d2 < 2; ++d2)
                            Using[d1] += Using_ref[d2] * sing.base_[d2][d1];
                    }

                    position_type spts = pos - pos_i * h;

                    auto bfunc = fem::P1_integration(spts, h);
                    auto ielem = fem::get_element(pos_i);

                    for (std::size_t je = 0; je < bfunc.size(); ++je)
                    {
                        auto u = sol.at(ielem[je]);
                        for (std::size_t d = 0; d < 2; ++d)
                            u[d] += Using[d] * bfunc[je] * theta / N / sing.H2_;
                    }
                }
            }
            PetscFunctionReturn(0);
        }

      std::cout<<"add singularity in fluid...\n";

      auto union_box_func = geometry::union_box<int, Dimensions>;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.rhs, 0, false);
      ierr = sol.fill(0.);CHKERRQ(ierr);

      auto boxp = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 1);
      std::array<double, Dimensions> hp;
      for(std::size_t d=0; d<Dimensions; ++d)
        hp[d] = 2*ctx.problem.ctx->h[d];

      auto solp = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.rhs, 1, false);
      ierr = solp.fill(0.);CHKERRQ(ierr);

      //Loop on particles couples
      for (std::size_t ipart=0; ipart<ctx.particles.size()-1; ++ipart)
      {
        auto p1 = ctx.particles[ipart];
        for (std::size_t jpart=ipart+1; jpart<ctx.particles.size(); ++jpart)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

          using shape_type = typename decltype(p1)::shape_type;

          // Velocity
          auto sing = singularity<shape_type, Dimensions>(p1, p2, h[0]);
          if (sing.is_singularity_)
          {
            auto pbox = sing.get_box(h);
            if (geometry::intersect(box, pbox))
            {
                auto pos = ctx.surf_points[ipart_1][isurf].second +
                           ctx.surf_points[ipart_1][isurf].first * h;
                if (geometry::point_inside(box, pos))
                {
                    auto pos_ref_part = sing.get_pos_in_part_ref(pos);
                    if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                    {
                        auto Using = sing.get_u_sing(pos);
                        for (std::size_t d = 0; d < Dimensions; ++d)
                            g[ipart_1][isurf][d] += Using[d];
                    }
                }
            }
          }

          // Pressure
          //auto singp = singularity<shape_type, Dimensions>(p1, p2, hp[0]);
          if (sing.is_singularity_)
          {
            auto pboxp = sing.get_box(hp);
            if (geometry::intersect(boxp, pboxp))
            {
              auto new_box = geometry::box_inside(boxp, pboxp);
              ierr = computesingularST_pressure(sing, p1, p2, solp, new_box, hp);CHKERRQ(ierr);
            }
          }
        }
      }

      ierr = sol.local_to_global(ADD_VALUES);CHKERRQ(ierr);
      ierr = solp.local_to_global(ADD_VALUES);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "add_singularity_to_last_sol"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode add_singularity_to_last_sol(Ctx& ctx, Vec vsol)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto union_box_func = geometry::union_box<int, Dimensions>;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, vsol, 0, false);

      ierr = sol.global_to_local(INSERT_VALUES);CHKERRQ(ierr);

      //Loop on particles couples
      for (std::size_t ipart=0; ipart<ctx.particles.size()-1; ++ipart)
      {
        auto p1 = ctx.particles[ipart];
        for (std::size_t jpart=ipart+1; jpart<ctx.particles.size(); ++jpart)
        {
          auto p2 = ctx.particles[jpart];

          using shape_type = typename decltype(p1)::shape_type;
          auto sing = singularity<shape_type, Dimensions>(p1, p2, h[0]);

          if (sing.is_singularity_)
          {
            auto pbox = sing.get_box(h);
            if (geometry::intersect(box, pbox))
            {
                auto pos = ctx.surf_points[ipart_2][isurf].second +
                           ctx.surf_points[ipart_2][isurf].first * h;
                if (geometry::point_inside(box, pos))
                {
                    auto pos_ref_part = sing.get_pos_in_part_ref(pos);
                    if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                    {
                        auto Using = sing.get_u_sing(pos);
                        for (std::size_t d = 0; d < Dimensions; ++d)
                            g[ipart_2][isurf][d] += Using[d];
                    }
                }
            }

            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "addsingularity"
        template<typename Shape>
        PetscErrorCode
        addsingularity(singularity<Shape, 2> sing, particle<Shape> const &p1,
                       particle<Shape> const &p2, petsc::petsc_vec<2> &sol,
                       geometry::box<int, 2> box,
                       std::array<double, 2> const &h)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            const int Dimensions = 2;
            using position_type = geometry::position<double, Dimensions>;
            using position_type_i = geometry::position<int, Dimensions>;

    #undef __FUNCT__
    #define __FUNCT__ "add_singularity_to_ureg"
    template<std::size_t Dimensions, class Particles>
    PetscErrorCode add_singularity_to_ureg(DM dm, std::array<double, Dimensions> h, Vec vsol, const Particles& particles)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;
      using position_type = geometry::position<double, Dimensions>;
      using position_type_i = geometry::position<int, Dimensions>;
      std::array<double, Dimensions> hp = {2*h[0], 2*h[1]};

      auto union_box_func = geometry::union_box<int, Dimensions>;

      auto box = fem::get_DM_bounds<Dimensions>(dm, 0);
      auto box_p = fem::get_DM_bounds<Dimensions>(dm, 1);

      auto sol_u = petsc::petsc_vec<Dimensions>(dm, vsol, 0, false);
      auto sol_p = petsc::petsc_vec<Dimensions>(dm, vsol, 1, false);

      //Loop on particles couples
      for (std::size_t ipart=0; ipart<particles.size()-1; ++ipart)
      {
        auto p1 = particles[ipart];
        for (std::size_t jpart=ipart+1; jpart<particles.size(); ++jpart)
        {
          auto p2 = particles[jpart];

          using shape_type = typename decltype(p1)::shape_type;
          auto sing = singularity<shape_type, Dimensions>(p1, p2, h[0]);
          auto sing_p = singularity<shape_type, Dimensions>(p1, p2, hp[0]);

          if (sing.is_singularity_)
          {
            auto pbox = sing.get_box(h);
            auto pbox_p = sing_p.get_box(hp);
            
            // VELOCITY
            if (geometry::intersect(box, pbox))
            {
              auto new_box = geometry::box_inside(box, pbox);
              for (std::size_t j=new_box.bottom_left[1]; j<new_box.upper_right[1]; ++j)
              {
                for (std::size_t i=new_box.bottom_left[0]; i<new_box.upper_right[0]; ++i)
                {
                  position_type_i pts_i = {i, j};
                  position_type pts = {i*h[0], j*h[1]};
                  if (!p1.contains(pts) and !p2.contains(pts))
                  {
                    auto u_sing = sing.get_u_sing(pts);
                    auto u = sol_u.at_g(pts_i);
                    //std::cout << "u " << u[0] << " " << u[1] << "\n";
                    for (std::size_t d=0; d<Dimensions; ++d)
                    {
                      u[d] += u_sing[d];
                    }
                    //std::cout << "u new" << u[0] << " " << u[1] << "\n";
                  }
                  // // Inside P1
                  // else if (p1.contains(pts))
                  // {
                  //   auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                  //   if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                  //   {
                  //     position_type pts_polar = {sqrt( (pts[0]- p1.center_[0])*(pts[0]- p1.center_[0]) + (pts[1]- p1.center_[1])*(pts[1]- p1.center_[1]) ), atan2(pts[1]- p1.center_[1], pts[0]- p1.center_[0])};
                  //     position_type pts_border = {p1.center_[0] + p1.shape_factors_[0]*cos(pts_polar[1]), p1.center_[1] + p1.shape_factors_[0]*sin(pts_polar[1])};

                  //     auto Using = sing.get_u_sing(pts_border);
                  //     auto chi = cafes::singularity::alphaTrunc(pts_polar[0]/p1.shape_factors_[0]);
                  //     auto u = sol_u.at_g(pts_i);
                  //     for (std::size_t d=0; d<Dimensions; ++d)
                  //       u[d] += Using[d]*chi;
                  //   }
                  // }
                  // // Inside P2
                  // else if (p2.contains(pts))
                  // {
                  //   auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                  //   if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                  //   {
                  //     position_type pts_polar = {sqrt( (pts[0]- p2.center_[0])*(pts[0]- p2.center_[0]) + (pts[1]- p2.center_[1])*(pts[1]- p2.center_[1]) ), atan2(pts[1]- p2.center_[1], pts[0]- p2.center_[0])};
                  //     position_type pts_border = {p2.center_[0] + p2.shape_factors_[0]*cos(pts_polar[1]), p2.center_[1] + p2.shape_factors_[0]*sin(pts_polar[1])};

                  //     auto Using = sing.get_u_sing(pts_border);
                  //     auto chi = cafes::singularity::alphaTrunc(pts_polar[0]/p2.shape_factors_[0]);
                  //     auto u = sol_u.at_g(pts_i);
                  //     for (std::size_t d=0; d<Dimensions; ++d)
                  //       u[d] += Using[d]*chi;
                  //   }
                  // }
                  // else if (p1.contains(pts))
                  // {
                  //   auto u = sol_u.at_g(pts_i);
                  //   u[0] += p1.velocity_[0];
                  //   u[1] += p1.velocity_[1];
                  // }
                }
              }
            }

            // PRESSURE
            if (geometry::intersect(box_p, pbox_p))
            {
              auto new_box = geometry::box_inside(box_p, pbox_p);
              for (std::size_t j=new_box.bottom_left[1]; j<new_box.upper_right[1]; ++j)
              {
                for (std::size_t i=new_box.bottom_left[0]; i<new_box.upper_right[0]; ++i)
                {
                  position_type_i pts_i = {i, j};
                  position_type pts = {i*hp[0], j*hp[1]};
                  if (!p1.contains(pts) and !p2.contains(pts))
                  {
                    auto p_sing = sing_p.get_p_sing(pts);
                    auto p = sol_p.at_g(pts_i);
                    p[0] += p_sing;
                  }
                }
              }
            }

          }
        }
      }
      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "compute_singular_forces"
    template<typename Shape>
    auto compute_singular_forces(singularity<Shape, 2> sing, std::size_t N)
    {
      using position_type = geometry::position<double, 2>;

            for (std::size_t j = box.bottom_left[1]; j < box.upper_right[1];
                 ++j)
            {
                for (std::size_t i = box.bottom_left[0]; i < box.upper_right[0];
                     ++i)
                {
                    position_type_i pts_i = {i, j};
                    auto ielem = fem::get_element(pts_i);

                    for (std::size_t js = 0; js < sing.scale; ++js)
                    {
                        for (std::size_t is = 0; is < sing.scale; ++is)
                        {
                            position_type pts = {i * h[0] + is * hs[0],
                                                 j * h[1] + js * hs[1]};

                            if (!p1.contains(pts) && !p2.contains(pts))
                            {
                                auto pos_ref_part =
                                    sing.get_pos_in_part_ref(pts);

                                if (std::abs(pos_ref_part[1]) <=
                                    sing.cutoff_dist_)
                                {
                                    position_type pts_loc = {is * hs[0],
                                                             js * hs[1]};

                                    auto bfunc =
                                        fem::P1_integration(pts_loc, h);

                                    auto Using = sing.get_u_sing(pts);

                                    for (std::size_t je = 0; je < bfunc.size();
                                         ++je)
                                    {
                                        auto u = sol.at(ielem[je]);

                                        for (std::size_t d = 0; d < Dimensions;
                                             ++d)
                                            u[d] += coef * Using[d] * bfunc[je];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // const int Dimensions = 2;
            // using position_type = geometry::position<Dimensions, double>;
            // using position_type_i = geometry::position<Dimensions, int>;

            // for(std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j)
            // {
            //   for(std::size_t i=box.bottom_left[0]; i<box.upper_right[0];
            //   ++i)
            //   {
            //     std::array<int, 2> pts_i = {{static_cast<int>(i),
            //     static_cast<int>(j)}}; position_type pts{i*h[0], j*h[1]};

            //     if (!p1.contains(pts) && !p2.contains(pts))
            //     {
            //       auto pos_ref_part = sing.get_pos_in_part_ref(pts);

            //       if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
            //       {
            //         auto Using = sing.get_u_sing(pts);
            //         auto u = sol.at(pts_i);
            //         for (std::size_t d=0; d<Dimensions; ++d)
            //           u[d] += Using[d];
            //       }
            //     }
            //   }
            // }
            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "add_singularity_in_fluid"
        template<std::size_t Dimensions, typename Ctx>
        PetscErrorCode add_singularity_in_fluid(Ctx &ctx)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            std::cout << "add singularity in fluid...\n";

            auto union_box_func = geometry::union_box<int, Dimensions>;

            auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
            auto &h = ctx.problem.ctx->h;

            auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm,
                                                    ctx.problem.rhs, 0, false);
            ierr = sol.fill(0.);
            CHKERRQ(ierr);

            // auto boxp = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm,
            // 1); std::array<double, Dimensions> hp; for(std::size_t d=0;
            // d<Dimensions; ++d)
            //   hp[d] = 2*ctx.problem.ctx->h[d];

            // auto solp = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm,
            // ctx.problem.rhs, 1, false); ierr = solp.fill(0.);CHKERRQ(ierr);

            // Loop on particles couples
            for (std::size_t ipart = 0; ipart < ctx.particles.size() - 1;
                 ++ipart)
            {
                auto p1 = ctx.particles[ipart];
                for (std::size_t jpart = ipart + 1;
                     jpart < ctx.particles.size(); ++jpart)
                {
                    auto p2 = ctx.particles[jpart];

                    using shape_type = typename decltype(p1)::shape_type;
                    auto sing =
                        singularity<shape_type, Dimensions>(p1, p2, h[0]);

                    if (sing.is_singularity_)
                    {
                        auto pbox = sing.get_box(h);
                        if (geometry::intersect(box, pbox))
                        {
                            auto new_box = geometry::box_inside(box, pbox);
                            ierr = computesingularST(sing, p1, p2, sol, new_box,
                                                     h);
                            CHKERRQ(ierr);
                        }

                        // auto pboxp = sing.get_box(hp);
                        // if (geometry::intersect(boxp, pboxp))
                        // {
                        //   auto new_box = geometry::box_inside(boxp, pboxp);
                        //   ierr = computesingularST_pressure(sing, p1, p2,
                        //   solp, new_box, hp);CHKERRQ(ierr);
                        // }
                    }
                }
            }

            ierr = sol.local_to_global(ADD_VALUES);
            CHKERRQ(ierr);
            // ierr = solp.local_to_global(ADD_VALUES);CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "add_singularity_to_last_sol"
        template<std::size_t Dimensions, typename Ctx>
        PetscErrorCode add_singularity_to_last_sol(Ctx &ctx, Vec vsol)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            auto union_box_func = geometry::union_box<int, Dimensions>;

            auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
            auto &h = ctx.problem.ctx->h;

            auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, vsol,
                                                    0, false);

            ierr = sol.global_to_local(INSERT_VALUES);
            CHKERRQ(ierr);

            // Loop on particles couples
            for (std::size_t ipart = 0; ipart < ctx.particles.size() - 1;
                 ++ipart)
            {
                auto p1 = ctx.particles[ipart];
                for (std::size_t jpart = ipart + 1;
                     jpart < ctx.particles.size(); ++jpart)
                {
                    auto p2 = ctx.particles[jpart];

                    using shape_type = typename decltype(p1)::shape_type;
                    auto sing =
                        singularity<shape_type, Dimensions>(p1, p2, h[0]);

                    if (sing.is_singularity_)
                    {
                        auto pbox = sing.get_box(h);
                        if (geometry::intersect(box, pbox))
                        {
                            auto new_box = geometry::box_inside(box, pbox);
                            ierr =
                                addsingularity(sing, p1, p2, sol, new_box, h);
                            CHKERRQ(ierr);
                        }
                    }
                }
            }

            ierr = sol.local_to_global(INSERT_VALUES);
            CHKERRQ(ierr);

            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "compute_singular_forces"
        template<typename Shape>
        auto compute_singular_forces(singularity<Shape, 2> sing, std::size_t N)
        {
            using position_type = geometry::position<double, 2>;

            double mu = 1.;
            double theta = std::asin(sing.cutoff_dist_ * sing.H1_);

            physics::force<2> force{};

            for (std::size_t i = -N; i < N; ++i)
            {
                double cos_t = std::cos(i * theta / N);
                double sin_t = std::sin(i * theta / N);
                std::array<std::array<double, 2>, 2> sigma;

                position_type pos{(cos_t - 1.) / sing.H1_, sin_t / sing.H1_};
                position_type normal{cos_t, sin_t};

                auto gradUsing = sing.get_grad_u_sing_ref(pos);
                auto psing = sing.get_p_sing(pos);

                sigma[0][0] = 2 * mu * gradUsing[0][0] - psing;
                sigma[0][1] = mu * (gradUsing[0][1] + gradUsing[1][0]);
                sigma[1][0] = sigma[0][1];
                sigma[1][1] = 2 * mu * gradUsing[1][1] - psing;

                for (std::size_t d1 = 0; d1 < 2; ++d1)
                    for (std::size_t d2 = 0; d2 < 2; ++d2)
                        force[d1] +=
                            theta / N / sing.H1_ * sigma[d1][d2] * normal[d2];
            }
            return force;
        }

#undef __FUNCT__
#define __FUNCT__ "compute_singular_forces_on_part1"
        template<typename Shape>
        auto compute_singular_forces_on_part1(singularity<Shape, 2> sing,
                                              std::size_t N)
        {
            using position_type = geometry::position<double, 2>;

            double mu = 1.;
            double theta = std::asin(sing.cutoff_dist_ * sing.H1_);

            physics::force<2> force{};

            for (std::size_t i = -N; i < N; ++i)
            {
                double cos_t = std::cos(i * theta / N);
                double sin_t = std::sin(i * theta / N);
                std::array<std::array<double, 2>, 2> sigma;

                position_type pos{(cos_t - 1.) / sing.H1_, sin_t / sing.H1_};
                position_type normal{cos_t, sin_t};

                auto gradUsing = sing.get_grad_u_sing_ref(pos);
                auto psing = sing.get_p_sing(pos);

                sigma[0][0] = 2 * mu * gradUsing[0][0] - psing;
                sigma[0][1] = mu * (gradUsing[0][1] + gradUsing[1][0]);
                sigma[1][0] = sigma[0][1];
                sigma[1][1] = 2 * mu * gradUsing[1][1] - psing;

                for (std::size_t d1 = 0; d1 < 2; ++d1)
                    for (std::size_t d2 = 0; d2 < 2; ++d2)
                        force[d1] +=
                            theta / N / sing.H1_ * sigma[d1][d2] * normal[d2];
            }
            return force;
        }

#undef __FUNCT__
#define __FUNCT__ "compute_singular_forces_on_part2"
        template<typename Shape>
        auto compute_singular_forces_on_part2(singularity<Shape, 2> sing,
                                              std::size_t N)
        {
            using position_type = geometry::position<double, 2>;

            double mu = 1.;
            double theta = std::asin(sing.cutoff_dist_ * sing.H1_);

            physics::force<2> force{};

            for (std::size_t i = -N; i < N; ++i)
            {
                double cos_t = std::cos(i * theta / N);
                double sin_t = std::sin(i * theta / N);
                std::array<std::array<double, 2>, 2> sigma;

                position_type pos{(cos_t - 1.) / sing.H1_, sin_t / sing.H1_};
                position_type normal{cos_t, sin_t};

                auto gradUsing = sing.get_grad_u_sing_ref(pos);
                auto psing = sing.get_p_sing(pos);

                sigma[0][0] = 2 * mu * gradUsing[0][0] - psing;
                sigma[0][1] = mu * (gradUsing[0][1] + gradUsing[1][0]);
                sigma[1][0] = sigma[0][1];
                sigma[1][1] = 2 * mu * gradUsing[1][1] - psing;

                for (std::size_t d1 = 0; d1 < 2; ++d1)
                    for (std::size_t d2 = 0; d2 < 2; ++d2)
                        force[d1] +=
                            theta / N / sing.H1_ * sigma[d1][d2] * normal[d2];
            }
            return force;
        }

#undef __FUNCT__
#define __FUNCT__ "compute_singular_forces_on_part1"
        template<typename Shape>
        auto compute_singular_forces_on_part1(singularity<Shape, 3> sing,
                                              std::size_t N)
        {
            using position_type = geometry::position<double, 3>;

            double mu = 1.;
            double theta = std::asin(sing.cutoff_dist_ * sing.H1_);
            double A = 2 * M_PI * theta / N / N / sing.H1_ / sing.H1_;

            physics::force<3> force{};
            physics::force<3> force_ref{};

            for (std::size_t i = 1; i <= N; ++i)
            {
                double cosi_t = std::cos(i * theta / N);
                double sini_t = std::sin(i * theta / N);
                double coef = A * sini_t;
                for (std::size_t j = 1; j <= N; ++j)
                {
                    double cosj_t = std::cos(2 * M_PI * j / N);
                    double sinj_t = std::sin(2 * M_PI * j / N);
                    std::array<std::array<double, 3>, 3> sigma;

                    position_type pos{sini_t * cosj_t / sing.H1_,
                                      sini_t * sinj_t / sing.H1_,
                                      (cosi_t - 1.) / sing.H1_};
                    position_type normal{sini_t * cosj_t, sini_t * sinj_t,
                                         cosi_t};

                    auto gradUsing = sing.get_grad_u_sing_ref(pos);
                    auto psing = sing.get_p_sing_ref(pos);

                    sigma[0][0] = 2 * mu * gradUsing[0][0] - psing;
                    sigma[0][1] = 2 * mu * gradUsing[0][1];
                    sigma[0][2] = mu * (gradUsing[0][2] + gradUsing[2][0]);

                    sigma[1][0] = sigma[0][1];
                    sigma[1][1] = 2 * mu * gradUsing[1][1] - psing;
                    sigma[1][2] = mu * (gradUsing[1][2] + gradUsing[2][1]);

                    sigma[2][0] = sigma[0][2];
                    sigma[2][1] = sigma[1][2];
                    sigma[2][2] = 2 * mu * gradUsing[2][2] - psing;

                    for (std::size_t d1 = 0; d1 < 3; ++d1)
                        for (std::size_t d2 = 0; d2 < 3; ++d2)
                            force_ref[d1] += coef * sigma[d1][d2] * normal[d2];

                    for (std::size_t d1 = 0; d1 < 3; ++d1)
                        for (std::size_t d2 = 0; d2 < 3; ++d2)
                            force[d1] += sing.base_[d2][d1] * force_ref[d2];
                }
            }
            return force;
        }

#undef __FUNCT__
#define __FUNCT__ "compute_singular_forces_on_part2"
        template<typename Shape>
        auto compute_singular_forces_on_part2(singularity<Shape, 3> sing,
                                              std::size_t N)
        {
            using position_type = geometry::position<double, 3>;

            double mu = 1.;
            double theta = std::asin(sing.cutoff_dist_ * sing.H2_);
            double A = 2 * M_PI * theta / N / N / sing.H2_ / sing.H2_;

            physics::force<3> force{};
            physics::force<3> force_ref{};

            for (std::size_t i = 1; i <= N; ++i)
            {
                double cosi_t = std::cos(i * theta / N);
                double sini_t = std::sin(i * theta / N);
                double coef = A * sini_t;
                for (std::size_t j = 1; j <= N; ++j)
                {
                    double cosj_t = std::cos(2 * M_PI * j / N);
                    double sinj_t = std::sin(2 * M_PI * j / N);
                    std::array<std::array<double, 3>, 3> sigma;

                    position_type pos{
                        sini_t * cosj_t / sing.H2_, sini_t * sinj_t / sing.H2_,
                        sing.contact_length_ + (1. - cosi_t) / sing.H2_};
                    position_type normal{sini_t * cosj_t, sini_t * sinj_t,
                                         -cosi_t};

                    auto gradUsing = sing.get_grad_u_sing_ref(pos);
                    auto psing = sing.get_p_sing_ref(pos);

                    sigma[0][0] = 2 * mu * gradUsing[0][0] - psing;
                    sigma[0][1] = 2 * mu * gradUsing[0][1];
                    sigma[0][2] = mu * (gradUsing[0][2] + gradUsing[2][0]);

                    sigma[1][0] = sigma[0][1];
                    sigma[1][1] = 2 * mu * gradUsing[1][1] - psing;
                    sigma[1][2] = mu * (gradUsing[1][2] + gradUsing[2][1]);

                    sigma[2][0] = sigma[0][2];
                    sigma[2][1] = sigma[1][2];
                    sigma[2][2] = 2 * mu * gradUsing[2][2] - psing;

                    for (std::size_t d1 = 0; d1 < 3; ++d1)
                        for (std::size_t d2 = 0; d2 < 3; ++d2)
                            force_ref[d1] += coef * sigma[d1][d2] * normal[d2];

                    for (std::size_t d1 = 0; d1 < 3; ++d1)
                        for (std::size_t d2 = 0; d2 < 3; ++d2)
                            force[d1] += sing.base_[d2][d1] * force_ref[d2];
                }
            }
            return force;
        }

#undef __FUNCT__
#define __FUNCT__ "compute_singular_forces"
        template<std::size_t Dimensions, typename Ctx>
        PetscErrorCode compute_singular_forces(Ctx &ctx, std::size_t N = 100)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            auto &h = ctx.problem.ctx->h;

            for (std::size_t ipart = 0; ipart < ctx.particles.size(); ++ipart)
                ctx.particles[ipart].force_.fill(0.);

            // Loop on particles couples
            for (std::size_t ipart = 0; ipart < ctx.particles.size() - 1;
                 ++ipart)
            {
                auto p1 = ctx.particles[ipart];
                for (std::size_t jpart = ipart + 1;
                     jpart < ctx.particles.size(); ++jpart)
                {
                    auto p2 = ctx.particles[jpart];

      vtkDoubleArray* pressure_sing = vtkDoubleArray::New();
      pressure_sing->SetName("pressure_sing");
      // pressure_sing[0] = vtkDoubleArray::New(); 
      // pressure_sing[0]->SetName("pressure_sing");
      // pressure_sing[1] = vtkDoubleArray::New(); 
      // pressure_sing[1]->SetName("div_velocity_sing");

      vtkDoubleArray* div_velocity_sing = vtkDoubleArray::New();
      div_velocity_sing->SetName("div_velocity_sing");

                    if (sing.is_singularity_)
                    {
                        ctx.particles[ipart].force_ -=
                            compute_singular_forces_on_part1(sing, N);
                        ctx.particles[jpart].force_ -=
                            compute_singular_forces_on_part2(sing, N);
                    }
                }
            }
            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "compute_singular_forces"
        template<std::size_t Dimensions, typename Shape>
        PetscErrorCode
        compute_singular_forces(std::vector<particle<Shape>> const &particles,
                                std::vector<physics::force<Dimensions>> &forces,
                                std::array<double, Dimensions> const &h,
                                std::size_t N = 100)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            for (std::size_t ipart = 0; ipart < particles.size(); ++ipart)
                forces[ipart].fill(0.);

            // Loop on particles couples
            for (std::size_t ipart = 0; ipart < particles.size() - 1; ++ipart)
            {
                auto p1 = particles[ipart];
                for (std::size_t jpart = ipart + 1; jpart < particles.size();
                     ++jpart)
                {
                    auto p2 = particles[jpart];

                    // using shape_type = typename decltype(p1)::shape_type;
                    auto sing = singularity<Shape, Dimensions>(p1, p2, h[0]);

                if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                {
                  add_sing = true;
                  position_type pts_loc = {is*hs[0], js*hs[1]};
                          
                  auto Using = sing.get_u_sing(pts);
                  auto gradUsing = sing.get_grad_u_sing(pts);
                  auto psing = sing.get_p_sing(pts);
                  auto divUsing =  sing.get_divu_sing(pts);
                  // Add points to vtk + singular value to vtk
                  velocity_sing->InsertNextTuple3(Using[0], Using[1], 0.);
                  gradx_velocity_sing->InsertNextTuple3(gradUsing[0][0], gradUsing[1][0], 0.);
                  grady_velocity_sing->InsertNextTuple3(gradUsing[0][1], gradUsing[1][1], 0.);
                  //pressure_sing->InsertNextValue(divUsing);
                  pressure_sing->InsertNextValue(psing);
                  //pressure_sing[1]->InsertNextValue(divUsing);
                }
              }

              if (!add_sing)
              {
                velocity_sing->InsertNextTuple3(0., 0., 0.);
                gradx_velocity_sing->InsertNextTuple3(0., 0., 0.);
                grady_velocity_sing->InsertNextTuple3(0., 0., 0.);
                pressure_sing->InsertNextValue(0.);
                //pressure_sing[1]->InsertNextValue(0.);
              }
            }
            PetscFunctionReturn(0);
        }
      }

      singDataSet->SetPoints(singPoints);
      singDataSet->GetPointData()->AddArray(velocity_sing);
      singDataSet->GetPointData()->AddArray(gradx_velocity_sing);
      singDataSet->GetPointData()->AddArray(grady_velocity_sing);
      singDataSet->GetPointData()->SetScalars(pressure_sing);
      //singDataSet->GetPointData()->SetScalars(div_velocity_sing);

      vtkXMLStructuredGridWriter* singDataWriter = vtkXMLStructuredGridWriter::New();
      std::stringstream oc;

      oc << path << "/" << filename << "_sing_" << rank << "_" << ".vts";//a changer
      singDataWriter->SetFileName(oc.str().data());
      //dataWriter->SetDataModeToAscii();
#if VTK_MAJOR_VERSION <= 5
      singDataWriter->SetInput(singDataSet);
#else
      singDataWriter->SetInputData(singDataSet);
#endif
      singDataWriter->Write();

      singPoints->Delete();
      velocity_sing->Delete();
      gradx_velocity_sing->Delete();
      grady_velocity_sing->Delete();
      singDataSet->Delete();
      singDataWriter->Delete();

      pressure_sing->Delete();
      //div_velocity_sing->Delete();

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "save_singularity_chi"
    template<std::size_t Dimensions, typename Shape>
    PetscErrorCode save_singularity_chi(const char* path,
                                    const char* filename,
                                    singularity<Shape, Dimensions> sing, 
                                    particle<Shape> const& p1, particle<Shape> const& p2,
                                    geometry::box<int, 2> box,
                                    std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      int rank;
      MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

      using position_type = geometry::position<double, Dimensions>;
      using position_type_i = geometry::position<int, Dimensions>;

      vtkStructuredGrid* singDataSet;

      singDataSet = vtkStructuredGrid::New();

      singDataSet->SetExtent(0, box.length(0)*sing.scale-1, 0, box.length(1)*sing.scale-1, 0, 0);

      vtkPoints* singPoints = vtkPoints::New();

      vtkDoubleArray* pressure_sing = vtkDoubleArray::New();
      pressure_sing->SetName("pressure_sing");

      std::array<double, Dimensions> hs = {{h[0]/sing.scale, h[1]/sing.scale}};
      double coef = 1./(sing.scale*sing.scale);

      for(std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j)
      {
        for(std::size_t js=0; js<sing.scale; ++js)
        {
          for(std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i)
          {
            position_type_i pts_i = {i, j};

            for(std::size_t is=0; is<sing.scale; ++is){
              position_type pts = {i*h[0] + is*hs[0], j*h[1] + js*hs[1]};
              singPoints->InsertNextPoint(pts[0], pts[1], 0.);

              bool add_sing = false;

              if (!p1.contains(pts) && !p2.contains(pts))
              {
                auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                if (std::abs(pos_ref_part[1]) <= sing.cutoff_dist_)
                {
                  add_sing = true;
                  position_type pts_loc = {is*hs[0], js*hs[1]};
                          
                  auto psing = sing.get_p_sing(pts);
                  // Add points to vtk + singular value to vtk
                  pressure_sing->InsertNextValue(psing);
                }
              }
              else if (p1.contains(pts))
              {
                position_type pts_polar = {sqrt( (pts[0]- p1.center_[0])*(pts[0]- p1.center_[0]) + (pts[1]- p1.center_[1])*(pts[1]- p1.center_[1]) ), atan2(pts[1]- p1.center_[1], pts[0]- p1.center_[0])};
                position_type pts_border = {p1.center_[0] + p1.shape_factors_[0]*cos(pts_polar[1]), p1.center_[1] + p1.shape_factors_[0]*sin(pts_polar[1])};
                auto pos_ref_part = sing.get_pos_in_part_ref(pts_border);
                add_sing = true;
                auto psing = sing.get_p_sing(pts_border)*cafes::singularity::alphaTrunc(pts_polar[0]/p1.shape_factors_[0]);
                // Add points to vtk + singular value to vtk
                pressure_sing->InsertNextValue(psing);
              }
              else if (p2.contains(pts))
              {
                position_type pts_polar = {sqrt( (pts[0]- p2.center_[0])*(pts[0]- p2.center_[0]) + (pts[1]- p2.center_[1])*(pts[1]- p2.center_[1]) ), atan2(pts[1]- p2.center_[1], pts[0]- p2.center_[0])};
                position_type pts_border = {p2.center_[0] + p2.shape_factors_[0]*cos(pts_polar[1]), p2.center_[1] + p2.shape_factors_[0]*sin(pts_polar[1])};
                auto pos_ref_part = sing.get_pos_in_part_ref(pts_border);
                add_sing = true;
                auto psing = sing.get_p_sing(pts_border)*cafes::singularity::alphaTrunc(pts_polar[0]/p2.shape_factors_[0]);
                // Add points to vtk + singular value to vtk
                pressure_sing->InsertNextValue(psing);
              }
              

              if (!add_sing)
              {
                pressure_sing->InsertNextValue(0.);
              }
            }
          }
        }
      }

      singDataSet->SetPoints(singPoints);
      singDataSet->GetPointData()->SetScalars(pressure_sing);

      vtkXMLStructuredGridWriter* singDataWriter = vtkXMLStructuredGridWriter::New();
      std::stringstream oc;

      oc << path << "/" << filename << "_sing_" << 0 << "_" << 1 << ".vts";//a changer
      singDataWriter->SetFileName(oc.str().data());
      //dataWriter->SetDataModeToAscii();
#if VTK_MAJOR_VERSION <= 5
      singDataWriter->SetInput(singDataSet);
#else
      singDataWriter->SetInputData(singDataSet);
#endif
      singDataWriter->Write();

      singPoints->Delete();
      singDataSet->Delete();
      singDataWriter->Delete();

      pressure_sing->Delete();

      PetscFunctionReturn(0);
    }


    #undef __FUNCT__
    #define __FUNCT__ "save_singularity"
    template<typename Shape, std::size_t Dimensions>
    PetscErrorCode save_singularity(const char* path,
                                    const char* filename,
                                    singularity<Shape, Dimensions> sing, 
                                    particle<Shape> const& p1, particle<Shape> const& p2,
                                    geometry::box<int, 3> box,
                                    std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      using position_type = geometry::position<double, Dimensions>;
      using position_type_i = geometry::position<int, Dimensions>;

      vtkStructuredGrid* singDataSet;

      singDataSet = vtkStructuredGrid::New();

      singDataSet->SetExtent(0, box.length(0)*sing.scale-1, 0, box.length(1)*sing.scale-1, 0, box.length(2)*sing.scale-1);

      vtkPoints* singPoints = vtkPoints::New();

#undef __FUNCT__
#define __FUNCT__ "add_singularity_to_surf"
        template<std::size_t Dimensions, typename Ctx>
        PetscErrorCode add_singularity_to_surf(
            Ctx &ctx,
            std::vector<std::vector<std::array<double, Dimensions>>> &g)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            std::cout << "add singularity to surf...\n";

            auto union_box_func = geometry::union_box<int, Dimensions>;

            auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
            auto &h = ctx.problem.ctx->h;

            // Loop on particles couples
            for (std::size_t ipart = 0; ipart < ctx.particles.size() - 1;
                 ++ipart)
            {
                auto p1 = ctx.particles[ipart];
                for (std::size_t jpart = ipart + 1;
                     jpart < ctx.particles.size(); ++jpart)
                {
                    auto p2 = ctx.particles[jpart];

                    using shape_type = typename decltype(p1)::shape_type;
                    auto sing =
                        singularity<shape_type, Dimensions>(p1, p2, h[0]);

                    if (sing.is_singularity_)
                    {
                        auto pbox = sing.get_box(h);
                        if (geometry::intersect(box, pbox))
                        {
                            auto new_box = geometry::box_inside(box, pbox);
                            geometry::box<double, Dimensions> new_box_d{
                                new_box.bottom_left, new_box.upper_right};
                            new_box_d.bottom_left *= h;
                            new_box_d.upper_right *= h;
                            ierr = computesingularBC(ctx, sing, ipart, jpart,
                                                     new_box_d, g);
                            CHKERRQ(ierr);
                        }
                    }
                }
            }

            PetscFunctionReturn(0);
        }

#undef __FUNCT__
#define __FUNCT__ "add_singularity_to_surf"
        template<std::size_t Dimensions, typename Ctx>
        PetscErrorCode
        add_singularity_to_surf(Ctx &ctx, petsc::petsc_vec<Dimensions> &sol)
        {
            PetscErrorCode ierr;
            PetscFunctionBeginUser;

            std::cout << "add singularity to surf new...\n";

            auto union_box_func = geometry::union_box<int, Dimensions>;

            auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
            auto &h = ctx.problem.ctx->h;

            // Loop on particles couples
            for (std::size_t ipart = 0; ipart < ctx.particles.size() - 1;
                 ++ipart)
            {
                auto p1 = ctx.particles[ipart];
                for (std::size_t jpart = ipart + 1;
                     jpart < ctx.particles.size(); ++jpart)
                {
                    auto p2 = ctx.particles[jpart];

                    using shape_type = typename decltype(p1)::shape_type;
                    auto sing =
                        singularity<shape_type, Dimensions>(p1, p2, h[0]);

                    // uncomment for 2D and implement for 3D
                    // if (sing.is_singularity_)
                    // {
                    //   ierr = computesingularBC(sing, sol, h,
                    //   box);CHKERRQ(ierr);
                    // }
                }
            }

            PetscFunctionReturn(0);
        }

    #undef __FUNCT__
    #define __FUNCT__ "save_singularity_chi"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode save_singularity_chi(const char* path,
            const char* filename,
            Ctx& ctx)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto union_box_func = geometry::union_box<int, Dimensions>;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      //Loop on particles couples
      for (std::size_t ipart=0; ipart<ctx.particles.size()-1; ++ipart)
      {
        auto p1 = ctx.particles[ipart];
        for (std::size_t jpart=ipart+1; jpart<ctx.particles.size(); ++jpart)
        {
          auto p2 = ctx.particles[jpart];

          using shape_type = typename decltype(p1)::shape_type;
          auto sing = singularity<shape_type, Dimensions>(p1, p2, h[0]);

          if (sing.is_singularity_)
          {
            // auto pbox = union_box_func({geometry::floor((p1.center_ - sing.cutoff_dist_)/h), 
            //                             geometry::ceil((p1.center_ + sing.cutoff_dist_)/h)},
            //                            {geometry::floor((p2.center_ - sing.cutoff_dist_)/h), 
            //                             geometry::ceil((p2.center_ + sing.cutoff_dist_)/h)});
            auto pbox = sing.get_box(h);

            if (geometry::intersect(box, pbox))
            {
              auto new_box = geometry::box_inside(box, pbox);

              ierr = save_singularity_chi(path, filename, sing, p1, p2, new_box, h);CHKERRQ(ierr);
            }
          }
        }
      }

      PetscFunctionReturn(0);
    }

  }
}

#endif
