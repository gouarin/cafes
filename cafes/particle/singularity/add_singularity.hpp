#ifndef CAFES_PARTICLE_SINGULARITY_ADD_SINGULARITY_HPP_INCLUDED
#define CAFES_PARTICLE_SINGULARITY_ADD_SINGULARITY_HPP_INCLUDED

#include <particle/particle.hpp>
#include <particle/singularity/singularity.hpp>
#include <particle/geometry/box.hpp>
#include <petsc/vec.hpp>

#include <petsc.h>
#include <iostream>
#include <sstream>

#include "vtkDoubleArray.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkStructuredGrid.h"
#include "vtkXMLStructuredGridWriter.h"

namespace cafes
{
  namespace singularity
  {
    #undef __FUNCT__
    #define __FUNCT__ "computesingularST"
    template<typename Shape>
    PetscErrorCode computesingularST(singularity<2, Shape> sing, 
                                     particle<Shape> const& p1, particle<Shape> const& p2,
                                     petsc::petsc_vec<2>& sol,
                                     geometry::box<2, int> box,
                                     std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 2;
      using position_type = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;

      std::array<double, Dimensions> hs = {{h[0]/sing.scale, h[1]/sing.scale}};
      double coef = 1./(sing.scale*sing.scale);

      for(std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j)
      {
        for(std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i)
        {
          position_type_i pts_i = {i, j};
          auto ielem = fem::get_element(pts_i);

          for(std::size_t js=0; js<sing.scale; ++js)
          {
            for(std::size_t is=0; is<sing.scale; ++is){
              position_type pts = {i*h[0] + is*hs[0], j*h[1] + js*hs[1]};

              if (!p1.contains(pts) && !p2.contains(pts))
              {
                auto pos_ref_part = sing.get_pos_in_part_ref(pts);

                if (abs(pos_ref_part[1]) < sing.cutoff_dist_)
                {
                  position_type pts_loc = {is*hs[0], js*hs[1]};
                  
                  auto bfunc = fem::P1_integration_grad(pts_loc, h);

                  auto gradUsing = sing.get_grad_u_sing(pts);
                  auto psing = sing.get_p_sing(pts);
                  
                  for (std::size_t j=0; j<bfunc.size(); ++j){
                    auto u = sol.at(ielem[j]);

                    for (std::size_t d1=0; d1<Dimensions; ++d1)
                    {
                      for (std::size_t d2=0; d2<Dimensions; ++d2)
                        u[d1] -= coef*gradUsing[d1][d2]*bfunc[j][d2];
                      u[d1] += coef*psing*bfunc[j][d1];
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
    #define __FUNCT__ "computesingularST"
    template<typename Shape>
    PetscErrorCode computesingularST(singularity<3, Shape> sing, 
                                     particle<Shape> const& p1, particle<Shape> const& p2,
                                     petsc::petsc_vec<3>& sol,
                                     geometry::box<3, int> box,
                                     std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      const int Dimensions = 3;
      using position_type = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;

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
              for(std::size_t js=0; js<sing.scale; ++js)
              {
                for(std::size_t is=0; is<sing.scale; ++is){
                  position_type pts = {i*h[0] + is*hs[0], 
                                       j*h[1] + js*hs[1],
                                       k*h[2] + ks*hs[2] };

                  if (!p1.contains(pts) && !p2.contains(pts))
                  {
                    position_type pts_loc = {is*hs[0], 
                                             js*hs[1],
                                             ks*hs[2]};
                    
                    auto bfunc = fem::P1_integration_grad(pts_loc, h);

                    auto gradUsing = sing.get_grad_u_sing(pts);
                    auto psing = sing.get_p_sing(pts);
                    
                    for (std::size_t j=0; j<bfunc.size(); ++j){
                      auto u = sol.at(ielem[j]);

                      for (std::size_t d1=0; d1<Dimensions; ++d1){
                        for (std::size_t d2=0; d2<Dimensions; ++d2)
                          u[d1] -= coef*gradUsing[d1][d2]*bfunc[j][d2];
                        u[d1] += coef*psing*bfunc[j][d1];
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
    template<std::size_t Dimensions, typename Shape, typename Ctx>
    PetscErrorCode computesingularBC(Ctx& ctx,
                                     singularity<Dimensions, Shape> sing, 
                                     std::size_t ipart_1, std::size_t ipart_2,
                                     std::vector<std::vector<std::array<double, Dimensions>>>& g
                                    )
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      for(std::size_t isurf=0; isurf<ctx.surf_points[ipart_1].size(); ++isurf)
      {
        auto pos = ctx.surf_points[ipart_1][isurf].second;
        auto pos_ref_part = sing.get_pos_in_part_ref(pos);
        if (abs(pos_ref_part[1]) < sing.cutoff_dist_)
        {
          auto Using = sing.get_u_sing(pos);
          for (std::size_t d=0; d<Dimensions; ++d)
            g[ipart_1][isurf][d] -= Using[d];
        }
      }

      for(std::size_t isurf=0; isurf<ctx.surf_points[ipart_2].size(); ++isurf)
      {
        auto pos = ctx.surf_points[ipart_2][isurf].second;
        auto pos_ref_part = sing.get_pos_in_part_ref(pos);
        if (abs(pos_ref_part[1]) < sing.cutoff_dist_)
        {
          auto Using = sing.get_u_sing(pos);
          for (std::size_t d=0; d<Dimensions; ++d)
            g[ipart_2][isurf][d] -= Using[d];
        }
      }

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "add_singularity_in_fluid"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode add_singularity_in_fluid(Ctx& ctx)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto union_box_func = geometry::union_box<Dimensions, int>;

      auto box = fem::get_DM_bounds<Dimensions>(ctx.problem.ctx->dm, 0);
      auto& h = ctx.problem.ctx->h;

      auto sol = petsc::petsc_vec<Dimensions>(ctx.problem.ctx->dm, ctx.problem.rhs, 0, false);
      ierr = sol.fill(0.);CHKERRQ(ierr);

      //Loop on particles couples
      for (std::size_t ipart=0; ipart<ctx.particles.size()-1; ++ipart)
      {
        auto p1 = ctx.particles[ipart];
        for (std::size_t jpart=ipart+1; jpart<ctx.particles.size(); ++jpart)
        {
          auto p2 = ctx.particles[jpart];

          using shape_type = typename decltype(p1)::shape_type;
          auto sing = singularity<Dimensions, shape_type>(p1, p2, h[0]);

          if (sing.is_singularity_)
          {
            auto p1_bl = (p1.center_ - sing.cutoff_dist_)/h;
            auto p1_ur = (p1.center_ + sing.cutoff_dist_)/h;

            auto p2_bl = (p2.center_ - sing.cutoff_dist_)/h;
            auto p2_ur = (p2.center_ + sing.cutoff_dist_)/h;

            std::transform(p1_bl.begin(), p1_bl.end(), p1_bl.begin(), [](auto& p){return std::floor(p);});
            std::transform(p1_ur.begin(), p1_ur.end(), p1_ur.begin(), [](auto& p){return std::ceil(p);});
            std::transform(p2_bl.begin(), p2_bl.end(), p2_bl.begin(), [](auto& p){return std::floor(p);});
            std::transform(p2_ur.begin(), p2_ur.end(), p2_ur.begin(), [](auto& p){return std::ceil(p);});

            auto pbox = union_box_func({p1_bl, p1_ur}, {p2_bl, p2_ur});

            // auto pbox = union_box_func({(p1.center_ - sing.cutoff_dist_)/h, (p1.center_ + sing.cutoff_dist_)/h},
            //                            {(p2.center_ - sing.cutoff_dist_)/h, (p2.center_ + sing.cutoff_dist_)/h});

            if (geometry::intersect(box, pbox))
            {
              auto new_box = geometry::box_inside(box, pbox);
              ierr = computesingularST(sing, p1, p2, sol, new_box, h);CHKERRQ(ierr);
            }
          }
        }
      }

      ierr = sol.local_to_global(ADD_VALUES);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "add_singularity_to_surf"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode add_singularity_to_surf(Ctx& ctx, std::vector<std::vector<std::array<double, Dimensions>>>& g)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto union_box_func = geometry::union_box<Dimensions, int>;

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
          auto sing = singularity<Dimensions, shape_type>(p1, p2, h[0]);

          if (sing.is_singularity_)
          {
            auto p1_bl = (p1.center_ - sing.cutoff_dist_)/h;
            auto p1_ur = (p1.center_ + sing.cutoff_dist_)/h;

            auto p2_bl = (p2.center_ - sing.cutoff_dist_)/h;
            auto p2_ur = (p2.center_ + sing.cutoff_dist_)/h;

            std::transform(p1_bl.begin(), p1_bl.end(), p1_bl.begin(), [](auto& p){return std::floor(p);});
            std::transform(p1_ur.begin(), p1_ur.end(), p1_ur.begin(), [](auto& p){return std::ceil(p);});
            std::transform(p2_bl.begin(), p2_bl.end(), p2_bl.begin(), [](auto& p){return std::floor(p);});
            std::transform(p2_ur.begin(), p2_ur.end(), p2_ur.begin(), [](auto& p){return std::ceil(p);});

            auto pbox = union_box_func({p1_bl, p1_ur}, {p2_bl, p2_ur});
            // auto pbox = union_box_func({(p1.center_ - sing.cutoff_dist_)/h, (p1.center_ + sing.cutoff_dist_)/h},
            //                           {(p2.center_ - sing.cutoff_dist_)/h, (p2.center_ + sing.cutoff_dist_)/h});

            if (geometry::intersect(box, pbox))
            {
              ierr = computesingularBC(ctx, sing, ipart, jpart, g);CHKERRQ(ierr);
            }
          }
        }
      }

      PetscFunctionReturn(0);
    }


    #undef __FUNCT__
    #define __FUNCT__ "save_singularity"
    template<std::size_t Dimensions, typename Shape>
    PetscErrorCode save_singularity(const char* path,
                                    const char* filename,
                                    singularity<Dimensions, Shape> sing, 
                                    particle<Shape> const& p1, particle<Shape> const& p2,
                                    geometry::box<2, int> box,
                                    std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      using position_type = geometry::position<Dimensions, double>;
      using position_type_i = geometry::position<Dimensions, int>;

      vtkStructuredGrid* singDataSet;

      singDataSet = vtkStructuredGrid::New();

      singDataSet->SetExtent(0, box.length(0)*sing.scale-1, 0, box.length(1)*sing.scale-1, 0, 0);

      vtkPoints* singPoints = vtkPoints::New();

      vtkDoubleArray* velocity_sing = vtkDoubleArray::New();
      velocity_sing->SetNumberOfComponents(3);
      velocity_sing->SetName("velocity_sing");

      vtkDoubleArray* gradx_velocity_sing = vtkDoubleArray::New();
      gradx_velocity_sing->SetNumberOfComponents(3);
      gradx_velocity_sing->SetName("gradx_velocity_sing");

      vtkDoubleArray* grady_velocity_sing = vtkDoubleArray::New();
      grady_velocity_sing->SetNumberOfComponents(3);
      grady_velocity_sing->SetName("grady_velocity_sing");

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

                if (abs(pos_ref_part[1]) < sing.cutoff_dist_)
                {
                  add_sing = true;
                  position_type pts_loc = {is*hs[0], js*hs[1]};
                  
                  auto Using = sing.get_u_sing(pts);
                  auto gradUsing = sing.get_grad_u_sing(pts);
                  auto psing = sing.get_p_sing(pts);
                  // Add points to vtk + singular value to vtk
                  velocity_sing->InsertNextTuple3(Using[0], Using[1], 0.);
                  gradx_velocity_sing->InsertNextTuple3(gradUsing[0][0], gradUsing[0][1], 0.);
                  grady_velocity_sing->InsertNextTuple3(gradUsing[1][0], gradUsing[1][1], 0.);
                  pressure_sing->InsertNextValue(psing);
                }
              }

              if (!add_sing)
              {
                velocity_sing->InsertNextTuple3(0., 0., 0.);
                gradx_velocity_sing->InsertNextTuple3(0., 0., 0.);
                grady_velocity_sing->InsertNextTuple3(0., 0., 0.);
                pressure_sing->InsertNextValue(0.);
              }
            }
          }
        }
      }

      singDataSet->SetPoints(singPoints);
      singDataSet->GetPointData()->AddArray(velocity_sing);
      singDataSet->GetPointData()->AddArray(gradx_velocity_sing);
      singDataSet->GetPointData()->AddArray(grady_velocity_sing);
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
      velocity_sing->Delete();
      gradx_velocity_sing->Delete();
      grady_velocity_sing->Delete();
      singDataSet->Delete();
      singDataWriter->Delete();

      pressure_sing->Delete();

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "save_singularity"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode save_singularity(const char* path,
                                    const char* filename,
                                    Ctx& ctx)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto union_box_func = geometry::union_box<Dimensions, int>;

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
          auto sing = singularity<Dimensions, shape_type>(p1, p2, h[0]);

          if (sing.is_singularity_)
          {
            auto p1_bl = (p1.center_ - sing.cutoff_dist_)/h;
            auto p1_ur = (p1.center_ + sing.cutoff_dist_)/h;

            auto p2_bl = (p2.center_ - sing.cutoff_dist_)/h;
            auto p2_ur = (p2.center_ + sing.cutoff_dist_)/h;

            std::transform(p1_bl.begin(), p1_bl.end(), p1_bl.begin(), [](auto& p){return std::floor(p);});
            std::transform(p1_ur.begin(), p1_ur.end(), p1_ur.begin(), [](auto& p){return std::ceil(p);});
            std::transform(p2_bl.begin(), p2_bl.end(), p2_bl.begin(), [](auto& p){return std::floor(p);});
            std::transform(p2_ur.begin(), p2_ur.end(), p2_ur.begin(), [](auto& p){return std::ceil(p);});

            auto pbox = union_box_func({p1_bl, p1_ur}, {p2_bl, p2_ur});

            // auto pbox = union_box_func({(p1.center_ - sing.cutoff_dist_)/h, (p1.center_ + sing.cutoff_dist_)/h},
            //                            {(p2.center_ - sing.cutoff_dist_)/h, (p2.center_ + sing.cutoff_dist_)/h});

            if (geometry::intersect(box, pbox))
            {
              auto new_box = geometry::box_inside(box, pbox);
              ierr = save_singularity(path, filename, sing, p1, p2, new_box, h);CHKERRQ(ierr);
            }
          }
        }
      }

      PetscFunctionReturn(0);
    }

  }
}

#endif