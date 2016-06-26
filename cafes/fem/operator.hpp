#ifndef CAFES_FEM_OPERATOR_HPP_INCLUDED
#define CAFES_FEM_OPERATOR_HPP_INCLUDED

#include <fem/matElem.hpp>
#include <fem/mesh.hpp>
#include <fem/quadrature.hpp>
#include <particle/geometry/position.hpp>
#include <petsc/vec.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <petsc.h>

namespace cafes
{
  namespace fem
  {

    template<typename Box, typename Function, typename Position>
    void iterate_impl(Box const& b, Function&& f, Position& p, std::integral_constant<std::size_t,0> const&)
    {
      std::forward<Function>(f)(p);
    }

    template<typename Box, typename Function, typename Position, typename Index>
    void iterate_impl(Box const& b, Function&& f, Position& p, Index const&)
    {
      static constexpr std::size_t n = Index::value-1;

      for( p[n] = b.bottom_left[n]; p[n] < b.upper_right[n]; ++p[n] )
      {
        iterate_impl(b, std::forward<Function>(f), p, std::integral_constant<std::size_t,n>{});
      }
    }

    template<typename Box, typename Function>
    void iterate(Box const& b, Function&& f)
    {
      typename Box::position_type pos;
      iterate_impl(b, f, pos, std::integral_constant<std::size_t,Box::dimensions>{});
    }


    auto const kernel_1 = [](auto const& x, auto& y, auto const& matelem){
      auto const kernel_pos = [&](auto const& pos){
        auto const ielem = get_element(pos);
        auto const nbasis = ielem.size();

        for(std::size_t k1=0; k1<nbasis; ++k1){
          auto uy = y.at(ielem[k1]);
          for(std::size_t k2=0; k2<nbasis; ++k2){
            auto ux = x.at(ielem[k2]);
            for(std::size_t d=0; d<x.dof_; ++d){
                uy[d] += ux[d]*matelem[k1*nbasis+k2];
            }
          }
        }
      };
      return kernel_pos;
    };

    auto const kernel_11 = [](auto const& x, auto& y, auto const& matelem){
      auto const kernel_pos = [&](auto const& pos){
        auto const ielem = get_element(pos);
        auto const nbasis = ielem.size();

        for(std::size_t k1=0; k1<nbasis; ++k1)
        {
          auto uy = y.at(ielem[k1]);
          for(std::size_t k2=0; k2<nbasis; ++k2)
          {
            auto ux = x.at(ielem[k2]);
            for(std::size_t d1=0; d1<x.dof_; ++d1)
            {
              for(std::size_t d2=0; d2<x.dof_; ++d2)
              {
                uy[d1] += ux[d2]*matelem[k1*nbasis+k2][d1*x.dof_ + d2];
              }
            }
          }
        }
      };
      return kernel_pos;
    };

    auto const kernel_2 = [](auto& x, auto const& matelem){
      auto const kernel_pos = [&](auto const& pos){
        auto const ielem = get_element(pos);
        auto const nbasis = ielem.size();

        for(std::size_t ie=0; ie<nbasis; ++ie){
          auto u = x.at(ielem[ie]);
          for(std::size_t d=0; d<x.dof_; ++d){
            u[d] += matelem[ie*nbasis+ie];
          }
        }
      };
      return kernel_pos;
    };

    auto const kernel_3 = [](auto const& x1, auto const& x2, auto& y1, auto& y2, auto const& matelem){
      auto const kernel_pos = [&](auto const& pos){
        auto const ielem_p = get_element(pos);
        auto const ielem_v = get_element_4Q1(pos);
        auto const nbasis_p = ielem_p.size();
        auto const nbasis_v = ielem_v.size();

        for(std::size_t ie_v=0; ie_v<nbasis_v; ++ie_v){
          auto ux = x1.at(ielem_v[ie_v]);
          auto uy = y1.at(ielem_v[ie_v]);
          for(std::size_t ie_p=0; ie_p<nbasis_p; ++ie_p){
            auto uxp = x2.at(ielem_p[ie_p]);
            auto uyp = y2.at(ielem_p[ie_p]);
            for(std::size_t d=0; d<x1.dof_; ++d){
              uyp[0] += ux[d]*matelem[ie_p*nbasis_v+ie_v][d];
              uy[d] -= uxp[0]*matelem[ie_p*nbasis_v+ie_v][d];
            }
          }
        }
      };
      return kernel_pos;
    };

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_mult"
    template<typename MatElem, typename Function, std::size_t Dimensions>
    PetscErrorCode diag_block_mult(petsc::petsc_vec<Dimensions>& x, petsc::petsc_vec<Dimensions>& y,
                                   MatElem const& matelem, Function&& kernel)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto box = get_DM_bounds<Dimensions>(x.dm_);

      iterate(box, kernel(x, y, matelem));

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_mult_diag"
    template<typename MatElem, typename Function, std::size_t Dimensions>
    PetscErrorCode diag_block_mult_diag(petsc::petsc_vec<Dimensions>& x,
                                        MatElem const& matelem, Function&& kernel)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto box = get_DM_bounds<Dimensions>(x.dm_);

      iterate(box, kernel(x, matelem));

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "offset_block_mult"
    template<typename MatElem, typename Function, std::size_t Dimensions>
    PetscErrorCode offset_block_mult(petsc::petsc_vec<Dimensions> const& x1,
                                     petsc::petsc_vec<Dimensions> const& x2,
                                     petsc::petsc_vec<Dimensions>& y1,
                                     petsc::petsc_vec<Dimensions>& y2,
                                     MatElem const& matelem, Function&& kernel)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto box = get_DM_bounds<Dimensions>(x2.dm_);

      iterate(box, kernel(x1, x2, y1, y2, matelem));

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "laplacian_mult"
    template<std::size_t Dimensions>
    PetscErrorCode laplacian_mult(petsc::petsc_vec<Dimensions>& x, 
                                  petsc::petsc_vec<Dimensions>& y, 
                                  std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = diag_block_mult(x, y, getMatElemLaplacian(h), kernel_1);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "laplacian_mult_diag"
    template<std::size_t Dimensions>
    PetscErrorCode laplacian_mult_diag(petsc::petsc_vec<Dimensions>& x, std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = diag_block_mult_diag(x, getMatElemLaplacian(h), kernel_2);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "mass_mult"
    template<std::size_t Dimensions>
    PetscErrorCode mass_mult(petsc::petsc_vec<Dimensions>& x, petsc::petsc_vec<Dimensions>& y, std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = diag_block_mult(x, y, getMatElemMass(h), kernel_1);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "mass_mult_diag"
    template<std::size_t Dimensions>
    PetscErrorCode mass_mult_diag(petsc::petsc_vec<Dimensions>& x, std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = diag_block_mult_diag(x, getMatElemMass(h), kernel_2);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "strain_tensor_mult"
    template<std::size_t Dimensions>
    PetscErrorCode strain_tensor_mult(petsc::petsc_vec<Dimensions>& x, petsc::petsc_vec<Dimensions>& y, std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = diag_block_mult(x, y, getMatElemStrainTensor(h), kernel_11);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "B_and_BT_mult"
    template<std::size_t Dimensions>
    PetscErrorCode B_and_BT_mult(petsc::petsc_vec<Dimensions> const& x1,
                                 petsc::petsc_vec<Dimensions> const& x2,
                                 petsc::petsc_vec<Dimensions>& y1,
                                 petsc::petsc_vec<Dimensions>& y2,
                                 std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = offset_block_mult(x1, x2, y1, y2, getMatElemPressure(h), kernel_3);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }
  }
}
#endif