#ifndef CAFES_FEM_OPERATOR_HPP_INCLUDED
#define CAFES_FEM_OPERATOR_HPP_INCLUDED

#include <fem/matElem.hpp>
#include <fem/mesh.hpp>
#include <array>
#include <iostream>
#include <petsc.h>

namespace cafes
{
  namespace fem
  {

    /////////////////////////////////////////////////////////////////////////////
    //
    //  2D operators
    //
    /////////////////////////////////////////////////////////////////////////////
    #undef __FUNCT__
    #define __FUNCT__ "diag_block_mult"
    PetscErrorCode diag_block_mult(DM dm, Vec x, Vec y, std::array<double, 2> h, int dof,
                                   PetscErrorCode(*matelem)(PetscReal*, std::array<double, 2> const&))
    {
      PetscErrorCode ierr;
      PetscScalar tmp[4][2];
      PetscScalar val[4][2];
      PetscReal MatElem[16];
      PetscScalar ***px, ***py;

      PetscFunctionBeginUser;

      ierr = (*matelem)(MatElem, h);

      auto box = get_DM_bounds<2>(dm);

      ierr = DMDAVecGetArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(dm, y, &py);CHKERRQ(ierr);

      for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
        for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
          for(std::size_t k1=0; k1<4; k1++){
            std::size_t l1 = i + ind2d[k1][0];
            std::size_t p1 = j + ind2d[k1][1];
            for(std::size_t k2=0; k2<dof; ++k2){
              tmp[k1][k2] = px[p1][l1][k2];
              val[k1][k2] = 0.;
            }
          }

          for(std::size_t k1=0; k1<4; ++k1){
            for(std::size_t k2=0; k2<4; ++k2){
                for(std::size_t k3=0; k3<dof; ++k3){
                    val[k1][k3] += tmp[k2][k3]*MatElem[k1*4+k2];
                }
            }
          }

          for(std::size_t k1=0; k1<4; ++k1){
            std::size_t l1 = i + ind2d[k1][0];
            std::size_t p1 = j + ind2d[k1][1];
            for(std::size_t k2=0; k2<dof; ++k2)
              py[p1][l1][k2] += val[k1][k2];
          }
        }
      }
      ierr = DMDAVecRestoreArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dm, y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_mult_diag"
    PetscErrorCode diag_block_mult_diag(DM dm, Vec x, std::array<double, 2> h, int dof,
                                        PetscErrorCode(*matelem)(PetscReal*, std::array<double, 2> const&))
    {
      PetscErrorCode ierr;
      PetscReal MatElem[16];
      PetscScalar ***px;

      PetscFunctionBeginUser;

      ierr = (*matelem)(MatElem, h);

      auto box = get_DM_bounds<2>(dm);

      ierr = DMDAVecGetArrayDOF(dm, x, &px);CHKERRQ(ierr);

      for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
        for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
          for(std::size_t k1=0; k1<4; k1++){
            std::size_t l1 = i + ind2d[k1][0];
            std::size_t p1 = j + ind2d[k1][1];
            for(std::size_t k2=0; k2<dof; ++k2){
                px[p1][l1][k2] += MatElem[k1*4+k1];
            }
          }
        }
      }
      ierr = DMDAVecRestoreArrayDOF(dm, x, &px);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "laplacian_mult"
    PetscErrorCode laplacian_mult(DM dm, Vec x, Vec y, std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult(dm, x, y, h, info.dof, getMatElemLaplacian2d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "laplacian_mult_diag"
    PetscErrorCode laplacian_mult_diag(DM dm, Vec x, std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult_diag(dm, x, h, info.dof, getMatElemLaplacian2d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "mass_mult"
    PetscErrorCode mass_mult(DM dm, Vec x, Vec y, std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult(dm, x, y, h, info.dof, getMatElemMass2d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "mass_mult_diag"
    PetscErrorCode mass_mult_diag(DM dm, Vec x, std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult_diag(dm, x, h, info.dof, getMatElemMass2d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "strain_tensor_mult"
    PetscErrorCode strain_tensor_mult(DM dm, Vec x, Vec y, std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscScalar tmp[4][2];
      PetscScalar val[4][2];
      PetscReal MatElem[16][4];
      PetscScalar ***px, ***py;

      PetscFunctionBeginUser;
      ierr = getMatElemStrainTensor2d(MatElem, h);

      auto box = get_DM_bounds<2>(dm);

      ierr = DMDAVecGetArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(dm, y, &py);CHKERRQ(ierr);

      for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
        for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
          for(std::size_t k1=0; k1<4; k1++){
            std::size_t l1 = i + ind2d[k1][0];
            std::size_t p1 = j + ind2d[k1][1];
            for(std::size_t k2=0; k2<2; ++k2){
              tmp[k1][k2] = px[p1][l1][k2];
              val[k1][k2] = 0.;
            }
          }

          for(std::size_t k1=0; k1<4; ++k1){
            for(std::size_t k2=0; k2<4; ++k2){
              val[k1][0] += tmp[k2][0]*MatElem[k1*4+k2][0] 
                          + tmp[k2][1]*MatElem[k1*4+k2][1];

              val[k1][1] += tmp[k2][0]*MatElem[k1*4+k2][2] 
                          + tmp[k2][1]*MatElem[k1*4+k2][3];
            }
          }

          for(std::size_t k1=0; k1<4; ++k1){
            std::size_t l1 = i + ind2d[k1][0];
            std::size_t p1 = j + ind2d[k1][1];
            for(std::size_t k2=0; k2<2; ++k2)
              py[p1][l1][k2] += val[k1][k2];
          }
        }
      }
      ierr = DMDAVecRestoreArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dm, y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "B_and_BT_mult"
    PetscErrorCode B_and_BT_mult(DM dm, Vec xu, Vec xp, Vec yu, Vec yp, std::array<double, 2> const& h)
    {
      PetscErrorCode ierr;
      PetscScalar tmp[4];
      PetscScalar val[4][2];
      PetscReal MatElemOnu[36], MatElemOnv[36];
      DM dav, dap;
      PetscScalar ***pxu, ***pyu;
      PetscScalar **pxp, **pyp;

      PetscFunctionBeginUser;

      ierr = DMCompositeGetEntries(dm, &dav, &dap);CHKERRQ(ierr);
      auto box = get_DM_bounds<2>(dap);

      ierr = DMDAVecGetArrayDOFRead(dav, xu, &pxu);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(dav, yu, &pyu);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayRead(dap, xp, &pxp);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(dap, yp, &pyp);CHKERRQ(ierr);

      ierr = getMatElemPressure2d(MatElemOnu, MatElemOnv, h);

      for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
        for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
          for(std::size_t k1=0; k1<4; ++k1){
            std::size_t l1 = i + ind2d[k1][0];
            std::size_t p1 = j + ind2d[k1][1];
            tmp[k1] = pxp[p1][l1];
            val[k1][0] = 0.;
            val[k1][1] = 0.;
          }

          for(std::size_t k1=0; k1<9; ++k1){
            std::size_t l1 = 2*i + indpu2d[k1][0];
            std::size_t p1 = 2*j + indpu2d[k1][1];

            auto tmp1v = pxu[p1][l1][0];
            auto tmp2v = pxu[p1][l1][1];

            auto val1 = 0.;
            auto val2 = 0.;

            for(std::size_t k2=0; k2<4; ++k2){
              val1 -= tmp[k2]*MatElemOnu[k2*9+k1];
              val2 -= tmp[k2]*MatElemOnv[k2*9+k1];

              val[k2][0] += tmp1v*MatElemOnu[k2*9+k1];
              val[k2][1] += tmp2v*MatElemOnv[k2*9+k1];
            }

            pyu[p1][l1][0] += val1;
            pyu[p1][l1][1] += val2;
          }

          for(std::size_t k1=0; k1<4; ++k1){
            std::size_t l1 = i + ind2d[k1][0];
            std::size_t p1 = j + ind2d[k1][1];
            pyp[p1][l1] += val[k1][0] + val[k1][1];
          }
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dav, xu, &pxu);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dav, yu, &pyu);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayRead(dap, xp, &pxp);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(dap, yp, &pyp);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    /////////////////////////////////////////////////////////////////////////////
    //
    //  3D operators
    //
    /////////////////////////////////////////////////////////////////////////////
    #undef __FUNCT__
    #define __FUNCT__ "diag_block_mult"
    PetscErrorCode diag_block_mult(DM dm, Vec x, Vec y, std::array<double, 3> const& h, int dof,
                                   PetscErrorCode(*matelem)(PetscReal*, std::array<double, 3> const&))
    {
      PetscErrorCode ierr;
      PetscInt nbasis = 8;
      PetscScalar tmp[nbasis][3];
      PetscScalar val[nbasis][3];
      PetscReal MatElem[nbasis*nbasis];
      PetscScalar ****px, ****py;

      PetscFunctionBeginUser;

      ierr = (*matelem)(MatElem, h);

      auto box = get_DM_bounds<3>(dm);

      ierr = DMDAVecGetArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(dm, y, &py);CHKERRQ(ierr);

      for (std::size_t k=box.bottom_left[2]; k<box.upper_right[2]; ++k){
          for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
            for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
              for(std::size_t k1=0; k1<nbasis; k1++){
                std::size_t l1 = i + ind3d[k1][0];
                std::size_t p1 = j + ind3d[k1][1];
                std::size_t r1 = k + ind3d[k1][2];
                for(std::size_t k2=0; k2<dof; ++k2){
                  tmp[k1][k2] = px[r1][p1][l1][k2];
                  val[k1][k2] = 0.;
                }
              }

              for(std::size_t k1=0; k1<nbasis; ++k1){
                for(std::size_t k2=0; k2<nbasis; ++k2){
                    for(std::size_t k3=0; k3<dof; ++k3){
                        val[k1][k3] += tmp[k2][k3]*MatElem[k1*nbasis+k2];
                    }
                }
              }

              for(std::size_t k1=0; k1<nbasis; ++k1){
                std::size_t l1 = i + ind3d[k1][0];
                std::size_t p1 = j + ind3d[k1][1];
                std::size_t r1 = k + ind3d[k1][2];
                for(std::size_t k2=0; k2<dof; ++k2)
                  py[r1][p1][l1][k2] += val[k1][k2];
              }
            }
          }
      }
      ierr = DMDAVecRestoreArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dm, y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_mult_diag"
    PetscErrorCode diag_block_mult_diag(DM dm, Vec x, std::array<double, 3> const& h, int dof,
                                        PetscErrorCode(*matelem)(PetscReal*, std::array<double, 3> const&))
    {
      PetscErrorCode ierr;
      PetscInt nbasis = 8;
      PetscReal MatElem[nbasis*nbasis];
      PetscScalar ****px;

      PetscFunctionBeginUser;

      ierr = (*matelem)(MatElem, h);

      auto box = get_DM_bounds<3>(dm);

      ierr = DMDAVecGetArrayDOF(dm, x, &px);CHKERRQ(ierr);

      for (std::size_t k=box.bottom_left[2]; k<box.upper_right[2]; ++k){
          for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
            for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
              for(std::size_t k1=0; k1<nbasis; k1++){
                std::size_t l1 = i + ind3d[k1][0];
                std::size_t p1 = j + ind3d[k1][1];
                std::size_t r1 = k + ind3d[k1][2];
                for(std::size_t k2=0; k2<dof; ++k2){
                    px[r1][p1][l1][k2] += MatElem[k1*nbasis+k1];
                }
              }
            }
          }
      }
      ierr = DMDAVecRestoreArrayDOF(dm, x, &px);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "laplacian_mult"
    PetscErrorCode laplacian_mult(DM dm, Vec x, Vec y, std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult(dm, x, y, h, info.dof, getMatElemLaplacian3d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "laplacian_mult_diag"
    PetscErrorCode laplacian_mult_diag(DM dm, Vec x, std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult_diag(dm, x, h, info.dof, getMatElemLaplacian3d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "mass_mult"
    PetscErrorCode mass_mult(DM dm, Vec x, Vec y, std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult(dm, x, y, h, info.dof, getMatElemMass3d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "mass_mult_diag"
    PetscErrorCode mass_mult_diag(DM dm, Vec x, std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      DMDALocalInfo info;
      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);
      ierr = diag_block_mult_diag(dm, x, h, info.dof, getMatElemMass3d);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "strain_tensor_mult"
    PetscErrorCode strain_tensor_mult(DM dm, Vec x, Vec y, std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscInt nbasisfunc = 8, nbasis2 = nbasisfunc*nbasisfunc;
      PetscScalar tmp[nbasisfunc][3];
      PetscScalar val[nbasisfunc][3];
      PetscScalar ****px, ****py;
      PetscReal MatElem[nbasis2][9];

      PetscFunctionBeginUser;

      ierr = getMatElemStrainTensor3d(MatElem, h);

      auto box = get_DM_bounds<3>(dm);

      ierr = DMDAVecGetArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(dm, y, &py);CHKERRQ(ierr);

      for (std::size_t k=box.bottom_left[2]; k<box.upper_right[2]; ++k){
        for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
          for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
            for(std::size_t k1=0; k1<nbasisfunc; ++k1){
              std::size_t l1 = i + ind3d[k1][0];
              std::size_t p1 = j + ind3d[k1][1];
              std::size_t r1 = k + ind3d[k1][2];
              for(std::size_t k2=0; k2<3; ++k2){
                tmp[k1][k2] = px[r1][p1][l1][k2];
                val[k1][k2] = 0.;
              }
            }

            for(std::size_t k1=0; k1<nbasisfunc; ++k1){
              for(std::size_t k2=0; k2<nbasisfunc; ++k2){
                val[k1][0] += tmp[k2][0]*MatElem[k1*nbasisfunc+k2][0] 
                            + tmp[k2][1]*MatElem[k1*nbasisfunc+k2][1] 
                            + tmp[k2][2]*MatElem[k1*nbasisfunc+k2][2];
                
                val[k1][1] += tmp[k2][0]*MatElem[k1*nbasisfunc+k2][3] 
                            + tmp[k2][1]*MatElem[k1*nbasisfunc+k2][4] 
                            + tmp[k2][2]*MatElem[k1*nbasisfunc+k2][5];
                
                val[k1][2] += tmp[k2][0]*MatElem[k1*nbasisfunc+k2][6] 
                            + tmp[k2][1]*MatElem[k1*nbasisfunc+k2][7] 
                            + tmp[k2][2]*MatElem[k1*nbasisfunc+k2][8];
              }
            }

            for(std::size_t k1=0; k1<nbasisfunc; ++k1){
              std::size_t l1 = i + ind3d[k1][0];
              std::size_t p1 = j + ind3d[k1][1];
              std::size_t r1 = k + ind3d[k1][2];
              for(std::size_t k2=0; k2<3; ++k2)
                py[r1][p1][l1][k2] += val[k1][k2];
            }
          }
        }
      }
      ierr = DMDAVecRestoreArrayDOFRead(dm, x, &px);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dm, y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "B_and_BT_mult"
    PetscErrorCode B_and_BT_mult(DM dm, Vec xu, Vec xp, Vec yu, Vec yp, std::array<double, 3> const& h)
    {
      PetscErrorCode ierr;
      PetscScalar tmpu[27][3], tmpp[8], vv[8];
      PetscReal MatElemOnu[216], MatElemOnv[216], MatElemOnw[216];
      DM dav, dap;
      DMDALocalInfo infop;
      PetscScalar ****pxu, ****pyu;
      PetscScalar ***pxp, ***pyp;

      PetscFunctionBeginUser;

      ierr = DMCompositeGetEntries(dm, &dav, &dap);CHKERRQ(ierr);

      ierr = getMatElemPressure3d(MatElemOnu, MatElemOnv, MatElemOnw, h);

      auto box = get_DM_bounds<3>(dap);


      ierr = DMDAVecGetArrayDOFRead(dav, xu, &pxu);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayDOF(dav, yu, &pyu);CHKERRQ(ierr);
      ierr = DMDAVecGetArrayRead(dap, xp, &pxp);CHKERRQ(ierr);
      ierr = DMDAVecGetArray(dap, yp, &pyp);CHKERRQ(ierr);

      for (std::size_t k=box.bottom_left[2]; k<box.upper_right[2]; ++k){
        for (std::size_t j=box.bottom_left[1]; j<box.upper_right[1]; ++j){
          for (std::size_t i=box.bottom_left[0]; i<box.upper_right[0]; ++i){
            for(std::size_t k1=0; k1<27; ++k1){
              std::size_t l1 = 2*i + indpu3d[k1][0];
              std::size_t p1 = 2*j + indpu3d[k1][1];
              std::size_t r1 = 2*k + indpu3d[k1][2];
              for(std::size_t k2=0; k2<3; ++k2)
                tmpu[k1][k2] = pxu[r1][p1][l1][k2];
            }

            for(std::size_t k1=0; k1<8; k1++){
              std::size_t l1 = i + ind3d[k1][0];
              std::size_t p1 = j + ind3d[k1][1];
              std::size_t r1 = k + ind3d[k1][2];
              tmpp[k1] = pxp[r1][p1][l1];
              vv[k1] = 0.;
            }

            for(std::size_t k1=0; k1<27; k1++){
              std::size_t l1 = 2*i + indpu3d[k1][0];
              std::size_t p1 = 2*j + indpu3d[k1][1];
              std::size_t r1 = 2*k + indpu3d[k1][2];
              auto val1 = 0;
              auto val2 = 0;
              auto val3 = 0;

              for(std::size_t k2=0; k2<8; k2++){
                val1 -= tmpp[k2]*MatElemOnu[k2*27+k1];
                val2 -= tmpp[k2]*MatElemOnv[k2*27+k1];
                val3 -= tmpp[k2]*MatElemOnw[k2*27+k1];

                vv[k2] += tmpu[k1][0]*MatElemOnu[k2*27+k1]
                         +tmpu[k1][1]*MatElemOnv[k2*27+k1]
                         +tmpu[k1][2]*MatElemOnw[k2*27+k1];
              }

              pyu[r1][p1][l1][0] += val1;
              pyu[r1][p1][l1][1] += val2;
              pyu[r1][p1][l1][2] += val3;
            }

            for(std::size_t k1=0; k1<8; k1++){
              std::size_t l1 = i + ind3d[k1][0];
              std::size_t p1 = j + ind3d[k1][1];
              std::size_t r1 = k + ind3d[k1][2];
              pyp[r1][p1][l1] += vv[k1];
            }
          }
        }
      }

      ierr = DMDAVecRestoreArrayDOFRead(dav, xu, &pxu);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayDOF(dav, yu, &pyu);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArrayRead(dap, xp, &pxp);CHKERRQ(ierr);
      ierr = DMDAVecRestoreArray(dap, yp, &pyp);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  }
}
#endif