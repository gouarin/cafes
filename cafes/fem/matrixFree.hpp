#ifndef CAFES_FEM_MATRIXFREE_HPP_INCLUDED
#define CAFES_FEM_MATRIXFREE_HPP_INCLUDED

#include <fem/operator.hpp>
#include <iostream>
#include <petsc.h>

namespace cafes
{
  namespace fem
  {
    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix"    
    template<typename CTX>
    typename std::enable_if<CTX::ndm_type::value==1, PetscErrorCode>::type
    diag_block_matrix(Mat A, Vec x, Vec y)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);

      Vec xlocal, ylocal;
      ierr = DMGetLocalVector(ctx->dm, &xlocal);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(ctx->dm, x, INSERT_VALUES, xlocal);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(ctx->dm, x, INSERT_VALUES, xlocal);CHKERRQ(ierr);

      ierr = VecSet(y, 0.);CHKERRQ(ierr);
      ierr = DMGetLocalVector(ctx->dm, &ylocal);CHKERRQ(ierr);
      ierr = VecSet(ylocal, 0.);CHKERRQ(ierr);

      ierr = ctx->apply(ctx->dm, xlocal, ylocal, ctx->h);CHKERRQ(ierr);

      ierr = DMLocalToGlobalBegin(ctx->dm, ylocal, ADD_VALUES, y);CHKERRQ(ierr);
      ierr = DMLocalToGlobalEnd(ctx->dm, ylocal, ADD_VALUES, y);CHKERRQ(ierr);

      ierr = DMRestoreLocalVector(ctx->dm, &xlocal);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(ctx->dm, &ylocal);CHKERRQ(ierr);

      if(ctx->set_bc_){
         ierr = SetDirichletOnVec(ctx->dm, ctx->bc_, x, y);CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix"    
    template<typename CTX>
    typename std::enable_if<(CTX::ndm_type::value>1), PetscErrorCode>::type
    diag_block_matrix(Mat A, Vec x, Vec y)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);

      int ndm;
      ierr = DMCompositeGetNumberDM(ctx->dm, &ndm);CHKERRQ(ierr);

      DM dm[ndm];
      ierr = DMCompositeGetEntriesArray(ctx->dm, dm);CHKERRQ(ierr);
        
      Vec xdm[ndm], ydm[ndm];
      ierr = DMCompositeGetAccessArray(ctx->dm, x, ndm, NULL, xdm);CHKERRQ(ierr);
      ierr = DMCompositeGetAccessArray(ctx->dm, y, ndm, NULL, ydm);CHKERRQ(ierr);

      ierr = VecSet(y, 0.);CHKERRQ(ierr);

      for(std::size_t i=0; i<ndm; ++i){
        Vec xlocal, ylocal;
        ierr = DMGetLocalVector(dm[i], &xlocal);CHKERRQ(ierr);
        ierr = DMGlobalToLocalBegin(dm[i], xdm[i], INSERT_VALUES, xlocal);CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(dm[i], xdm[i], INSERT_VALUES, xlocal);CHKERRQ(ierr);
   
        ierr = DMGetLocalVector(dm[i], &ylocal);CHKERRQ(ierr);
        ierr = VecSet(ylocal, 0.);CHKERRQ(ierr);

        ierr = ctx->apply(dm[i], xlocal, ylocal, ctx->h[i]);CHKERRQ(ierr);

        ierr = DMLocalToGlobalBegin(dm[i], ylocal, ADD_VALUES, ydm[i]);CHKERRQ(ierr);
        ierr = DMLocalToGlobalEnd(dm[i], ylocal, ADD_VALUES, ydm[i]);CHKERRQ(ierr);

        ierr = DMRestoreLocalVector(dm[i], &xlocal);CHKERRQ(ierr);
        ierr = DMRestoreLocalVector(dm[i], &ylocal);CHKERRQ(ierr);
      }

      ierr = DMCompositeRestoreAccessArray(ctx->dm, x, ndm, NULL, xdm);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccessArray(ctx->dm, y, ndm, NULL, ydm);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix_diag"
    template<typename CTX>
    typename std::enable_if<CTX::ndm_type::value==1, PetscErrorCode>::type
    diag_block_matrix_diag(Mat A, Vec x)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);

      ierr = VecSet(x, 0.);CHKERRQ(ierr);

      Vec xlocal;
      ierr = DMGetLocalVector(ctx->dm, &xlocal);CHKERRQ(ierr);
      ierr = VecSet(xlocal, 0.);CHKERRQ(ierr);

      ierr = ctx->apply_diag(ctx->dm, xlocal, ctx->h);CHKERRQ(ierr);

      ierr = DMLocalToGlobalBegin(ctx->dm, xlocal, ADD_VALUES, x);CHKERRQ(ierr);
      ierr = DMLocalToGlobalEnd(ctx->dm, xlocal, ADD_VALUES, x);CHKERRQ(ierr);

      ierr = DMRestoreLocalVector(ctx->dm, &xlocal);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix_diag"
    template<typename CTX>
    typename std::enable_if<(CTX::ndm_type::value>1), PetscErrorCode>::type
    diag_block_matrix_diag(Mat A, Vec x)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);

      ierr = VecSet(x, 0.);CHKERRQ(ierr);

      int ndm;
      ierr = DMCompositeGetNumberDM(ctx->dm, &ndm);CHKERRQ(ierr);

      DM dm[ndm];
      ierr = DMCompositeGetEntriesArray(ctx->dm, dm);CHKERRQ(ierr);
        
      Vec xdm[ndm];
      ierr = DMCompositeGetAccessArray(ctx->dm, x, ndm, NULL, xdm);CHKERRQ(ierr);

      for(std::size_t i=0; i<ndm; ++i){
        Vec xlocal;
        ierr = DMGetLocalVector(dm[i], &xlocal);CHKERRQ(ierr);
        ierr = VecSet(xlocal, 0.);CHKERRQ(ierr);
   
        ierr = ctx->apply_diag(dm[i], xlocal, ctx->h[i]);CHKERRQ(ierr);

        ierr = DMLocalToGlobalBegin(dm[i], xlocal, ADD_VALUES, xdm[i]);CHKERRQ(ierr);
        ierr = DMLocalToGlobalEnd(dm[i], xlocal, ADD_VALUES, xdm[i]);CHKERRQ(ierr);

        ierr = DMRestoreLocalVector(dm[i], &xlocal);CHKERRQ(ierr);
      }

      ierr = DMCompositeRestoreAccessArray(ctx->dm, x, ndm, NULL, xdm);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "stokes_matrix"
    template<typename CTX>
    PetscErrorCode stokes_matrix(Mat A, Vec x, Vec y)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      DM dav, dap;
      Vec xu, xp, xulocal, xplocal;
      Vec yu, yp, yulocal, yplocal;

      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, (void**) &ctx);CHKERRQ(ierr);

      ierr = DMCompositeGetEntries(ctx->dm, &dav, &dap);CHKERRQ(ierr);

      ierr = DMCompositeGetAccess(ctx->dm, x, &xu, &xp);CHKERRQ(ierr);

      ierr = DMGetLocalVector(dav, &xulocal);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(dav, xu, INSERT_VALUES, xulocal);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dav, xu, INSERT_VALUES, xulocal);CHKERRQ(ierr);

      ierr = DMGetLocalVector(dap, &xplocal);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(dap, xp, INSERT_VALUES, xplocal);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dap, xp, INSERT_VALUES, xplocal);CHKERRQ(ierr);

      ierr = VecSet(y, 0.);CHKERRQ(ierr);
      ierr = DMGetLocalVector(dav, &yulocal);CHKERRQ(ierr);
      ierr = DMGetLocalVector(dap, &yplocal);CHKERRQ(ierr);
      ierr = VecSet(yulocal, 0.);CHKERRQ(ierr);
      ierr = VecSet(yplocal, 0.);CHKERRQ(ierr);

      ierr = ctx->apply(dav, xulocal, yulocal, ctx->h);CHKERRQ(ierr);
      ierr = B_and_BT_mult(ctx->dm, xulocal, xplocal, yulocal, yplocal, ctx->h);CHKERRQ(ierr);

      ierr = DMCompositeGetAccess(ctx->dm, y, &yu, &yp);CHKERRQ(ierr);

      ierr = DMLocalToGlobalBegin(dav, yulocal, ADD_VALUES, yu);CHKERRQ(ierr);
      ierr = DMLocalToGlobalEnd(dav, yulocal, ADD_VALUES, yu);CHKERRQ(ierr);

      ierr = DMLocalToGlobalBegin(dap, yplocal, ADD_VALUES, yp);CHKERRQ(ierr);
      ierr = DMLocalToGlobalEnd(dap, yplocal, ADD_VALUES, yp);CHKERRQ(ierr);

      if (ctx->set_bc_){
        ierr = SetDirichletOnVec(dav, ctx->bc_, xu, yu);CHKERRQ(ierr);
      }

      ierr = DMCompositeRestoreAccess(ctx->dm, x, &xu, &xp);CHKERRQ(ierr);
      ierr = DMCompositeRestoreAccess(ctx->dm, y, &yu, &yp);CHKERRQ(ierr);

      ierr = DMRestoreLocalVector(dav, &xulocal);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dav, &yulocal);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dap, &xplocal);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dap, &yplocal);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    template<typename CTX>
    Mat make_matrix(CTX *ctx, 
                    PetscErrorCode(*method)(Mat, Vec, Vec) = diag_block_matrix<CTX>, 
                    PetscErrorCode(*method_diag)(Mat, Vec) = diag_block_matrix_diag<CTX>)
    {
      Mat A;
      int localsize, totalsize;
    
      get_DM_sizes(ctx->dm, localsize, totalsize);
      MatCreateShell(PETSC_COMM_WORLD, localsize, localsize, totalsize, totalsize, ctx, &A);
      MatShellSetOperation(A, MATOP_MULT, (void(*)())(*method));
      if (ctx->apply_diag)
        MatShellSetOperation(A, MATOP_GET_DIAGONAL, (void(*)())(*method_diag));

      MatSetFromOptions(A);

      return A;
    }
  }
}
#endif