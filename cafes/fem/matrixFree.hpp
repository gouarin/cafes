#ifndef CAFES_FEM_MATRIXFREE_HPP_INCLUDED
#define CAFES_FEM_MATRIXFREE_HPP_INCLUDED

#include <fem/operator.hpp>
#include <petsc/vec.hpp>
#include <iostream>
#include <array>
#include <petsc.h>

namespace cafes
{
  namespace fem
  {
    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix"
    template<typename CTX, typename Dimensions = typename CTX::dimension_type>
    typename std::enable_if<CTX::ndm_type::value==1, PetscErrorCode>::type
    diag_block_matrix(Mat A, Vec x, Vec y)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);
      ierr = VecSet(y, 0.);CHKERRQ(ierr);

      auto xpetsc = petsc::petsc_vec<Dimensions::value>(ctx->dm, x, 0);
      auto ypetsc = petsc::petsc_vec<Dimensions::value>(ctx->dm, y, 0, false);

      ierr = xpetsc.global_to_local(INSERT_VALUES);CHKERRQ(ierr);
      ierr = ypetsc.fill(0.);CHKERRQ(ierr);

      ierr = ctx->apply(xpetsc, ypetsc, ctx->h);CHKERRQ(ierr);

      ierr = ypetsc.local_to_global(ADD_VALUES);CHKERRQ(ierr);

      if(ctx->set_bc_){
         ierr = SetDirichletOnVec(xpetsc, ypetsc, ctx->bc_);CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix"
    template<typename CTX, typename Dimensions = typename CTX::dimension_type>
    typename std::enable_if<(CTX::ndm_type::value>1), PetscErrorCode>::type
    diag_block_matrix(Mat A, Vec x, Vec y)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);
      ierr = VecSet(y, 0.);CHKERRQ(ierr);

      int ndm;
      ierr = DMCompositeGetNumberDM(ctx->dm, &ndm);CHKERRQ(ierr);

      for(std::size_t i=0; i<ndm; ++i){
        auto xpetsc = petsc::petsc_vec<Dimensions::value>(ctx->dm, x, i);
        auto ypetsc = petsc::petsc_vec<Dimensions::value>(ctx->dm, y, i, false);

        ierr = xpetsc.global_to_local(INSERT_VALUES);CHKERRQ(ierr);
        ierr = ypetsc.fill(0.);CHKERRQ(ierr);

        ierr = ctx->apply(xpetsc, ypetsc, ctx->h[i]);CHKERRQ(ierr);

        ierr = ypetsc.local_to_global(ADD_VALUES);CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }


    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix_diag"
    template<typename CTX, typename Dimensions = typename CTX::dimension_type>
    typename std::enable_if<CTX::ndm_type::value==1, PetscErrorCode>::type
    diag_block_matrix_diag(Mat A, Vec x)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);

      auto xpetsc = petsc::petsc_vec<Dimensions::value>(ctx->dm, x, 0, false);
      ierr = xpetsc.fill(0.);CHKERRQ(ierr);

      ierr = ctx->apply_diag(xpetsc, ctx->h);CHKERRQ(ierr);

      ierr = xpetsc.fill_global(0.);CHKERRQ(ierr);
      ierr = xpetsc.local_to_global(ADD_VALUES);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "diag_block_matrix_diag"
    template<typename CTX, typename Dimensions = typename CTX::dimension_type>
    typename std::enable_if<(CTX::ndm_type::value>1), PetscErrorCode>::type
    diag_block_matrix_diag(Mat A, Vec x)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);

      int ndm;
      ierr = DMCompositeGetNumberDM(ctx->dm, &ndm);CHKERRQ(ierr);

      for(std::size_t i=0; i<ndm; ++i){
        auto xpetsc = petsc::petsc_vec<Dimensions::value>(ctx->dm, x, i, false);
        ierr = xpetsc.fill(0.);CHKERRQ(ierr);

        ierr = ctx->apply_diag(xpetsc, ctx->h[i]);CHKERRQ(ierr);

        ierr = xpetsc.fill_global(0.);CHKERRQ(ierr);
        ierr = xpetsc.local_to_global(ADD_VALUES);CHKERRQ(ierr);
      }

      PetscFunctionReturn(0);
    }

    #undef __FUNCT__
    #define __FUNCT__ "stokes_matrix"
    template<typename CTX, typename Dimensions = typename CTX::dimension_type>
    PetscErrorCode stokes_matrix(Mat A, Vec x, Vec y)
    {
      CTX *ctx;
      PetscErrorCode ierr;
      DM dav, dap;
      PetscFunctionBeginUser;

      ierr = MatShellGetContext(A, (void**) &ctx);CHKERRQ(ierr);
      ierr = VecSet(y, 0.);CHKERRQ(ierr);

      using petsc_type = petsc::petsc_vec<Dimensions::value>;

      std::array<petsc_type, 2> xpetsc{petsc_type(ctx->dm, x, 0),
                                       petsc_type(ctx->dm, x, 1)};

      std::array<petsc_type, 2> ypetsc{petsc_type(ctx->dm, y, 0, false),
                                       petsc_type(ctx->dm, y, 1, false)};

      for(std::size_t i=0; i<xpetsc.size(); ++i)
      {
        ierr = xpetsc[i].global_to_local(INSERT_VALUES);CHKERRQ(ierr);
        ierr = ypetsc[i].fill(0.);CHKERRQ(ierr);
      }

      ierr = ctx->apply(xpetsc[0], ypetsc[0], ctx->h);CHKERRQ(ierr);
      ierr = B_and_BT_mult(xpetsc[0], xpetsc[1], ypetsc[0], ypetsc[1], ctx->h);CHKERRQ(ierr);

      for(std::size_t i=0; i<xpetsc.size(); ++i)
      {
        ierr = ypetsc[i].local_to_global(ADD_VALUES);CHKERRQ(ierr);
      }

      if (ctx->set_bc_){
        ierr = SetDirichletOnVec(xpetsc[0], ypetsc[0], ctx->bc_);CHKERRQ(ierr);
      }

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

    #undef __FUNCT__
    #define __FUNCT__ "apply_mass_matrix"
    template<std::size_t Dimensions>
    PetscErrorCode apply_mass_matrix(DM dm, Vec x, Vec y, std::array<double, Dimensions> const& h)
    {
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      auto hp{h};
      std::for_each(hp.begin(), hp.end(), [](auto x){x*=2;});

      using ctx = problem::context<Dimensions, 2>;
      ctx s{dm, {{h, hp }}, mass_mult};
      Mat mass = make_matrix<ctx>(&s, diag_block_matrix<ctx>);

      ierr = MatMult(mass, x, y);CHKERRQ(ierr);
      ierr = MatDestroy(&mass);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

  }
}
#endif