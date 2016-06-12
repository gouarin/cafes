#ifndef PARTICLE_PROBLEM_laplacian_HPP_INCLUDED
#define PARTICLE_PROBLEM_laplacian_HPP_INCLUDED

#include <problem/problem.hpp>
#include <problem/context.hpp>
#include <problem/options.hpp>
#include <fem/matrixFree.hpp>
#include <fem/rhs.hpp>
#include <fem/mesh.hpp>
#include <petsc.h>
#include <iostream>
#include <algorithm>

namespace cafes
{
  namespace problem
  {

    //!
    //! Laplacian class 
    //! 
    template<std::size_t Dimensions>
    struct laplacian : public Problem<Dimensions>{
      //using parent = Problem<Dimensions>::Problem;
      using Ctx = context<Dimensions>;
      fem::rhs_conditions<Dimensions> rhsc_;
      options<Dimensions> opt{};

      Ctx *ctx; //!< Context for the PETSc solver
      Vec sol;  //!< The solution of Laplacian problem
      Vec rhs;  //!< The RHS of Laplacian problem
      Mat A;    //!< The matrix of Laplacian problem
      Mat P;    //!< The preconditioner of Laplacian problem
      KSP ksp;  //!< The solver of Laplacian problem

      laplacian(fem::dirichlet_conditions<Dimensions> bc, fem::rhs_conditions<Dimensions> rhsc={nullptr})
      {
        opt.process_options();

        DM mesh;
        fem::createLaplacianMesh<Dimensions>(mesh, opt.mx, opt.xperiod);

        DMCreateGlobalVector(mesh, &sol);
        VecDuplicate(sol, &rhs);
        VecSet(rhs, 0.);

        std::array<double, Dimensions> h;
        for(std::size_t i = 0; i<Dimensions; ++i)
          h[i] = opt.lx[i]/(opt.mx[i]-1);

        // set Laplacian matrix
        PetscErrorCode(*method)(DM, Vec, Vec, std::array<double, Dimensions> const&);
        if (opt.strain_tensor){
          method = fem::strain_tensor_mult;
        }
        else
          method = fem::laplacian_mult;

        // fix this to avoid raw pointer !!
        rhsc_ = rhsc;
        ctx = new Ctx{mesh, h, method, fem::laplacian_mult_diag};
        ctx->set_dirichlet_bc(bc);
        A = fem::make_matrix<Ctx>(ctx);
        MatSetDM(A, mesh);
        MatSetFromOptions(A);
        DMSetMatType(mesh, MATSHELL);
      }

      #undef __FUNCT__
      #define __FUNCT__ "setup_RHS"
      virtual PetscErrorCode setup_RHS() override 
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;

        if (rhsc_.has_condition()){
          ierr = fem::set_rhs<Dimensions, 1>(ctx->dm, rhs, rhsc_, ctx->h);CHKERRQ(ierr);
        }
        
        ierr = SetDirichletOnRHS(ctx->dm, ctx->bc_, rhs, ctx->h);CHKERRQ(ierr);
        
        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "setup_KSP"
      virtual PetscErrorCode setup_KSP() override
      {
        PetscErrorCode ierr;
        PC pc;
        PetscFunctionBeginUser;

        ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
        ierr = KSPSetDM(ksp, ctx->dm);CHKERRQ(ierr);
        ierr = KSPSetDMActive(ksp, PETSC_FALSE);CHKERRQ(ierr);
        ierr = KSPSetOptionsPrefix(ksp, "laplacian_");CHKERRQ(ierr);

        ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);

        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPSetUp(ksp);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "solve"
      virtual PetscErrorCode solve() override
      {
        PetscErrorCode ierr;
        PC pc;
        PetscFunctionBeginUser;

        ierr = KSPSolve(ksp, rhs, sol);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }
    };
  } 

  template<std::size_t Dimensions, std::size_t Ndof>
  problem::laplacian<Dimensions> make_laplacian(fem::dirichlet_conditions<Dimensions> const& dc, fem::rhs_conditions<Ndof> const& rhs)
  {
    return {dc, rhs}; 
  }

}

#endif