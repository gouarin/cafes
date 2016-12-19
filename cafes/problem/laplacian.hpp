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

#ifndef PARTICLE_PROBLEM_laplacian_HPP_INCLUDED
#define PARTICLE_PROBLEM_laplacian_HPP_INCLUDED

#include <problem/problem.hpp>
#include <problem/context.hpp>
#include <problem/options.hpp>
#include <fem/matrixFree.hpp>
#include <problem/stokes.hpp>

#include <fem/rhs.hpp>
#include <fem/mesh.hpp>
#include <petsc.h>
#include <iostream>
#include <algorithm>

namespace cafes
{
  namespace problem
  {

    #undef __FUNCT__
    #define __FUNCT__ "setMGSolver"
    PetscErrorCode setMGSolver(KSP ksp, DM dm){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;
      PC pc;
      PetscInt MGlevels;
      DMDALocalInfo  info;

      PetscFunctionBeginUser;

      ierr = DMDAGetLocalInfo(dm, &info);CHKERRQ(ierr);

      /* Set MG solver on velocity field*/
      ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);

      ierr = KSPSetType(ksp, KSPFGMRES);CHKERRQ(ierr);      
      ierr = PCSetType(pc, PCMG);CHKERRQ(ierr);
      
      auto i = (info.mx<info.my)? info.mx: info.my;
      i = (info.mz == 1 || i<info.mz)? i: info.mz;

      MGlevels = 1;
      while(i > 8){
        i >>= 1;
        MGlevels++;
      }

      ierr = PCMGSetLevels(pc, MGlevels, PETSC_NULL);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }

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
        PetscErrorCode(*method)(petsc::petsc_vec<Dimensions>&,
                                petsc::petsc_vec<Dimensions>&,
                                std::array<double, Dimensions> const&);
        method = fem::laplacian_mult;

        // fix this to avoid raw pointer !!
        rhsc_ = rhsc;
        ctx = new Ctx{mesh, h, method, fem::diag_laplacian_mult};
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
        
        auto petsc_rhs = petsc::petsc_vec<Dimensions>(ctx->dm, rhs, 0);
        ierr = SetDirichletOnRHS(petsc_rhs, ctx->bc_, ctx->h);CHKERRQ(ierr);
        
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
        //ierr = KSPSetType(ksp, KSPGCR);CHKERRQ(ierr);
        ierr = KSPSetDMActive(ksp, PETSC_FALSE);CHKERRQ(ierr);
        ierr = KSPSetOptionsPrefix(ksp, "laplacian_");CHKERRQ(ierr);

        ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);

        ierr = setMGSolver(ksp, ctx->dm);CHKERRQ(ierr);

        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
        ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);

        PetscBool same;
        ierr = PetscObjectTypeCompare((PetscObject)pc, PCMG, &same);CHKERRQ(ierr);

        // if MG is set
        if (same) {
          std::cout << "set MG !!!!!\n";

          PetscErrorCode(*method)(petsc::petsc_vec<Dimensions>&,
                                  petsc::petsc_vec<Dimensions>&,
                                  std::array<double, Dimensions> const&);
          
          method = fem::laplacian_mult;

          PetscInt MGlevels;
          KSP smoother;
          ierr = PCMGGetLevels(pc, &MGlevels);CHKERRQ(ierr);

          for(std::size_t i=0; i<MGlevels; ++i){
            auto mg_h = ctx->h;
            std::for_each(mg_h.begin(), mg_h.end(), [&](auto& x){x*=(1<<(MGlevels-1-i));});

            auto *mg_ctx = new Ctx{ctx->dm, mg_h, method, fem::diag_laplacian_mult};
            mg_ctx->set_dirichlet_bc(ctx->bc_);

            PC pcsmoother;
            ierr = PCMGGetSmoother(pc, i, &smoother);CHKERRQ(ierr);

            ierr = KSPSetComputeOperators(smoother, createLevelMatrices<Ctx>, (void *) mg_ctx);CHKERRQ(ierr);
            ierr = KSPSetType(smoother, KSPCG);CHKERRQ(ierr);
            ierr = KSPGetPC(smoother, &pcsmoother);CHKERRQ(ierr);
            ierr = PCSetType(pcsmoother, PCJACOBI);CHKERRQ(ierr);
          }

          ierr = PCSetUp(pc);CHKERRQ(ierr);
        }
        ierr = KSPSetUp(ksp);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "solve"
      virtual PetscErrorCode solve() override
      {
        PetscErrorCode ierr;
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

  template<std::size_t Dimensions>
  problem::laplacian<Dimensions> make_laplacian(fem::dirichlet_conditions<Dimensions> const& dc)
  {
    return {dc}; 
  }

}

#endif