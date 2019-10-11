#include <cafes.hpp>
#include <fem/assembling.hpp>
#include <petsc.h>

void ones(const PetscReal x[], PetscScalar *u){
  *u = 1;
}

void zeros(const PetscReal x[], PetscScalar *u){
  *u = 0;
}

void quadratic_u_2d(const PetscReal x[], PetscScalar *u){
  *u = x[0]*x[0] + x[1]*x[1];
}

void quadratic_v_2d(const PetscReal x[], PetscScalar *u){
  *u = 2*x[0]*x[0] - 2*x[0]*x[1];
}

void three(const PetscReal x[], PetscScalar *u){
  *u = 3.;
}

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  std::size_t const dim = 2;

  ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

//   auto bc = cafes::make_bc<dim>({ {{quadratic_u_2d , quadratic_v_2d}} // left
//                                 , {{quadratic_u_2d , quadratic_v_2d}} // right
//                                 , {{quadratic_u_2d , quadratic_v_2d}} // bottom
//                                 , {{quadratic_u_2d , quadratic_v_2d}} // top
//                                 });
//   auto rhs = cafes::make_rhs<dim>({{ three, three}});

  auto bc = cafes::make_bc<dim>({ {{zeros , zeros}} // left
                                , {{zeros , zeros}} // right
                                , {{zeros , zeros}} // bottom
                                , {{ones , zeros}} // top
                                });

  auto rhs = cafes::make_rhs<dim>({{ zeros, zeros}});
  
  auto st = cafes::make_stokes<dim>(bc, rhs);

  Mat A;
  PreallocateMat(st.ctx, st.opt, MATAIJ, &A, false);

  DM dav, dap;
  ierr = DMCompositeGetEntries(st.ctx->dm, &dav, &dap);CHKERRQ(ierr);

  MatSetDM(A, st.ctx->dm);
  MatSetFromOptions(A);

  cafes::laplacian_assembling(dav, A, st.ctx->h);
  cafes::B_BT_assembling(st.ctx->dm, A, st.ctx->h);

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  cafes::fem::SetDirichletOnMat(A, bc);

  st.setup_RHS();

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  KSPSetDM(ksp, st.ctx->dm);CHKERRQ(ierr);
  KSPSetDMActive(ksp, PETSC_FALSE);CHKERRQ(ierr);
  KSPSetOptionsPrefix(ksp, "stokes_");CHKERRQ(ierr);
  KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
  KSPSetFromOptions(ksp);
  
  KSPSolve(ksp, st.rhs, st.sol);

  ierr = cafes::io::save_VTK("Resultats", "test", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}