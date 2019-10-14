// ./preallocate -nx 10 -ny 10 -stokes_ksp_type gmres -stokes_pc_type none -stokes_ksp_monitor

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

  Mat A, P;
  PreallocateMat(st.ctx, st.opt, MATAIJ, &A, false);
  PreallocateMat(st.ctx, st.opt, MATAIJ, &P, true);

  DM dav, dap;
  ierr = DMCompositeGetEntries(st.ctx->dm, &dav, &dap);CHKERRQ(ierr);

  MatSetDM(A, st.ctx->dm);
  MatSetDM(P, st.ctx->dm);
  MatSetFromOptions(A);
  MatSetFromOptions(P);

  cafes::laplacian_assembling(dav, A, st.ctx->h);
  cafes::B_BT_assembling(st.ctx->dm, A, st.ctx->h);
  cafes::laplacian_assembling(dav, P, st.ctx->h);
  std::array<double, 2> hp; hp[0] = 2*st.ctx->h[0]; hp[1] = 2*st.ctx->h[1];
  cafes::mass_assembling(st.ctx->dm, P, hp);

  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);

  //ierr = MatView(P, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  cafes::fem::SetDirichletOnMat(A, bc);
  cafes::fem::SetDirichletOnMat(P, bc);

  st.setup_RHS();

  KSP ksp;
  KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
  KSPSetDM(ksp, st.ctx->dm);CHKERRQ(ierr);
  KSPSetDMActive(ksp, PETSC_FALSE);CHKERRQ(ierr);
  KSPSetOptionsPrefix(ksp, "stokes_");CHKERRQ(ierr);
  KSPSetOperators(ksp, A, P);CHKERRQ(ierr);
  KSPSetFromOptions(ksp);
  
  KSPSolve(ksp, st.rhs, st.sol);

  ierr = cafes::io::save_VTK("Resultats", "test", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}