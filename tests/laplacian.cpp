#include <cafes.hpp>
#include <petsc.h>

void quadratic_u_2d(const PetscReal x[], PetscScalar *u){
  *u = x[0]*x[0] + x[1]*x[1];
}

void quadratic_v_2d(const PetscReal x[], PetscScalar *u){
  *u = 2*x[0]*x[0] - 2*x[0]*x[1];
}

void three(const PetscReal x[], PetscScalar *u){
  *u = 3.;
}

void ones(const PetscReal x[], PetscScalar *u){
  *u = 1.;
}

void zeros(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  std::size_t const dim = 2;

  ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

  // auto bc = cafes::make_bc<dim>({ {{quadratic_u_2d}} // left
  //                               , {{quadratic_u_2d}} // right
  //                               , {{quadratic_u_2d}} // bottom
  //                               , {{quadratic_u_2d}} // top
  //                               });

  // auto rhs = cafes::make_rhs<dim>({{three}});

  auto bc = cafes::make_bc<dim>({ {{ones}} // left
                                , {{zeros}} // right
                                , {{zeros}} // bottom
                                , {{zeros}} // top
                                });

  auto rhs = cafes::make_rhs<dim>({{zeros}});
  //auto rhs = cafes::make_rhs<dim>({{nullptr}}); // doesn't work but it should !!!
  auto lap = cafes::make_laplacian<dim>(bc, rhs);
  //auto lap = cafes::make_laplacian<dim>(bc);

  lap.setup_RHS();
  lap.setup_KSP();
  lap.solve();

  ierr = cafes::io::save_VTK_laplacian("Resultats", "rhs", lap.rhs, lap.ctx->dm, lap.ctx->h);CHKERRQ(ierr);
  ierr = cafes::io::save_VTK_laplacian("Resultats", "test", lap.sol, lap.ctx->dm, lap.ctx->h);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}