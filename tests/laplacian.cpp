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

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  std::size_t const dim = 2;

  ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

  auto bc = cafes::make_bc<dim>({ {{quadratic_u_2d}} // left
                                , {{quadratic_u_2d}} // right
                                , {{quadratic_u_2d}} // bottom
                                , {{quadratic_u_2d}} // top
                                });

  auto rhs = cafes::make_rhs<dim>({{ three}});
  
  auto lap = cafes::make_laplacian<dim>(bc, rhs);

  lap.setup_RHS();
  lap.setup_KSP();
  lap.solve();

  ierr = cafes::io::save_VTK("Resultats", "test", lap.sol, lap.ctx->dm, lap.ctx->h);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}