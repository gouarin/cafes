#include <cafes.hpp>
#include <petsc.h>

void quadratic_u_3d(const PetscReal x[], PetscScalar *u){
  *u = x[0]*x[0] + x[1]*x[1];
}

void quadratic_v_3d(const PetscReal x[], PetscScalar *u){
  *u = x[1]*x[1] + x[2]*x[2];
}

void quadratic_w_3d(const PetscReal x[], PetscScalar *u){
  *u = x[0]*x[0] + x[1]*x[1] - 2*(x[0] + x[1])*x[2];
}


void three(const PetscReal x[], PetscScalar *u){
  *u = 3.;
}

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  std::size_t const dim = 3;

  ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

  auto bc = cafes::make_bc<dim>({ {{quadratic_u_3d , quadratic_v_3d, quadratic_w_3d}} // left
                                , {{quadratic_u_3d , quadratic_v_3d, quadratic_w_3d}} // right
                                , {{quadratic_u_3d , quadratic_v_3d, quadratic_w_3d}} // bottom
                                , {{quadratic_u_3d , quadratic_v_3d, quadratic_w_3d}}   // top
                                , {{quadratic_u_3d , quadratic_v_3d, quadratic_w_3d}}  // front
                                , {{quadratic_u_3d , quadratic_v_3d, quadratic_w_3d}}  // back
                                });

  auto rhs = cafes::make_rhs<dim>({{ three, three, three }});
  
  auto st = cafes::make_stokes<dim>(bc, rhs);

  st.setup_RHS();
  st.setup_KSP();
  st.solve();

  ierr = cafes::io::save_VTK("Resultats", "test", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}