#include <cafes.hpp>
#include <petsc.h>

void zeros(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

void ones(const PetscReal x[], PetscScalar *u){
  *u = 1.;
}

void ones_m(const PetscReal x[], PetscScalar *u){
  *u = -1.;
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 3;

    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({ {{zeros , zeros, zeros}} // left
                                  , {{zeros , zeros, zeros}} // right
                                  , {{ones_m, zeros, zeros}} // bottom
                                  , {{ones  , zeros, zeros}}   // top
                                  , {{zeros , zeros, zeros}}  // front
                                  , {{zeros , zeros, zeros}}  // back
                                  });

    auto rhs = cafes::make_rhs<dim>({{ zeros, zeros, zeros }});
    
    auto st = cafes::make_stokes<dim>(bc, rhs);

    double R1   = 0.05;

    auto se1 = cafes::make_sphere( {.5, .5, .5}, R1, 0);

    std::vector<cafes::particle<decltype(se1)>> pt{ cafes::make_particle(se1, {0., 0., 0.}, 0.) };

    auto s = cafes::make_NtoD(pt, st, {.1, .1});

    ierr = s.create_Mat_and_Vec();CHKERRQ(ierr);
    ierr = s.setup_RHS();CHKERRQ(ierr);
    ierr = s.setup_KSP();CHKERRQ(ierr);
    ierr = s.solve();CHKERRQ(ierr);

    ierr = cafes::io::save_hdf5("Resultats", "test_avec_sing", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
    ierr = cafes::io::saveParticles("Resultats", "test_avec_sing", pt);CHKERRQ(ierr);

    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}
