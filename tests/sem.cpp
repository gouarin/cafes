#include <cafes.hpp>
#include <petsc.h>

void zeros(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

void ones(const PetscReal x[], PetscScalar *u){
  *u = 1.;
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 2;

    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({ {{zeros , zeros}} 
                                  , {{zeros , zeros}}
                                  , {{zeros, zeros}}
                                  , {{zeros, zeros}}
                                  });
    auto rhs = cafes::make_rhs<dim>({{ zeros, zeros }});
    auto st = cafes::make_stokes<dim>(bc, rhs);

    auto se1 = cafes::make_ellipsoid<2>( {.3,.5}, {.1,.2}, 1. );
    auto se2 = cafes::make_ellipsoid<2>( {.7,.5}, {.1,.1}, 2. );
    std::vector<cafes::particle<decltype(se1)>> pt{ cafes::make_particle(se1, {1000.,0.}, 1.),
                                                   cafes::make_particle(se2, {1000.,0.}, 1.)
                                                 };

    auto s = cafes::make_SEM(pt, st, .1);
    ierr = s.create_Mat_and_Vec();CHKERRQ(ierr);
    ierr = st.setup_KSP();CHKERRQ(ierr);
    ierr = s.setup_RHS();CHKERRQ(ierr);
    ierr = s.setup_KSP();CHKERRQ(ierr);
    ierr = s.solve();CHKERRQ(ierr);

    ierr = cafes::io::save_VTK("Resultats", "test", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
    // PetscViewer      viewer;
    // ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "solution.h5", FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    // ierr = VecView(st.sol, viewer);CHKERRQ(ierr);
    // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}