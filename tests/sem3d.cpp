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
    std::size_t const dim = 3;

    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({ {{zeros, zeros, zeros}} 
                                  , {{zeros, zeros, zeros}}
                                  , {{zeros, zeros, zeros}}
                                  , {{zeros, zeros, zeros}}
                                  , {{zeros, zeros, zeros}}
                                  , {{zeros, zeros, zeros}}
                                  });
    auto rhs = cafes::make_rhs<dim>({{ ones, zeros, zeros }});
    auto st = cafes::make_stokes<dim>(bc, rhs);

    auto se = cafes::make_ellipsoid<dim>( {.5,.5, .5}, {.1,.1, .1}, 1., 1.);
    std::vector<cafes::particle<decltype(se)>> pt{ cafes::make_particle(se, {0.,150., 1.}, 0.25)
                                                 
                                                 };

    auto s = cafes::make_SEM(pt, st, {{.1, .1}});
    ierr = s.create_Mat_and_Vec();CHKERRQ(ierr);
    s.setup_RHS();
    s.setup_KSP();
    s.solve();

    // PetscViewer      viewer;
    // ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "solution.h5", FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    // ierr = VecView(st.sol, viewer);CHKERRQ(ierr);
    // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}