#include <problem/stokes.hpp>
#include <fem/bc.hpp>

void zeros(const PetscReal x[], PetscScalar *u){
  *u = 0.;
}

void ones(const PetscReal x[], PetscScalar *u){
  *u = 1.;
}

auto make_stokes(std::integral_constant<std::size_t, 2> const&){
  cafes::fem::rhs_conditions<2>       rhs{{ ones, zeros }};
  cafes::fem::dirichlet_conditions<2> dc = { {{zeros , zeros}} 
                                           , {{zeros , zeros}}
                                           , {{zeros, zeros}}
                                           , {{zeros, zeros}}
                                           };

  cafes::problem::stokes<2> st{dc, rhs}; 
  return st;
}

auto make_stokes(std::integral_constant<std::size_t, 3> const&){
  cafes::fem::rhs_conditions<3>       rhs{{ zeros, zeros, zeros }};
  cafes::fem::dirichlet_conditions<3> dc = { {{ones , zeros, zeros}} 
                                           , {{ones  , zeros, zeros}}
                                           , {{zeros, zeros, zeros}}
                                           , {{zeros, zeros, zeros}}
                                           };

  cafes::problem::stokes<3> st{dc, rhs}; 
  return st;
}

template<std::size_t Dimensions>
auto make_stokes(){
 return make_stokes(std::integral_constant<std::size_t, Dimensions>{});
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 2;

    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

    auto st = make_stokes<dim>();
    st.setup_RHS();
    st.setup_KSP();
    st.solve();

    // PetscViewer      viewer;
    // ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, "solution.h5", FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    // ierr = VecView(st.sol, viewer);CHKERRQ(ierr);
    // ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}