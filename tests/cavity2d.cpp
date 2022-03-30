#include <cafes.hpp>
#include <petsc.h>

void zeros(const PetscReal x[], PetscScalar *u)
{
    *u = 0.;
}

void ones(const PetscReal x[], PetscScalar *u)
{
    *u = 1.;
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 2;

    ierr = PetscInitialize(&argc, &argv, (char *)0, (char *)0);CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({
        {{zeros, zeros}}, // left
        {{zeros, zeros}}, // right
        {{zeros, zeros}}, // bottom
        {{ones, zeros}}, // top
    });

    auto rhs = cafes::make_rhs<dim>({{zeros, zeros}});
    auto st = cafes::make_stokes<dim>(bc, rhs);

    st.setup_RHS();
    st.setup_KSP();
    st.solve();

    int const mx = st.opt.mx[0] - 1;
    std::string stout0 = "cavity_";
    stout0.append(std::to_string(mx));
    ierr = cafes::io::save_hdf5("Resultats", stout0.c_str(), st.sol, st.ctx->dm, st.ctx->h);
    CHKERRQ(ierr);
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;
}