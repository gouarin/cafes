#include <cafes.hpp>
#include <petsc.h>

void quadratic_u_2d(const PetscReal x[], PetscScalar *u)
{
    *u = x[0] * x[0] + x[1] * x[1];
}

void quadratic_v_2d(const PetscReal x[], PetscScalar *u)
{
    *u = 2 * x[0] * x[0] - 2 * x[0] * x[1];
}

void three(const PetscReal x[], PetscScalar *u)
{
    *u = 3.;
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 2;

    ierr = PetscInitialize(&argc, &argv, (char *)0, (char *)0);
    CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({
        {{quadratic_u_2d, quadratic_v_2d}}, // left
        {{quadratic_u_2d, quadratic_v_2d}}, // right
        {{quadratic_u_2d, quadratic_v_2d}}, // bottom
        {{quadratic_u_2d, quadratic_v_2d}}  // top
    });

    auto rhs = cafes::make_rhs<dim>({{three, three}});
    auto st = cafes::make_stokes<dim>(bc, rhs);

    st.setup_RHS();
    st.setup_KSP();
    st.solve();

    int const mx = st.opt.mx[0] - 1;
    std::string stout0 = "stokes_analytic_";
    stout0.append(std::to_string(mx));

    ierr = cafes::io::save_hdf5("Resultats", stout0.c_str(), st.sol, st.ctx->dm, st.ctx->h);
    CHKERRQ(ierr);
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;
}