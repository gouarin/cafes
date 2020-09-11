#include <cafes.hpp>
#include <petsc.h>

#include <particle/singularity/add_singularity.hpp>

void zeros(const PetscReal x[], PetscScalar *u)
{
    *u = 0.;
}

void ones(const PetscReal x[], PetscScalar *u)
{
    *u = 1.;
}

void ones_m(const PetscReal x[], PetscScalar *u)
{
    *u = -1.;
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 2;

    ierr = PetscInitialize(&argc, &argv, (char *)0, (char *)0);
    CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({
        {{zeros, zeros}} // gauche
        ,
        {{zeros, zeros}} // droite
        ,
        {{zeros, zeros}} // bas
        ,
        {{zeros, zeros}} // haut
    });

    auto rhs = cafes::make_rhs<dim>({{zeros, zeros}});

    auto st = cafes::make_stokes<dim>(bc, rhs);
    int const mx = st.opt.mx[0] - 1;

    // auto c = cafes::myparticle({.5,.5}, {1.,0.}, 0.1);
    // c.discretize_surface(10);

    double R1 = .1;
    double distance = .15;

    auto se1 = cafes::make_circle({.5 - .5*(distance+R1), .5}, R1, 0);
    auto se2 = cafes::make_circle({.5 + .5*(distance+R1), .5}, R1, 0);
    std::vector<cafes::particle<decltype(se1)>> pt{
        cafes::make_particle_with_velocity(se1, {1., 0.}, 0.),
        cafes::make_particle_with_velocity(se2, {-1., 0.}, 0.)
    };

    auto s = cafes::make_DtoN(pt, st, st.ctx->h[0]);

    ierr = s.create_Mat_and_Vec();
    CHKERRQ(ierr);

    ierr = s.setup_RHS();
    CHKERRQ(ierr);
    ierr = s.setup_KSP();
    CHKERRQ(ierr);
    ierr = s.solve();
    CHKERRQ(ierr);

    std::string stout = "two_parts_solution_";
    stout.append(std::to_string(mx));
    stout.append("_distance_");
    stout.append(std::to_string(distance));
    const char * stw = stout.c_str();
    ierr = cafes::io::save_hdf5("Resultats", stw, st.sol, st.ctx->dm,
                              st.ctx->h);
    CHKERRQ(ierr);

    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;
}
