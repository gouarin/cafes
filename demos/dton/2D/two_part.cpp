#include <cafes.hpp>
#include <petsc.h>

#include <particle/singularity/add_singularity.hpp>

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
    std::size_t const dim = 2;

    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({ {{zeros , zeros}} // gauche
                                  , {{zeros , zeros}} // droite
                                  , {{zeros, zeros}} // bas
                                  , {{zeros, zeros}}   // haut
                                  });

    auto rhs = cafes::make_rhs<dim>({{ zeros, zeros }});
    
    auto st = cafes::make_stokes<dim>(bc, rhs);

    double R1   = .1;
    double R2   = .1;
    double dist = .1/6;

    auto se1 = cafes::make_circle( { .5-R1-dist/2, .5}, R1, 0);
    auto se2 = cafes::make_circle( { .5+R2+dist/2, .5}, R2, 0);
    std::vector<cafes::particle<decltype(se1)>> pt{ cafes::make_particle(se1, { 1., 0.}, {0.,0.,0.}),
	                                            cafes::make_particle(se2, {-1., 0.}, {0.,0.,0.})};

    auto s = cafes::make_DtoN(pt, st, .1);
    
    ierr = s.create_Mat_and_Vec();CHKERRQ(ierr);
    ierr = cafes::singularity::save_singularity<dim, decltype(*(s.ctx))>("Resultats", "two_part", *(s.ctx));CHKERRQ(ierr);

    ierr = s.setup_RHS();CHKERRQ(ierr);
    ierr = s.setup_KSP();CHKERRQ(ierr);
    ierr = s.solve();CHKERRQ(ierr);

    //ierr = cafes::io::save_VTK("Resultats", "two_part", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
    ierr = cafes::io::saveParticles("Resultats", "two_part", pt);CHKERRQ(ierr);

    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}
