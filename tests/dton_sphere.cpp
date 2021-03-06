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

    // double radius = .1;
    // double z = 4*radius;
    // double y = 4*radius;
    // double x = 4*radius;

    double R1   = 0.1;
    double R2   = .08;
    double dist = .08/12;


    auto se1 = cafes::make_sphere( {.5-R1-dist, .5, .5}, R1, 0);
    auto se2 = cafes::make_sphere( {.5+R2+dist, .5, .5}, R2, 0);
    // std::vector<cafes::particle<decltype(se)>> pt;
    // while(z<.9){
    //   y = 4*radius;
    //   while(y<.9){
    //     x = 4*radius;
    //     while(x<.9){
    //       se = cafes::make_sphere( { x, y, z}, radius, 0);
    //       pt.push_back(cafes::make_particle(se, {1., 0.}, {0.,0.,0.}));
    //       x += 4*radius;
    //     }
    //     y += 4*radius;
    //   }
    //   z += 4*radius;
    // }

    // std::cout << pt.size() << "\n";
    std::vector<cafes::particle<decltype(se1)>> pt{ cafes::make_particle(se1, {1., 0., 0.}, {0.,0.,0.}),
                                                    cafes::make_particle(se2, {-1., 0., 0.}, {0.,0.,0.})};

    auto s = cafes::make_DtoN(pt, st, {.1, .1});
    
    ierr = s.create_Mat_and_Vec();CHKERRQ(ierr);
    ierr = s.setup_RHS();CHKERRQ(ierr);
    ierr = s.setup_KSP();CHKERRQ(ierr);

    ierr = cafes::singularity::compute_singular_forces<dim, decltype(*(s.ctx))>(*(s.ctx), 100);CHKERRQ(ierr);
    std::cout << "force on particle 0 -> " << s.ctx->particles[0].force_ << "\n";

    //ierr = cafes::singularity::save_singularity<dim, decltype(*(s.ctx))>("Resultats", "test_avec_tronc", *(s.ctx));CHKERRQ(ierr);
    ierr = s.solve();CHKERRQ(ierr);


    ierr = cafes::io::save_hdf5("Resultats", "test_avec_sing", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
    ierr = cafes::io::saveParticles("Resultats", "test_avec_sing", pt);CHKERRQ(ierr);
    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}
