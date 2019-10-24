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
    std::size_t const dim = 2;

    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({ {{zeros , zeros}} 
                                  , {{zeros , zeros}}
                                  , {{zeros, zeros}}
                                  , {{zeros, zeros}}
                                  });
    auto rhs = cafes::make_rhs<dim>({{ zeros, zeros }});
    auto st = cafes::make_stokes<dim>(bc, rhs);

    // auto se1 = cafes::make_ellipsoid<2>( {.3,.5}, {.1,.2}, 2., M_PI/2);
    // auto se3 = cafes::make_ellipsoid<2>( { .5, .5}, {.1,.1}, 2., 0);
    // auto se2 = cafes::make_ellipsoid<2>( {.7,.5}, {.1,.1}, 1., 0);

    double r = 0.05;
    auto se1 = cafes::make_circle( {.5 - r - r/4, .5}, r);
    auto se2 = cafes::make_circle( {.5 + r + r/4, .5}, r);
    std::vector<cafes::particle<decltype(se1)>> pt{ cafes::make_particle_with_force(se1, {1., 0}, 1.),
                                                    cafes::make_particle_with_force(se2, {1., 0}, 1.)
                                                 };

    auto s = cafes::make_SEM(pt, st, st.ctx->h[0]);
    ierr = s.create_Mat_and_Vec();CHKERRQ(ierr);
    ierr = s.setup_RHS();CHKERRQ(ierr);
    ierr = s.setup_KSP();CHKERRQ(ierr);
    ierr = s.solve();CHKERRQ(ierr);

    s.get_new_velocities();
    std::cout << s.parts_[0].velocity_[0] << " " << s.parts_[0].velocity_[1] << "\n";
    std::cout << s.parts_[1].velocity_[0] << " " << s.parts_[1].velocity_[1] << "\n";
    ierr = cafes::io::save_VTK("Resultats", "test", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);

    // Refinement and interpolation test
    std::array<int, 2> refine = {2,2};
    Vec sol_refine;
    DM dm_refine;
    std::array<double, 2> h_refine;
    h_refine[0] = st.ctx->h[0]/refine[0];
    h_refine[1] = st.ctx->h[1]/refine[1];
    std::cout << h_refine[0] << " " << h_refine[1] << std::endl;

    ierr = cafes::posttraitement::linear_interpolation(st.ctx->dm, st.sol, dm_refine, sol_refine, refine, st.ctx->h);CHKERRQ(ierr);
    ierr = cafes::io::save_VTK("Resultats", "test_refine", sol_refine, dm_refine, h_refine);CHKERRQ(ierr);
    
    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}