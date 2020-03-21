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
                                  , {{zeros , zeros, zeros}} // bottom
                                  , {{zeros , zeros, zeros}}   // top
                                  , {{zeros , zeros, zeros}}  // front
                                  , {{zeros , zeros, zeros}}  // back
                                  });

    auto rhs = cafes::make_rhs<dim>({{ zeros, zeros, zeros }});
    
    auto st = cafes::make_stokes<dim>(bc, rhs);

    double R1   = 0.1;

    auto se1 = cafes::make_sphere( {.5, .5, .5}, R1, 0);

    std::vector<cafes::particle<decltype(se1)>> pt{ cafes::make_particle(se1, {.1, 0., 0.}, {0.,0.,0.}) };

    auto s = cafes::make_DtoN(pt, st, {.1, .1});
    
    ierr = s.create_Mat_and_Vec();CHKERRQ(ierr);
    ierr = s.setup_RHS();CHKERRQ(ierr);
    ierr = s.setup_KSP();CHKERRQ(ierr);

    ierr = s.solve();CHKERRQ(ierr);


    std::vector<double> mean;
    std::vector<double> cross_prod;

    mean.resize(pt.size()*dim);
    cross_prod.resize(pt.size()*dim);

    int num = 0;
    PetscScalar const *pg;
    ierr = VecGetArrayRead(s.sol, &pg);CHKERRQ(ierr);

    auto box = cafes::fem::get_DM_bounds<dim>(st.ctx->dm, 0);
    auto& h = st.ctx->h;
    for(auto& p: s.parts_){
      auto pbox = p.bounding_box(st.ctx->h);
      if (cafes::geometry::intersect(box, pbox)){
        auto new_box = cafes::geometry::overlap_box(box, pbox);
        auto pts = cafes::find_fluid_points_insides(p, new_box, st.ctx->h);
        for(auto& ind: pts){
          std::array<double, dim> g;
          for(std::size_t d=0; d<dim; ++d){
            mean[d] += pg[num];
            g[d] = pg[num++];
          }
          cafes::geometry::position<dim, double> pos { ind[0]*h[0] - p.center_[0], ind[1]*h[1] - p.center_[1], ind[2]*h[2] - p.center_[2] };
          cross_prod[0] += pos[1]*g[2] - pos[2]*g[1];
          cross_prod[1] += pos[2]*g[0] - pos[0]*g[2];
          cross_prod[2] += pos[0]*g[1] - pos[1]*g[0];
        }       
      }
    }
    ierr = VecRestoreArrayRead(s.sol, &pg);CHKERRQ(ierr);

    for(std::size_t d=0; d<dim; ++d){
      mean[d] *= s.parts_[0].volume()/s.num_[0];
      cross_prod[d] *= s.parts_[0].volume()/s.num_[0];
    }

    std::cout << "force: ";
    for(std::size_t d=0; d<dim; ++d)
      std::cout << mean[d] << " ";
    std::cout << "\n";
    std::cout << "torque: ";
    for(std::size_t d=0; d<dim; ++d)
      std::cout << cross_prod[d] << " ";
    std::cout << "\n";

    ierr = cafes::io::save_hdf5("Resultats", "dton_one_part", st.sol, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);
    ierr = cafes::io::saveParticles("Resultats", "dton_one_part", pt);CHKERRQ(ierr);

    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}
