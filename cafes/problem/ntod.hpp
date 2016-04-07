#ifndef CAFES_PROBLEM_NTOD_HPP_INCLUDED
#define CAFES_PROBLEM_NTOD_HPP_INCLUDED

#include <problem/dton.hpp>

namespace cafes
{
  namespace problem
  {

    #undef __FUNCT__
    #define __FUNCT__ "NtoD_matrix"
    template<std::size_t Dimensions, typename Ctx>
    PetscErrorCode NtoD_matrix(Mat A, Vec x, Vec y){
      PetscErrorCode ierr;
      PetscFunctionBeginUser;

      std::size_t torque_size = (Dimensions==2)?1:3;

      Ctx *ctx;
      ierr = MatShellGetContext(A, (void**) &ctx);CHKERRQ(ierr);

      auto box = fem::get_DM_bounds<Dimensions>(ctx->problem.ctx->problem.ctx->dm, 0);
      auto& h = ctx->problem.ctx->problem.ctx->h;

      // set velocity and angular velocity on particles
      std::size_t num = 0;

      if (ctx->compute_rhs)
      {
        ctx->problem.ctx->compute_rhs = true;
        ctx->problem.ctx->add_rigid_motion = false;
        ctx->problem.ctx->compute_singularity = true;
      }
      else
      {
        ctx->problem.ctx->compute_rhs = false;
        ctx->problem.ctx->add_rigid_motion = true;
        ctx->problem.ctx->compute_singularity = true;

        // set velocity and angular velocity on particles
        PetscScalar const *px;
        ierr = VecGetArrayRead(x, &px);CHKERRQ(ierr);

        std::size_t num_print = 0;
        for (std::size_t ipart=0; ipart<ctx->particles.size(); ++ipart)
        {
          auto p = ctx->particles[ipart];
          auto pbox = p.bounding_box(h);
          if (geometry::intersect(box, pbox)){
            std::cout << "particle " << ipart << ":\n";
            std::cout << "        velocity: ";
            for(std::size_t d=0; d<Dimensions; ++d)
              std::cout << px[num_print++] << " ";
            std::cout << "\n";
            std::cout << "    ang_velocity: ";
            for(std::size_t d=0; d<Dimensions; ++d)
              std::cout << px[num_print++] << " ";
            std::cout << "\n";

            for(std::size_t d=0; d<Dimensions; ++d)
              ctx->particles[ipart].velocity_[d] = px[num++];
            for(std::size_t d=0; d<torque_size; ++d)
              ctx->particles[ipart].angular_velocity_[d] = px[num++];
          }
        }
        ierr = VecRestoreArrayRead(x, &px);CHKERRQ(ierr);
      }

      ierr = ctx->problem.setup_RHS();CHKERRQ(ierr);

      // solve DtoN with these velocities
      ctx->problem.ctx->compute_rhs = false;
      ctx->problem.ctx->add_rigid_motion = false;
      ctx->problem.ctx->compute_singularity = false;

      ierr = ctx->problem.solve();CHKERRQ(ierr);

      if (ctx->compute_rhs)
      {
        ctx->problem.ctx->compute_rhs = true;
        ctx->problem.ctx->add_rigid_motion = false;
        ctx->problem.ctx->compute_singularity = true;
      }
      else
      {
        ctx->problem.ctx->compute_rhs = false;
        ctx->problem.ctx->add_rigid_motion = true;
        ctx->problem.ctx->compute_singularity = true;
      }
      
      //ierr = ctx->problem.solve_last_problem();CHKERRQ(ierr);
      
      if (ctx->compute_singularity)
      {
        ierr = singularity::compute_singular_forces(ctx->particles, ctx->forces, h, 100);CHKERRQ(ierr);
      }

      std::vector<double> mean, mean_local;
      std::vector<double> cross_prod, cross_prod_local;

      mean_local.resize(ctx->particles.size()*Dimensions);
      cross_prod_local.resize(ctx->particles.size()*torque_size);

      std::fill(mean_local.begin(), mean_local.end(), 0.);
      std::fill(cross_prod_local.begin(), cross_prod_local.end(), 0.);

      mean.resize(ctx->particles.size()*Dimensions);
      cross_prod.resize(ctx->particles.size()*torque_size);
      std::fill(mean.begin(), mean.end(), 0.);
      std::fill(cross_prod.begin(), cross_prod.end(), 0.);

      num = 0;
      PetscScalar const *pg;
      ierr = VecGetArrayRead(ctx->problem.sol, &pg);CHKERRQ(ierr);

      for(std::size_t ipart=0; ipart<ctx->particles.size(); ++ipart)
      {
        auto p = ctx->particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          auto new_box = geometry::overlap_box(box, pbox);
          auto pts = find_fluid_points_insides(p, new_box, h);
          for(auto& ind: pts){
            std::array<double, Dimensions> g;
            for(std::size_t d=0; d<Dimensions; ++d){
              mean_local[ipart*Dimensions + d] += pg[num];
              g[d] = pg[num++];
            }
            cafes::geometry::position<Dimensions, double> pos { ind[0]*h[0] - p.center_[0], ind[1]*h[1] - p.center_[1], ind[2]*h[2] - p.center_[2] };
            cross_prod_local[ipart*Dimensions    ] += pos[1]*g[2] - pos[2]*g[1];
            cross_prod_local[ipart*Dimensions + 1] += pos[2]*g[0] - pos[0]*g[2];
            cross_prod_local[ipart*Dimensions + 2] += pos[0]*g[1] - pos[1]*g[0];
          }       
        }
        for(std::size_t d=0; d<Dimensions; ++d){
          mean_local[ipart*Dimensions + d] *= p.volume()/ctx->num[ipart];
          cross_prod_local[ipart*Dimensions + d] *= p.volume()/ctx->num[ipart];
        }

      }
      ierr = VecRestoreArrayRead(ctx->problem.sol, &pg);CHKERRQ(ierr);

      MPI_Allreduce(mean_local.data(), mean.data(), mean.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
      MPI_Allreduce(cross_prod_local.data(), cross_prod.data(), cross_prod.size(), MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

      // set y with the forces and the torques computed by DtoN
      num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for (std::size_t ipart=0; ipart<ctx->particles.size(); ++ipart)
      {
        auto p = ctx->particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox)){
          if (ctx->compute_singularity){
            
            std::cout << "particle " << ipart << ":\n";
            std::cout << "    mean: ";
            for(std::size_t d=0; d<Dimensions; ++d)
              std::cout << mean[ipart*Dimensions + d] << " ";
            std::cout << "\n";
            std::cout << "    sing: ";
            for(std::size_t d=0; d<Dimensions; ++d)
              std::cout << ctx->forces[ipart][d] << " ";
            std::cout << "\n";

            for(std::size_t d=0; d<Dimensions; ++d)
              py[num++] = ctx->forces[ipart][d] + mean[ipart*Dimensions + d];
            for(std::size_t d=0; d<torque_size; ++d)
              //py[num++] = ctx->particles[ipart].torque_[d] + vol*cross_prod[ipart*Dimensions + d];
              py[num++] = cross_prod[ipart*Dimensions + d];
          }
          else
          {
            for(std::size_t d=0; d<Dimensions; ++d)
              py[num++] = mean[ipart*Dimensions + d];
            for(std::size_t d=0; d<torque_size; ++d)
              py[num++] = cross_prod[ipart*Dimensions + d];
          }
        }
      }
      ierr = VecRestoreArray(y, &py);CHKERRQ(ierr);

      PetscFunctionReturn(0);
    }


    template<typename Shape, std::size_t Dimensions, typename Problem_type>
    struct NtoD : public Problem<Dimensions>
    {
      using force_type   = physics::force<Dimensions>;

      DtoN<Shape, Dimensions, Problem_type> dton_;

      std::vector<force_type> forces_;

      Vec sol;
      Vec rhs;
      Mat A;
      KSP ksp;

      using Ctx = particle_context<Dimensions, Shape, typename problem::DtoN<Shape, Dimensions, Problem_type> >;
      Ctx *ctx;

      using dpart_type = typename std::conditional<Dimensions == 2, 
                                  double, 
                                  std::array<double, 2>>::type;

      NtoD(std::vector<particle<Shape>>& parts, Problem_type& p, dpart_type dpart):
      dton_{parts, p, dpart}
      {
        dton_.create_Mat_and_Vec();
        dton_.setup_KSP();
        dton_.default_flags_ = false;
        forces_.resize(dton_.parts_.size());
      }

      #undef __FUNCT__
      #define __FUNCT__ "create_Mat_and_Vec"
      PetscErrorCode create_Mat_and_Vec()
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;

        ctx = new Ctx{dton_,
                      dton_.parts_,
                      dton_.surf_points_,
                      dton_.radial_vec_,
                      dton_.nb_surf_points_,
                      dton_.num_,
                      forces_,
                      dton_.scale_,
                      false,
                      false,
                      false};

        auto box = fem::get_DM_bounds<Dimensions>(ctx->problem.ctx->problem.ctx->dm, 0);
        auto& h = ctx->problem.ctx->problem.ctx->h;

        std::size_t size = 0;
        for(auto& p: ctx->particles){
          auto pbox = p.bounding_box(h);
          if (geometry::intersect(box, pbox)){
            size++;
          }
        }

        std::size_t force_size = Dimensions;
        std::size_t torque_size = (Dimensions == 2)?1:3;
        std::size_t local_size = size*(force_size + torque_size);

        ierr = MatCreateShell(PETSC_COMM_WORLD, local_size, local_size, PETSC_DECIDE, PETSC_DECIDE, ctx, &A);CHKERRQ(ierr);
        ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)(void))NtoD_matrix<Dimensions, Ctx>);CHKERRQ(ierr);

        ierr = MatCreateVecs(A, &sol, &rhs);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "setup_RHS"
      virtual PetscErrorCode setup_RHS() override 
      {
        PetscErrorCode ierr;
        PetscFunctionBeginUser;

        ctx->compute_rhs = true;
        ctx->compute_singularity = true;

        ierr = MatMult(A, sol, rhs);CHKERRQ(ierr);
        ierr = VecScale(rhs, -1.);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "setup_KSP"
      virtual PetscErrorCode setup_KSP() override
      {
        PetscErrorCode ierr;
        PC             pc;
        PetscFunctionBegin;

        ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);CHKERRQ(ierr);
        ierr = KSPSetOptionsPrefix(ksp, "NtoD_");CHKERRQ(ierr);

        ierr = KSPSetOperators(ksp, A, A);CHKERRQ(ierr);
        ierr = KSPGetPC(ksp, &pc);CHKERRQ(ierr);
        ierr = PCSetType(pc, PCNONE);CHKERRQ(ierr);
        ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }

      #undef __FUNCT__
      #define __FUNCT__ "solve"
      virtual PetscErrorCode solve() override
      {
        PetscErrorCode ierr;
        PetscFunctionBegin;

        ctx->compute_rhs = false;

        ierr = KSPSolve(ksp, rhs, sol);CHKERRQ(ierr);

        auto box = fem::get_DM_bounds<Dimensions>(ctx->problem.ctx->problem.ctx->dm, 0);
        auto& h = ctx->problem.ctx->problem.ctx->h;
        std::size_t num = 0, num_print=0;
        std::size_t torque_size = (Dimensions == 2)?1:3;

        PetscScalar *psol;
        ierr = VecGetArray(sol, &psol);CHKERRQ(ierr);

        for (std::size_t ipart=0; ipart<ctx->particles.size(); ++ipart)
        {
          auto p = ctx->particles[ipart];
          auto pbox = p.bounding_box(h);
          if (geometry::intersect(box, pbox)){
            for(std::size_t d=0; d<Dimensions; ++d)
              ctx->particles[ipart].velocity_[d] = psol[num++];
            for(std::size_t d=0; d<torque_size; ++d)
              ctx->particles[ipart].angular_velocity_[d] = psol[num++];

            std::cout << "particle " << ipart << ":\n";
            std::cout << "    velocity: ";
            for(std::size_t d=0; d<Dimensions; ++d)
              std::cout << psol[num_print++] << " ";
            std::cout << "\n    angular_velocity: ";
            for(std::size_t d=0; d<torque_size; ++d)
              std::cout << psol[num_print++] << " ";
            std::cout << "\n";
          }
        }
        ierr = VecRestoreArray(sol, &psol);CHKERRQ(ierr);

        // solve the problem with the right control
        dton_.default_flags_ = true;
        ierr = dton_.setup_RHS();CHKERRQ(ierr);
        ierr = dton_.solve();CHKERRQ(ierr);

        PetscFunctionReturn(0);
      }      
    };
  }

  template<typename PL, typename Problem_type, typename Dimensions = typename PL::value_type::dimension_type> 
  auto make_NtoD(PL& pt, Problem_type& p, 
                 typename std::conditional<Dimensions::value == 2, double, 
                                           std::array<double, 2>>::type const& dpart)
  {
    using s_t = typename PL::value_type::shape_type;
    return problem::NtoD<s_t, Dimensions::value, Problem_type>{pt, p, dpart};
  }

}
#endif