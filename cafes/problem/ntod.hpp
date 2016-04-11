#ifndef CAFES_PROBLEM_NTOD_HPP_INCLUDED
#define CAFES_PROBLEM_NTOD_HPP_INCLUDED

#include <problem/dton.hpp>
#include <particle/forces_torques.hpp>

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

      ierr = VecSet(ctx->problem.rhs, 0.);CHKERRQ(ierr);
      ierr = VecSet(ctx->problem.sol, 0.);CHKERRQ(ierr);

      auto box = fem::get_DM_bounds<Dimensions>(ctx->problem.ctx->problem.ctx->dm, 0);
      auto& h = ctx->problem.ctx->problem.ctx->h;

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
          if (geometry::intersect(box, pbox))
          {
            for(std::size_t d=0; d<Dimensions; ++d)
              ctx->particles[ipart].velocity_[d] = px[num++];
            if (Dimensions == 2)
              ctx->particles[ipart].angular_velocity_[2] = px[num++];
            else
              for(std::size_t d=0; d<Dimensions; ++d)
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

      std::vector<double> forces, torques;

      forces.resize(ctx->particles.size()*Dimensions);
      torques.resize(ctx->particles.size()*torque_size);

      std::fill(forces.begin(), forces.end(), 0.);
      std::fill(torques.begin(), torques.end(), 0.);

      ierr = forces_torques_with_control(ctx->particles,
                                         ctx->problem.sol,
                                         box,
                                         forces,
                                         torques,
                                         ctx->num,
                                         h,
                                         ctx->compute_singularity);CHKERRQ(ierr);

      // set y with the forces and the torques computed by DtoN
      num = 0;
      PetscScalar *py;
      ierr = VecGetArray(y, &py);CHKERRQ(ierr);

      for (std::size_t ipart=0; ipart<ctx->particles.size(); ++ipart)
      {
        auto p = ctx->particles[ipart];
        auto pbox = p.bounding_box(h);
        if (geometry::intersect(box, pbox))
        {
          for(std::size_t d=0; d<Dimensions; ++d)
            py[num++] = forces[ipart*Dimensions + d];
          for(std::size_t d=0; d<torque_size; ++d)
            py[num++] = torques[ipart*Dimensions + d];
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
      Vec stokes_sol_save;
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
        VecDuplicate(dton_.problem_.sol, &stokes_sol_save);
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

        // dton_.ctx->compute_rhs = true;
        // dton_.ctx->add_rigid_motion = false;
        // dton_.ctx->compute_singularity = true;

        // ierr = dton_.solve_last_problem();CHKERRQ(ierr);

        // // save Stokes sol
        // ierr = VecCopy(dton_.problem_.sol, stokes_sol_save);CHKERRQ(ierr);
        // ierr = cafes::io::save_VTK("Resultats", "stokes_sol_rhs", dton_.problem_.sol, dton_.problem_.ctx->dm, dton_.problem_.ctx->h);CHKERRQ(ierr);

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
        ctx->compute_singularity = true;

        ierr = KSPSolve(ksp, rhs, sol);CHKERRQ(ierr);

        auto box = fem::get_DM_bounds<Dimensions>(ctx->problem.ctx->problem.ctx->dm, 0);
        auto& h = ctx->problem.ctx->problem.ctx->h;
        std::size_t num = 0, num_print=0;
        std::size_t torque_size = (Dimensions == 2)?1:3;

        PetscScalar *psol;
        ierr = VecGetArray(sol, &psol);CHKERRQ(ierr);
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

        for (std::size_t ipart=0; ipart<ctx->particles.size(); ++ipart)
        {
          auto p = ctx->particles[ipart];
          auto pbox = p.bounding_box(h);
          if (geometry::intersect(box, pbox)){
            for(std::size_t d=0; d<Dimensions; ++d)
              ctx->particles[ipart].velocity_[d] = psol[num++];
            for(std::size_t d=0; d<torque_size; ++d)
              ctx->particles[ipart].angular_velocity_[d] = psol[num++];

            std::cout << "[" << rank << "] " << "particle " << ipart << ":\n";
            std::cout << "[" << rank << "] " << "    velocity: ";
            for(std::size_t d=0; d<Dimensions; ++d)
              std::cout << psol[num_print++] << " ";
            std::cout << "\n" << "[" << rank << "] " <<  "    angular_velocity: ";
            for(std::size_t d=0; d<torque_size; ++d)
              std::cout << psol[num_print++] << " ";
            std::cout << "\n";
          }
        }
        ierr = VecRestoreArray(sol, &psol);CHKERRQ(ierr);

        // solve the problem with the right control
        dton_.default_flags_ = true;
        ierr = VecSet(dton_.sol, 0.);CHKERRQ(ierr);
        ierr = VecSet(dton_.rhs, 0.);CHKERRQ(ierr);
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