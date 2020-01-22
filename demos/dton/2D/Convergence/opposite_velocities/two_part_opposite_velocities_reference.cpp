#include <cafes.hpp>
#include <petsc.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <particle/singularity/add_singularity.hpp>

void zeros(const PetscReal x[], PetscScalar *u)
{
    *u = 0.;
}

void ones(const PetscReal x[], PetscScalar *u)
{
    *u = 10.;
}

void ones_m(const PetscReal x[], PetscScalar *u)
{
    *u = -1.;
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 2;
    std::string saverep = "Reference/";
    const char * srep = saverep.c_str();


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
    

    double R1 = .1;
    double R2 = .1;
    double dist = .1/5;

    auto se1 = cafes::make_circle({.5 - R1 - dist / 2, .5}, R1, 0);
    auto se2 = cafes::make_circle({.5 + R2 + dist / 2, .5}, R2, 0);

    std::vector<cafes::particle<decltype(se1)>> pt{
        cafes::make_particle_with_velocity(se1, {1., 0.}, 0.),
        cafes::make_particle_with_velocity(se2, {-1., 0.}, 0.)};

    auto s = cafes::make_DtoN(pt, st, .1);

    ierr = s.create_Mat_and_Vec();
    CHKERRQ(ierr);

    ierr = s.setup_RHS();
    CHKERRQ(ierr);
    ierr = s.setup_KSP();
    CHKERRQ(ierr);
    ierr = s.solve();
    CHKERRQ(ierr);

    // COARSE DMDA INFO
    DM davCoarse, dapCoarse;
    DMDALocalInfo infopCoarse, infovCoarse;
    ierr = DMCompositeGetEntries(st.ctx->dm, &davCoarse, &dapCoarse);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(davCoarse, &infovCoarse);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dapCoarse, &infopCoarse);CHKERRQ(ierr);
    std::cout << "Size Coarse Mesh.   Pressure : " << infopCoarse.mx << " " << infopCoarse.my << " Velocity : " << infovCoarse.mx << " " << infovCoarse.my  << std::endl;

    std::string stout0 = "two_part_rhs_";
    stout0.append(std::to_string(mx));
    const char * stw0 = stout0.c_str();
    ierr = cafes::io::save_VTK(srep, stw0, s.sol_rhs, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);

    stout0 = "two_part_tmp_";
    stout0.append(std::to_string(mx));
    stw0 = stout0.c_str();
    ierr = cafes::io::save_VTK(srep, stw0, s.sol_tmp, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);

    stout0 = "two_part_reg_";
    stout0.append(std::to_string(mx));
    stw0 = stout0.c_str();
    ierr = cafes::io::save_VTK(srep, stw0, s.sol_reg, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);
    
    // DMDA INFO
    DM dav, dap;
    DMDALocalInfo infop, infov;

    ierr = DMCompositeGetEntries(st.ctx->dm, &dav, &dap);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dav, &infov);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);

    // ZEROS IN PARTICLES (REG SOLUTION)
    auto solu = cafes::petsc::petsc_vec<dim>(st.ctx->dm, s.sol_reg, 0, false);
    auto solp = cafes::petsc::petsc_vec<dim>(st.ctx->dm, s.sol_reg, 1, false);

    for (std::size_t j =0; j< infov.my; j++ )
    {
        for (std::size_t i = 0; i <infov.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto pts = cafes::geometry::position<double, dim>{i*st.ctx->h[0], j*st.ctx->h[1]};
            auto usol = solu.at_g(pos);
            if (se1.contains(pts) or se2.contains(pts))
            {
                usol[0] = 0.;
                usol[1] = 0.;
            }
        }
    }

    for (std::size_t j =0; j<infop.my; j++ )
    {
        for (std::size_t i = 0; i <infop.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto pts = cafes::geometry::position<double, dim>{2*i*st.ctx->h[0], 2*j*st.ctx->h[1]};
            auto psol = solp.at_g(pos);
            if (se1.contains(pts) or se2.contains(pts))
            {
                psol[0] = 0.;
            }
        }
    }

    std::string stout1 = "two_part_reg_zerosInParts";
    stout1.append(std::to_string(mx));
    const char * stw1 = stout1.c_str();
    ierr = cafes::io::save_VTK(srep, stw1, s.sol_reg, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);

    // ADD SINGULARITY TO UREG
    cafes::singularity::add_singularity_to_ureg(st.ctx->dm, st.ctx->h, s.sol_reg, pt);

    stout1 = "two_part_total_";
    stout1.append(std::to_string(mx));
    stw1 = stout1.c_str();
    ierr = cafes::io::save_VTK(srep, stw1, s.sol_reg, st.ctx->dm, st.ctx->h);CHKERRQ(ierr);

    // SAVE TO READ BACK

    auto soluref = cafes::petsc::petsc_vec<dim>(st.ctx->dm, s.sol_reg, 0, false);
    auto solpref = cafes::petsc::petsc_vec<dim>(st.ctx->dm, s.sol_reg, 1, false);

    std::ofstream outrefv, outrefp;
    std::string savevel = saverep + "reference_sem_velocity_" + std::to_string(mx) + ".txt";
    std::string savepres = saverep + "reference_sem_pressure_" + std::to_string(mx) + ".txt";
    const char * svel = savevel.c_str();
    const char * spres = savepres.c_str();
    outrefv.open(svel, std::ofstream::out);
    outrefp.open(spres, std::ofstream::out);

     for (std::size_t j =0; j< infov.my; j++ )
    {
        for (std::size_t i = 0; i <infov.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto pts = cafes::geometry::position<double, dim>{i*st.ctx->h[0], j*st.ctx->h[1]};
            auto usol = soluref.at_g(pos);
            outrefv << usol[0] << " " << usol[1] << "\n";
        }
    }

    for (std::size_t j =0; j<infop.my; j++ )
    {
        for (std::size_t i = 0; i <infop.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto pts = cafes::geometry::position<double, dim>{2*i*st.ctx->h[0], 2*j*st.ctx->h[1]};
            auto psol = solpref.at_g(pos);
            outrefp << psol[0] << "\n";
            
        }
    }
    outrefv.close();
    outrefp.close();
    
    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;
}
