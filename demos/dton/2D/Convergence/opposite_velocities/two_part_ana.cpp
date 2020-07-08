#include <cafes.hpp>
#include <petsc.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <particle/singularity/add_singularity.hpp>

double C_u = 0., C_v = 0.;
double center_x = .5, center_y = .5;
double R = 0.1;

void quadratic_u_2d(const PetscReal x[], PetscScalar *u)
{
    double xx = x[0]-center_x;
    double yy = x[1]-center_y;
    *u = -yy*xx*xx + yy*(R*R-yy*yy) + C_u;
}

void quadratic_v_2d(const PetscReal x[], PetscScalar *u)
{
    double xx = x[0]-center_x;
    double yy = x[1]-center_y;
    *u = xx*yy*yy - xx*(R*R-xx*xx) + C_v;
}

void zeros(const PetscReal x[], PetscScalar *u)
{
    *u = 0.;
}

void fx(const PetscReal x[], PetscScalar *u)
{
    double xx = x[0]-center_x;
    double yy = x[1]-center_y;
    *u = 8*yy;
}

void fy(const PetscReal x[], PetscScalar *u)
{
    double xx = x[0]-center_x;
    double yy = x[1]-center_y;
    *u = -8*xx;
}

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    std::size_t const dim = 2;
    int const nref = 2048;
    std::string saverep = "Resultats";
    std::string refrep = "Reference";
    const char * srep = saverep.c_str();
    const char * rrep = refrep.c_str();

    ierr = PetscInitialize(&argc, &argv, (char *)0, (char *)0);
    CHKERRQ(ierr);

    auto bc = cafes::make_bc<dim>({
        {{quadratic_u_2d, quadratic_v_2d}}, // left
        {{quadratic_u_2d, quadratic_v_2d}}, // right
        {{quadratic_u_2d, quadratic_v_2d}}, // bottom
        {{quadratic_u_2d, quadratic_v_2d}}  // top
    });

    auto rhs = cafes::make_rhs<dim>({{fx, fy}});
    auto st = cafes::make_stokes<dim>(bc, rhs);

    int const mx = st.opt.mx[0] - 1;
    
    auto se = cafes::make_circle({center_x, center_y}, R, 0);

    std::vector<cafes::particle<decltype(se)>> pt{
        cafes::make_particle_with_velocity(se, {C_u, C_v}, 0.)};

    auto s = cafes::make_DtoN(pt, st, .01);

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

    std::string stout0 = "two_part_ana_rhs_";
    stout0.append(std::to_string(mx));
    const char * stw0 = stout0.c_str();
    ierr = cafes::io::save_hdf5(srep, stw0, s.sol_rhs, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);

    stout0 = "two_part_ana_tmp_";
    stout0.append(std::to_string(mx));
    stw0 = stout0.c_str();
    ierr = cafes::io::save_hdf5(srep, stw0, s.sol_tmp, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);

    {
        // ZEROS IN PARTICLES (REFINED SOLUTION)
        auto bbox = cafes::fem::get_DM_bounds<dim>(st.ctx->dm, 0);
        auto bboxp = cafes::fem::get_DM_bounds<dim>(st.ctx->dm, 1);
        auto solu = cafes::petsc::petsc_vec<dim>(st.ctx->dm, s.sol_reg, 0, false);
        auto solp = cafes::petsc::petsc_vec<dim>(st.ctx->dm, s.sol_reg, 1, false);

        for(std::size_t j=bbox.bottom_left[1]; j<bbox.upper_right[1]; ++j)
        {
            for(std::size_t i=bbox.bottom_left[0]; i<bbox.upper_right[0]; ++i)
            {
                auto pos = cafes::geometry::position<int, dim>{i,j};
                auto pts = cafes::geometry::position<double, dim>{i*st.ctx->h[0], j*st.ctx->h[0]};
                auto usol = solu.at_g(pos);
                for(std::size_t ipart=0; ipart<pt.size(); ++ipart)
                {
                    auto part = pt[ipart];
                    if (part.contains(pts))
                    {
                        usol[0] = part.velocity_[0];
                        usol[1] = part.velocity_[1];
                    }
                }
            }
        }
    }

    stout0 = "two_part_ana_reg_";
    stout0.append(std::to_string(mx));
    stw0 = stout0.c_str();
    ierr = cafes::io::save_hdf5(srep, stw0, s.sol_reg, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);

    stout0 = "two_part_ana_";
    stout0.append(std::to_string(mx));
    stw0 = stout0.c_str();
    ierr = cafes::io::save_hdf5(srep, stw0, st.sol, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);

    // REFINEMENT AND INTERPOLATION
    int const fine = nref/mx;
    std::array<int, dim> refine = {fine,fine};
    std::cout << "Refinment factors : " << refine[0] << " " << refine[1] << std::endl;
    Vec sol_refine;
    DM dm_refine;
    std::array<double, dim> h_refine;
    h_refine[0] = st.ctx->h[0]/refine[0];
    h_refine[1] = st.ctx->h[1]/refine[1];
    std::cout << "Interpolating to reference mesh..." << std::endl;
    ierr = cafes::posttraitement::linear_interpolation(st.ctx->dm, s.sol_reg, dm_refine, sol_refine, refine, st.ctx->h);CHKERRQ(ierr);
    std::cout << "Done." << std::endl;
    
    // // DMDA INFO
    // DM dav, dap;
    // DMDALocalInfo infop, infov;

    // ierr = DMCompositeGetEntries(dm_refine, &dav, &dap);CHKERRQ(ierr);
    // ierr = DMDAGetLocalInfo(dav, &infov);CHKERRQ(ierr);
    // ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);
    
    // ZEROS IN PARTICLES (REFINED SOLUTION)
    auto bbox = cafes::fem::get_DM_bounds<dim>(dm_refine, 0);
    auto bboxp = cafes::fem::get_DM_bounds<dim>(dm_refine, 1);
    auto solu = cafes::petsc::petsc_vec<dim>(dm_refine, sol_refine, 0, false);
    auto solp = cafes::petsc::petsc_vec<dim>(dm_refine, sol_refine, 1, false);

    for(std::size_t j=bbox.bottom_left[1]; j<bbox.upper_right[1]; ++j)
    {
        for(std::size_t i=bbox.bottom_left[0]; i<bbox.upper_right[0]; ++i)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto pts = cafes::geometry::position<double, dim>{i*h_refine[0], j*h_refine[1]};
            auto usol = solu.at_g(pos);
            for(std::size_t ipart=0; ipart<pt.size(); ++ipart)
            {
                auto part = pt[ipart];
                if (part.contains(pts))
                {
                    usol[0] = part.velocity_[0];
                    usol[1] = part.velocity_[1];
                }
            }
        }
    }

    // for(std::size_t j=bboxp.bottom_left[1]; j<bboxp.upper_right[1]; ++j)
    // {
    //     for(std::size_t i=bboxp.bottom_left[0]; i<bboxp.upper_right[0]; ++i)
    //     {
    //         auto pos = cafes::geometry::position<int, dim>{i,j};
    //         auto pts = cafes::geometry::position<double, dim>{2*i*h_refine[0], 2*j*h_refine[1]};
    //         auto psol = solp.at_g(pos);
    //         if (se1.contains(pts) or se2.contains(pts))
    //         {
    //             psol[0] = 0.;
    //         }
    //     }
    // }

    stout0 = "two_part_ana_refined_";
    stout0.append(std::to_string(mx));
    stw0 = stout0.c_str();
    ierr = cafes::io::save_hdf5(srep, stw0, sol_refine, dm_refine, h_refine);CHKERRQ(ierr);

    // ADD SINGULARITY TO UREG REFINED
    cafes::singularity::add_singularity_to_ureg(dm_refine, h_refine, sol_refine, pt);

    stout0 = "two_part_ana_total_";
    stout0.append(std::to_string(mx));
    stw0 = stout0.c_str();
    ierr = cafes::io::save_hdf5(srep, stw0, sol_refine, dm_refine, h_refine);CHKERRQ(ierr);
}
