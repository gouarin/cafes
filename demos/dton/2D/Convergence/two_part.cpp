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
    int const nref = 600;

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
    

    double R1 = .05;
    double R2 = .05;
    double dist = .05/5;

    auto se1 = cafes::make_circle({.5 - R1 - dist / 2, .5}, R1, 0);
    auto se2 = cafes::make_circle({.5 + R2 + dist / 2, .5}, R2, 0);
    // std::vector<cafes::particle<decltype(se1)>> pt{ cafes::make_particle(se1,
    // { 1., 0.}, {0.,0.,0.}),
    //                                           cafes::make_particle(se2, {-1.,
    //                                           0.}, {0.,0.,0.})};
    std::vector<cafes::particle<decltype(se1)>> pt{
        cafes::make_particle_with_velocity(se1, {1., 0.}, 0.),
        cafes::make_particle_with_velocity(se2, {1., 0.}, 0.)};

    auto s = cafes::make_DtoN(pt, st, .1);

    ierr = s.create_Mat_and_Vec();
    CHKERRQ(ierr);
    // ierr = cafes::singularity::save_singularity<dim,
    // decltype(*(s.ctx))>("Resultats", "two_part", *(s.ctx));CHKERRQ(ierr);

    ierr = s.setup_RHS();
    CHKERRQ(ierr);
    ierr = s.setup_KSP();
    CHKERRQ(ierr);
    ierr = s.solve();
    CHKERRQ(ierr);

    DM davCoarse, dapCoarse;
    DMDALocalInfo infopCoarse, infovCoarse;
    ierr = DMCompositeGetEntries(st.ctx->dm, &davCoarse, &dapCoarse);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(davCoarse, &infovCoarse);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dapCoarse, &infopCoarse);CHKERRQ(ierr);
    cout << "Size Coarse Mesh.   Pressure : " << infopCoarse.mx << " " << infopCoarse.my << " Velocity : " << infovCoarse.mx << " " << infovCoarse.my  << "\n";


    ierr = cafes::io::save_VTK("Resultats", "two_part", st.sol, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);
    std::string stout0 = "two_part_";
    stout0.append(std::to_string(mx));
    const char * stw0 = stout0.c_str();
    ierr = cafes::io::save_VTK("Resultats", stw0, s.sol_tmp, st.ctx->dm,
                               st.ctx->h);
    CHKERRQ(ierr);
    // ierr = cafes::io::saveParticles("Resultats", "two_part",
    // pt);CHKERRQ(ierr);

    // Refinement and interpolation test
    int const fine = nref/mx;
    std::array<int, dim> refine = {fine,fine};
    cout << "Refinment factors : " << refine[0] << " " << refine[1] << "\n";
    Vec sol_refine;
    DM dm_refine;
    std::array<double, dim> h_refine;
    h_refine[0] = st.ctx->h[0]/refine[0];
    h_refine[1] = st.ctx->h[1]/refine[1];
    //std::cout << h_refine[0] << " " << h_refine[1] << std::endl;

    ierr = cafes::posttraitement::linear_interpolation(st.ctx->dm, s.sol_tmp, dm_refine, sol_refine, refine, st.ctx->h);CHKERRQ(ierr);
    
    std::string stout1 = "two_part_refined_";
    stout1.append(std::to_string(mx));
    const char * stw1 = stout1.c_str();
    ierr = cafes::io::save_VTK("Resultats", stw1, sol_refine, dm_refine, h_refine);CHKERRQ(ierr);

    // READ REFERENCE SOLUTION
    Vec solref;
    DM dav, dap;
    DMDALocalInfo infop, infov;
    ifstream file;
    std::stringstream op, ov;
    double tmp;

    ierr = DMCompositeGetEntries(dm_refine, &dav, &dap);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dav, &infov);CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(dap, &infop);CHKERRQ(ierr);
    cout << "Size Reference Mesh.   Pressure : " << infop.mx << " " << infop.my << " Velocity : " << infov.mx << " " << infov.my  << "\n";

    ierr = DMCreateGlobalVector(dm_refine, &solref);CHKERRQ(ierr);
    auto solrefu = cafes::petsc::petsc_vec<dim>(dm_refine, solref, 0);
    auto solrefp = cafes::petsc::petsc_vec<dim>(dm_refine, solref, 1);
    ierr = solrefu.global_to_local(INSERT_VALUES);CHKERRQ(ierr);
    ierr = solrefp.global_to_local(INSERT_VALUES);CHKERRQ(ierr);

    // Fill pressure
    op << "Reference/Rsur5_nref120_sameVelocity" << "/" <<  "pressure.txt";
    file.open(op.str());
    file >> tmp; file >> tmp; file >> tmp;
    for (std::size_t j =0; j<infop.my; j++ )
    {
        for (std::size_t i = 0; i < infop.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto u = solrefp.at(pos);
            file >> u[0];
            //cout << u[0] << endl;
        }
    }
    file.close();

    // Fill velocity
    ov << "Reference/Rsur5_nref120_sameVelocity" << "/" <<  "velocity.txt";
    file.open(ov.str());
    file >> tmp; file >> tmp; file >> tmp;
    for (std::size_t j =0; j<infov.my; j++ )
    {
        for (std::size_t i = 0; i < infov.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto u = solrefu.at(pos);
            for (std::size_t d=0; d<dim; d++)
            {
                file >> u[d];
            }
        }
    }
    file.close();

    ierr = solrefu.local_to_global(INSERT_VALUES);CHKERRQ(ierr);
    ierr = solrefp.local_to_global(INSERT_VALUES);CHKERRQ(ierr);

    ierr = cafes::io::save_VTK("Resultats", "reference_solution_256", solref, dm_refine, h_refine);CHKERRQ(ierr);

    // COMPUTE ERROR
    Vec error;
    ierr = DMCreateGlobalVector(dm_refine, &error);CHKERRQ(ierr);
    auto erroru = cafes::petsc::petsc_vec<dim>(dm_refine, error, 0, false);
    auto solu = cafes::petsc::petsc_vec<dim>(dm_refine, sol_refine, 0, true);
    auto errorp = cafes::petsc::petsc_vec<dim>(dm_refine, error, 1, false);
    auto solp = cafes::petsc::petsc_vec<dim>(dm_refine, sol_refine, 1, true);

    for (std::size_t j =0; j< infov.my; j++ )
    {
        for (std::size_t i = 0; i <infov.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto pts = cafes::geometry::position<double, dim>{i*h_refine[0], j*h_refine[1]};
            auto u = erroru.at_g(pos);
            auto uref = solrefu.at_g(pos);
            auto usol = solu.at_g(pos);
            if (!se1.contains(pts) && !se2.contains(pts))
            {
                u[0] = uref[0] - usol[0];
                u[1] = uref[1] - usol[1];
            }
            else
            {
                u[0] = 0.;
                u[1] = 0.;
            }
        }
    }

    for (std::size_t j =0; j<infop.my; j++ )
    {
        for (std::size_t i = 0; i <infop.mx; i++)
        {
            auto pos = cafes::geometry::position<int, dim>{i,j};
            auto pts = cafes::geometry::position<double, dim>{2*i*h_refine[0], 2*j*h_refine[1]};
            auto u = errorp.at_g(pos);
            auto uref = solrefp.at_g(pos);
            auto usol = solp.at_g(pos);
            if (!se1.contains(pts) && !se2.contains(pts))
            {
                u[0] = uref[0] - usol[0];
            }
            else
            {
                u[0] = 0.;
            }
        }
    }

    std::string stout2 = "error_";
    stout2.append(std::to_string(mx));
    const char * stw2 = stout2.c_str();
    ierr = cafes::io::save_VTK("Resultats", stw2, error, dm_refine, h_refine);CHKERRQ(ierr);

    // COMPUTE L2 AND H1 ERRORS
    
    // INTEGRATION KERNEL
    auto const kernel = [](auto& errorL2, const auto& error, const auto& h)
    {
    auto const kernel_pos = [&](auto const &pos)
    {
        auto ielem = cafes::fem::get_element(pos);
        for (std::size_t k=0; k<ielem.size(); k++)
        {
            auto u = error.at(ielem[k]);
            for (std::size_t d=0; d<dim; d++)
            {
                errorL2 += .25*u[d]*u[d]*h[0]*h[1];
            }
        }
    };
    return kernel_pos;
    };

    // VELOCITY INTEGRATION
    double erroruL2=0.;
    auto box = cafes::fem::get_DM_bounds<dim>(dav);
    auto errortest = cafes::petsc::petsc_vec<dim>(dm_refine, error, 0, true);
    cafes::algorithm::iterate(box, kernel(erroruL2, errortest, h_refine));
    erroruL2 = std::sqrt(erroruL2);
    cout << "L2 error velocity : " << erroruL2 << "\n";

    // PRESSURE INTEGRATION
    double errorpL2=0.;
    auto boxp = cafes::fem::get_DM_bounds<dim>(dap);
    auto errortestp = cafes::petsc::petsc_vec<dim>(dm_refine, error, 1, true);
    cafes::algorithm::iterate(boxp, kernel(errorpL2, errortestp, h_refine));
    errorpL2 = std::sqrt(errorpL2);
    cout << "L2 error pressure : " << errorpL2 << "\n";

    std::string stout3 = "Resultats/L2errors_";
    std::string stmx = std::to_string(mx);
    stmx = std::string(3 - stmx.length(), '0').append(stmx);
    stout3.append(stmx);
    stout3.append(".txt");
    const char * stw3 = stout3.c_str();
    std::ofstream ofs (stw3, std::ofstream::out);
    ofs << st.ctx->h[0] << " " << erroruL2 << " " << errorpL2 <<  " " << s.kspiter << endl;
    ofs.close();

    // // OTHER WAY TO INTEGRATE
    // erroruL2 = 0;
    // for (std::size_t j=0; j< infov.my-1; j++)
    // {
    //     for (std::size_t i = 0; i < infov.mx-1; i++)
    //     {
    //         auto pos = cafes::geometry::position<int,dim>{i, j};
    //         auto ielem = cafes::fem::get_element(pos);
    //         for (std::size_t k=0; k<ielem.size(); k++)
    //         {
    //             auto u = erroru.at(ielem[k]);
    //             for (std::size_t d=0; d<dim; d++)
    //             {
    //                 erroruL2 += .25*u[d]*u[d]*h_refine[0]*h_refine[1];
    //             }
    //         }
            
    //     }
    // }
    // erroruL2 = std::sqrt(erroruL2);
    // cout << "L2 error velocity : " << erroruL2 << "\n";

    // ierr = cafes::io::read_solution_TXT("./Reference/Rsur5_nref25", solref, dmref, npoints); //256 cells in x and y in pressure

    ierr = PetscFinalize();
    CHKERRQ(ierr);

    return 0;
}
