#include <fem/matrixFree.hpp>
#include <fem/mesh.hpp>
#include <petsc.h>
#include <array>
#include <type_traits>
#include <iostream>
#include <chrono>

template<std::size_t Dimensions, std::size_t Ndm=1>
struct context{
    template<std::size_t N> using int_ = std::integral_constant<std::size_t, N>;
    using ndm_type = int_<Ndm>;
    DM dm;
    using array1d = std::array<double, Dimensions>;
    using array2d = std::array<std::array<double, Dimensions>, Ndm>;
    typename std::conditional<Ndm == 1, array1d, array2d>::type const h;
    PetscErrorCode(*apply)(DM, Vec, Vec, std::array<double, Dimensions> const&);
};

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    ierr = PetscInitialize(&argc, &argv,  (char *)0, (char *)0);CHKERRQ(ierr);

    DM mesh;
    ierr = cafes::fem::createMesh<2>(mesh, {{1000, 1000}}, {{PETSC_FALSE, PETSC_FALSE}});CHKERRQ(ierr);

    DM dav, dap;
    ierr = DMCompositeGetEntries(mesh, &dav, &dap);CHKERRQ(ierr);

    Vec x, y;

    // STOKES
    ierr = DMCreateGlobalVector(mesh, &x);CHKERRQ(ierr);
    ierr = VecDuplicate(x, &y);CHKERRQ(ierr);
    
    using ctx = context<2>;
    ctx s{mesh, {{.1, .1}}, cafes::fem::laplacian_mult};
    Mat A = cafes::fem::make_matrix<ctx>(s, cafes::fem::stokes_matrix<ctx>);

    // // Autre
    // ierr = DMCreateGlobalVector(dav, &x);CHKERRQ(ierr);
    // ierr = VecDuplicate(x, &y);CHKERRQ(ierr);

    // using ctx = context<2>;
    // ctx s{ dav, {{.1, .1}}, cafes::fem::mass_mult};
    // Mat A = cafes::fem::make_matrix<ctx>(s, cafes::fem::diag_block_matrix<ctx>);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    for(size_t i=0; i<10; ++i){
        ierr = MatMult(A, x, y);CHKERRQ(ierr); 
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

    ierr = PetscFinalize();CHKERRQ(ierr);

    return 0;
}