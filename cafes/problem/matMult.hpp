#ifndef PARTICLE_PROBLEM_MATMULT_HPP_INCLUDED
#define PARTICLE_PROBLEM_MATMULT_HPP_INCLUDED

#include <petsc.h>
#include <iostream>

struct context
{
    int a, b;
};


namespace cafes
{
    namespace problem
    {
    
        template<std::size_t Dimensions>
        PetscErrorCode matMult_n(Mat A, Vec X, Vec y)
        {
            PetscErrorCode ierr;
            context *ctx;
            PetscFunctionBegin;

            ierr = MatShellGetContext(A, &ctx);CHKERRQ(ierr);
            std::cout << Dimensions << "\n";
            std::cout << ctx->a << ", " << ctx->b << "\n";

            PetscFunctionReturn(0);
        }

        using petsc_fn_t = void(*)();

        template<std::size_t Dimensions>
        petsc_fn_t matMult()
        {
          return petsc_fn_t(&matMult_n<Dimensions>);
        }
    }
}
#endif