// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, 
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
// IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
// NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
// OF SUCH DAMAGE.

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