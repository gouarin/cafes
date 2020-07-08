// Copyright (c) 2016, Loic Gouarin <loic.gouarin@math.u-psud.fr>
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef PARTICLE_PROBLEM_OPTIONS_HPP_INCLUDED
#define PARTICLE_PROBLEM_OPTIONS_HPP_INCLUDED

#include <array>
#include <petsc.h>

namespace cafes
{
    namespace problem
    {
        template<std::size_t Dimensions>
        struct options
        {
            std::array<PetscInt, Dimensions> mx;
            std::array<double, Dimensions> lx;
            std::array<PetscBool, Dimensions> xperiod;
            PetscInt order = 1;
            PetscBool strain_tensor = PETSC_FALSE;
            PetscBool assembling = PETSC_FALSE;
            PetscBool pmm = PETSC_FALSE;
            PetscBool compute_singularity = PETSC_FALSE;

            options()
            {
                mx.fill(17);
                lx.fill(1.);
                xperiod.fill(PETSC_FALSE);
            }

#undef __FUNCT__
#define __FUNCT__ "process_options_"
            PetscErrorCode
            process_options_(std::integral_constant<int, 2> const &)
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "",
                                         "Problem Options", "");
                CHKERRQ(ierr);
                ierr = PetscOptionsInt("-mx",
                                       "The number of points in x direction",
                                       "options.hpp", mx[0], &mx[0], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsInt("-my",
                                       "The number of points in y direction",
                                       "options.hpp", mx[1], &mx[1], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsReal("-Lx", "The size of x", "options.hpp",
                                        lx[0], &lx[0], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsReal("-Ly", "The size of y", "options.hpp",
                                        lx[1], &lx[1], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsInt("-order",
                                       "Discretization order for the velocity (1 or 2)",
                                       "options.hpp", order, &order, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-xperiod", "set periodic condition in x direction",
                    "options.hpp", xperiod[0], &xperiod[0], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-yperiod", "set periodic condition in y direction",
                    "options.hpp", xperiod[1], &xperiod[1], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-strain_tensor",
                    "Use strain tensor term instead of the Laplacian",
                    "options.hpp", strain_tensor, &strain_tensor, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-assembling", "Assemble the Stokes matrix", "options.hpp",
                    assembling, &assembling, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-PMM", "MG with Pressure Mass Matrix preconditionner",
                    "options.hpp", pmm, &pmm, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-compute_singularity", "Compute singularity in Stokes solver",
                    "options.hpp", compute_singularity, &compute_singularity, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsEnd();
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "process_options_"
            PetscErrorCode
            process_options_(std::integral_constant<int, 3> const &)
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "",
                                         "Problem Options", "");
                CHKERRQ(ierr);
                ierr = PetscOptionsInt("-mx",
                                       "The number of points in x direction",
                                       "options.hpp", mx[0], &mx[0], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsInt("-my",
                                       "The number of points in y direction",
                                       "options.hpp", mx[1], &mx[1], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsInt("-mz",
                                       "The number of points in z direction",
                                       "options.hpp", mx[2], &mx[2], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsReal("-Lx", "The size of x", "options.hpp",
                                        lx[0], &lx[0], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsReal("-Ly", "The size of y", "options.hpp",
                                        lx[1], &lx[1], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsReal("-Lz", "The size of z", "options.hpp",
                                        lx[2], &lx[2], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsInt("-order",
                                       "Discretization order for the velocity (1 or 2)",
                                       "options.hpp", order, &order, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-xperiod", "set periodic condition in x direction",
                    "options.hpp", xperiod[0], &xperiod[0], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-yperiod", "set periodic condition in y direction",
                    "options.hpp", xperiod[1], &xperiod[1], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-zperiod", "set periodic condition in z direction",
                    "options.hpp", xperiod[2], &xperiod[2], nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-strain_tensor",
                    "Use strain tensor term instead of the Laplacian",
                    "options.hpp", strain_tensor, &strain_tensor, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-assembling", "Assemble the Stokes matrix", "options.hpp",
                    assembling, &assembling, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsBool(
                    "-PMM", "MG with Pressure Mass Matrix preconditionner",
                    "options.hpp", pmm, &pmm, nullptr);
                CHKERRQ(ierr);
                ierr = PetscOptionsEnd();
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }

#undef __FUNCT__
#define __FUNCT__ "process_options"
            PetscErrorCode process_options()
            {
                PetscErrorCode ierr;
                PetscFunctionBeginUser;

                ierr = process_options_(std::integral_constant<int, Dimensions>{});
                CHKERRQ(ierr);

                PetscFunctionReturn(0);
            }
        };
    } // namespace problem
} // namespace cafes
#endif
