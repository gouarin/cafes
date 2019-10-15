#!/bin/bash

exec=../../../build/demos/two_part
nproc=1

mpirun -np $nproc $exec \
    -strain_tensor \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -dton_ksp_type gcr \
    -dton_ksp_monitor_true_residual \
    -dton_ksp_max_it 10 \
    -assembling \
    -dton_ksp_initial_guess_nonzero true \
    -mx 65 \
    -my 65 -help


# mpirun -np $nproc $exec \
#     -pmm \
#     -strain_tensor \
#     -stokes_ksp_rtol 1e-2 \
#     -dton_ksp_type fgmres \
#     -dton_ksp_monitor_true_residual \
#     -dton_ksp_max_it 10 \
#     -stokes_ksp_monitor \
#     -mx 65 \
#     -my 65