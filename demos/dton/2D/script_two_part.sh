#!/bin/bash

exec=../../../build/demos/two_part
nproc=4

mpirun -np $nproc $exec \
    -strain_tensor \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -dton_ksp_type gmres \
    -dton_ksp_monitor_true_residual \
    -dton_ksp_max_it 100 \
    -assembling \
    -dton_ksp_norm_type UNPRECONDITIONED \
    -mx 65 \
    -my 65


# mpirun -np $nproc $exec \
#     -pmm \
#     -strain_tensor \
#     -stokes_ksp_rtol 1e-2 \
#     -stokes_ksp_atol 1e-8 \
#     -dton_ksp_type gmres \
#     -dton_ksp_monitor_true_residual \
#     -dton_ksp_max_it 10 \
#     -stokes_ksp_monitor \
#     -dton_ksp_norm_type UNPRECONDITIONED \
#     -mx 65 \
#     -my 65
