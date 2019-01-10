#!/bin/bash

exec=../../../build/demos/two_part
nproc=1

mpirun -np $nproc $exec \
    -pmm \
    -strain_tensor \
    -stokes_ksp_rtol 1e-3 \
    -stokes_ksp_atol 1e-10 \
    -dton_ksp_type fgmres \
    -dton_ksp_monitor_true_residual \
    -dton_ksp_max_it 10 \
    -stokes_ksp_monitor \
    -mx 65 \
    -my 65