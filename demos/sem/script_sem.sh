#!/bin/bash

exec=../../build/demos/sem2d
nproc=1

mpirun -np $nproc $exec \
    -pmm \
    -strain_tensor \
    -stokes_ksp_rtol 1e-4 \
    -stokes_ksp_atol 1e-12 \
    -sem_ksp_rtol 1e-2 \
    -sem_ksp_type fgmres \
    -sem_ksp_monitor_true_residual \
    -sem_ksp_max_it 10 \
    -stokes_ksp_monitor \
    -mx 129 \
    -my 129