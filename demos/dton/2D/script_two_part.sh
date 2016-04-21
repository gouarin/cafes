#!/bin/bash

exec=../../../build/demos/two_part
nproc=1

mpiexec -n $nproc $exec \
    -pmm \
    -strain_tensor \
    -stokes_ksp_rtol 1e-6 \
    -dton_ksp_rtol 1e-2 \
    -dton_ksp_type gcr \
    -dton_ksp_monitor \
    -mx 65 \
    -my 65