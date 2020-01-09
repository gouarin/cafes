#!/bin/bash

exec=../../../../../build/demos/convergence_same_velocities_reference
nproc=1

for i in '401'
do
mpirun -np $nproc $exec \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -dton_ksp_type lgmres \
    -dton_ksp_max_it 500 \
    -dton_ksp_rtol 1e-5 \
    -dton_ksp_monitor_true_residual \
    -assembling \
    -mx $i \
    -my $i \

done