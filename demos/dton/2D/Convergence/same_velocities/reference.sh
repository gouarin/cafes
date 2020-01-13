#!/bin/bash

exec=../../../../../build/demos/convergence_same_velocities_reference
nproc=1

for i in '257'
do
mpirun -np $nproc $exec \
    -assembling \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -stokes_ksp_monitor_true_residual \
    -dton_ksp_type lgmres \
    -dton_ksp_max_it 500 \
    -dton_ksp_rtol 1e-10 \
    -dton_ksp_monitor_true_residual \
    -mx $i \
    -my $i \
    # 
    # -PMM \
    #-stokes_ksp_rtol 1e-6 \

done