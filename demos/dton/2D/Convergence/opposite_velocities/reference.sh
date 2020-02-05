#!/bin/bash

exec=../../../../../build/demos/convergence_opposite_velocities_reference
nproc=2

for i in '129'
do
mpirun -np $nproc $exec \
    -assembling \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -stokes_ksp_monitor_true_residual \
    -stokes_pc_factor_mat_solver_type superlu_dist \
    -dton_ksp_type lgmres \
    -dton_ksp_max_it 10000 \
    -dton_ksp_rtol 1e-5 \
    -dton_ksp_monitor_true_residual \
    -mx $i \
    -my $i \
    # 
    # -PMM \
    #-stokes_ksp_rtol 1e-6 \

done
