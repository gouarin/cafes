#!/bin/bash

exec=../../../../../build/demos/convergence_opposite_velocities_reference
nx=2049

$exec \
-assembling \
-stokes_ksp_type preonly \
-stokes_pc_type lu \
-stokes_pc_factor_mat_solver_type mumps \
-mat_mumps_icntl_2 1 \
-mat_mumps_icntl_4 2 \
-dton_ksp_type lgmres \
-dton_ksp_max_it 10000 \
-dton_ksp_rtol 1e-6 \
-dton_ksp_monitor_true_residual \
-mx $nx \
-my $nx \
