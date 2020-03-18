#!/bin/bash

exec=../../../../../build/demos/convergence
nproc=1

for i in '33' '65' '129' '257' '513' '1025'
do
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
    -mx $i \
    -my $i \
    -compute_singularity
done

    # -assembling \
    # -stokes_ksp_type preonly \
    # -stokes_pc_type lu \
    # -dton_ksp_type lgmres \
    # -dton_ksp_max_it 10 \
    # -dton_ksp_rtol 1e-10 \
    # -dton_ksp_monitor_true_residual \
    # -mx $i \
    # -my $i
    #-stokes_ksp_type preonly \
    #-stokes_pc_type lu \
    #-stokes_pc_factor_mat_solver_type mumps \
    #-mat_superlu_dist_replacetinypivot \
    #-dton_ksp_initial_guess_nonzero true \
    #-assembling \
    #-dton_ksp_type gmres \
    #-dton_ksp_max_it 5 \
    #-dton_ksp_rtol 1e-5 \
    #-dton_ksp_monitor_true_residual \
    #-mx $i \
    #-my $i \
    #-stokes_pc_factor_mat_solver_type mumps \
    #-start_in_debugger lldb \
    # -xperiod \
    # -yperiod \
    # -stokes_ksp_type preonly \
    # -stokes_pc_type lu \
    #-xperiod \
    #-dton_ksp_initial_guess_nonzero true \
    #-dton_ksp_norm_type UNPRECONDITIONED \
    # -strain_tensor \
    #-on_error_attach_debugger lldb \
    #-dton_ksp_rtol 1e-4 \
    # -assembling \
    # -stokes_ksp_type preonly \
    # -stokes_pc_type lu \

    
# done

#cat Resultats/L2errors* > Resultats/convergence.txt
#rm -f Resultats/L2errors*


# mpirun -np $nproc $exec \
#     -pmm \
#     -stokes_ksp_rtol 1e-2 \
#     -stokes_ksp_atol 1e-8 \
#     -dton_ksp_type gmres \
#     -dton_ksp_monitor_true_residual \
#     -dton_ksp_max_it 5000 \
#     -stokes_ksp_monitor \
#     -dton_ksp_norm_type UNPRECONDITIONED \
#     -mx 65 \
#     -my 65 \
#     #-strain_tensor \
