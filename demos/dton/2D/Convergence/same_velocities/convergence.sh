#!/bin/bash

exec=../../../../../build/demos/convergence_same_velocities
nproc=1

for i in  '9' '17' '33' '65' '129' #'11' '21' '26' '51' '101' #'201' #'126' '201' '251' #'101' '126' '201' '251'
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
    #-xperiod \
    #-dton_ksp_initial_guess_nonzero true \
    #-dton_ksp_norm_type UNPRECONDITIONED \
    # -strain_tensor \
    #-on_error_attach_debugger lldb \
    #-dton_ksp_rtol 1e-4 \

    
done
#-dton_ksp_atol 5.e-4 \

cat Results/L2errors* > Results/convergence.txt
rm -f Results/L2errors*


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
