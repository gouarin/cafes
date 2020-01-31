#!/bin/bash

exec=../../../../../build/demos/convergence
nproc=2

for i in '65' #'11' '21' '26' '51' '101' '126' '201' '251' #'101' '126' '201' '251'
do
mpirun -np $nproc -use-hwthread-cpus $exec \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -assembling \
    -strain_tensor \
    -dton_ksp_type gmres \
    -dton_ksp_max_it 5000 \
    -dton_ksp_rtol 1e-5 \
    -dton_ksp_monitor_true_residual \
    -mx $i \
    -my $i \
    -on_error_attach_debugger lldb \
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

    
done

# cat Resultats/L2errors* > Resultats/convergence.txt
# rm -f Resultats/L2errors*


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
