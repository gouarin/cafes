#!/bin/bash

exec=../../../../build/demos/convergence
nproc=1

for i in '16' '21' '26' '31' '51' '61' '76' '101' '151'
do
mpirun -np $nproc $exec \
    -strain_tensor \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -dton_ksp_type minres \
    -dton_ksp_monitor_true_residual \
    -dton_ksp_max_it 100 \
    -assembling \
    -mx $i \
    -my $i
    #-dton_ksp_norm_type UNPRECONDITIONED \
done
#-dton_ksp_atol 5.e-4 \

cat Resultats/L2errors* > Resultats/convergence.txt
rm -f Resultats/L2errors*


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
