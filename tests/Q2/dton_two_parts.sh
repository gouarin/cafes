#!/bin/bash

exec=../../build/tests/dton_two_parts
nproc=1

for i in '21' '25' '31' '35' '41' '45' '51' '55' '61'
do
$exec \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -dton_ksp_monitor_true_residual \
    -dton_ksp_rtol 1e-5 \
    -dton_ksp_type gmres \
    -assembling \
    -order 2 \
    -mx $i \
    -my $i 
done
