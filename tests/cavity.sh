#!/bin/bash

#exec=../build/tests/stokes_analytic_2d
exec=../build/tests/cavity2d
nproc=1

# for i in '257' '513' '1025'
# for i in '33' '65' '129' '257'
# do
#     $exec \
#     -pmm \
#     -stokes_ksp_monitor \
#     -mx $i \
#     -my $i
# done

#for i in '33' '65' '129'
for i in '33'
do
    $exec \
    -assembling \
    -stokes_ksp_type preonly \
    -stokes_pc_type lu \
    -mx $i \
    -my $i
done
