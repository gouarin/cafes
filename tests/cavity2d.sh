#!/bin/sh

mpirun -np 1 ./tests/cavity2d \
       -PMM \
       -mx 513 -my 513 \
       -stokes_ksp_monitor  \
       -stokes_fieldsplit_0_ksp_type gcr \
       -stokes_fieldsplit_0_ksp_rtol 1e-3 \
       -stokes_fieldsplit_0_pc_type mg \
       -stokes_fieldsplit_0_pc_mg_levels 9 \
       -stokes_fieldsplit_0_mg_levels_ksp_type richardson \
       -stokes_fieldsplit_0_mg_levels_pc_type jacobi \
       -stokes_fieldsplit_1_ksp_type cg \
       -stokes_fieldsplit_1_pc_type jacobi \
       -stokes_fieldsplit_1_ksp_rtol 1e-2 \
       -assembling \
       # -stokes_fieldsplit_0_mg_coarse_pc_type telescope \
       # -stokes_fieldsplit_0_mg_coarse_pc_telescope_reduction_factor 8 \
       # -stokes_fieldsplit_0_mg_coarse_telescope_pc_type mg \
       # -stokes_fieldsplit_0_mg_coarse_telescope_pc_mg_levels 3 \


# mpirun -np 8 ./tests/cavity2d \
#        -PMM \
#        -mx 257 -my 257 \
#        -stokes_ksp_monitor  \
#        -stokes_fieldsplit_0_ksp_type gcr \
#        -stokes_fieldsplit_0_pc_type mg \
#        -stokes_fieldsplit_0_pc_mg_levels 4 \
#        -stokes_fieldsplit_0_mg_coarse_pc_type telescope \
#        -stokes_fieldsplit_0_mg_coarse_pc_telescope_reduction_factor 5 \
#        -stokes_fieldsplit_0_mg_coarse_telescope_pc_type mg \
#        -stokes_fieldsplit_0_mg_coarse_telescope_pc_mg_levels 2 \
#        -stokes_fieldsplit_0_mg_levels_ksp_type cg \
#        -stokes_fieldsplit_0_mg_levels_pc_type jacobi \
#        -stokes_fieldsplit_0_mg_levels_ksp_rtol 1e-3 \
#        -stokes_fieldsplit_1_ksp_type fgmres \
#        -stokes_fieldsplit_1_pc_type jacobi \
#        -stokes_fieldsplit_1_ksp_rtol 1e-4 \
#        -assembling \
#        # -stokes_fieldsplit_0_ksp_monitor \
#        # -stokes_ksp_view
#        # -stokes_fieldsplit_1_ksp_monitor \

# mpirun -np 4 ./tests/cavity2d \
#        -PMM \
#        -mx 65 -my 65 \
#        -stokes_ksp_monitor  \
#        -stokes_fieldsplit_0_ksp_type gcr \
#        -stokes_fieldsplit_0_pc_type mg \
#        -stokes_fieldsplit_0_pc_mg_levels 2 \
#        -stokes_fieldsplit_0_mg_coarse_pc_type telescope \
#        -stokes_fieldsplit_0_mg_coarse_pc_telescope_reduction_factor 4 \
#        -stokes_fieldsplit_0_mg_coarse_telescope_pc_type mg \
#        -stokes_fieldsplit_0_mg_coarse_telescope_pc_mg_levels 2 \
#        -stokes_fieldsplit_1_ksp_type fgmres \
#        -stokes_fieldsplit_1_pc_type jacobi \
#        -assembling \
#        -stokes_fieldsplit_0_ksp_monitor \
#        -stokes_fieldsplit_1_ksp_monitor \
#        # -help
#        # -stokes_fieldsplit_0_mg_levels_1_ksp_type cg \
#        # -stokes_fieldsplit_0_mg_levels_1_pc_type jacobi \
#        # -stokes_ksp_view