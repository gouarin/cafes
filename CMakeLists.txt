cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

include_directories(cafes external_packages/qhull)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

FIND_PACKAGE(PETSc)
FIND_PACKAGE(MPI)
FIND_PACKAGE(VTK REQUIRED)
#Find_package(HDF5)

INCLUDE (${VTK_USE_FILE})
#include_directories(${HDF5_C_INCLUDE_DIR} ${PETSC_INCLUDE_CONF} ${PETSC_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
include_directories(${PETSC_INCLUDE_CONF} ${PETSC_INCLUDE_DIR} ${MPI_INCLUDE_PATH})

ADD_SUBDIRECTORY(external_packages/qhull)
ADD_SUBDIRECTORY(tests)
