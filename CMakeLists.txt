cmake_minimum_required(VERSION 3.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

SET(CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_DEBUG} --coverage")
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE "${CMAKE_EXE_LINKER_FLAGS_DEBUG} --coverage")
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} --coverage")

include_directories(cafes external_packages/qhull)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(EXT_PROJECTS_DIR ./external_packages)

FIND_PACKAGE(PETSc)
FIND_PACKAGE(MPI)
FIND_PACKAGE(VTK REQUIRED)
#Find_package(HDF5)

INCLUDE (${VTK_USE_FILE})
#include_directories(${HDF5_C_INCLUDE_DIR} ${PETSC_INCLUDE_CONF} ${PETSC_INCLUDE_DIR} ${MPI_INCLUDE_PATH})
include_directories(${PETSC_INCLUDE_CONF} ${PETSC_INCLUDE_DIR} ${MPI_INCLUDE_PATH})

# Includes Catch in the project:
add_subdirectory(${EXT_PROJECTS_DIR}/catch)
include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
enable_testing(true)  # Enables unit-testing.

ADD_SUBDIRECTORY(${EXT_PROJECTS_DIR}/qhull)
ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(demos)
ADD_SUBDIRECTORY(cafes/particle/geometry/tests)

