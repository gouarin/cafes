cmake_minimum_required(VERSION 3.0)

PROJECT(cafes LANGUAGES CXX C)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/" ${CMAKE_MODULE_PATH})

include_directories(cafes external_packages/qhull)

set (CMAKE_CXX_STANDARD 14)

include(FindPkgConfig)
pkg_check_modules(PETSC PETSc)

FIND_PACKAGE(MPI)
Find_package(HDF5)

include_directories(${PETSC_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS})
link_directories(${PETSC_LIBRARY_DIRS})

ADD_SUBDIRECTORY(external_packages/qhull)
ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(demos)
