#!/bin/sh
set -ex

CURRENT_DIR=`pwd`

VERSION=$1
PETSC_DIR=$CACHE_DIRECTORY/petsc-$VERSION
PETSC_ARCH=arch-linux2-c-opt

if [ ! -d $PETSC_DIR ]; then
  mkdir -p $PETSC_DIR
  wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-$1.tar.gz
  tar xzf petsc-lite-$1.tar.gz --strip-components=1 -C $PETSC_DIR
  cd $PETSC_DIR
  ./configure --with-cc=mpicc --with-cxx=mpicxx --with-debugging=0 --with-fc=0
  make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
  cd $CURRENT_DIR
fi