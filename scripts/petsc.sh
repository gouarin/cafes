#!/bin/sh
set -ex

CURRENT_DIR=`pwd`

#VERSION="3.6.4 3.7.2"
VERSION="3.7.2"
PETSC_ARCH=arch-linux2-c-opt

for v in $VERSION
do
  PETSC_DIR=$CACHE_DIRECTORY/petsc-$v
  echo "build petsc $v"
  if [ ! -d $PETSC_DIR ]; then
    mkdir -p $PETSC_DIR
    curl -O http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-$v.tar.gz
    tar xzf petsc-lite-$v.tar.gz --strip-components=1 -C $PETSC_DIR
    cd $PETSC_DIR
    ./configure --with-cc=mpicc --with-cxx=mpicxx --with-debugging=0 --with-fc=0
    make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
  fi
done
cd $CURRENT_DIR
