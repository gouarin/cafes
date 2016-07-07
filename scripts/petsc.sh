#!/bin/sh
set -ex

CURRENT_DIR=`pwd`
PETSC_DIR=$HOME/petsc
PETSC_ARCH=arch-linux2-c-opt

echo $CURRENT_DIR
echo $1

mkdir -p $PETSC_DIR
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-$1.tar.gz
tar xzf petsc-lite-$1.tar.gz --strip-components=1 -C $PETSC_DIR
cd $PETSC_DIR
./configure --with-cc=mpicc --with-cxx=mpicxx --with-debugging=0 --with-dynamic-libraries --with-fc=0 --with-shared-libraries
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSC_ARCH all
cd $CURRENT_DIR