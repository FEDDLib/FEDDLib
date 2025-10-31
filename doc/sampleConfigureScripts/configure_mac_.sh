#!/bin/bash
BUILD_TYPE=DEBUG

MPI_C_COMPILER=`which mpicc`
MPI_CXX_COMPILER=`which mpicxx`

TRILINOS_DIR=/Users/chris/c-code/MyTrilinosSecond/installed
SOURCE_DIR=/Users/chris/c-code/fedd/feddlib

rm -rf CMake*
cmake \
    -D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
    -D CMAKE_C_COMPILER:PATH=${MPI_C_COMPILER} \
    -D CMAKE_CXX_COMPILER:PATH=${MPI_CXX_COMPILER} \
    -D CMAKE_CXX_FLAGS:STRING="-D FROSCH_Epetra64 -D HAVE_EXPLICIT_INSTANTIATION -Wno-deprecated -Wno-sign-compare -Wno-unused-variable -Wno-pedantic" \
    -D CMAKE_CXX_STANDARD_LIBRARIES:STRING="/opt/local/lib/gcc7/libgfortran.dylib /opt/local/lib/openmpi-gcc7/libmpi_mpifh.dylib" \
    -D CMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -D MPI_BIN_DIR="/opt/local/bin/" \
    -D MPI_EXEC:FILEPATH="mpirun" \
    -D MPI_EXEC_PRE_NUMPROCS_FLAGS="--oversubscribe" \
    -D FEDDLib_ENABLE_ALL_PACKAGES:BOOL=ON \
    -D FEDDLib_ENABLE_TESTS:BOOL=ON \
    -D TPL_FIND_SHARED_LIBS:BOOL=ON \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_Trilinos:BOOL=ON \
    -D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
    -D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib \
${SOURCE_DIR}
