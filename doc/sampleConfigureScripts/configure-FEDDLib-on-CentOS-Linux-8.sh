#!/bin/bash

BUILD_TYPE=DEBUG

MPI_C_COMPILER=mpicc
MPI_CXX_COMPILER=mpicxx
MPI_FORTRAN_COMPILER=mpifort

BASE_FEDDLIB=$HOME/dev/feddlib/Jascha

TRILINOS_DIR=$HOME/dev/trilinos/feddlib/install
SOURCE_DIR=$BASE_FEDDLIB/src
INSTALL_DIR=$BASE_FEDDLIB/install
BUILD_DIR=$BASE_FEDDLIB/build
#ACEGEN_INTERFACE_DIR=$HOME/dev/feddlib/interface2/install

spack load hdf5 # hdf5@1.14.6

# Note on MUMPS vs HDF5 vs Z library
# With MUMPS, zlib is loaded, which shadows the system libraries that are used for HDF5 if that is compiled independently. Thus, I have built HDF5 via Spack as well to avoid inconsistencies.

# Note on HDF5
# It requires libdl and libz. We instruct the linker to include them via -ldl and -lz below.

### Compiler flags
# Macros are set s.t. "#ifdef MACRO" can be used in the code.
MACROS="-D HAVE_EXPLICIT_INSTANTIATION"
# FLAGS is always used
# FLAGS_RELEASE is used only if build type is RELEASE
# FLAGS_DEBUG is used only if build type is DEBUG
FLAGS="$MACROS -Wno-cpp -Wno-deprecated-declarations -ldl -lz" # -Wno-deprecated -Wno-sign-compare -Wno-unused-variable -fpermissive"
#FLAGS = "$MACROS"
FLAGS_RELEASE="-O3 -march=native"
#FLAGS_RELEASE=""
FLAGS_DEBUG="-O0 -ggdb -debug all -traceback -ftrapuv"
#FLAGS_DEBUG=""

rm -rf $BUILD_DIR/CMake*

#--trace-expand
cmake \
-D CMAKE_EXPORT_COMPILE_COMMANDS:BOOL=True \
-D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
-D CMAKE_C_COMPILER=$MPI_C_COMPILER \
-D CMAKE_CXX_COMPILER=$MPI_CXX_COMPILER \
-D CMAKE_Fortran_COMPILER=$MPI_FORTRAN_COMPILER \
-D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
-D CMAKE_C_FLAGS:STRING="$FLAGS" \
-D CMAKE_CXX_FLAGS:STRING="$FLAGS" \
-D CMAKE_C_FLAGS_RELEASE:STRING="$FLAGS_RELEASE" \
-D CMAKE_CXX_FLAGS_RELEASE:STRING="$FLAGS_RELEASE" \
-D CMAKE_C_FLAGS_DEBUG:STRING="$FLAGS_DEBUG" \
-D CMAKE_CXX_FLAGS_DEBUG:STRING="$FLAGS_DEBUG" \
-D CMAKE_CXX_STANDARD=17 \
-D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
-D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
-D FEDDLib_ENABLE_ALL_PACKAGES:BOOL=ON \
-D FEDDLib_ENABLE_TESTS:BOOL=ON \
-D TPL_FIND_SHARED_LIBS:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D TPL_ENABLE_Trilinos:BOOL=ON \
-D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
-D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib64 \
-D TPL_ENABLE_AceGENInterface:BOOL=OFF \
-S ${SOURCE_DIR} -B ${BUILD_DIR}
# > trace.txt 2>&1

# MPI_EXEC_MAX_NUMPROCS: Tests requiring more than the specified number of processes are excluded from CTest.

#-D TPL_AceGENInterface_LIBRARIES:STRING="$ACEGEN_INTERFACE_DIR/lib/libinterface2.a;$ACEGEN_INTERFACE_DIR/lib/libaceutility.a" \
#-D TPL_AceGENInterface_INCLUDE_DIRS:STRING="$ACEGEN_INTERFACE_DIR/include" \


