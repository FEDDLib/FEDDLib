#!/bin/bash
# FEDDLib on macOS / Apple silicon (M4) against the Trilinos 16.1.0 from
# configure-Trilinos-on-macOS-AppleSilicon.sh (Homebrew Open MPI, Apple
# Accelerate via veclibfort, Open MPI-built spack HDF5), with the Interface2
# AceGen interface, as used for the local builds of 2026-09.
#
# Usage: run from an empty build directory, then `ninja` (and `ninja install`).
# Options (environment variables):
#   ACEGEN=OFF   build without the Interface2 AceGen interface (the tests that
#                need it are skipped)
#   ASAN=ON      AddressSanitizer build. Only CMAKE_CXX_FLAGS gets the flags:
#                setting them in the linker flags breaks CMake's mpifort
#                check. Run with ASAN_OPTIONS=detect_container_overflow=0,
#                since the uninstrumented Trilinos shares libc++ containers
#                (false container-overflow reports in Teuchos XML parsing).
#   BUILD_TYPE, TRILINOS_DIR, INSTALL_DIR, INTERFACE_DIR, HDF5_DIR, ZLIB_DIR
#
# Interface2 (https://github.com/nrsharan/Interface2) is built with its
# config_scripts/do-config-interface2-macos-accelerate.sh (installs to
# ~/opt/interface2). Without MKL, replace "pardisomkl" by "klu2" in the
# parameter files of the tests/examples.

BUILD_TYPE=${BUILD_TYPE:-RelWithDebInfo}

SOURCE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
TRILINOS_DIR=${TRILINOS_DIR:-$HOME/opt/trilinos-16.1.0}
INSTALL_DIR=${INSTALL_DIR:-$HOME/opt/FEDDLib}
INTERFACE_DIR=${INTERFACE_DIR:-$HOME/opt/interface2}

export PATH=$HOME/dev/spack/bin:$PATH
HDF5_DIR=${HDF5_DIR:-$(spack location -i hdf5+mpi ^openmpi)}
# the zlib the spack HDF5 was built against
ZLIB_DIR=${ZLIB_DIR:-$(spack location -i zlib-ng)}

CXX_FLAGS="-D HAVE_EXPLICIT_INSTANTIATION -Wno-deprecated -Wno-sign-compare -Wno-unused-variable -Wno-write-strings"
if [ "${ASAN:-OFF}" = ON ]; then
  CXX_FLAGS="$CXX_FLAGS -fsanitize=address -fno-omit-frame-pointer -g"
fi

if [ "${ACEGEN:-ON}" = ON ]; then
  ACEGEN_ARGS=(-D TPL_ENABLE_AceGENInterface:BOOL=ON
               -D TPL_AceGENInterface_LIBRARIES:STRING="$INTERFACE_DIR/lib/libinterface2.a;"
               -D TPL_AceGENInterface_INCLUDE_DIRS:STRING="$INTERFACE_DIR/include")
else
  ACEGEN_ARGS=(-D TPL_ENABLE_AceGENInterface:BOOL=OFF)
fi

rm -rf CMake*
cmake \
-G Ninja \
-D CMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
-D CMAKE_C_COMPILER=mpicc \
-D CMAKE_CXX_COMPILER=mpicxx \
-D CMAKE_Fortran_COMPILER=mpifort \
-D CMAKE_INSTALL_PREFIX:STRING=${INSTALL_DIR} \
-D CMAKE_CXX_FLAGS:STRING="${CXX_FLAGS}" \
-D CMAKE_CXX_STANDARD:STRING=17 \
-D MPI_EXEC_MAX_NUMPROCS:STRING=6 \
-D FEDDLib_ENABLE_ALL_PACKAGES:BOOL=ON \
-D FEDDLib_ENABLE_TESTS:BOOL=ON \
-D TPL_FIND_SHARED_LIBS:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D TPL_ENABLE_Trilinos:BOOL=ON \
-D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
-D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib \
-D TPL_HDF5_INCLUDE_DIRS:PATH=$HDF5_DIR/include \
-D TPL_HDF5_LIBRARY_DIRS:PATH=$HDF5_DIR/lib \
-D TPL_Z_INCLUDE_DIRS:PATH=$ZLIB_DIR/include \
-D TPL_Z_LIBRARY_DIRS:PATH=$ZLIB_DIR/lib \
"${ACEGEN_ARGS[@]}" \
${SOURCE_DIR}
