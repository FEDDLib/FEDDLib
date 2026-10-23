#!/bin/bash -l
# Trilinos 16.1.0 for FEDDLib on macOS / Apple silicon (M4), as used for the
# local install ~/opt/trilinos-16.1.0 (2026-09). It is the macOS variant of
# the cluster's Trilinos configuration: same package set and preprocessor
# flags, with these substitutions for what the Mac provides:
#   - compilers: Homebrew Open MPI wrappers (Apple clang, gfortran);
#   - BLAS/LAPACK: Apple Accelerate via spack's veclibfort shim (no MKL on
#     Apple silicon; do not use spack's netlib-lapack, see
#     install-dependencies-on-macOS-AppleSilicon-via-spack.sh);
#   - no PARDISO_MKL (use Amesos2's KLU2, "klu2", instead of "pardisomkl");
#   - HDF5 and ParMETIS are spack builds against the same Homebrew Open MPI;
#   - tests off (disk space).
# Tpetra supports a single global ordinal type per build: this one uses
# long long (FEDDLib).
# Rythmos no longer exists in 16.1.0; PackagesList.cmake lists it as allowed
# missing, so the enable is kept only for parity with the cluster script.
#
# Usage: run from an empty build directory, then `ninja install`.
# Dependencies: install-dependencies-on-macOS-AppleSilicon-via-spack.sh.
TYPE=RELEASE

BASE_DIR=${TRILINOS_SOURCE_DIR:-$HOME/dev/Trilinos/Trilinos-16.1.0} # trilinos-release-16-1-0
INSTALL_DIR=${TRILINOS_INSTALL_DIR:-$HOME/opt/trilinos-16.1.0}

export PATH=$HOME/dev/spack/bin:$PATH

MPI_BIN_DIR=$(dirname $(which mpiexec))

HDF5=$(spack location -i hdf5+mpi ^openmpi)
BOOST=$(spack location -i boost)
PARMETIS=$(spack location -i parmetis ^openmpi)
METIS=$(spack location -i metis)
VECLIBFORT=$(spack location -i veclibfort)

rm -rf CMake*

cmake \
    -G Ninja \
    -D CMAKE_BUILD_TYPE:STRING=$TYPE \
    -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
    -D CMAKE_C_FLAGS:STRING="-O3 -Wno-format-security -DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DNDEBUG  -Wno-deprecated" \
    -D CMAKE_CXX_FLAGS:STRING="-O3 -Wno-invalid-specialization -Wno-format-security -DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DNDEBUG  -Wno-deprecated" \
    -D CMAKE_Fortran_FLAGS:STRING="-O2" \
    -D MPI_BIN_DIR=$MPI_BIN_DIR \
    -D CMAKE_C_COMPILER=mpicc \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_Fortran_COMPILER=mpifort \
    -D CMAKE_CXX_STANDARD:STRING=17 \
    -D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
    -D Trilinos_ENABLE_Fortran:BOOL=ON \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D Trilinos_ENABLE_Amesos:BOOL=ON \
    -D Trilinos_ENABLE_Amesos2:BOOL=ON \
    -D Trilinos_ENABLE_Anasazi:BOOL=ON \
    -D Trilinos_ENABLE_AztecOO:BOOL=ON \
    -D Trilinos_ENABLE_Belos:BOOL=ON \
    -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
    -D Trilinos_ENABLE_Isorropia:BOOL=ON \
    -D Trilinos_ENABLE_Kokkos:BOOL=ON \
    -D Trilinos_ENABLE_ML:BOOL=ON \
    -D Trilinos_ENABLE_MueLu:BOOL=ON \
    -D Trilinos_ENABLE_NOX:BOOL=ON \
    -D Trilinos_ENABLE_OpenMP:BOOL=OFF \
    -D Trilinos_ENABLE_Rythmos:BOOL=ON \
    -D Trilinos_ENABLE_ShyLU_DD:BOOL=ON \
    -D Trilinos_ENABLE_ShyLU_DDCore:BOOL=ON \
    -D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
    -D Trilinos_ENABLE_Teuchos:BOOL=ON \
    -D Trilinos_ENABLE_Thyra:BOOL=ON \
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \
    -D Trilinos_ENABLE_Zoltan:BOOL=ON \
    -D Trilinos_ENABLE_Zoltan2:BOOL=ON \
    -D TPL_FIND_SHARED_LIBS:BOOL=ON \
    -D TPL_ENABLE_DLlib:BOOL=OFF \
    -D TPL_ENABLE_Pthread:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D TPL_ENABLE_LAPACK:BOOL=ON \
    -D TPL_ENABLE_Matio:BOOL=OFF \
    -D TPL_ENABLE_METIS:BOOL=ON \
    -D TPL_ENABLE_BLAS:BOOL=ON \
    -D TPL_BLAS_LIBRARIES:STRING="$VECLIBFORT/lib/libvecLibFort.dylib" \
    -D TPL_LAPACK_LIBRARIES:STRING="$VECLIBFORT/lib/libvecLibFort.dylib" \
    -D METIS_LIBRARY_DIRS:PATH=$METIS/lib \
    -D METIS_INCLUDE_DIRS:PATH=$METIS/include \
    -D TPL_ENABLE_ParMETIS:BOOL=ON \
    -D ParMETIS_LIBRARY_DIRS:PATH="$PARMETIS/lib;$METIS/lib"  \
    -D ParMETIS_INCLUDE_DIRS:PATH="$PARMETIS/include;$METIS/include" \
    -D Amesos_ENABLE_PARDISO_MKL:BOOL=OFF \
    -D Amesos2_ENABLE_PARDISO_MKL:BOOL=OFF \
    -D TPL_ENABLE_PARDISO_MKL:BOOL=OFF \
    -D TPL_ENABLE_MKL:BOOL=OFF \
    -D Boost_INCLUDE_DIRS:PATH=$BOOST/include \
    -D EpetraExt_USING_HDF5:BOOL=ON \
    -D TPL_ENABLE_HDF5:BOOL=ON \
    -D HDF5_LIBRARY_DIRS:PATH=$HDF5/lib \
    -D HDF5_INCLUDE_DIRS:PATH=$HDF5/include \
    -D Trilinos_ENABLE_TESTS:BOOL=OFF \
    -D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
    -D NOX_ENABLE_TESTS:BOOL=OFF \
    -D Trilinos_ENABLE_CONFIGURE_TIMING=ON \
    -D Trilinos_ENABLE_PACKAGE_CONFIGURE_TIMING=ON \
$BASE_DIR
