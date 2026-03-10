#!/bin/bash

#TYPE=DEBUG
TYPE=RELEASE

BASE_DIR=$HOME/dev/trilinos/feddlib/src
INSTALL_DIR=$HOME/dev/trilinos/feddlib/install
#HDF5=$HOME/opt/hdf5-1.14.1-2
BOOST=$HOME/opt/boost_1_89_0
PARMETIS=$HOME/opt/graph_partitioning/parmetis-repo-install # parmetis@4.0.3
METIS=$HOME/opt/graph_partitioning/metis-repo-install # metis@5.2.1
#UMFPACK=$HOME/opt/umfpack/UMFPACK-5.3.0
OPENBLAS=$HOME/opt/OpenBLAS-0.3.30/lib/libopenblas.a
BLAS=$OPENBLAS
LAPACK=$OPENBLAS
#spack load hdf5 # hdf5@1.14.6

# Note 2025-11-08: To use MUMPS with the METIS version from the repository, the Spack configure scripts need to be updated, since METIS now depends on GKlib. To avoid this hassle, use an older version of METIS.
spack load mumps # mumps@5.8.1

# Note on MUMPS vs HDF5 vs Z library
# With MUMPS, zlib-ng is loaded, which shadows the system libraries that are used for HDF5 if that is compiled independently. Thus, I have built HDF5 via Spack as well to avoid inconsistencies.

# The settings
# -D Tpetra_INST_INT_LONG_LONG:BOOL=OFF \
# -D Tpetra_INST_INT_INT:BOOL=ON \
# are for MUMPS. By default the top setting would be on.

rm -rf CMake*

cmake \
-D CMAKE_BUILD_TYPE:STRING=$TYPE \
-D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
-D CMAKE_C_FLAGS:STRING="-D H5_HAVE_PARALLEL -D SHYLU_NODEBASKER" \
-D CMAKE_CXX_FLAGS:STRING="-D H5_HAVE_PARALLEL -D SHYLU_NODEBASKER" \
-D CMAKE_C_COMPILER=mpicc \
-D CMAKE_CXX_COMPILER=mpicxx \
-D CMAKE_Fortran_COMPILER=mpifort \
-D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
-D CMAKE_CXX_STANDARD:STRING=17 \
-D BUILD_SHARED_LIBS:BOOL=OFF \
-D Trilinos_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_EXAMPLES:BOOL=OFF \
-D Trilinos_ENABLE_Fortran:BOOL=ON \
-D Trilinos_ENABLE_OpenMP:BOOL=OFF \
-D Trilinos_ENABLE_Kokkos:BOOL=ON \
-D Trilinos_ENABLE_Teuchos:BOOL=ON \
-D Trilinos_ENABLE_Epetra:BOOL=OFF \
-D Trilinos_ENABLE_EpetraExt:BOOL=OFF \
-D Trilinos_ENABLE_PyTrilinos:BOOL=OFF \
-D Trilinos_ENABLE_Xpetra:BOOL=ON \
-D Trilinos_ENABLE_Tpetra:BOOL=ON \
-D Tpetra_INST_INT_LONG_LONG:BOOL=OFF \
-D Tpetra_INST_INT_INT:BOOL=ON \
-D Trilinos_ENABLE_AztecOO:BOOL=OFF \
-D Trilinos_ENABLE_Belos:BOOL=ON \
-D Trilinos_ENABLE_Anasazi:BOOL=ON \
-D Trilinos_ENABLE_Amesos:BOOL=OFF \
-D Trilinos_ENABLE_Amesos2:BOOL=ON \
-D Amesos2_ENABLE_MUMPS:BOOL=ON \
-D Amesos2_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_NOX:BOOL=ON \
-D NOX_ENABLE_TESTS:BOOL=ON \
-D Trilinos_ENABLE_Zoltan:BOOL=ON \
-D Zoltan_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_Zoltan2:BOOL=ON \
-D Trilinos_ENABLE_Ifpack2:BOOL=ON \
-D Trilinos_ENABLE_Thyra:BOOL=ON \
-D Thyra_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_ShyLU:BOOL=OFF \
-D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \
-D Stratimikos_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_MueLu:BOOL=OFF \
-D Trilinos_ENABLE_ML:BOOL=OFF \
-D Trilinos_ENABLE_Sacado:Bool=OFF \
-D Sacado_ENABLE_TESTS:BOOL=OFF \
-D Trilinos_ENABLE_Panzer:BOOL=OFF \
-D Trilinos_ENABLE_Intrepid2:BOOL=OFF \
-D Trilinos_ENABLE_Shards:BOOL=OFF \
-D Trilinos_ENABLE_STK:BOOL=OFF \
-D Trilinos_ENABLE_Galeri:BOOL=ON \
-D Trilinos_ENABLE_Teko:BOOL=ON \
-D TPL_FIND_SHARED_LIBS:BOOL=ON \
-D TPL_ENABLE_DLlib:BOOL=OFF \
-D TPL_ENABLE_Pthread:BOOL=OFF \
-D TPL_ENABLE_MPI:BOOL=ON \
-D TPL_ENABLE_LAPACK:BOOL=ON \
-D TPL_ENABLE_BLAS:BOOL=ON \
-D TPL_ENABLE_Matio:BOOL=OFF \
-D TPL_ENABLE_METIS:BOOL=ON \
-D METIS_LIBRARY_DIRS:PATH="$METIS/lib" \
-D METIS_INCLUDE_DIRS:PATH="$METIS/include" \
-D TPL_ENABLE_ParMETIS:BOOL=ON \
-D ParMETIS_LIBRARY_DIRS:PATH="$PARMETIS/lib" \
-D ParMETIS_INCLUDE_DIRS:PATH="$PARMETIS/include" \
-D TPL_ENABLE_UMFPACK:BOOL=OFF \
-D TPL_BLAS_LIBRARIES:STRING="$BLAS;-lgfortran;-lm" \
-D TPL_LAPACK_LIBRARIES:STRING="$LAPACK;-lgfortran;-lm" \
-D TPL_ENABLE_Boost:BOOL=ON \
-D Boost_INCLUDE_DIRS:PATH=$BOOST/include \
-D TPL_ENABLE_HDF5:BOOL=OFF \
-D Trilinos_EXTRA_LINK_FLAGS:STRING="" \
-D Trilinos_ENABLE_CONFIGURE_TIMING=ON \
-D Trilinos_ENABLE_PACKAGE_CONFIGURE_TIMING=ON \
$BASE_DIR

#-D HDF5_LIBRARY_DIRS:PATH=$HDF5/lib \
#-D HDF5_INCLUDE_DIRS:PATH=$HDF5/include \
#-D TPL_UMFPACK_LIBRARIES:STRING="$UMFPACK/lib/libumfpack.a;$UMFPACK/lib/libamd.a" \
#-D UMFPACK_INCLUDE_DIRS:PATH=$UMFPACK/include \
