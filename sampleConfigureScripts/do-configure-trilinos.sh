#!/bin/bash

TYPE=RELEASE

BASE_DIR=/home/lea/Software/TRILINOS/source/Trilinos/ #/src/mytrilinos
INSTALL_DIR=/home/lea/Software/TRILINOS/install/ #/opt/mytrilinos

#HDF5=/opt/gnu/hdf5-1.8.19/
HDF5=/home/lea/Software/TPL/hdf5_1.10.4/src/
#BOOST
BOOST=/home/lea/Software/TPL/boost_1_79_0/
#PARMETIS
PARMETIS=/home/lea/Software/TPL/parmetis-4.0.3/
#UMFPACK
UMFPACK=/home/lea/Software/TPL/UMFPACK/UMFPACK/
#BLAS
BLAS= /home/lea/Software/TPL/BLAS-3.10.0/
#MUMPS
MUMPS=/home/lea/Software/TPL/MUMPS/MUMPS_4.7.3/

rm -rf CMake*

cmake \
    -D CMAKE_BUILD_TYPE:STRING=$TYPE \
    -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
    -D CMAKE_C_FLAGS:STRING="-Wno-format-security -DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DNDEBUG" \
    -D CMAKE_CXX_FLAGS:STRING="-Wno-format-security -DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DNDEBUG" \
    -D CMAKE_CXX_STANDARD:STRING=14 \
    -D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
    -D Trilinos_ENABLE_Fortran:BOOL=ON \
    -D PYTHON_EXECUTABLE:STRING="/usr/bin/python3" \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D BUILD_SHARED_LIBS:BOOL=OFF \
    -D Trilinos_ENABLE_Fortran:BOOL=ON \
    -D Trilinos_ENABLE_OpenMP:BOOL=OFF \
    -D Trilinos_ENABLE_Teuchos:BOOL=ON \
    -D Trilinos_ENABLE_Epetra:BOOL=ON \
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
    -D Trilinos_ENABLE_AztecOO:BOOL=ON \
    -D Trilinos_ENABLE_Belos:BOOL=ON \
    -D Trilinos_ENABLE_Anasazi:BOOL=ON \
    -D Trilinos_ENABLE_Amesos:BOOL=ON \
    -D Trilinos_ENABLE_Amesos2:BOOL=ON \
    -D Trilinos_ENABLE_NOX:BOOL=ON \
    -D Trilinos_ENABLE_Zoltan:BOOL=ON \
    -D Trilinos_ENABLE_Zoltan2:BOOL=ON \
    -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
    -D Trilinos_ENABLE_ML:BOOL=ON \
    -D Trilinos_ENABLE_Rythmos:BOOL=OFF \
    -D Trilinos_ENABLE_Thyra:BOOL=ON \
    -D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
    -D Trilinos_ENABLE_MueLu:BOOL=OFF \
    -D Tpetra_INST_INT_INT:BOOL=OFF \
    -D Tpetra_INST_INT_LONG_LONG:BOOL=ON \
    -D TPL_FIND_SHARED_LIBS:BOOL=ON \
    -D TPL_ENABLE_DLlib:BOOL=OFF \
    -D TPL_ENABLE_Pthread:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D TPL_ENABLE_LAPACK:BOOL=ON \
    -D LAPACK_INCLUDE_DIRS=$BLAS \
    -D LAPACK_LIBRARIES:STRING="/home/lea/Software/TPL/OpenBLAS-0.3.17/libopenblas.so" \
    -D Trilinos_ENABLE_OpenMP:BOOL=OFF \
    -D TPL_ENABLE_Matio:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_METIS:BOOL=ON \
    -D TPL_ENABLE_BLAS:BOOL=ON \
    -D BLAS_INCLUDE_DIRS=$BLAS \
    -D BLAS_LIBRARIES:STRING="/home/lea/Software/TPL/OpenBLAS-0.3.17/libopenblas.so" \
    -D METIS_LIBRARY_DIRS:PATH=$PARMETIS/build/Linux-x86_64/libmetis \
    -D METIS_INCLUDE_DIRS:PATH=$PARMETIS/metis/include \
    -D TPL_ENABLE_ParMETIS:BOOL=ON \
    -D ParMETIS_LIBRARY_DIRS:PATH=$PARMETIS/build/Linux-x86_64/libparmetis  \
    -D ParMETIS_INCLUDE_DIRS:PATH=$PARMETIS/include \
    -D TPL_ENABLE_UMFPACK:BOOL=ON \
    -D TPL_UMFPACK_LIBRARIES:STRING="/home/lea/Software/TPL/UMFPACK/UMFPACK/Lib/libumfpack.a;/home/lea/Software/TPL/UMFPACK/AMD/Lib/libamd.a" \
    -D UMFPACK_INCLUDE_DIRS:PATH=$UMFPACK/Include \
    -D Amesos_ENABLE_MUMPS:BOOL=OFF \
    -D Amesos2_ENABLE_MUMPS:BOOL=ON \
    -D TPL_ENABLE_MUMPS:BOOL=ON \
    -D MUMPS_LIBRARY_DIRS:PATH=$MUMPS/lib \
    -D MUMPS_INCLUDE_DIRS:PATH=$MUMPS/include \
    -D MUMPS_LIBRARY_NAMES:STRING="dmumps;pord" \
    -D Amesos_ENABLE_SCALAPACK:BOOL=ON \
    -D SCALAPACK_INCLUDE_DIRS:FILEPATH="/home/lea/Software/TPL/scalapack-2.2.0/SRC" \
	-D SCALAPACK_LIBRARY_DIRS:FILEPATH="/home/lea/Software/TPL/scalapack-2.2.0" \
	-D SCALAPACK_LIBRARY_NAMES:STRING="scalapack" \
	-D Amesos_ENABLE_BLACS:BOOL=ON \
	-D BLACS_INCLUDE_DIRS:FILEPATH="/home/lea/Software/TPL/BLACS/SRC/MPI" \
	-D BLACS_LIBRARY_DIRS:FILEPATH="/home/lea/Software/TPL/BLACS/SRC/MPI" \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D Boost_INCLUDE_DIRS:PATH=$BOOST \
    -D EpetraExt_USING_HDF5:BOOL=ON \
    -D TPL_ENABLE_HDF5:BOOL=ON \
    -D HDF5_LIBRARY_DIRS:PATH=$HDF5/.libs \
    -D HDF5_INCLUDE_DIRS:PATH=$HDF5/ \
    -D Trilinos_ENABLE_CONFIGURE_TIMING=ON \
	-D Trilinos_ENABLE_PACKAGE_CONFIGURE_TIMING=ON \
$BASE_DIR
