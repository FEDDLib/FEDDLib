#!/bin/bash                                                                                                                  

#! Change paths in base and install dir to fit your setup !#
TYPE=DEBUG
BASE_DIR=/home/lsassmannshausen/Trilinos/source/Trilinos/ # Trilinos source directory
INSTALL_DIR=/home/lsassmannshausen/Trilinos/install/ # Trilinos install directory

rm -rf CMake*

cmake \
    -D CMAKE_BUILD_TYPE:STRING=$TYPE \
    -D CMAKE_INSTALL_PREFIX:STRING=$INSTALL_DIR \
    -D CMAKE_C_FLAGS:STRING="-O2 -Wno-format-security -DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DNDEBUG -fPIC -DFROSCH_TIMER_DETAILS=2" \
    -D CMAKE_CXX_FLAGS:STRING="-O2 -Wno-format-security -DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DNDEBUG -fPIC -DFROSCH_TIMER_DETAILS=2" \
    -D CMAKE_CXX_STANDARD:STRING=17 \
    -D CMAKE_C_COMPILER=mpicc \
    -D CMAKE_CXX_COMPILER=mpicxx \
    -D CMAKE_Fortran_COMPILER=mpifort \
    -D MPI_EXEC_MAX_NUMPROCS:STRING=8 \
    -D Trilinos_ENABLE_Fortran:BOOL=ON \
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
    -D Trilinos_ENABLE_Rythmos:BOOL=ON \
    -D Trilinos_ENABLE_Thyra:BOOL=ON \
    -D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
    -D Trilinos_ENABLE_MueLu:BOOL=OFF \
    -D TPL_FIND_SHARED_LIBS:BOOL=ON \
    -D TPL_ENABLE_DLlib:BOOL=OFF \
    -D TPL_ENABLE_Pthread:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D TPL_ENABLE_LAPACK:BOOL=ON \
    -D Tpetra_INST_INT_INT:BOOL=ON \
    -D Tpetra_INST_INT_LONG_LONG:BOOL=OFF \
    -D TPL_ENABLE_MPI:BOOL=ON \
    -D TPL_ENABLE_METIS:BOOL=ON \
    -D TPL_ENABLE_BLAS:BOOL=ON \
    -D TPL_ENABLE_ParMETIS:BOOL=ON \
    -D TPL_ParMETIS_INCLUDE:STRING="${METISROOT}/include/metis.h;${PARMETISROOT}/include/parmetis.h" \
    -D Amesos_ENABLE_MUMPS:BOOL=OFF \
    -D Amesos2_ENABLE_MUMPS:BOOL=ON \
    -D TPL_ENABLE_Boost:BOOL=ON \
    -D Boost_INCLUDE_DIRS:PATH=$BOOST \
    -D EpetraExt_USING_HDF5:BOOL=ON \
    -D TPL_ENABLE_HDF5:BOOL=ON \
    -D HDF5_LIBRARY_DIRS:PATH=$HDF5/.libs \
    -D HDF5_INCLUDE_DIRS:PATH=$HDF5/ \
    -D Trilinos_ENABLE_TESTS:BOOL=ON \
    -D NOX_ENABLE_ThyraCore:BOOL=ON \
    -D NOX_ENABLE_TESTS:BOOL=ON \
    -D Trilinos_ENABLE_CONFIGURE_TIMING=ON \
    -D Trilinos_ENABLE_PACKAGE_CONFIGURE_TIMING=ON \
$BASE_DIR

