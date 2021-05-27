################################################################################
#
# Set up env on a CEE RHEL6 system for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=GNU
fi

if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  unset ATDM_CONFIG_KOKKOS_ARCH
else
  echo
  echo "***"
  echo "*** ERROR: Specifying KOKKOS_ARCH is not supported on RHEL6 ATDM builds"
  echo "*** remove '$ATDM_CONFIG_KOKKOS_ARCH' from JOB_NAME=$JOB_NAME"
  echo "***"
  return
fi

echo "Using CEE RHEL6 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON

# Get ATDM_CONFIG_NUM_CORES_ON_MACHINE for this machine
source $ATDM_SCRIPT_DIR/utils/get_num_cores_on_machine.sh

if [ "$ATDM_CONFIG_NUM_CORES_ON_MACHINE" -gt "16" ] ; then
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=16
  # NOTE: We get links crashing if we try to use to many processes.  ToDo: We
  # should limit the number of processes that ninja uses to link instead of
  # reducing the overrall parallel build level like this.
else
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=$ATDM_CONFIG_NUM_CORES_ON_MACHINE
fi

export ATDM_CONFIG_BUILD_COUNT=$ATDM_CONFIG_MAX_NUM_CORES_TO_USE
# NOTE: Use as many build processes and there are cores by default.

module purge

if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
fi
# NOTE: Above, we use 1/2 as many executors as

if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
  module load sparc-dev/gcc
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/SRC
  export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib/libsuperlu_dist_4.2.a
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
  module load sparc-dev/intel-17.0.1_intelmpi-5.1.2
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  export ATDM_CONFIG_MPI_EXEC=mpirun
  export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=-np
  export ATDM_CONFIG_OPENMP_FORTRAN_FLAGS=-fopenmp
  export ATDM_CONFIG_OPENMP_FORTRAN_LIB_NAMES=gomp
  export ATDM_CONFIG_OPENMP_GOMP_LIBRARY=-lgomp
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
  export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib/libsuperlu_dist.a
elif [ "$ATDM_CONFIG_COMPILER" == "CLANG" ]; then
  module load sparc-dev/clang-5.0.1_openmpi-1.10.2
  #export OMPI_CXX=`which icpc`
  #export OMPI_CC=`which icc`
  #export OMPI_FC=`which ifort`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/SRC
  export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib/libsuperlu_dist_4.2.a
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
  return
fi

# Use updated Ninja and CMake
module load atdm-env
module load atdm-cmake/3.11.1
module load atdm-ninja_fortran/1.7.2

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_BINUTILS_LIBS="/usr/lib64/libbfd.so;/usr/lib64/libiberty.a"
# NOTE: Above, we have to explicitly set the libs to use libbdf.so instead of
# libbdf.a because the former works and the latter does not and TriBITS is set
# up to only find static libs by default!

#export ATDM_CONFIG_BLAS_LIBS="-L${CBLAS_ROOT}/mkl/lib/intel64;-L${CBLAS_ROOT}/lib/intel64;-lmkl_intel_lp64;-lmkl_intel_thread;-lmkl_core;-liomp5"
#export ATDM_CONFIG_LAPACK_LIBS="-L${CBLAS_ROOT}/mkl/lib/intel64"

# NOTE: The above does not work.  For some reason, the library 'iomp5' can't
# be found at runtime.  Instead, you have to explicitly list out the library
# files in order as shown below.  Very sad.

atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${CBLAS_ROOT}/mkl/lib/intel64 .so \
  mkl_intel_lp64 mkl_intel_thread mkl_core

atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${CBLAS_ROOT}/lib/intel64 .so \
  iomp5

export ATDM_CONFIG_LAPACK_LIBS=${ATDM_CONFIG_BLAS_LIBS}

# NOTE: HDF5_ROOT and NETCDF_ROOT should already be set in env from above
# module loads!

export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
