#!/bin/bash
# Third-party libraries for Trilinos 16.1.0 and FEDDLib on macOS / Apple
# silicon (M4, macOS 26), as used for the local builds of 2026-09 (versions
# from `spack find`). The spack checkout is ~/dev/spack.
#
# Toolchain from Homebrew and Xcode:
#   xcode-select --install                  # Apple clang 21
#   brew install open-mpi gcc cmake ninja   # Open MPI 5.0.9, gcc/gfortran 15.2, cmake 4.2, ninja 1.13
#
# Pitfalls on Apple silicon:
#   - Register Homebrew's Open MPI in spack and build every MPI dependency
#     against it (^openmpi). HDF5/ParMETIS built against spack's MPICH, mixed
#     with Homebrew Open MPI, fail at run time with misleading messages
#     ("MPICH ... before initializing", "HDF5: infinite loop closing library").
#   - BLAS/LAPACK: use Apple Accelerate through the veclibfort shim. spack's
#     netlib-lapack/OpenBLAS for target=m4 compile their Fortran with
#     -march=armv9.2-a+sme2, which GCC treats as including SVE; the M4 cannot
#     execute SVE instructions (EXC_BAD_INSTRUCTION, e.g. in dnrm2_ inside
#     FROSch's coarse QR).
#   - There is no MKL on Apple silicon, so no PARDISO: use Amesos2's "klu2"
#     wherever the parameter files say "pardisomkl".

export PATH=$HOME/dev/spack/bin:$PATH

spack external find openmpi                 # Homebrew Open MPI 5.0.9 (prefix /opt/homebrew)

spack install veclibfort                    # 0.4.3
spack install metis                         # 5.1.0
spack install parmetis ^openmpi             # 4.0.3
spack install hdf5+mpi+hl ^openmpi          # 1.14.6 (pulls in zlib-ng 2.2.4)
spack install boost                         # 1.88.0 (Trilinos only needs the headers)
