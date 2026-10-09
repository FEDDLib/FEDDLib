# Building FEDDLib on macOS / Apple silicon

The scripts used for the local FEDDLib stack on an Apple M4 (macOS 26, 2026-09),
in build order:

1. `install-dependencies-on-macOS-AppleSilicon-via-spack.sh`: Homebrew toolchain
   (Apple clang, gfortran, Open MPI) and the spack-built third-party libraries
   (veclibfort, METIS, ParMETIS, HDF5, Boost).
2. `configure-Trilinos-on-macOS-AppleSilicon.sh`: Trilinos 16.1.0 (release tag
   `trilinos-release-16-1-0`) with global ordinal `long long`, installed to
   `~/opt/trilinos-16.1.0`. FEDDLib needs Trilinos 16.1.0: later versions drop
   deprecated packages it uses.
3. Interface2 (optional, for the AceGen elements): its
   `config_scripts/do-config-interface2-macos-accelerate.sh`, installed to
   `~/opt/interface2`.
4. `configure-FEDDLib-on-macOS-AppleSilicon.sh`: FEDDLib itself (`ACEGEN=OFF`
   without Interface2, `ASAN=ON` for an AddressSanitizer build).

Every configure script runs from an empty build directory and states its
options and paths at the top.

## Differences from a cluster build

- BLAS/LAPACK come from Apple Accelerate (via veclibfort) instead of MKL, so
  there is no PARDISO: parameter files have to use Amesos2's `klu2` instead of
  `pardisomkl`. With KLU2 as FROSch's overlapping solver the artery examples
  converge; replacing it by an incomplete factorization (Ifpack2 RILUK) makes the
  linear solves stall.
- All MPI libraries must come from the same MPI (here Homebrew Open MPI);
  spack's default MPICH must not be mixed in.
- Memory: the artery examples need about 7 GB on 6 processes; with 16 GB
  some process counts are killed for lack of memory.
