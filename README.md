# FEDDLib

This is a C++ based library for finite element problems that are foremost solved with domain decomposition methods. FEDDLib is short for ***F**inite **E**lement and **D**omain **D**ecomposition **Lib**rary*. The FEDDLib is based in large parts on the open source library [Trilinos](https://trilinos.github.io/) and uses [CMake](https://cmake.org/) with [TriBITS](https://tribits.org/) ([repository](https://github.com/TriBITSPub/TriBITS) on GitHub) as its build system. A Doxygen documentation is hosted at [feddlib.github.io/FEDDLib](https://feddlib.github.io/FEDDLib/).

## Installation Prerequisites

In order to install the FEDDLib, you need to have an installed Trilinos version that can be linked against. The FEDDLib repository contains a link to a specific commit of the Trilinos repository as a git submodule. This commit serves as a reference, a version of Trilinos that the current state of the FEDDLib (its branch *develop*) has been tested against.

### Trilinos

Help for an installation of Trilinos can be found in the `doc` folder.

Trilinos contains its own direct solver KLU2 in Amesos2. As a result, you are not required to install any additional direct solvers (e.g., [PardisoMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html), [UMFPACK](https://github.com/PetterS/SuiteSparse/tree/master/UMFPACK), [MUMPS](https://mumps-solver.org/index.php)). The following Trilinos dependencies need to be installed: [METIS](https://github.com/KarypisLab/METIS), [ParMETIS](https://github.com/KarypisLab/ParMETIS), [BOOST](https://www.boost.org/), [BLAS](https://github.com/OpenMathLib/OpenBLAS).

The most important CMake configuration flags for Trilinos are

```
    -D Trilinos_ENABLE_Teuchos:BOOL=ON \        # Smart pointers and other tools
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \         # Linear algebra backend (matrices, vectors etc.)
    -D Trilinos_ENABLE_Belos:BOOL=ON \          # Iterative linear solver
    -D Trilinos_ENABLE_Amesos2:BOOL=ON \        # Direct solver
    -D Trilinos_ENABLE_NOX:BOOL=ON \            # Nonlinear solver
    -D Trilinos_ENABLE_Zoltan2:BOOL=ON \        # Graph partitioning
    -D Trilinos_ENABLE_Thyra:BOOL=ON \          # Abstract linear algebra (linear operator etc.)
    -D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \ # Preconditioner (Overlapping Schwarz)
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \    # Abstract interface to preconditioners etc.
```

### FEDDLib

To install the FEDDLib, [HDF5](https://www.hdfgroup.org/solutions/hdf5/) needs to be installed as a direct dependency. Note that HDF5 itself depends on the compression library z. MUMPS also depends on [zlib](https://zlib.net/) (or [zlib-ng](https://github.com/zlib-ng/zlib-ng)). There is a risk of and a risk in mixing different z libraries, depending on how the dependencies are installed.

HDF5 is a required dependency of the FEDDLib, which means that it does not have to be specifically enabled via CMake flags.

```
    -D TPL_ENABLE_Trilinos:BOOL=ON \
    -D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
    -D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib \ 
```

Check your Trilinos installation folder if the library was placed in the `lib` or the `lib64` folder and make the corresponding adjustment above (`$TRILINOS_DIR/lib` vs. `$TRILINOS_DIR/lib64`).

### AceGEN Interface

The FEDDLib has an optional dependency on an interface to [AceGEN](https://www.wolfram.com/products/applications/acegen/). To enable its use and installation, you can use the following CMake flags and adapt them to your setup.

```
    -D TPL_ENABLE_AceGENInterface:BOOL=ON \
    -D TPL_AceGENInterface_LIBRARIES:STRING=" "\
    -D TPL_AceGENInterface_INCLUDE_DIRS:STRING=" "\
```


## Simulation Setup from a User Perspective

1. **Mesh**
   * As a starting point, a discretization of the physical domain is needed. FEDDLib provides functionalities to read meshes that are written in the INRIA mesh format (`.mesh`) that is used by Medit. For mesh generation, one can use, for example, [Gmsh](https://gmsh.info/) and export to an `.mesh` file directly. A post-processing step is required, however, to set the flags of vertices. Alternatively, Gmsh's internal msh format can be converted to `.mesh` files (see `meshes/Msh2Mesh/` for a MATLAB and Python script). Various meshes are included in the directory `meshes`.

2. **Simulation**
   * The FEDDLib contains assembly routines for specific model problems such as the Poisson equation, linear and nonlinear elasticity, a Stokes problem, the (transient) Navier-Stokes equations, or fluid-structure interaction problems; see `feddlib/problems/tests` and `feddlib/problems/example`. The user can select one of these problem types and perform simulations for their computational domain adjusting the parameter files accordingly (see `doc/parameterFiles` for a documentation on parameters).
   * To change the boundary conditions of a model problem, e.g., `feddlib/problems/tests/stokes/`, one has to modify the corresponding main file and set the desired boundary condition, e.g., the no-slip condition, corresponding to a specific boundary flag in the mesh. The code must then be rebuilt.
   * The simulation can then be run from the command line via `mpirun -np N ./problem_name.exe`. In the case of parallel computations, the partitioning into N subdomains corresponding to N MPI ranks is done automatically.

3. **Visualization**
   * For post-processing, the FEDDLib offers an exporter for HDF5 (`.h5`) files with [XDMF](https://xdmf.org/) (`.xmf`) files (the XDMF files serve as descriptors for the compressed HDF5 files, which allow complex data structures). To visualize soluttions, [ParaView](https://www.paraview.org/) can be used by opening `.xmf` files.
 
