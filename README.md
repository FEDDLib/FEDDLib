# FEDDLib
This is a C++ based library. FEDDLib is short for 'Finite Element and Domain Decomposition Library'. The FEDDLib is based in large parts on the open source Library Trilinos (trilinos.github.io)


REQUIEREMENTS:
  -> In order to install the FEDDLib you are requiered to have a running Trilinos version to link with. 
  
Include and Library:
-D TPL_ENABLE_Trilinos:BOOL=ON \
-D Trilinos_INCLUDE_DIRS:PATH=$TRILINOS_DIR/include \
-D Trilinos_LIBRARY_DIRS:PATH=$TRILINOS_DIR/lib \ # in newer Trilinos Versions (after Jan 2024 you will need  $TRILINOS_DIR/lib_64)

Trilinos:
Basic requierements for Trilinos build are the following packages:
    -D Trilinos_ENABLE_Teuchos:BOOL=ON \  # Smart Pointer
    -D Trilinos_ENABLE_Epetra:BOOL=ON \   # Linear Algebra
    -D Trilinos_ENABLE_Tpetra:BOOL=ON \   # Linear Algebra
    -D Trilinos_ENABLE_EpetraExt:BOOL=ON \ # Exporter
    -D Trilinos_ENABLE_Belos:BOOL=ON \   # Linear Iterative Solver
    -D Trilinos_ENABLE_Amesos2:BOOL=ON \ # Direct Solver
    -D Trilinos_ENABLE_NOX:BOOL=ON \ # Nonlinear Solver
    -D Trilinos_ENABLE_Zoltan2:BOOL=ON \ # Graph Partitioning
    -D Trilinos_ENABLE_Thyra:BOOL=ON \ # Preconditioner Interfacing
    -D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \ # Preconditioner
    -D Trilinos_ENABLE_Stratimikos:BOOL=ON \ # Preconditioner Interfacing
Trilinos contains its own direct solver KLU2 in Amesos2. As a result you are no reuqiered to install any additional direct solvers (i.e. pardisomkl,umfpack, mumps).
Additionally you will ne to install the packages: METIS, PARMETIS, BOOST, HDF5, BLAS

AceGEN Interface:
If you want to use the developed AceGEN interface to full extend you need to specify the Library and Include files.
-D TPL_ENABLE_AceGENInterface:BOOL=ON \
-D TPL_AceGENInterface_LIBRARIES:STRING=" "\
-D TPL_AceGENInterface_INCLUDE_DIRS:STRING=" "\
