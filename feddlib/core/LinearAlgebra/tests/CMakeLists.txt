### Map
# tribits_add_executable(
# 	map
# 	SOURCES map.cpp
# )

# tribits_add_test(
# 	map
#     NAME map
#     ARGS
#     NUM_MPI_PROCS 4
# )

### Map Tpetra
tribits_add_executable(
	map_tpetra
	SOURCES map_tpetra.cpp
)

tribits_add_test(
	map_tpetra
    NAME map_tpetra
    ARGS
    NUM_MPI_PROCS 4
)

### Matrix
tribits_add_executable(
	matrix
	SOURCES matrix.cpp
)

tribits_add_test(
	matrix
    NAME matrix
    ARGS
    NUM_MPI_PROCS 4
)


### Matrix
tribits_add_executable(
    blockMatrixPtrAccess
    SOURCES blockMatrixPtrAccess.cpp
)

tribits_add_test(
    blockMatrixPtrAccess
    NAME blockMatrixPtrAccess
    ARGS
    NUM_MPI_PROCS 2
)


### Matrix assembly
tribits_add_executable(
	matrix_assembly
	SOURCES matrix_assembly.cpp
)

tribits_add_test(
	matrix_assembly
    NAME matrix_assembly
    ARGS
    NUM_MPI_PROCS 2
)


### MultiVector
tribits_add_executable(
multivector
SOURCES multivector.cpp
)

tribits_add_test(
multivector
NAME multivector
ARGS
NUM_MPI_PROCS 4
)
