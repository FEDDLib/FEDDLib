INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### BlockViews
#TRIBITS_ADD_EXECUTABLE(
#	blockView
#	SOURCES blockView.cpp
#)

#TRIBITS_ADD_TEST(
#	blockView
#    NAME blockView
#    ARGS
#    NUM_MPI_PROCS 2
#)

### Map
TRIBITS_ADD_EXECUTABLE(
	map
	SOURCES map.cpp
)

TRIBITS_ADD_TEST(
	map
    NAME map
    ARGS
    NUM_MPI_PROCS 4
)

### Map Tpetra
TRIBITS_ADD_EXECUTABLE(
	map_tpetra
	SOURCES map_tpetra.cpp
)

TRIBITS_ADD_TEST(
	map_tpetra
    NAME map_tpetra
    ARGS
    NUM_MPI_PROCS 4
)

### Matrix
TRIBITS_ADD_EXECUTABLE(
	matrix
	SOURCES matrix.cpp
)

TRIBITS_ADD_TEST(
	matrix
    NAME matrix
    ARGS
    NUM_MPI_PROCS 4
)


### Matrix
TRIBITS_ADD_EXECUTABLE(
    blockMatrixPtrAccess
    SOURCES blockMatrixPtrAccess.cpp
)

TRIBITS_ADD_TEST(
    blockMatrixPtrAccess
    NAME blockMatrixPtrAccess
    ARGS
    NUM_MPI_PROCS 2
)


### Matrix assembly
TRIBITS_ADD_EXECUTABLE(
	matrix_assembly
	SOURCES matrix_assembly.cpp
)

TRIBITS_ADD_TEST(
	matrix_assembly
    NAME matrix_assembly
    ARGS
    NUM_MPI_PROCS 2
)


### MultiVector
TRIBITS_ADD_EXECUTABLE(
multivector
SOURCES multivector.cpp
)

TRIBITS_ADD_TEST(
multivector
NAME multivector
ARGS
NUM_MPI_PROCS 4
)

#TRIBITS_COPY_FILES_TO_BINARY_DIR(data_matrix
#	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#	SOURCE_FILES
#	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
#	DEST_FILES
#	EXEDEPS matrix
#)


### Consistent Partitioning
#TRIBITS_ADD_EXECUTABLE(
#consistentPartitioning
#SOURCES consistentPartitioning.cpp
#)

#TRIBITS_ADD_TEST(
#consistentPartitioning
#NAME consistentPartitioning
#ARGS
#NUM_MPI_PROCS 2
#)
