INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### neumannBC
TRIBITS_ADD_EXECUTABLE(
neumannBC
SOURCES neumannBC.cpp
)

TRIBITS_ADD_TEST(
neumannBC
NAME neumannBC
ARGS
NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_bc_2dSquare
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS neumannBC
)

### Write and Read File HDF5
TRIBITS_ADD_EXECUTABLE(
writeRead
SOURCES write_read.cpp
)

TRIBITS_ADD_TEST(
writeRead
NAME writeRead
ARGS
NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_writeRead
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS writeRead
)