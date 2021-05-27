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
