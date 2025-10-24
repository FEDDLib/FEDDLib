### neumannBC
tribits_add_executable(
neumannBC
SOURCES neumannBC.cpp
)

tribits_add_test(
neumannBC
NAME neumannBC
ARGS
NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_bc_2dSquare
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS neumannBC
)

### Write and Read File HDF5
tribits_add_executable(
writeRead
SOURCES write_read.cpp
)

tribits_add_test(
writeRead
NAME writeRead
ARGS
NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_writeRead
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS writeRead
)
