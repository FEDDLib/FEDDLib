### Domain
tribits_add_executable(
domain
SOURCES domain.cpp
)

tribits_add_test(
domain
NAME domain
ARGS
NUM_MPI_PROCS 4
)


### Elements
tribits_add_executable(
elements
SOURCES elements.cpp
)

tribits_add_test(
elements
NAME elements
ARGS
NUM_MPI_PROCS 4
)

### FE
tribits_add_executable(
fe
SOURCES fe.cpp
)

tribits_add_test(
fe
NAME fe
ARGS
NUM_MPI_PROCS 4
)

### QuadRules
tribits_add_executable(
quadRules
SOURCES quadRules.cpp
)

tribits_add_test(
quadRules
NAME quadRules
ARGS
NUM_MPI_PROCS 1
)

### TestQuadRules
tribits_add_executable(
testQuadRules
SOURCES testQuadRules.cpp
)

tribits_add_test(
testQuadRules
NAME testQuadRules
ARGS
NUM_MPI_PROCS 1
)


### SurfaceIntegral
tribits_add_executable(
surfaceIntegral
SOURCES surfaceIntegral.cpp
)

tribits_add_test(
surfaceIntegral
NAME surfaceIntegral
ARGS
NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_mesh_2dSquare
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS surfaceIntegral
)
