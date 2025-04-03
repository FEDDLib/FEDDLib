INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### Domain
TRIBITS_ADD_EXECUTABLE(
domain
SOURCES domain.cpp
)

TRIBITS_ADD_TEST(
domain
NAME domain
ARGS
NUM_MPI_PROCS 4
)


### Elements
TRIBITS_ADD_EXECUTABLE(
elements
SOURCES elements.cpp
)

TRIBITS_ADD_TEST(
elements
NAME elements
ARGS
NUM_MPI_PROCS 4
)

### FE
TRIBITS_ADD_EXECUTABLE(
fe
SOURCES fe.cpp
)

TRIBITS_ADD_TEST(
fe
NAME fe
ARGS
NUM_MPI_PROCS 4
)

### QuadRules
TRIBITS_ADD_EXECUTABLE(
quadRules
SOURCES quadRules.cpp
)

TRIBITS_ADD_TEST(
quadRules
NAME quadRules
ARGS
NUM_MPI_PROCS 1
)

### TestQuadRules
TRIBITS_ADD_EXECUTABLE(
testQuadRules
SOURCES testQuadRules.cpp
)

TRIBITS_ADD_TEST(
testQuadRules
NAME testQuadRules
ARGS
NUM_MPI_PROCS 1
)


### SurfaceIntegral
TRIBITS_ADD_EXECUTABLE(
surfaceIntegral
SOURCES surfaceIntegral.cpp
)

TRIBITS_ADD_TEST(
surfaceIntegral
NAME surfaceIntegral
ARGS
NUM_MPI_PROCS 2
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_mesh_2dSquare
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS surfaceIntegral
)
