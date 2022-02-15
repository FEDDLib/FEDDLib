INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### FE LAPLACE
TRIBITS_ADD_EXECUTABLE(
fe_test_laplace
SOURCES fe_test_laplace.cpp
)

TRIBITS_ADD_TEST(
fe_test_laplace
NAME fe_test_laplace
ARGS
NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fe_test_laplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemLaplace.xml 
	EXEDEPS fe_test_laplace
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fe_test_laplace
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES square.mesh
	EXEDEPS fe_test_laplace
)

### NAVIER STOKES
TRIBITS_ADD_EXECUTABLE(
fe_test_navierstokes
SOURCES fe_test_navierstokes.cpp
)

TRIBITS_ADD_TEST(
fe_test_navierstokes
NAME fe_test_navierstokes
ARGS
NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fe_test_navierstokes
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemNavierStokes.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemNavierStokes.xml 
	EXEDEPS fe_test_navierstokes
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fe_test_navierstokes
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES square.mesh
	EXEDEPS fe_test_navierstokes
)
