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
	SOURCE_FILES cube4.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cube4.mesh
	EXEDEPS fe_test_navierstokes
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fe_test_navierstokes2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES square.mesh
	EXEDEPS fe_test_navierstokes
)

### FE LINELAS
TRIBITS_ADD_EXECUTABLE(
fe_test_linElas
SOURCES fe_test_linElas.cpp
)

TRIBITS_ADD_TEST(
fe_test_linElas
NAME fe_test_linElas
ARGS
NUM_MPI_PROCS 1
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fe_test_linElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemLinElas.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemLinElas.xml 
	EXEDEPS fe_test_linElas
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fe_test_linElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES tetrahedron.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES tetrahedron.mesh
	EXEDEPS fe_test_linElas
)

### FE NONLINELAS
TRIBITS_ADD_EXECUTABLE(
fe_test_nonlinElas
SOURCES fe_test_nonLinElas.cpp
)

TRIBITS_ADD_TEST(
fe_test_nonlinElas
NAME fe_test_nonlinElas
ARGS
NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fe_test_nonlinElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemNonLinElas.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemNonLinElas.xml 
	EXEDEPS fe_test_nonlinElas
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_fe_test_nonlinElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cube.mesh
	EXEDEPS fe_test_nonlinElas
)

if(FEDD_HAVE_ACEGENINTERFACE)
	### Test for new interfacce
	### FE LAPLACE
	TRIBITS_ADD_EXECUTABLE(
	interface_test
	SOURCES test.cpp
	)

	TRIBITS_ADD_TEST(
	interface_test
	NAME interface_test
	ARGS
	NUM_MPI_PROCS 1
	)

endif()
