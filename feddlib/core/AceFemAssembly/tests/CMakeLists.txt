include(AddSubdirectories)

### FE LAPLACE
tribits_add_executable(
fe_test_laplace
SOURCES fe_test_laplace.cpp
)

tribits_add_test(
fe_test_laplace
NAME fe_test_laplace
ARGS
NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_fe_test_laplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemLaplace.xml 
	EXEDEPS fe_test_laplace
)

tribits_copy_files_to_binary_dir(mesh_fe_test_laplace
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES square.mesh
	EXEDEPS fe_test_laplace
)

### NAVIER STOKES
tribits_add_executable(
fe_test_navierstokes
SOURCES fe_test_navierstokes.cpp
)

tribits_add_test(
fe_test_navierstokes
NAME fe_test_navierstokes
ARGS
NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_fe_test_navierstokes
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemNavierStokes.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemNavierStokes.xml 
	EXEDEPS fe_test_navierstokes
)

tribits_copy_files_to_binary_dir(mesh_fe_test_navierstokes
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube4.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cube4.mesh
	EXEDEPS fe_test_navierstokes
)

tribits_copy_files_to_binary_dir(mesh_fe_test_navierstokes2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES square.mesh
	EXEDEPS fe_test_navierstokes
)

### FE LINELAS
tribits_add_executable(
fe_test_linElas
SOURCES fe_test_linElas.cpp
)

tribits_add_test(
fe_test_linElas
NAME fe_test_linElas
ARGS
NUM_MPI_PROCS 1
)

tribits_copy_files_to_binary_dir(data_fe_test_linElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemLinElas.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemLinElas.xml 
	EXEDEPS fe_test_linElas
)

tribits_copy_files_to_binary_dir(mesh_fe_test_linElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES tetrahedron.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES tetrahedron.mesh
	EXEDEPS fe_test_linElas
)

### FE NONLINELAS
tribits_add_executable(
fe_test_nonlinElas
SOURCES fe_test_nonLinElas.cpp
)

tribits_add_test(
fe_test_nonlinElas
NAME fe_test_nonlinElas
ARGS
NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_fe_test_nonlinElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblemNonLinElas.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblemNonLinElas.xml 
	EXEDEPS fe_test_nonlinElas
)

tribits_copy_files_to_binary_dir(mesh_fe_test_nonlinElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cube.mesh
	EXEDEPS fe_test_nonlinElas
)

if(FEDD_HAVE_ACEGENINTERFACE)
	### Test for new interfacce
	### FE LAPLACE
	tribits_add_executable(
	interface_test
	SOURCES test.cpp
	)

	tribits_add_test(
	interface_test
	NAME interface_test
	ARGS
	NUM_MPI_PROCS 1
	)

endif()
