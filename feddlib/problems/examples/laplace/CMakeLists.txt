tribits_add_executable(
	laplace
	SOURCES main.cpp
)

tribits_add_test(
	laplace
    NAME laplace
    ARGS
    NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_laplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml 
	EXEDEPS laplace
)

tribits_copy_files_to_binary_dir(mesh_laplace
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS laplace
)
