tribits_add_executable(
	steadyLinElasAssFE
	SOURCES main.cpp
)

tribits_add_test(
	steadyLinElasAssFE
	NAME steadyLinElasAssFE
	ARGS
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_steadyLinElasAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyLinElasAssFE
)

tribits_copy_files_to_binary_dir(mesh_steadyLinElasAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube_h_1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_h_1.mesh
	EXEDEPS steadyLinElasAssFE
)
