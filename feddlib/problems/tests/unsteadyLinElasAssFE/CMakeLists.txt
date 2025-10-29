tribits_add_executable(
	unsteadyLinElasAssFE
	SOURCES main.cpp
)

tribits_add_test(
	unsteadyLinElasAssFE
	NAME unsteadyLinElasAssFE
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_unsteadyLinElasAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyLinElasAssFE
)

tribits_copy_files_to_binary_dir(mesh_unsteadyLinElasAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS unsteadyLinElasAssFE
)
