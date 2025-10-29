tribits_add_executable(
	unsteadyLinElas
	SOURCES main.cpp
)

tribits_add_test(
	unsteadyLinElas
	NAME unsteadyLinElas
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_unsteadyLinElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyLinElas
)

tribits_copy_files_to_binary_dir(mesh_unsteadyLinElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_solid_h002.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_solid_h002.mesh
	EXEDEPS unsteadyLinElas
)
