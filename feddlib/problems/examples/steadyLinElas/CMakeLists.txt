tribits_add_executable(
	steadyLinElas
	SOURCES main.cpp
)

tribits_add_test(
	steadyLinElas
	NAME steadyLinElas
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_steadyLinElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyLinElas
)

tribits_copy_files_to_binary_dir(mesh_steadyLinElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_solid_h002.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_solid_h002.mesh
	EXEDEPS steadyLinElas
)
