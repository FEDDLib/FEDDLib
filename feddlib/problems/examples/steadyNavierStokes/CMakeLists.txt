tribits_add_executable(
	steadyNavierStokes
	SOURCES main.cpp
)

tribits_add_test(
	steadyNavierStokes
	NAME steadyNavierStokes
	ARGS
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_steadyNavierStokes
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	EXEDEPS steadyNavierStokes
)

tribits_copy_files_to_binary_dir(mesh_steadyNavierStokes
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_fluid_h004.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_fluid_h004.mesh
	EXEDEPS steadyNavierStokes
)