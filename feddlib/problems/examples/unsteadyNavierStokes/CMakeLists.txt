tribits_add_executable(
	unsteadyNavierStokes
	SOURCES main.cpp
)

tribits_add_test(
	unsteadyNavierStokes
	NAME unsteadyNavierStokes
	ARGS
	#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_unsteadyNavierStokes
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersPrecBlock.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersPrecBlock.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS unsteadyNavierStokes
)

tribits_copy_files_to_binary_dir(mesh_unsteadyNavierStokes
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_fluid_h004.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES dfg_fsi_fluid_h004.mesh
	EXEDEPS unsteadyNavierStokes
)
