tribits_add_executable(
	steadyNavierStokesAssFE
	SOURCES main.cpp
)

tribits_add_test(
	steadyNavierStokesAssFE
	NAME steadyNavierStokesAssFE
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_steadyNavierStokesAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyNavierStokesAssFE
)

tribits_copy_files_to_binary_dir(mesh_steadyNavierStokesAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_fluid_h008.mesh DFG3DCylinder_1k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_fluid_h008.mesh DFG3DCylinder_1k.mesh
	EXEDEPS steadyNavierStokesAssFE
)




