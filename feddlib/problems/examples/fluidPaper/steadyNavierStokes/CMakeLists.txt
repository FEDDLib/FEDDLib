TRIBITS_ADD_EXECUTABLE(
	steadyNavierStokesPaper
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyNavierStokesPaper
	NAME steadyNavierStokesPaper
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyNavierStokesPaper
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyNavierStokesPaper
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNavierStokesPaper
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_fluid_h004.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_fluid_h004.mesh
	EXEDEPS steadyNavierStokesPaper
)
