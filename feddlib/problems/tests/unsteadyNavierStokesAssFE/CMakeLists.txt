TRIBITS_ADD_EXECUTABLE(
	unsteadyNavierStokesAssFE
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyNavierStokesAssFE
	NAME unsteadyNavierStokesAssFE
	ARGS
	#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNavierStokesAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS unsteadyNavierStokesAssFE
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNavierStokesAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_fluid_h008.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES dfg_fsi_fluid_h008.mesh
	EXEDEPS unsteadyNavierStokesAssFE
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNavierStokesAssFE3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES DFG3DCylinder_1k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES DFG3DCylinder_1k.mesh
	EXEDEPS unsteadyNavierStokesAssFE
)

