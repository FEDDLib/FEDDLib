TRIBITS_ADD_EXECUTABLE(
	unsteadyNavierStokes_Restart
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyNavierStokes_Restart
	NAME unsteadyNavierStokes_Restart_4
	ARGS
	#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
	unsteadyNavierStokes_Restart
	NAME unsteadyNavierStokes_Restart_6
	ARGS
	#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 6
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNavierStokes_Restart
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml Solutionu.h5 Solutionp.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml Solutionu.h5 Solutionp.h5
	EXEDEPS unsteadyNavierStokes_Restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNavierStokes_Restart
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES BFS3dCC.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES BFS3dCC.mesh
	EXEDEPS unsteadyNavierStokes_Restart
)


