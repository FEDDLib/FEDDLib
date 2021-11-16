TRIBITS_ADD_EXECUTABLE(
	unsteadyNavierStokesAdaptive
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyNavierStokesAdaptive
	NAME unsteadyNavierStokesAdaptive
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNavierStokesAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS unsteadyNavierStokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNavierStokesAdaptiveBFS
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES bfs.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES bfs.mesh
	EXEDEPS unsteadyNavierStokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNavierStokesAdaptiveBFS3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES bfs3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES bfs3D.mesh
	EXEDEPS unsteadyNavierStokesAdaptive
)


