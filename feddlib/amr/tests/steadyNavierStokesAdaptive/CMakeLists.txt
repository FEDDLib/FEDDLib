TRIBITS_ADD_EXECUTABLE(
	steadyNavierStokesAdaptive
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyNavierStokesAdaptive
	NAME steadyNavierStokesAdaptive
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyNavierStokesAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyNavierStokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNavierStokesAdaptiveBFS
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES bfs.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES bfs.mesh
	EXEDEPS steadyNavierStokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNavierStokesAdaptiveBFS3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES bfs3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES bfs3D.mesh
	EXEDEPS steadyNavierStokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNavierStokesAdaptiveLDC
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES ldc2D2.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES ldc2D2.mesh
	EXEDEPS steadyNavierStokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNavierStokesAdaptiveLDC3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES ldc3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES ldc3D.mesh
	EXEDEPS steadyNavierStokesAdaptive
)
