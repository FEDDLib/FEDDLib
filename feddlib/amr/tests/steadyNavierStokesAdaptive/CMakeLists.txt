tribits_add_executable(
	steadyNavierStokesAdaptive
	SOURCES main.cpp
)

tribits_add_test(
	steadyNavierStokesAdaptive
	NAME steadyNavierStokesAdaptive
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_steadyNavierStokesAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyNavierStokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_steadyNavierStokesAdaptiveBFS
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES bfs.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES bfs.mesh
	EXEDEPS steadyNavierStokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_steadyNavierStokesAdaptiveBFS3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES bfs3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES bfs3D.mesh
	EXEDEPS steadyNavierStokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_steadyNavierStokesAdaptiveLDC
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES ldc2D2.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES ldc2D2.mesh
	EXEDEPS steadyNavierStokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_steadyNavierStokesAdaptiveLDC3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES ldc3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES ldc3D.mesh
	EXEDEPS steadyNavierStokesAdaptive
)
