tribits_add_executable(
	stokesAdaptive
	SOURCES main.cpp
)

tribits_add_test(
	stokesAdaptive
	NAME stokesAdaptive
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

tribits_copy_files_to_binary_dir(data_stokesAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	EXEDEPS stokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_stokesAdaptive3D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES BFS3dC.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES BFS3dC.mesh
	EXEDEPS stokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_stokesAdaptive3D2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES ldc3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES ldc3D.mesh
	EXEDEPS stokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_stokesAdaptive2D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES Verfuerth2.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES Verfuerth2.mesh
	EXEDEPS stokesAdaptive
)

tribits_copy_files_to_binary_dir(mesh_stokesAdaptive2D2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES circle.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES circle.mesh
	EXEDEPS stokesAdaptive
)
