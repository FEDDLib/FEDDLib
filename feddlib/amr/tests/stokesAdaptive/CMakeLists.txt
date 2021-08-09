TRIBITS_ADD_EXECUTABLE(
	stokesAdaptive
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	stokesAdaptive
	NAME stokesAdaptive
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_stokesAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	EXEDEPS stokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_stokesAdaptive3D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES BFS3dC.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES BFS3dC.mesh
	EXEDEPS stokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_stokesAdaptive3D2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES ldc3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES ldc3D.mesh
	EXEDEPS stokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_stokesAdaptive2D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES Verfuerth2.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES Verfuerth2.mesh
	EXEDEPS stokesAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_stokesAdaptive2D2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES circle.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES circle.mesh
	EXEDEPS stokesAdaptive
)
