TRIBITS_ADD_EXECUTABLE(
	laplaceAdaptive
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	laplaceAdaptive
    NAME laplaceAdaptive2D
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 8
)

TRIBITS_ADD_TEST(
	laplaceAdaptive
    NAME laplaceAdaptive3D
    ARGS 
	--problemfile=parametersProblem3D.xml
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplaceAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml 
	EXEDEPS laplaceAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_laplaceAdaptive2D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES lShape4.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape4.mesh
	EXEDEPS laplaceAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_laplaceAdaptive2D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES lShape1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape1.mesh
	EXEDEPS laplaceAdaptive
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_laplaceAdaptive3D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES lShape3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape3D.mesh
	EXEDEPS laplaceAdaptive
)


