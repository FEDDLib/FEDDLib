TRIBITS_ADD_EXECUTABLE(
	laplaceAdaptiveInterpol
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	laplaceAdaptiveInterpol
    NAME laplaceAdaptiveInterpol2D
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 1
)

TRIBITS_ADD_TEST(
	laplaceAdaptiveInterpol
    NAME laplaceAdaptiveInterpol3D
    ARGS 
	--problemfile=parametersProblem3D.xml
    NUM_MPI_PROCS 1
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplaceAdaptiveInterpol
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml 
	EXEDEPS laplaceAdaptiveInterpol
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_laplaceAdaptiveInterpol2D
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES lShape4.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape4.mesh
	EXEDEPS laplaceAdaptiveInterpol
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_laplaceAdaptiveInterpol2D2
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES square_paper.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES square_paper.mesh
	EXEDEPS laplaceAdaptiveInterpol
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_laplaceAdaptiveInterpol2D3
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES lShape3D.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES lShape3D.mesh
	EXEDEPS laplaceAdaptiveInterpol
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_laplaceAdaptiveInterpol3D1
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES cube_14.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
        DEST_FILES cube_14.mesh
	EXEDEPS laplaceAdaptiveInterpol
)


