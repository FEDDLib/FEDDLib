TRIBITS_ADD_EXECUTABLE(
	laplaceBlocks
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	laplaceBlocks
    NAME laplaceBlocks
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplaceBlocks
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersPrec1.xml parametersPrec2.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersPrec1.xml parametersPrec2.xml parametersProblem.xml parametersSolver.xml 
	EXEDEPS laplaceBlocks
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_laplaceBlocks
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES square.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES square.mesh
    EXEDEPS laplaceBlocks
)