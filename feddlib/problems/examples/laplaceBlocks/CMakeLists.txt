tribits_add_executable(
	laplaceBlocks
	SOURCES main.cpp
)

tribits_add_test(
	laplaceBlocks
    NAME laplaceBlocks
    ARGS
    NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_laplaceBlocks
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersPrec1.xml parametersPrec2.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersPrec1.xml parametersPrec2.xml parametersProblem.xml parametersSolver.xml 
	EXEDEPS laplaceBlocks
)


tribits_copy_files_to_binary_dir(meshes_laplaceBlocks
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES square.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES square.mesh
    EXEDEPS laplaceBlocks
)