### mesh_nodeFlagsP2
tribits_add_executable(
    refinementFlag_uniform
    SOURCES refinementFlag_uniform.cpp
)

tribits_add_executable(
    refinementFlag_area
    SOURCES refinementFlag_area.cpp
)

tribits_add_test(
	refinementFlag_uniform
    NAME refinementFlagsTestUniform
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 4
)

tribits_add_test(
	refinementFlag_area
    NAME refinementFlagsTestArea
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_refinementFlag
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblem.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblem.xml
	EXEDEPS refinementFlag_uniform
)

tribits_copy_files_to_binary_dir(mesh_refinementFlag
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES meshNodeTestP23D_2.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES meshNodeTestP23D_2.mesh
    EXEDEPS refinementFlag_uniform
)

### 
