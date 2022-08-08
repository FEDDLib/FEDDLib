### mesh_nodeFlagsP2
TRIBITS_ADD_EXECUTABLE(
    refinementFlag_uniform
    SOURCES refinementFlag_uniform.cpp
)

TRIBITS_ADD_EXECUTABLE(
    refinementFlag_area
    SOURCES refinementFlag_area.cpp
)

TRIBITS_ADD_TEST(
	refinementFlag_uniform
    NAME refinementFlagsTestUniform
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
	refinementFlag_area
    NAME refinementFlagsTestArea
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_refinementFlag
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblem.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblem.xml
	EXEDEPS refinementFlag_uniform
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_refinementFlag
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES meshNodeTestP23D_2.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES meshNodeTestP23D_2.mesh
    EXEDEPS refinementFlag_uniform
)

### 
