TRIBITS_ADD_EXECUTABLE(
	steadyLinElasAssFE
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyLinElasAssFE
	NAME steadyLinElasAssFE
	ARGS
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyLinElasAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyLinElasAssFE
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyLinElasAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube_h_1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_h_1.mesh
	EXEDEPS steadyLinElasAssFE
)
