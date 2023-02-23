TRIBITS_ADD_EXECUTABLE(
	steadyNonLinElasAssFE
	SOURCES main.cpp
)

#TRIBITS_ADD_TEST(
#	steadyNonLinElasAssFE
#	NAME steadyNonLinElasAssFE
#	ARGS
#	NUM_MPI_PROCS 4
#)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyNonLinElasAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyNonLinElasAssFE
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNonLinElasAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube_0_1.mesh cube_0_5.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_0_1.mesh cube_0_5.mesh
	EXEDEPS steadyNonLinElasAssFE
)
