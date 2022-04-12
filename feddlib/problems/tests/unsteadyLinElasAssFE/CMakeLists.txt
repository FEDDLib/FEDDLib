TRIBITS_ADD_EXECUTABLE(
	unsteadyLinElasAssFE
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyLinElasAssFE
	NAME unsteadyLinElasAssFE
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyLinElasAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyLinElasAssFE
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyLinElasAssFE
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
	EXEDEPS unsteadyLinElasAssFE
)
