TRIBITS_ADD_EXECUTABLE(
	unsteadyNonLinElasticity
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyNonLinElasticity
	NAME unsteadyNonLinElasticity
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNonLinElasticity
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyNonLinElasticity
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNonLinElasticity
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES square_solid.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES square_solid.mesh
    EXEDEPS unsteadyNonLinElasticity
)
