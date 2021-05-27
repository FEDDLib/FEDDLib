TRIBITS_ADD_EXECUTABLE(
	unsteadyLinElas
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyLinElas
	NAME unsteadyLinElas
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 2
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyLinElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyLinElas
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyLinElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES dfg_fsi_solid_h002.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_solid_h002.mesh
	EXEDEPS unsteadyLinElas
)
