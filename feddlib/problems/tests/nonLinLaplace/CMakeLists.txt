TRIBITS_ADD_EXECUTABLE(
  nonLinLaplace
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
  nonLinLaplace
  NAME nonLinLaplace
	ARGS
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_nonLinLaplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
  EXEDEPS nonLinLaplace
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_nonLinLaplace
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES square.mesh
  EXEDEPS nonLinLaplace
)
