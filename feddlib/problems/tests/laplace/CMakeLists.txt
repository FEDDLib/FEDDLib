TRIBITS_ADD_EXECUTABLE(
	laplace
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	laplace
    NAME laplace
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplace
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml 
	EXEDEPS laplace
)
