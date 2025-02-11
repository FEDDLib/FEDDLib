IF(${MPI_EXEC_MAX_NUMPROCS} LESS 8)
    MESSAGE("  The performace test steadyLinElas_Perf requires MPI_EXEC_MAX_NUMPROCS >= 8. The test will not be performed.")
ENDIF()

TRIBITS_ADD_EXECUTABLE(
	steadyLinElas_Perf
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyLinElas_Perf
	NAME steadyLinElas_Perf
	ARGS
	NUM_MPI_PROCS 8
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyLinElas_Perf
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyLinElas_Perf
)
