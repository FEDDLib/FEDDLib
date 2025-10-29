if(${MPI_EXEC_MAX_NUMPROCS} LESS 8)
    message("  The performace test steadyLinElas_Perf requires MPI_EXEC_MAX_NUMPROCS >= 8. The test will not be performed.")
endif()

tribits_add_executable(
	steadyLinElas_Perf
	SOURCES main.cpp
)

tribits_add_test(
	steadyLinElas_Perf
	NAME steadyLinElas_Perf
	ARGS
	NUM_MPI_PROCS 8
)

tribits_copy_files_to_binary_dir(data_steadyLinElas_Perf
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyLinElas_Perf
)
