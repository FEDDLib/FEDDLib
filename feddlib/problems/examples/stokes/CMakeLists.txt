tribits_add_executable(
	stokes
	SOURCES main.cpp
)

tribits_add_test(
	stokes
	NAME stokes
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_stokes
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	EXEDEPS stokes
)
