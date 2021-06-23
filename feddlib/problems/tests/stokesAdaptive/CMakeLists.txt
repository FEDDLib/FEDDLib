TRIBITS_ADD_EXECUTABLE(
	stokesAdaptive
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	stokesAdaptive
	NAME stokesAdaptive
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_stokesAdaptive
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersPrecBlock.xml
	EXEDEPS stokesAdaptive
)
