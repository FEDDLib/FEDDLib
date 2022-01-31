INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### FE
TRIBITS_ADD_EXECUTABLE(
fe_test
SOURCES fe_test.cpp
)

TRIBITS_ADD_TEST(
fe_test
NAME fe_test
ARGS
NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fe_test
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblem.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblem.xml 
	EXEDEPS fe_test
)
