INCLUDE(TribitsCopyFilesToBinaryDir)
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(AddSubdirectories)

### neumannBC
TRIBITS_ADD_EXECUTABLE(
	parameterList
	SOURCES parameterList.cpp
)

TRIBITS_ADD_TEST(
	parameterList
	NAME parameterList
	ARGS
	NUM_MPI_PROCS 1
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_parameterList
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parameterlist.yaml sublist1.yaml sublist2.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parameterlist.yaml sublist1.yaml sublist2.xml
	EXEDEPS parameterList
)

