### neumannBC
tribits_add_executable(
	parameterList
	SOURCES parameterList.cpp
)

tribits_add_test(
	parameterList
	NAME parameterList
	ARGS
	NUM_MPI_PROCS 1
)

tribits_copy_files_to_binary_dir(data_parameterList
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parameterlist.yaml sublist1.yaml sublist2.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parameterlist.yaml sublist1.yaml sublist2.xml
	EXEDEPS parameterList
)

