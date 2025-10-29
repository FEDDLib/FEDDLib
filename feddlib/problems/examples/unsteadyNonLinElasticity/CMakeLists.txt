tribits_add_executable(
	unsteadyNonLinElasticity
	SOURCES main.cpp
)

tribits_add_test(
	unsteadyNonLinElasticity
	NAME unsteadyNonLinElasticity
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_unsteadyNonLinElasticity
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyNonLinElasticity
)

tribits_copy_files_to_binary_dir(mesh_unsteadyNonLinElasticity
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES square_solid.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES square_solid.mesh
    EXEDEPS unsteadyNonLinElasticity
)
