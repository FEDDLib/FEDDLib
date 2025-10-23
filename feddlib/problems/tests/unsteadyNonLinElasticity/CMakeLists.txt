tribits_add_executable(
	unsteadyNonLinElasticityAssFE
	SOURCES main.cpp
)

tribits_add_test(
	unsteadyNonLinElasticityAssFE
	NAME unsteadyNonLinElasticityAssFE
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_unsteadyNonLinElasticityAssFE
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS unsteadyNonLinElasticityAssFE
)

tribits_copy_files_to_binary_dir(mesh_unsteadyNonLinElasticityAssFE
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES cube_h_1.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES cube_h_1.mesh
    EXEDEPS unsteadyNonLinElasticityAssFE
)
