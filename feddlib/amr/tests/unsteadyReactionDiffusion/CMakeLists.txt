tribits_add_executable(
	unsteadyReactionDiffusion
	SOURCES main.cpp
)

tribits_add_test(
	unsteadyReactionDiffusion
	NAME unsteadyReactionDiffusion
    ARGS
    NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_steadyReactionDiffusion
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml 
	EXEDEPS unsteadyReactionDiffusion
)

tribits_copy_files_to_binary_dir(mesh_unsteadyReactionDiffusion
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube_0_5.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_0_5.mesh
	EXEDEPS unsteadyReactionDiffusion
)

