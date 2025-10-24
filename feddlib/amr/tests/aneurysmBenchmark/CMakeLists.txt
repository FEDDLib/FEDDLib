tribits_add_executable(
	aneurysmBenchmark
	SOURCES main.cpp
)

tribits_add_test(
	aneurysmBenchmark
	NAME aneurysmBenchmark
	ARGS
	NUM_MPI_PROCS 16
)

tribits_copy_files_to_binary_dir(data_aneurysmBenchmark
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersProblem_Laplace.xml parametersSolver_Laplace.xml parametersPrec_Laplace.xml  
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersProblem_Laplace.xml parametersSolver_Laplace.xml parametersPrec_Laplace.xml 
	EXEDEPS aneurysmBenchmark
)

tribits_copy_files_to_binary_dir(mesh_aneurysmBenchmark
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/Adaptive
	SOURCE_FILES aneurysmBenchmark.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES aneurysmBenchmark.mesh
	EXEDEPS aneurysmBenchmark
)
