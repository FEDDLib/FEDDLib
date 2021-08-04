TRIBITS_ADD_EXECUTABLE(
	aneurysmBenchmark
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	aneurysmBenchmark
	NAME aneurysmBenchmark
	ARGS
	NUM_MPI_PROCS 16
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_aneurysmBenchmark
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersProblem_Laplace.xml parametersSolver_Laplace.xml parametersPrec_Laplace.xml  parametersPrec_Stokes.xml parametersProblem_Stokes.xml parametersSolver_Stokes.xml parametersTeko_Stokes.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml parametersProblem_Laplace.xml parametersSolver_Laplace.xml parametersPrec_Laplace.xml parametersPrec_Stokes.xml parametersProblem_Stokes.xml parametersSolver_Stokes.xml parametersTeko_Stokes.xml
	EXEDEPS aneurysmBenchmark
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_aneurysmBenchmark
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES aneurysmBenchmark.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES aneurysmBenchmark.mesh
	EXEDEPS aneurysmBenchmark
)
