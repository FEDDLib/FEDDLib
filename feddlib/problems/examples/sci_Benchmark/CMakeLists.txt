TRIBITS_ADD_EXECUTABLE(
	sci_Benchmark
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    sci_Benchmark
    NAME sci_Benchmark
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sciBenchmark
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrec.xml parametersProblemSCI.xml parametersSolverSCI.xml parametersPrecStructure.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrec.xml parametersProblemSCI.xml parametersSolverSCI.xml parametersPrecStructure.xml
	EXEDEPS sci_Benchmark
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciBenchmark
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES structureBenchmark3.mesh 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES structureBenchmark3.mesh 
	EXEDEPS sci_Benchmark
)
