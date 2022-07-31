TRIBITS_ADD_EXECUTABLE(
	sci_linear
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    sci_linear
    NAME sci_linear
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sciTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrec.xml parametersProblemSCI.xml parametersSolverSCI.xml parametersPrecStructure.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrec.xml parametersProblemSCI.xml parametersSolverSCI.xml parametersPrecStructure.xml
	EXEDEPS sci_linear
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES cube_0_1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_0_1.mesh
	EXEDEPS sci_linear
)
