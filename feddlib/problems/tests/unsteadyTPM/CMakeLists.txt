TRIBITS_ADD_EXECUTABLE(
	unsteadyTPM
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyTPM
	NAME unsteadyTPM
	ARGS
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadytpm
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersProblemSmall.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersProblemSmall.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS unsteadyTPM
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadytpm
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES tpm_square_h01.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES tpm_square_h01.mesh
    EXEDEPS unsteadyTPM
)
