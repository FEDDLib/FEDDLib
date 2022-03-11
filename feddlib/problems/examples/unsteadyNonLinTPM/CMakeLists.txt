TRIBITS_ADD_EXECUTABLE(
	unsteadyNonLinTPM
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	unsteadyNonLinTPM
	NAME unsteadyNonLinTPM
	ARGS
	NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
    unsteadyNonLinTPM
    NAME unsteadyNonLinTPM3D
    ARGS "--problemfile=parametersProblem3D.xml --precfile=parametersPrec3D.xml"
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadynonlintpm
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersPrec3D.xml parametersProblem.xml parametersProblem3D.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersPrec3D.xml parametersProblem.xml parametersProblem3D.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS unsteadyNonLinTPM
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadynonlintpm
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES tpm_square_h01.mesh bucket3D.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES tpm_square_h01.mesh bucket3D.mesh
    EXEDEPS unsteadyNonLinTPM
)
