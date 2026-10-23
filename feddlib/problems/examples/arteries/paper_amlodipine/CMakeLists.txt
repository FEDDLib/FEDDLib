TRIBITS_ADD_EXECUTABLE(
	paper_amlodipine
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_paper_amlodipine
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml
	EXEDEPS paper_amlodipine
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_paper_amlodipine
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES Artery_dan_SCI.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES Artery_dan_SCI.mesh
	EXEDEPS paper_amlodipine
)
