TRIBITS_ADD_EXECUTABLE(
	steadyGeneralizedNonNewtonianAssFE_PowerLaw 
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyGeneralizedNonNewtonianAssFE_PowerLaw  
	NAME steadyGeneralizedNonNewtonianAssFE_PowerLaw  
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyGeneralizedNonNewtonianAssFE_PowerLaw   
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyGeneralizedNonNewtonianAssFE_PowerLaw  
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyGeneralizedNonNewtonianAssFE_PowerLaw   
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES rectangle_200.mesh 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES rectangle_200.mesh 
	EXEDEPS steadyGeneralizedNonNewtonianAssFE_PowerLaw  
)




