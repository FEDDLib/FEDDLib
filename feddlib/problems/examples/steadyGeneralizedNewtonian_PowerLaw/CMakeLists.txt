TRIBITS_ADD_EXECUTABLE(
	steadyGeneralizedNewtonianAssFE_PowerLaw 
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyGeneralizedNewtonianAssFE_PowerLaw  
	NAME steadyGeneralizedNewtonianAssFE_PowerLaw  
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyGeneralizedNewtonianAssFE_PowerLaw   
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyGeneralizedNewtonianAssFE_PowerLaw  
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyGeneralizedNewtonianAssFE_PowerLaw   
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES rectangle_200.mesh 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES rectangle_200.mesh 
	EXEDEPS steadyGeneralizedNewtonianAssFE_PowerLaw  
)




