tribits_add_executable(
	steadyGeneralizedNewtonianAssFE_PowerLaw 
	SOURCES main.cpp
)

tribits_add_test(
	steadyGeneralizedNewtonianAssFE_PowerLaw  
	NAME steadyGeneralizedNewtonianAssFE_PowerLaw  
	ARGS
#--precfile=parametersPrecMumps.xml
	NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_steadyGeneralizedNewtonianAssFE_PowerLaw   
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml parametersProblem.xml parametersSolver.xml parametersTeko.xml
	EXEDEPS steadyGeneralizedNewtonianAssFE_PowerLaw  
)

tribits_copy_files_to_binary_dir(mesh_steadyGeneralizedNewtonianAssFE_PowerLaw   
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES rectangle_200.mesh 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES rectangle_200.mesh 
	EXEDEPS steadyGeneralizedNewtonianAssFE_PowerLaw  
)




