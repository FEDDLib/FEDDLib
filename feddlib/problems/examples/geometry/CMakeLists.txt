tribits_add_executable(
	geometry
    SOURCES main.cpp
)

tribits_add_test(
    geometry
    NAME geometry
    ARGS
    NUM_MPI_PROCS 4
)

tribits_copy_files_to_binary_dir(data_geometryTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS geometry
)

tribits_copy_files_to_binary_dir(meshes_geometryTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
	EXEDEPS geometry
)
