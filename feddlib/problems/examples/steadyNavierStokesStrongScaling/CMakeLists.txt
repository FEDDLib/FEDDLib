TRIBITS_ADD_EXECUTABLE(
	steadyNavierStokes_scaling
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyNS
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersProblem.xml parametersSolver.xml parametersPrec.xml parametersPrecTeko.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersProblem.xml parametersSolver.xml parametersPrec.xml parametersPrecTeko.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	EXEDEPS steadyNavierStokes_scaling
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_steadyNS
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh fsi_benchmark_curved_fluid_5k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh fsi_benchmark_curved_fluid_5k.mesh
	EXEDEPS steadyNavierStokes_scaling
)
