TRIBITS_ADD_EXECUTABLE(
	fsi_arteryBenchmark
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsiArteryBenchmark
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersProblemFSI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml  parametersProblemFSI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	EXEDEPS fsi_arteryBenchmark
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsiArteryBenchmark
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fluidBenchmark3.mesh structureBenchmark3.mesh 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fluidBenchmark3.mesh structureBenchmark3.mesh 
	EXEDEPS fsi_arteryBenchmark
)
