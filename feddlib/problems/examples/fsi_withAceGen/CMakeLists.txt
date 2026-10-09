TRIBITS_ADD_EXECUTABLE(
	fsi_withAceGen
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsi_withAceGen
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersPrecChem.xml parametersProblemFSI.xml parametersProblemFSCI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSCI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersPrecChem.xml parametersProblemFSI.xml parametersProblemFSCI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSCI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml
	EXEDEPS fsi_withAceGen
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsi_withAceGen
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fluid_safe.mesh solid_safe.mesh fsi_tube_fluid_1_5_cm.mesh fsi_tube_solid_1_5_cm.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fluid_safe.mesh solid_safe.mesh fsi_tube_fluid_1_5_cm.mesh fsi_tube_solid_1_5_cm.mesh
	EXEDEPS fsi_withAceGen
)
