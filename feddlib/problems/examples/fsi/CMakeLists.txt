TRIBITS_ADD_EXECUTABLE(
	fsi
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    fsi
    NAME fsi
    ARGS
    NUM_MPI_PROCS 2
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsiTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersProblemFSI.xml parametersSolverFSI.xml parametersPrecFluidBlock.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml laplace_parabolic_parabolic_fsi_fluid_length_0_5_mm_P1.h5 laplace_parabolic_parabolic_fsi_fluid_length_0_5_mm_P2.h5 laplace_parabolic_parabolic_fsi_fluid_2mm_P1.h5 laplace_parabolic_parabolic_fsi_fluid_2mm_P2.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml  parametersProblemFSI.xml parametersSolverFSI.xml parametersPrecFluidBlock.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml laplace_parabolic_parabolic_fsi_fluid_length_0_5_mm_P1.h5 laplace_parabolic_parabolic_fsi_fluid_length_0_5_mm_P2.h5 laplace_parabolic_parabolic_fsi_fluid_2mm_P1.h5 laplace_parabolic_parabolic_fsi_fluid_2mm_P2.h5
	EXEDEPS fsi
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsiTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh fsi_fluid_length_0_5_mm.mesh fsi_solid_length_0_5_mm.mesh fsi_fluid_2mm.mesh fsi_solid_2mm.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh fsi_fluid_length_0_5_mm.mesh fsi_solid_length_0_5_mm.mesh fsi_fluid_2mm.mesh fsi_solid_2mm.mesh
	EXEDEPS fsi
)
