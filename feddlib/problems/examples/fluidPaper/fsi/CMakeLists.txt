TRIBITS_ADD_EXECUTABLE(
	fsiPaper
    SOURCES main.cpp
)

#TRIBITS_ADD_TEST(
#    fsiPaper
#    NAME fsiPaper
#    ARGS
#    NUM_MPI_PROCS 2
#)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsiTestPaper
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersProblemFSI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml parametersProblemFSITeko_SIMPLE_R.xml parametersProblemFSITeko_SIMPLE_R_CB.xml parametersPrecFluidTeko_SIMPLE_R.xml parametersPrecFluidTeko_SIMPLE_R_CB.xml sbatch240.sh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml  parametersProblemFSI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml parametersProblemFSITeko_SIMPLE_R.xml parametersProblemFSITeko_SIMPLE_R_CB.xml parametersPrecFluidTeko_SIMPLE_R.xml parametersPrecFluidTeko_SIMPLE_R_CB.xml sbatch240.sh
	EXEDEPS fsiPaper
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsiTestPaper
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES richter_3D_fluid.mesh richter_3D_solid.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES richter_3D_fluid.mesh richter_3D_solid.mesh
	EXEDEPS fsiPaper
)
