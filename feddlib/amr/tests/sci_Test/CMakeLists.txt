TRIBITS_ADD_EXECUTABLE(
	sci_Test
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    sci_Test
    NAME sci_Test
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sciscaling
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrec.xml parametersPrecTeko.xml parametersPrecBlock.xml parametersProblemSCI.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml sbatch216Mono.sh sbatch1000Mono.sh sbatch4096Mono.sh parametersPrec_RGDSW.xml #parametersProblemSCI_RGDSW.xml parametersProblemSCI_GDSW.xml parametersPrec_RGDSW_CB.xml  parametersPrec_GDSW_CB.xml parametersProblemSCI_RGDSW_CB.xml parametersProblemSCI_GDSW_CB.xml parametersPrec_RGDSW_CB_CM.xml  parametersPrec_GDSW_CB_CM.xml parametersProblemSCI_RGDSW_CB_CM.xml parametersProblemSCI_GDSW_CB_CM.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrec.xml parametersPrecTeko.xml parametersPrecBlock.xml parametersProblemSCI.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml sbatch216Mono.sh sbatch1000Mono.sh sbatch4096Mono.sh parametersPrec_RGDSW.xml  #parametersProblemSCI_RGDSW.xml parametersProblemSCI_GDSW.xml parametersPrec_RGDSW_CB.xml  parametersPrec_GDSW_CB.xml parametersProblemSCI_RGDSW_CB.xml parametersProblemSCI_GDSW_CB.xml parametersPrec_RGDSW_CB_CM.xml  parametersPrec_GDSW_CB_CM.xml parametersProblemSCI_RGDSW_CB_CM.xml parametersProblemSCI_GDSW_CB_CM.xml
	EXEDEPS sci_Test
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciscaling
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES cube_h_0_05.mesh cube_h_0_1.mesh cube_h_1.mesh artery_h_025.mesh artery_h_0_1.mesh  artery_h_0_8_L_0_25.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_h_0_05.mesh cube_h_0_1.mesh cube_h_1.mesh artery_h_025.mesh artery_h_0_1.mesh  artery_h_0_8_L_0_25.mesh
	EXEDEPS sci_Test 
)
