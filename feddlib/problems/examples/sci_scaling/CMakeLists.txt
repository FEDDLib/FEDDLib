TRIBITS_ADD_EXECUTABLE(
	sci_scaling
    SOURCES main.cpp
)

# The case needs the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_TEST(
      sci_scaling
      NAME sci_scaling
      ARGS --problemfile=parametersProblemSCI_test.xml
      NUM_MPI_PROCS 4
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sciscaling
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemSCI_test.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml parametersProblemStructure_plaque_full.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemSCI_test.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml parametersProblemStructure_plaque_full.xml
	EXEDEPS sci_scaling
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciscaling
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES  cube_h_1.mesh artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh Artery_Small.mesh Artery_Long.mesh plaque_4k.mesh arterie_with_plaque_initma_media_adventitia.mesh plaque_short_6_k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_h_1.mesh  artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh Artery_Small.mesh Artery_Long.mesh plaque_4k.mesh arterie_with_plaque_initma_media_adventitia.mesh plaque_short_6_k.mesh
	EXEDEPS sci_scaling 
)
