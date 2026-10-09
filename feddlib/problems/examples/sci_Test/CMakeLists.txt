TRIBITS_ADD_EXECUTABLE(
	sci_Test
    SOURCES main.cpp
)

# The case needs the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_TEST(
      sci_Test
      NAME sci_Test
      ARGS --problemfile=parametersProblemSCI_test.xml
      NUM_MPI_PROCS 4
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sciTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemSCI_test.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemSCI_test.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml
	EXEDEPS sci_Test
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES plaque_short_6_k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES plaque_short_6_k.mesh
	EXEDEPS sci_Test
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_develop_sciTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES artery_h_0_8_L_0_25.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES artery_h_0_8_L_0_25.mesh
	EXEDEPS sci_Test
)
