TRIBITS_ADD_EXECUTABLE(
	simpleTest
    SOURCES main.cpp
)

# The case needs the AceGen elements (Interface2)
#   simpleTest:       one time step (simulationParameters_test.xml)
#   simpleTest_heavy: five time steps (simulationParameters_test_heavy.xml), about
#                     5 minutes on 4 cores, in a directory of its own; HEAVY:
#                     configure with -D FEDDLib_TEST_CATEGORIES:STRING=HEAVY to add it
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_TEST(
      simpleTest
      NAME simpleTest
      ARGS --simulationsParameters=simulationParameters_test.xml
      NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_ADVANCED_TEST(
      simpleTest_heavy
      OVERALL_WORKING_DIRECTORY TEST_NAME
      TEST_0 COPY_FILES_TO_TEST_DIR materialParameters.xml preconditionerParameters_Chemistry.xml
          preconditionerParameters_Structure.xml simulationParameters_test_heavy.xml solverParameters.xml
      TEST_1 COPY_FILES_TO_TEST_DIR fsi_solid_2mm.mesh
          SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
      TEST_2 EXEC simpleTest
          ARGS --simulationsParameters=simulationParameters_test_heavy.xml
          NUM_MPI_PROCS 4
      OVERALL_NUM_MPI_PROCS 4
      CATEGORIES HEAVY
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_simpleTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml simulationParameters_test_heavy.xml solverParameters.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml simulationParameters_test_heavy.xml solverParameters.xml
	EXEDEPS simpleTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_simpleTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fsi_solid_2mm.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fsi_solid_2mm.mesh
	EXEDEPS simpleTest
)
