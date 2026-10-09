TRIBITS_ADD_EXECUTABLE(
	artery_dan_restart
    SOURCES main.cpp
)

# Restart tests, two runs of the same case each: the uninterrupted run
# (simulationParameters*.xml) writes checkpoints at t = 0.05 and 0.1 to
# restartCheckpoints/; the restart (simulationParameters*_restart.xml) continues
# from the checkpoint at t = 0.05, during the load ramp, and must reproduce the
# uninterrupted solution at t = 0.1 (restart during fiber reorientation is
# tested by sci_Cube_Restart). Each test runs in a directory of its own.
#   artery_restart:     a quarter of a short artery segment (artery_h_0_8_L_0_25_clean.mesh:
#                       artery_h_0_8_L_0_25.mesh without its two vertices that no element
#                       uses; the checkpoint files need node IDs without gaps)
#   artery_dan_restart: the artery_dan case (Artery_dan_SCI.mesh), about 20 minutes
#                       on 4 cores; HEAVY: configure with
#                       -D FEDDLib_TEST_CATEGORIES:STRING=HEAVY to add it
# The cases need the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_ADVANCED_TEST(
      artery_restart
      OVERALL_WORKING_DIRECTORY TEST_NAME
      TEST_0 COPY_FILES_TO_TEST_DIR materialParameters_dan.xml preconditionerParameters_Chemistry.xml
          preconditionerParameters_Structure.xml simulationParameters_short.xml
          simulationParameters_short_restart.xml solverParameters.xml
      TEST_1 COPY_FILES_TO_TEST_DIR artery_h_0_8_L_0_25_clean.mesh
          SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
      TEST_2 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters_short.xml
          NUM_MPI_PROCS 4
      TEST_3 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters_short_restart.xml
          NUM_MPI_PROCS 4
      OVERALL_NUM_MPI_PROCS 4
      FAIL_FAST
  )

  # artery_restart_adaptive: the same with adaptive time stepping and a repeated time step at
  # a checkpoint (see simulationParameters_short_adaptive.xml)
  TRIBITS_ADD_ADVANCED_TEST(
      artery_restart_adaptive
      OVERALL_WORKING_DIRECTORY TEST_NAME
      TEST_0 COPY_FILES_TO_TEST_DIR materialParameters_dan.xml preconditionerParameters_Chemistry.xml
          preconditionerParameters_Structure.xml simulationParameters_short_adaptive.xml
          simulationParameters_short_adaptive_restart.xml solverParameters.xml
      TEST_1 COPY_FILES_TO_TEST_DIR artery_h_0_8_L_0_25_clean.mesh
          SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
      TEST_2 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters_short_adaptive.xml
          NUM_MPI_PROCS 4
          PASS_REGULAR_EXPRESSION_ALL "did not converge.; repeating it with" "adaptive. repeated [1-9]"
          ALWAYS_FAIL_ON_NONZERO_RETURN
      TEST_3 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters_short_adaptive_restart.xml
          NUM_MPI_PROCS 4
      OVERALL_NUM_MPI_PROCS 4
      FAIL_FAST
  )

  TRIBITS_ADD_ADVANCED_TEST(
      artery_dan_restart
      OVERALL_WORKING_DIRECTORY TEST_NAME
      TEST_0 COPY_FILES_TO_TEST_DIR materialParameters_dan.xml preconditionerParameters_Chemistry.xml
          preconditionerParameters_Structure.xml simulationParameters.xml
          simulationParameters_restart.xml solverParameters.xml
      TEST_1 COPY_FILES_TO_TEST_DIR Artery_dan_SCI.mesh
          SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
      TEST_2 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters.xml
          NUM_MPI_PROCS 4
      TEST_3 EXEC artery_dan_restart
          ARGS --simulationsParameters=simulationParameters_restart.xml
          NUM_MPI_PROCS 4
      OVERALL_NUM_MPI_PROCS 4
      TIMEOUT 3600
      CATEGORIES HEAVY
      FAIL_FAST
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_artery_dan_restart
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml
	preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_restart.xml
	simulationParameters_short.xml simulationParameters_short_restart.xml simulationParameters_short_adaptive.xml
	simulationParameters_short_adaptive_restart.xml solverParameters.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml
	preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_restart.xml
	simulationParameters_short.xml simulationParameters_short_restart.xml simulationParameters_short_adaptive.xml
	simulationParameters_short_adaptive_restart.xml solverParameters.xml
	EXEDEPS artery_dan_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_artery_dan_restart
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES Artery_dan_SCI.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES Artery_dan_SCI.mesh
	EXEDEPS artery_dan_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_artery_restart
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES artery_h_0_8_L_0_25_clean.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES artery_h_0_8_L_0_25_clean.mesh
	EXEDEPS artery_dan_restart
)
