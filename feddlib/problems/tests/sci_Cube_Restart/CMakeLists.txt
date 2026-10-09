TRIBITS_ADD_EXECUTABLE(
	sci_Cube_Restart
    SOURCES main.cpp
)

# Restart test, two runs of the same case: the uninterrupted run
# (parametersProblemSCI.xml) writes checkpoints at t = 1.5 and 2.5 to
# restartCheckpoints/; the restart (parametersProblemSCI_restart.xml) continues
# from the checkpoint at t = 1.5, during fiber reorientation, and must reproduce
# the uninterrupted solution at t = 2.5.
# The structured cube needs n^3 ranks; 8 would exceed TriBITS' default
# MPI_EXEC_MAX_NUMPROCS of 4 (the test would be skipped), so it runs on 1.
# The case needs the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_ADVANCED_TEST(
      sci_Cube_Restart
      TEST_0 EXEC sci_Cube_Restart
          ARGS --problemfile=parametersProblemSCI.xml
          NUM_MPI_PROCS 1
      TEST_1 EXEC sci_Cube_Restart
          ARGS --problemfile=parametersProblemSCI_restart.xml
          NUM_MPI_PROCS 1
      OVERALL_NUM_MPI_PROCS 1
      FAIL_FAST
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sci_Cube_Restart
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemSCI_restart.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemSCI_restart.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml
	EXEDEPS sci_Cube_Restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sci_Cube_Restart
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES plaque_short_6_k.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES plaque_short_6_k.mesh
	EXEDEPS sci_Cube_Restart
)
