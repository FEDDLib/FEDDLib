TRIBITS_ADD_EXECUTABLE(
	ccb_hollow_cylinder
    SOURCES main.cpp
)

# The case needs the AceGen elements (Interface2)
if(TPL_ENABLE_AceGENInterface)
  TRIBITS_ADD_TEST(
      ccb_hollow_cylinder
      NAME ccb_hollow_cylinder
      ARGS --simulationsParameters=simulationParameters_test.xml
      NUM_MPI_PROCS 4
  )

  # Adaptive time stepping (DAESolverInTime::advanceInTimeSCI), in a directory of its own:
  # simulationParameters_adaptive.xml, in which no time step fails, gives the corner values of the
  # fixed time step size (simulationParameters_test.xml); in simulationParameters_adaptive_newton.xml
  # time steps the Newton iteration does not converge in are repeated with a smaller size.
  # simulationParameters_adaptive_element.xml: a time step the elements cannot compute their state
  # in is repeated, and accepted at the Minimum dt (Accept Element Failures At Minimum dt).
  TRIBITS_ADD_ADVANCED_TEST(
      ccb_hollow_cylinder_adaptive
      OVERALL_WORKING_DIRECTORY TEST_NAME
      TEST_0 COPY_FILES_TO_TEST_DIR materialParameters.xml preconditionerParameters_Chemistry.xml
          preconditionerParameters_Structure.xml solverParameters.xml simulationParameters_test.xml simulationParameters_adaptive.xml simulationParameters_adaptive_newton.xml simulationParameters_adaptive_element.xml
      TEST_1 COPY_FILES_TO_TEST_DIR hollow_cylinder_p1.mesh
          SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/ccb_hollow_cylinder
      TEST_2 EXEC ccb_hollow_cylinder
          ARGS --simulationsParameters=simulationParameters_test.xml
          NUM_MPI_PROCS 4
      TEST_3 CMND ${CMAKE_COMMAND}
          ARGS -E rename cornerValue.txt cornerValue_fixed.txt
      TEST_4 EXEC ccb_hollow_cylinder
          ARGS --simulationsParameters=simulationParameters_adaptive.xml
          NUM_MPI_PROCS 4
          PASS_REGULAR_EXPRESSION "adaptive. repeated 0 time step"
      TEST_5 CMND ${CMAKE_COMMAND}
          ARGS -E compare_files cornerValue.txt cornerValue_fixed.txt
      TEST_6 EXEC ccb_hollow_cylinder
          ARGS --simulationsParameters=simulationParameters_adaptive_newton.xml
          NUM_MPI_PROCS 4
          PASS_REGULAR_EXPRESSION_ALL "the Newton iteration did not converge.; repeating it with" "adaptive. repeated [1-9]"
      TEST_7 EXEC ccb_hollow_cylinder
          ARGS --simulationsParameters=simulationParameters_adaptive_element.xml
          NUM_MPI_PROCS 4
          PASS_REGULAR_EXPRESSION_ALL "an element could not compute its state.; repeating it with" "accepted the time step of 10 although" "adaptive. repeated [1-9]"
      OVERALL_NUM_MPI_PROCS 4
      FAIL_FAST
  )
endif()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_ccb_hollow_cylinder
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml simulationParameters_adaptive.xml simulationParameters_adaptive_newton.xml simulationParameters_adaptive_element.xml solverParameters.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters.xml preconditionerParameters_Chemistry.xml preconditionerParameters_Structure.xml simulationParameters.xml simulationParameters_test.xml simulationParameters_adaptive.xml simulationParameters_adaptive_newton.xml simulationParameters_adaptive_element.xml solverParameters.xml
	EXEDEPS ccb_hollow_cylinder
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_ccb_hollow_cylinder
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/ccb_hollow_cylinder
    SOURCE_FILES hollow_cylinder_p1.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES hollow_cylinder_p1.mesh
	EXEDEPS ccb_hollow_cylinder
)
