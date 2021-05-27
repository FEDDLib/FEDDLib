SET(Solver_HEADERS
    Solver/DAESolverInTime_decl.hpp
    Solver/DAESolverInTime_def.hpp
    Solver/DAESolverInTime.hpp
	Solver/LinearSolver_decl.hpp
	Solver/LinearSolver_def.hpp
	Solver/LinearSolver.hpp
    Solver/NonLinearSolver.hpp
    Solver/NonLinearSolver_decl.hpp
    Solver/NonLinearSolver_def.hpp
    Solver/Preconditioner.hpp
    Solver/Preconditioner_decl.hpp
    Solver/Preconditioner_def.hpp
    Solver/PreconditionerOperator.hpp
    Solver/PreconditionerOperator_decl.hpp
    Solver/PreconditionerOperator_def.hpp
    Solver/PrecBlock2x2.hpp
    Solver/PrecBlock2x2_decl.hpp
    Solver/PrecBlock2x2_def.hpp
    Solver/PrecOpFaCSI.hpp
    Solver/PrecOpFaCSI_decl.hpp
    Solver/PrecOpFaCSI_def.hpp
    Solver/TimeSteppingTools.hpp
CACHE INTERNAL "")

SET(Solver_SOURCES
	Solver/DAESolverInTime.cpp
    Solver/LinearSolver.cpp
    Solver/NonLinearSolver.cpp
    Solver/Preconditioner.cpp
    Solver/PreconditionerOperator.cpp
    Solver/PrecBlock2x2.cpp
    Solver/PrecOpFaCSI.cpp
	Solver/TimeSteppingTools.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Solver_HEADERS})
  STRING(REPLACE "Solver/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/problems/Solver"
  COMPONENT ${PACKAGE_NAME}
)
