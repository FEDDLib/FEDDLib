set(Solver_HEADERS
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

set(Solver_SOURCES
	Solver/DAESolverInTime.cpp
    Solver/LinearSolver.cpp
    Solver/NonLinearSolver.cpp
    Solver/Preconditioner.cpp
    Solver/PreconditionerOperator.cpp
    Solver/PrecBlock2x2.cpp
    Solver/PrecOpFaCSI.cpp
	Solver/TimeSteppingTools.cpp
CACHE INTERNAL "")

set(LOCAL_HEADERS)
foreach(INC ${Solver_HEADERS})
  string(REPLACE "Solver/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
