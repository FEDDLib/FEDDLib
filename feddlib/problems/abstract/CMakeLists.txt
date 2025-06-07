SET(problems_abstract_HEADERS
    #abstract/CoupledProblem_decl.hpp
    #abstract/CoupledProblem_def.hpp
    #abstract/CoupledProblem.hpp
	abstract/MinPrecProblem_decl.hpp
	abstract/MinPrecProblem_def.hpp
	abstract/MinPrecProblem.hpp
	abstract/NonLinearProblem_decl.hpp
	abstract/NonLinearProblem_def.hpp
	abstract/NonLinearProblem_def.hpp
	abstract/NonLinearProblem_def.hpp
	abstract/NonLinearProblem.hpp
	abstract/Problem_decl.hpp
	abstract/Problem_def.hpp
	abstract/Problem.hpp
	abstract/TimeProblem_decl.hpp
	abstract/TimeProblem_def.hpp
	abstract/TimeProblem.hpp
CACHE INTERNAL "")

SET(problems_abstract_SOURCES
    #abstract/CoupledProblem.cpp
	abstract/MinPrecProblem.cpp
	abstract/NonLinearProblem.cpp
	abstract/Problem.cpp
    abstract/TimeProblem.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Problems_HEADERS})
  STRING(REPLACE "abstract/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/problems/abstract"
  COMPONENT ${PACKAGE_NAME}
)
