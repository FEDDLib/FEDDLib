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

# Get Git commit hash (short)
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Create git_version.h file
set(GIT_VERSION_HEADER "${CMAKE_BINARY_DIR}/git_version.h")
file(WRITE ${GIT_VERSION_HEADER} "#pragma once\n")
file(APPEND ${GIT_VERSION_HEADER} "#define GIT_COMMIT_HASH \"${GIT_HASH}\"\n")
file(APPEND ${GIT_VERSION_HEADER} "#define GIT_BRANCH_NAME \"${GIT_BRANCH}\"\n")


INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/problems/abstract"
  COMPONENT ${PACKAGE_NAME}
)
