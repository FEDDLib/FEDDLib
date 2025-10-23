INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

tribits_package(problems)
SET(${PACKAGE_NAME}_VERSION "1.0.0")

#
# B) Set up package-specific options
#


FOREACH(TRILINOS_PACKAGE_NAME in ${Trilinos_PACKAGE_LIST})
    IF(${TRILINOS_PACKAGE_NAME} STREQUAL "Ifpack2")
        SET(FEDD_HAVE_IFPACK2 TRUE)
    ENDIF()
ENDFOREACH()

FOREACH(TRILINOS_PACKAGE_NAME in ${Trilinos_PACKAGE_LIST})
    IF(${TRILINOS_PACKAGE_NAME} STREQUAL "NOX")
        SET(FEDD_HAVE_NOX TRUE)
    ENDIF()
ENDFOREACH()

FOREACH(TRILINOS_PACKAGE_NAME in ${Trilinos_PACKAGE_LIST})
    IF(${TRILINOS_PACKAGE_NAME} STREQUAL "Teko")
        SET(FEDD_HAVE_TEKO TRUE)
    ENDIF()
ENDFOREACH()

tribits_configure_file(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS ${HEADERS}
${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

ADD_SUBDIRECTORY(abstract)
SET(HEADERS ${HEADERS} ${problems_abstract_HEADERS})
SET(SOURCES ${SOURCES} ${problems_abstract_SOURCES})

ADD_SUBDIRECTORY(Solver)
SET(HEADERS ${HEADERS} ${Solver_HEADERS})
SET(SOURCES ${SOURCES} ${Solver_SOURCES})

ADD_SUBDIRECTORY(specific)
SET(HEADERS ${HEADERS} ${problems_specific_HEADERS})
SET(SOURCES ${SOURCES} ${problems_specific_SOURCES})

# This is the list of files excluded from the tarball
tribits_exclude_files()

# Add library target
tribits_add_library(
    problems
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
)

tribits_add_test_directories(Solver/tests)
tribits_add_test_directories(tests)
#TRIBITS_ADD_TEST_DIRECTORIES(tests/UnitTest)
tribits_add_example_directories(examples)


#
# Do standard postprocessing
#

tribits_package_postprocess()

