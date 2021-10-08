INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

TRIBITS_PACKAGE(amr)
SET(${PACKAGE_NAME}_VERSION "1.0.0")

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#FOREACH(TRILINOS_PACKAGE_NAME in ${Trilinos_PACKAGE_LIST})
#    IF(${TRILINOS_PACKAGE_NAME} STREQUAL "Ifpack2")
#        SET(FEDD_HAVE_IFPACK2 TRUE)
#    ENDIF()
#ENDFOREACH()

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


TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS ${HEADERS}
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

SET(HEADERS ${HEADERS} 
	AdaptiveMeshRefinement.hpp
	AdaptiveMeshRefinement_decl.hpp
	AdaptiveMeshRefinement_def.hpp
	RefinementFactory.hpp
	RefinementFactory_decl.hpp
	RefinementFactory_def.hpp
	ErrorEstimation.hpp
	ErrorEstimation_decl.hpp
	ErrorEstimation_def.hpp
	ExporterParaViewAMR.hpp
	ExporterParaViewAMR_decl.hpp
	ExporterParaViewAMR_def.hpp
CACHE INTERNAL "")

SET(SOURCES ${SOURCES} 
	AdaptiveMeshRefinement.cpp
	RefinementFactory.cpp
	ErrorEstimation.cpp
	ExporterParaViewAMR.cpp
CACHE INTERNAL "")
	
#SET(LOCAL_HEADERS)
#FOREACH(INC ${AMR_HEADERS})
#  STRING(REPLACE "amr/ " "" LOCAL_INC ${INC})
#  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
#ENDFOREACH()

#INSTALL(
#  FILES ${LOCAL_HEADERS}
#  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/amr"
#  COMPONENT ${PACKAGE_NAME}
#)

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()
# This is the list of files excluded from the tarball
TRIBITS_EXCLUDE_FILES()

# Add library target
TRIBITS_ADD_LIBRARY(
amr
HEADERS ${HEADERS}
SOURCES ${SOURCES}
)

TRIBITS_ADD_TEST_DIRECTORIES(tests)


#TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()

