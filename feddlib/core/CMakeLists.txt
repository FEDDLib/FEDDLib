INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

tribits_package(core)
SET(${PACKAGE_NAME}_VERSION "1.0.0")

#
# B) Set up package-specific options
#


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


IF(HAVE_ACEGENINTERFACE)
  message (" Have ACEGEN Interface ")
  SET(FEDD_HAVE_ACEGENINTERFACE TRUE)
ENDIF()

tribits_configure_file(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS ${HEADERS}
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

SET(HEADERS ${HEADERS} FEDDCore.hpp)

ADD_SUBDIRECTORY(FE)
SET(HEADERS ${HEADERS} ${FE_HEADERS})
SET(SOURCES ${SOURCES} ${FE_SOURCES})

ADD_SUBDIRECTORY(General)
SET(HEADERS ${HEADERS} ${General_HEADERS})
SET(SOURCES ${SOURCES} ${General_SOURCES})

ADD_SUBDIRECTORY(IO)
SET(HEADERS ${HEADERS} ${IO_HEADERS})
SET(SOURCES ${SOURCES} ${IO_SOURCES})

ADD_SUBDIRECTORY(LinearAlgebra)
SET(HEADERS ${HEADERS} ${LinearAlgebra_HEADERS})
SET(SOURCES ${SOURCES} ${LinearAlgebra_SOURCES})

ADD_SUBDIRECTORY(Mesh)
SET(HEADERS ${HEADERS} ${Mesh_HEADERS})
SET(SOURCES ${SOURCES} ${Mesh_SOURCES})

ADD_SUBDIRECTORY(AceFemAssembly)
SET(HEADERS ${HEADERS} ${AceFemAssembly_HEADERS})
SET(SOURCES ${SOURCES} ${AceFemAssembly_SOURCES})

ADD_SUBDIRECTORY(AceFemAssembly/specific)
SET(HEADERS ${HEADERS} ${AceFemAssemblySpecific_HEADERS})
SET(SOURCES ${SOURCES} ${AceFemAssemblySpecific_SOURCES})

ADD_SUBDIRECTORY(AceFemAssembly/AceInterface)
SET(HEADERS ${HEADERS} ${AceFemAssemblyAceInterface_HEADERS})
SET(SOURCES ${SOURCES} ${AceFemAssemblyAceInterface_SOURCES})

#ADD_SUBDIRECTORY(Utils)
#SET(HEADERS ${HEADERS} ${Utils_HEADERS})
#SET(SOURCES ${SOURCES} ${Utils_SOURCES})

# This is the list of files excluded from the tarball
tribits_exclude_files()

# Add library target
tribits_add_library(
core
HEADERS ${HEADERS}
SOURCES ${SOURCES}
)

tribits_add_test_directories(FE/tests)
tribits_add_test_directories(General/tests)
tribits_add_test_directories(IO/tests)
tribits_add_test_directories(Mesh/tests)
tribits_add_test_directories(LinearAlgebra/tests)
tribits_add_test_directories(AceFemAssembly/tests)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# Do standard postprocessing
#

tribits_package_postprocess()

