#
# A) Define the package
# This is the TriBITS way of defining a package. It sets up all the necessary
# variables and infrastructure for the 'core' package.
#
tribits_package(core)
# Set the package version - this is used for packaging and versioning
SET(${PACKAGE_NAME}_VERSION "1.0.0")

#
# B) Set up package-specific options
# Here we check for specific Trilinos packages and set corresponding flags
#

tribits_add_debug_option()

tribits_add_show_deprecated_warnings_option()

# TPL detection and feature flags are handled in FindTPLTrilinos.cmake
# The FEDD_HAVE_* variables are defined there and available globally to all packages

# Generate a configuration header file with package-specific defines
# This creates core_config.h with all the FEDD_HAVE_* flags
tribits_configure_file(${PACKAGE_NAME}_config.h)

# Initialize empty lists for headers and sources
# These will be populated by the subdirectories
SET(HEADERS "")
SET(SOURCES "")

# Add the generated configuration header to the headers list
SET(HEADERS ${HEADERS}
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

# Add the main package header
SET(HEADERS ${HEADERS} FEDDCore.hpp)

# Add headers/sources from all subfolders
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

# This is the list of files excluded from the tarball
# TriBITS will exclude these files when creating distribution packages
tribits_exclude_files()

# Add library target
# This creates the actual library target 'core' with all collected headers and sources
# TriBITS will handle linking, installation, and other library properties
tribits_add_library(
core
HEADERS ${HEADERS}
SOURCES ${SOURCES}
)

# Add test directories
# Each subdirectory can have its own tests subdirectory
# TriBITS will automatically discover and build tests in these directories
tribits_add_test_directories(FE/tests)
tribits_add_test_directories(General/tests)
tribits_add_test_directories(IO/tests)
tribits_add_test_directories(Mesh/tests)
tribits_add_test_directories(LinearAlgebra/tests)
tribits_add_test_directories(AceFemAssembly/tests)

#
# Do standard postprocessing
# This performs final TriBITS package setup, dependency resolution, etc.
#
tribits_package_postprocess()
