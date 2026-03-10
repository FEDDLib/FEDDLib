#
# A) Define the package
# This is the TriBITS way of defining a package. It sets up all the necessary
# variables and infrastructure for the 'core' package.
#
tribits_package(core)
# Set the package version - this is used for packaging and versioning
set(${PACKAGE_NAME}_VERSION "1.0.0")

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
set(HEADERS "")
set(SOURCES "")

# Add the generated configuration header to the headers list
append_set(HEADERS
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

# Add the main package header
append_set(HEADERS FEDDCore.hpp)

# Add headers/sources from all subfolders
add_subdirectory(FE)
append_set(HEADERS ${FE_HEADERS})
append_set(SOURCES ${FE_SOURCES})

add_subdirectory(General)
append_set(HEADERS ${General_HEADERS})
append_set(SOURCES ${General_SOURCES})

add_subdirectory(IO)
append_set(HEADERS ${IO_HEADERS})
append_set(SOURCES ${IO_SOURCES})

add_subdirectory(LinearAlgebra)
append_set(HEADERS ${LinearAlgebra_HEADERS})
append_set(SOURCES ${LinearAlgebra_SOURCES})

add_subdirectory(Mesh)
append_set(HEADERS ${Mesh_HEADERS})
append_set(SOURCES ${Mesh_SOURCES})

add_subdirectory(AceFemAssembly)
append_set(HEADERS ${AceFemAssembly_HEADERS})
append_set(SOURCES ${AceFemAssembly_SOURCES})

add_subdirectory(AceFemAssembly/specific)
append_set(HEADERS ${AceFemAssemblySpecific_HEADERS})
append_set(SOURCES ${AceFemAssemblySpecific_SOURCES})

add_subdirectory(AceFemAssembly/AceInterface)
append_set(HEADERS ${AceFemAssemblyAceInterface_HEADERS})
append_set(SOURCES ${AceFemAssemblyAceInterface_SOURCES})

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

# Link Fortran runtime libraries
# The core library contains Fortran source files (from AceFemAssembly/AceInterface),
# so we need to ensure the Fortran runtime libraries are linked and propagated to
# consumers of the library. Using PUBLIC visibility ensures that when consumers
# link against libcore.a (or use imported targets), they automatically get the
# Fortran runtime dependencies.
if(CMAKE_Fortran_COMPILER)
    # CMake provides implicit link libraries detected from the Fortran compiler.
    # This includes libgfortran for GCC, libifcore for Intel, etc.
    # These are the libraries that the Fortran compiler implicitly links but
    # C++ compilers don't know about.
    if(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES)
        target_link_libraries(core
            PUBLIC
                ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
        )
    endif()

    # Also link standard Fortran libraries if available
    if(CMAKE_Fortran_STANDARD_LIBRARIES)
        target_link_libraries(core
            PUBLIC
                ${CMAKE_Fortran_STANDARD_LIBRARIES}
        )
    endif()
endif()

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
