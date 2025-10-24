#
# A) Define the project
#

# To be safe, define your minimum CMake version.  This may be newer than the
# min required by TriBITS.

CMAKE_MINIMUM_REQUIRED(VERSION 3.23.0 FATAL_ERROR)

INCLUDE(${CMAKE_SOURCE_DIR}/ProjectName.cmake)
PROJECT(${PROJECT_NAME} NONE)

# The version number.
SET(FEDDLIB_VERSION_MAJOR 1)
SET(FEDDLIB_VERSION_MINOR 0)

SET(${PROJECT_NAME}_ENABLE_Fortran ON CACHE BOOL "Enable Fortran compiler")

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#
# B) Pull in the TriBITS system and execute
#

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits/TriBITS.cmake)

# CMake requires this be in the top file and not in an include file :-(
CMAKE_MINIMUM_REQUIRED(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

# Add source dir to includes
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})

# Do all of the processing for this Tribits project
tribits_project()


FOREACH(ENABLED_PACKAGE ${FEDDLib_ENABLED_PACKAGES})
	STRING(TOUPPER ${ENABLED_PACKAGE} UPPER_PACKAGE)
	SET(FEDD_HAVE_${UPPER_PACKAGE} ON)
ENDFOREACH()

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/FEDDLib_Config.h.in
               ${CMAKE_BINARY_DIR}/FEDDLib_Config.h)

# Install the generated configuration header
INSTALL(FILES ${CMAKE_BINARY_DIR}/FEDDLib_Config.h DESTINATION include)

INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/feddlib/debug.areas
DESTINATION share/feddlib)
