#
# A) Define the project
#

# To be safe, define your minimum CMake version.  This may be newer than the
# min required by TriBITS.

cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

include(${CMAKE_SOURCE_DIR}/ProjectName.cmake)
project(${PROJECT_NAME} NONE)

# The version number.
set(FEDDLIB_VERSION_MAJOR 1)
set(FEDDLIB_VERSION_MINOR 0)

set(${PROJECT_NAME}_ENABLE_Fortran ON CACHE BOOL "Enable Fortran compiler")

# Set C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#
# B) Pull in the TriBITS system and execute
#

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits/TriBITS.cmake)

# CMake requires this be in the top file and not in an include file :-(
cmake_minimum_required(VERSION ${TRIBITS_CMAKE_MINIMUM_REQUIRED})

# Add source dir to includes
tribits_include_directories(${CMAKE_SOURCE_DIR})

# Do all of the processing for this Tribits project
tribits_project()


foreach(ENABLED_PACKAGE ${FEDDLib_ENABLED_PACKAGES})
	string(TOUPPER ${ENABLED_PACKAGE} UPPER_PACKAGE)
	set(FEDD_HAVE_${UPPER_PACKAGE} ON)
endforeach()

configure_file(${CMAKE_SOURCE_DIR}/FEDDLib_Config.h.in
               ${CMAKE_BINARY_DIR}/FEDDLib_Config.h)

# Install the generated configuration header
install(FILES ${CMAKE_BINARY_DIR}/FEDDLib_Config.h DESTINATION include)

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/feddlib/debug.areas
DESTINATION share/feddlib)
