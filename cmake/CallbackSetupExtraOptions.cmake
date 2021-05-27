
# We need to inject the XLib/cmake directory to find
# XLibCreateClientTemplateHeaders.cmake
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${XLib_SOURCE_DIR}/cmake")

MACRO(XLIB_DISABLE_PACKAGE_REQUIRING_CXX11  CXX11_PACKAGE_NAME_IN)
  IF ("${${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}}" STREQUAL "")
    MESSAGE(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}=OFF"
      " because ${PROJECT_NAME}_ENABLE_CXX11='${${PROJECT_NAME}_ENABLE_CXX11}'!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN} OFF)
  ELSEIF (${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN})
    MESSAGE( FATAL_ERROR
      "ERROR: Setting"
      " ${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}='${${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}}'"
      " is not consistent with "
      " ${PROJECT_NAME}_ENABLE_CXX11='${${PROJECT_NAME}_ENABLE_CXX11}'!"
      " ${CXX11_PACKAGE_NAME_IN} requires C++11 support!  Either don't"
      " enable the package ${CXX11_PACKAGE_NAME_IN} or enable support for C++11!")
  ELSE()
    # This package is already disabled which is just fine.
  ENDIF()
ENDMACRO()


MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)

  #MESSAGE("TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS got called!")

  SET(TPL_ENABLE_MPI OFF CACHE BOOL "Enable MPI support.")

  #
  # Set options for global enable/disable of float and complex
  #

  SET(XLib_ENABLE_FLOAT  OFF  CACHE  BOOL
    "Enable the float scalar type in all XLib packages by default.")

  SET(XLib_ENABLE_COMPLEX  OFF  CACHE  BOOL
    "Enable std::complex<T> scalar types in all XLib packages by default.")

  IF (XLib_ENABLE_COMPLEX  AND  Trilinos_ENABLE_FLOAT)
    SET(XLib_ENABLE_COMPLEX_FLOAT_DEFAULT  ON)
  ELSE()
    SET(XLib_ENABLE_COMPLEX_FLOAT_DEFAULT  OFF)
  ENDIF()
  SET(XLib_ENABLE_COMPLEX_FLOAT  ${Trilinos_ENABLE_COMPLEX_FLOAT_DEFAULT}
    CACHE  BOOL
    "Enable std::complex<float> scalar types in all XLib packages by default.")

  SET(XLib_ENABLE_COMPLEX_DOUBLE  ${Trilinos_ENABLE_COMPLEX}
    CACHE  BOOL
    "Enable std::complex<double> scalar types in all XLib packages by default.")

  OPTION(XLib_ENABLE_THREAD_SAFE
    "Enable thread safe code including RCP classes." OFF )

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
  IF (XLib_ENABLE_THREAD_SAFE AND NOT ${PROJECT_NAME}_ENABLE_CXX11)
    MESSAGE(FATAL_ERROR
      "You set XLib_ENABLE_THREAD_SAFE=ON, but ${PROJECT_NAME}' support"
      " for CXX11 is not enabled (${PROJECT_NAME}_ENABLE_CXX11=OFF)."
      "  This is not allowed.  Please enable ${PROJECT_NAME}_ENABLE_CXX11 in"
      " ${PROJECT_NAME} before attempting to enable XLib_ENABLE_THREAD_SAFE"
      " or leave XLib_ENABLE_THREAD_SAFE off.")
  ENDIF ()

  #
  # XLib Data Dir?  Is this still being used anywhere?
  #

  ADVANCED_SET(XLib_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path XLibData directory to find more tests and other stuff" )

  #
  # Put in disables based on various criteria
  #
    
  IF (NOT ${PROJECT_NAME}_ENABLE_CXX11)
    XLIB_DISABLE_PACKAGE_REQUIRING_CXX11("Kokkos")
    XLIB_DISABLE_PACKAGE_REQUIRING_CXX11("Tpetra")
  ENDIF()

  # Used by some XLib packages?
  SET(TRILINOS_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

ENDMACRO()
