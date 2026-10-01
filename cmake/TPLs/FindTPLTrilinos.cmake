include(TribitsTplDeclareLibraries)

if(Trilinos_LIBRARY_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} "${Trilinos_LIBRARY_DIRS}/cmake/Trilinos")
endif()
if (Trilinos_INCLUDE_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} ${Trilinos_INCLUDE_DIRS})
endif()

message("Trilinos_DIR: ${Trilinos_DIR}")

# Find TrilinosConfig.cmake and import it.
find_package (Trilinos NO_MODULE HINTS ${Trilinos_DIR})

# Stop CMake if Trilinos was not found.
if (NOT Trilinos_FOUND)
  message (FATAL_ERROR "Could not find Trilinos!")
endif ()

# if("${Trilinos_VERSION_MAJOR}" GREATER 10)
#   set (HAVE_TRILINOS_GT_10_6 TRUE)
#   message (STATUS "Using Trilinos > 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
# else()
#  if("${Trilinos_VERSION_MAJOR}" SMALLER 10)
#    message (STATUS "Using Trilinos <= 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
#  else()
#    if("${Trilinos_VERSION_MINOR}" GREATER 6)
#     set (HAVE_TRILINOS_GT_10_6 TRUE)
#     message (STATUS "Using Trilinos > 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
#    else()
#     message (STATUS "Using Trilinos <= 10.6 : " ${Trilinos_VERSION_MAJOR} "." ${Trilinos_VERSION_MINOR})
#    endif()
#  endif()
# endif ()

# Here it will be better just to raise a warning or have if(USE_TRILINOS_COMPILERS)
# Make sure to use same compilers and flags as Trilinos
if(NOT ${CMAKE_CXX_COMPILER} STREQUAL ${Trilinos_CXX_COMPILER})
  message(STATUS "the selected compiler differs from Trilinos CXX compiler")
  message(STATUS "CMAKE_CXX_COMPILER:    " ${CMAKE_CXX_COMPILER})
  message(STATUS "Trilinos_CXX_COMPILER: " ${Trilinos_CXX_COMPILER})
endif()
if(NOT ${CMAKE_C_COMPILER} STREQUAL ${Trilinos_C_COMPILER})
  message(STATUS "the selected compiler differs from Trilinos C compiler")
  message(STATUS "CMAKE_C_COMPILER:    " ${CMAKE_C_COMPILER})
  message(STATUS "Trilinos_C_COMPILER: " ${Trilinos_C_COMPILER})
endif()
if(NOT ${CMAKE_Fortran_COMPILER} STREQUAL ${Trilinos_Fortran_COMPILER})
  message(STATUS "the selected compiler differs from Trilinos Fortran compiler")
  message(STATUS "CMAKE_Fortran_COMPILER:    " ${CMAKE_C_COMPILER})
  message(STATUS "Trilinos_Fortran_COMPILER: " ${Trilinos_C_COMPILER})
endif()

# Optional Packages (to be moved outside with COMPONENTS ...)
list (APPEND FEDDLib_OPTIONAL_Trilinos_PKGS
  "NOX" "Teko" "Stratimikos" "ShyLU" "Zoltan" "Zoltan2" "MueLu" "Ifpack2" "Anasazi")

# Required packages (to be moved outside, like REQUIRED COMPONENTS ...)
list (APPEND FEDDLib_REQUIRED_Trilinos_PKGS
  "Amesos2" "Belos" "ShyLU_DDFROSch" "Stratimikos" "Teuchos" "Thyra" "Tpetra" "Xpetra")

# Start scanning Trilinos configuration
set (Trilinos_MISSING_REQUIRED_PACKAGE OFF) # Will be set to ON if a required package is missing
foreach (TYPE IN ITEMS "OPTIONAL" "REQUIRED")
  foreach (PKG IN LISTS FEDDLib_${TYPE}_Trilinos_PKGS)
    # Look for PKG
    list (FIND Trilinos_PACKAGE_LIST "${PKG}" PKG_FOUND)
    if (PKG_FOUND GREATER -1)
      if (TYPE STREQUAL "REQUIRED")
        message (STATUS "Trilinos :: required package     found:   ${PKG}")
      else ()
        message (STATUS "Trilinos :: optional package     found:   ${PKG}")
      endif ()
      list (APPEND FEDDLib_Trilinos_LIBRARIES "${${PKG}_LIBRARIES}")
      list (APPEND FEDDLib_Trilinos_TPL_INCLUDE_DIRS "${${PKG}_TPL_INCLUDE_DIRS}")
      list (APPEND FEDDLib_Trilinos_TPL_LIST "${${PKG}_TPL_LIST}")
      string (TOUPPER ${PKG} UPKG)
      set (${UPKG}_FOUND True)
      set (HAVE_TRILINOS_${UPKG} True)
    else ()
      if (TYPE STREQUAL "REQUIRED")
        message (SEND_ERROR "Trilinos :: required package NOT found:   ${PKG}")
        set (Trilinos_MISSING_REQUIRED_PACKAGE ON)
      else ()
        # Only give a status information about packages that are not available.
        # Replace STATUS with WARNING if you want something stronger.
        message (STATUS "Trilinos :: optional package NOT found:   ${PKG}   (Some tests might not compile properly.)")
      endif ()
    endif ()
  endforeach (PKG)
endforeach (TYPE)

# If there are any required packages missing: abort.
if (Trilinos_MISSING_REQUIRED_PACKAGE)
  message (FATAL_ERROR "Trilinos :: Required package(s) not found. See above for missing packages.")
endif ()

# Cleaning duplicates
list (REVERSE FEDDLib_Trilinos_TPL_LIST)
list (REMOVE_DUPLICATES FEDDLib_Trilinos_TPL_LIST)
list (REVERSE FEDDLib_Trilinos_TPL_LIST)
list (REVERSE FEDDLib_Trilinos_LIBRARIES)
list (REMOVE_DUPLICATES FEDDLib_Trilinos_LIBRARIES)
list (REVERSE FEDDLib_Trilinos_LIBRARIES)
list (REVERSE Trilinos_TPL_LIBRARIES)
list (REMOVE_DUPLICATES Trilinos_TPL_LIBRARIES)
list (REVERSE Trilinos_TPL_LIBRARIES)
set (FEDDLib_Trilinos_TPL_LIBRARIES ${Trilinos_TPL_LIBRARIES})

list (REMOVE_DUPLICATES FEDDLib_Trilinos_TPL_INCLUDE_DIRS)

list (APPEND FEDDLib_Trilinos_INCLUDE_DIRS
  ${Trilinos_INCLUDE_DIRS}
  ${FEDDLib_Trilinos_TPL_INCLUDE_DIRS})
# I think there's a better way to handle this ... CMake
# should take care of -L or -l or -rpath ...
set (FEDDLib_Trilinos_LIBS "-L${Trilinos_LIBRARY_DIRS}")
foreach (LIB IN LISTS FEDDLib_Trilinos_LIBRARIES)
  set (FEDDLib_Trilinos_LIBS "${FEDDLib_Trilinos_LIBS} -l${LIB}")
endforeach (LIB)
set (FEDDLib_Trilinos_LIBS ${FEDDLib_Trilinos_LIBS} ${FEDDLib_Trilinos_TPL_LIBRARIES})

# Detect available Trilinos packages and set FEDD_HAVE_* flags
# This runs after Trilinos is found, so Trilinos_PACKAGE_LIST is populated
foreach(TRILINOS_PACKAGE_NAME in ${Trilinos_PACKAGE_LIST})
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "Ifpack2")
        set(FEDD_HAVE_IFPACK2 TRUE)
        message(STATUS "FEDDLib:   FEDD_HAVE_IFPACK2   TRUE")
    endif()
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "NOX")
        set(FEDD_HAVE_NOX TRUE)
        message(STATUS "FEDDLib:   FEDD_HAVE_NOX       TRUE")
    endif()
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "Teko")
        set(FEDD_HAVE_TEKO TRUE)
        message(STATUS "FEDDLib:   FEDD_HAVE_TEKO      TRUE")
    endif()
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "Zoltan2")
        set(FEDD_HAVE_ZOLTAN2 TRUE)
        message(STATUS "FEDDLib:   FEDD_HAVE_ZOLTAN2   TRUE")
    endif()
endforeach()

# Uncomment the following to see all variables in the current scope.
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

# TPLs
set(TRILINOS_ENABLED_TPLS_LIST "")

# Define a function to check TPL status across packages.
# This will not differentiate between TPLs enabled only for individual packages, e.g., Amesos2_ENABLE_METIS.
# Important: It will check all variables for the searched TPL name, so it is not entirely robust.
function(check_tpl_status TPL_NAME)
    # The actual variable names depend on the internal Trilinos setup, 
    # but the pattern should be <PackageName>_ENABLE_<TPL_NAME>.

    # We find any variable that contains the searched name.
    get_cmake_property(_all_vars VARIABLES)
    foreach(_var ${_all_vars})
        if("${_var}" MATCHES "([A-Za-z0-9]+_ENABLE_${TPL_NAME})$")
            if(${_var}) # Check if the variable's value is TRUE/ON
                list(APPEND TRILINOS_ENABLED_TPLS_LIST "${TPL_NAME}")
                set(TRILINOS_ENABLED_TPLS_LIST ${TRILINOS_ENABLED_TPLS_LIST} PARENT_SCOPE)
                set(FOUND_TPL TRUE PARENT_SCOPE)
                return() # Found on occurence of the TPL name, so return.
            endif()
        endif()
    endforeach()
endfunction()

message(STATUS "Check all variables in scope for containing TPL names.")
foreach (TPL IN ITEMS "Boost" "BLAS" "LAPACK" "METIS" "ParMETIS" "KLU2" "UMFPACK" "SuperLU" "SuperLUDist" "MUMPS" "PARDISO_MKL" "HDF5")
  set(FOUND_TPL FALSE)
  check_tpl_status(${TPL})
  if (FOUND_TPL)
    string (TOUPPER ${TPL} UTPL)
    set (${UTPL}_IS_IN_TRILINOS True)
    message (STATUS "Trilinos :: TPL     found:   ${TPL}")
  else ()
    message (STATUS "Trilinos :: TPL NOT found:   ${TPL}")
  endif ()
endforeach (TPL)

# TRILINOS_ENABLED_TPLS_LIST contains all found TPLs.
message(STATUS "Summary of enabled TPLs: ${TRILINOS_ENABLED_TPLS_LIST}")

# Filling variables needed by the TriBITS system
# set (TPL_Trilinos_INCLUDE_DIRS ${FEDDLib_Trilinos_INCLUDE_DIRS})
# set (TPL_Trilinos_LIBRARY_DIRS Trilinos::all_selected_libs) 
# set (TPL_Trilinos_LIBRARIES Trilinos::all_selected_libs)

# This generates the TrilinosConfig.cmake file for installation together with FEDDLib libs and headers.
# Since Trilinos types are exposed direclty in the FEDDLib, Trilinos headers are a compile time dependency.
# It is possible that the generated TrilinosConfig.cmake file is not complete, as this use-case has not been tested to date [27.11.25]. 
# In this case the install location of Trilinos would need to be passed to the compiler using -I/path/to/Trilinos/install 
# or the commented code above needs to be fixed.
tribits_tpl_find_include_dirs_and_libraries(Trilinos)
message(STATUS "End of processing Trilinos installation.\n")
