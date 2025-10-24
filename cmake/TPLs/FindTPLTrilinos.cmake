include(TribitsTplDeclareLibraries)

# TRIBITS_TPL_DECLARE_LIBRARIES( Trilinos
#   REQUIRED_HEADERS Epetra_Comm.h
#   REQUIRED_LIBS_NAMES "epetra"
#   )

if(Trilinos_LIBRARY_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} "${Trilinos_LIBRARY_DIRS}/cmake/Trilinos")
endif()
if (Trilinos_INCLUDE_DIRS)
    set (Trilinos_DIR ${Trilinos_DIR} ${Trilinos_INCLUDE_DIRS})
endif()

message("Trilinos_DIR: ${Trilinos_DIR}")

# Here I am looking for TrilinosConfig.cmake and I will import it.
find_package (Trilinos NO_MODULE HINTS ${Trilinos_DIR})

# Stop cmake if Trilinos is not found.
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
  "NOX" "Thyra" "Rythmos" "Teko" "Stratimikos" "Isorropia" "ShyLU" "Zoltan2" "MueLu")

# Required packages (to be moved outside, like REQUIRED COMPONENTS ...)
list (APPEND FEDDLib_REQUIRED_Trilinos_PKGS
  "Belos" "Epetra" "EpetraExt" "ShyLU_DDFROSch" "Stratimikos" "Teko" "Teuchos" "Thyra" "Tpetra" "Xpetra")

# Start scanning Trilinos configuration
foreach (TYPE IN ITEMS "OPTIONAL" "REQUIRED")
  foreach (PKG IN LISTS FEDDLib_${TYPE}_Trilinos_PKGS)
    # Look for PKG
    list (FIND Trilinos_PACKAGE_LIST "${PKG}" PKG_FOUND)
    if (PKG_FOUND GREATER -1)
      # Found! Let's announce it!
      message (STATUS "Trilinos :: ${PKG} Found!")
      list (APPEND FEDDLib_Trilinos_LIBRARIES "${${PKG}_LIBRARIES}")
      list (APPEND FEDDLib_Trilinos_TPL_INCLUDE_DIRS "${${PKG}_TPL_INCLUDE_DIRS}")
      list (APPEND FEDDLib_Trilinos_TPL_LIST "${${PKG}_TPL_LIST}")
      string (TOUPPER ${PKG} UPKG)
      set (${UPKG}_FOUND True)
      set (HAVE_TRILINOS_${UPKG} True)
    else ()
      if (TYPE STREQUAL "REQUIRED")
        message (FATAL_ERROR "Trilinos :: ${PKG} NOT Found!")
      else ()
        message (WARNING "Trilinos :: ${PKG} NOT Found! Some test might not compile properly ...")
      endif ()
    endif ()
  endforeach (PKG)
endforeach (TYPE)

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

# TPLs
foreach (TPL IN ITEMS "ParMETIS" "Boost" "LAPACK" "BLAS" "UMFPACK" "SuperLU" "SuperLUDist" "HDF5")
    list (FIND FEDDLib_Trilinos_TPL_LIST ${TPL} TPL_FOUND)
  if (TPL_FOUND GREATER -1)
    string (TOUPPER ${TPL} UTPL)
    set (${UTPL}_IS_IN_TRILINOS True)
  endif()
endforeach (TPL)

# Filling variables needed by the TriBITS system
set (TPL_Trilinos_INCLUDE_DIRS ${FEDDLib_Trilinos_INCLUDE_DIRS})
set (TPL_Trilinos_LIBRARY_DIRS Trilinos::all_selected_libs) 
set (TPL_Trilinos_LIBRARIES Trilinos::all_selected_libs)

# Detect available Trilinos packages and set FEDD_HAVE_* flags
# This runs after Trilinos is found, so Trilinos_PACKAGE_LIST is populated
foreach(TRILINOS_PACKAGE_NAME in ${Trilinos_PACKAGE_LIST})
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "Ifpack2")
        set(FEDD_HAVE_IFPACK2 TRUE)
        message(STATUS "FEDDLib: Found Trilinos Ifpack2 package")
    endif()
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "NOX")
        set(FEDD_HAVE_NOX TRUE)
        message(STATUS "FEDDLib: Found Trilinos NOX package")
    endif()
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "Teko")
        set(FEDD_HAVE_TEKO TRUE)
        message(STATUS "FEDDLib: Found Trilinos Teko package")
    endif()
    if(${TRILINOS_PACKAGE_NAME} STREQUAL "Zoltan2")
        set(FEDD_HAVE_ZOLTAN2 TRUE)
        message(STATUS "FEDDLib: Found Trilinos Zoltan2 package")
    endif()
endforeach()

# Summary of detected features
message(STATUS "FEDDLib: Feature detection summary:")
message(STATUS "  - FEDD_HAVE_IFPACK2: ${FEDD_HAVE_IFPACK2}")
message(STATUS "  - FEDD_HAVE_NOX: ${FEDD_HAVE_NOX}")
message(STATUS "  - FEDD_HAVE_TEKO: ${FEDD_HAVE_TEKO}")



