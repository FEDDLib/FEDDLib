set(REQUIRED_TPLS MPI Trilinos HDF5 Z)

# The dl library is part of Apples core functions. This is also the case for, e.g., Ubuntu, where the functions are in libc. Only if libdl needs to be explicitly linked, it should be included here. We need it for HDF5 to avoid having the user pass "-ldl -lz" as a compiler flag.
enable_language(C)
include(CheckFunctionExists)
check_function_exists(dlopen DLOPEN_IN_LIBC)
if (UNIX AND NOT APPLE AND NOT DLOPEN_IN_LIBC)
  list(APPEND REQUIRED_TPLS DL)
endif()

tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES
  LIB_OPTIONAL_PACKAGES
  TEST_REQUIRED_PACKAGES
  TEST_OPTIONAL_PACKAGES
  LIB_REQUIRED_TPLS ${REQUIRED_TPLS}
  LIB_OPTIONAL_TPLS AceGENInterface
  TEST_REQUIRED_TPLS
  TEST_OPTIONAL_TPLS
)
