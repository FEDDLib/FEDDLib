#
# A) Define the package
#

tribits_package(problems)
set(${PACKAGE_NAME}_VERSION "1.0.0")

#
# B) Set up package-specific options
#

tribits_add_debug_option()

tribits_add_show_deprecated_warnings_option()

tribits_configure_file(${PACKAGE_NAME}_config.h)

set(HEADERS "")
set(SOURCES "")

append_set(HEADERS
${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

add_subdirectory(abstract)
append_set(HEADERS ${problems_abstract_HEADERS})
append_set(SOURCES ${problems_abstract_SOURCES})

add_subdirectory(Solver)
append_set(HEADERS ${Solver_HEADERS})
append_set(SOURCES ${Solver_SOURCES})

add_subdirectory(specific)
append_set(HEADERS ${problems_specific_HEADERS})
append_set(SOURCES ${problems_specific_SOURCES})

# This is the list of files excluded from the tarball
tribits_exclude_files()

# Add library target
tribits_add_library(
    problems
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
)

tribits_add_test_directories(Solver/tests)
tribits_add_test_directories(tests)
tribits_add_example_directories(examples)

#
# Do standard postprocessing
#

tribits_package_postprocess()
