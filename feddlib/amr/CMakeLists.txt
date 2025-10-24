#
# A) Define the package
#

tribits_package(amr)
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

append_set(HEADERS
	AdaptiveMeshRefinement.hpp
	AdaptiveMeshRefinement_decl.hpp
	AdaptiveMeshRefinement_def.hpp
	RefinementFactory.hpp
	RefinementFactory_decl.hpp
	RefinementFactory_def.hpp
	ErrorEstimation.hpp
	ErrorEstimation_decl.hpp
	ErrorEstimation_def.hpp
	ExporterParaViewAMR.hpp
	ExporterParaViewAMR_decl.hpp
	ExporterParaViewAMR_def.hpp
)

append_set(SOURCES
	AdaptiveMeshRefinement.cpp
	RefinementFactory.cpp
	ErrorEstimation.cpp
	ExporterParaViewAMR.cpp
)

# This is the list of files excluded from the tarball
tribits_exclude_files()

# Add library target
tribits_add_library(
amr
HEADERS ${HEADERS}
SOURCES ${SOURCES}
)

tribits_add_test_directories(tests)

#
# Do standard postprocessing
#

tribits_package_postprocess()
