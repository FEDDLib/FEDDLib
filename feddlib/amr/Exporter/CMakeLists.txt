SET(General_HEADERS
	Exporter/ExporterParaViewAMR.hpp
	Exporter/ExporterParaViewAMR_decl.hpp
	Exporter/ExporterParaViewAMR_def.hpp
CACHE INTERNAL "")

SET(General_SOURCES
    Exporter/ExporterParaViewAMR.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${General_HEADERS})
  STRING(REPLACE "Exporter/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/amr/Exporter"
  COMPONENT ${PACKAGE_NAME}
)
