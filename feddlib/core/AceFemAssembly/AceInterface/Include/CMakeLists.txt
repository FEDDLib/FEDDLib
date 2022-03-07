SET(AceFemAssemblyAceInterfaceInclude_HEADERS
	AceFemAssembly/AceInterface/Include/NeoHookQuadraticTets.hpp
CACHE INTERNAL "")

SET(AceFemAssemblyTestFE_SOURCES
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${AceFemAssemblyAceInterfaceInclude_HEADERS})
  STRING(REPLACE "AceInterface/Include/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/AceFemAssembly/AceInterface/Include"
  COMPONENT ${PACKAGE_NAME}
)
