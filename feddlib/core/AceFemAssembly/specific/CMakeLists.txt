SET(AceFemAssemblySpecific_HEADERS
	AceFemAssembly/specific/AssembleFEAceLaplace_decl.hpp
	AceFemAssembly/specific/AssembleFEAceLaplace_def.hpp
	AceFemAssembly/specific/AssembleFEAceLaplace.hpp
	AceFemAssembly/specific/AssembleFEAceNavierStokes_decl.hpp
	AceFemAssembly/specific/AssembleFEAceNavierStokes_def.hpp
	AceFemAssembly/specific/AssembleFEAceNavierStokes.hpp
CACHE INTERNAL "")

SET(AceFemAssemblySpecific_SOURCES
	AceFemAssembly/specific/AssembleFEAceLaplace.cpp
	AceFemAssembly/specific/AssembleFEAceNavierStokes.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${AceFemAssemblySpecific_HEADERS})
  STRING(REPLACE "specific/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/AceFemAssembly/specific"
  COMPONENT ${PACKAGE_NAME}
)
