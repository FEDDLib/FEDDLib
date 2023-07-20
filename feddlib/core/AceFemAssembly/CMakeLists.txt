SET(AceFemAssembly_HEADERS
	AceFemAssembly/AssembleFE_decl.hpp
	AceFemAssembly/AssembleFE_def.hpp
	AceFemAssembly/AssembleFE.hpp
	AceFemAssembly/AssembleFEBlock_decl.hpp
	AceFemAssembly/AssembleFEBlock_def.hpp
	AceFemAssembly/AssembleFEBlock.hpp
	AceFemAssembly/Helper.hpp
	AceFemAssembly/AssembleFEFactory_decl.hpp
	AceFemAssembly/AssembleFEFactory_def.hpp
	AceFemAssembly/AssembleFEFactory.hpp
CACHE INTERNAL "")

SET(AceFemAssembly_SOURCES
	AceFemAssembly/AssembleFE.cpp
	AceFemAssembly/AssembleFEBlock.cpp
	AceFemAssembly/Helper.cpp
	AceFemAssembly/AssembleFEFactory.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${AceFemAssembly_HEADERS})
  STRING(REPLACE "AceFemAssembly/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/AceFemAssembly"
  COMPONENT ${PACKAGE_NAME}
)


