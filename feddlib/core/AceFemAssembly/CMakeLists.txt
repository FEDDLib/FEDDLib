SET(AceFemAssembly_HEADERS
	AceFemAssembly/AssembleFE_decl.hpp
	AceFemAssembly/AssembleFE_def.hpp
	AceFemAssembly/AssembleFE.hpp
#	AceFemAssembly/specific/AssembleFEAceLaplace_decl.hpp
#	AceFemAssembly/specific/AssembleFEAceLaplace_def.hpp
#	AceFemAssembly/specific/AssembleFEAceLaplace.hpp
	AceFemAssembly/AssembleFEFactory_decl.hpp
	AceFemAssembly/AssembleFEFactory_def.hpp
	AceFemAssembly/AssembleFEFactory.hpp
#	AceFemAssembly/AssembleFEGeneral_decl.hpp
#	AceFemAssembly/AssembleFEGeneral_def.hpp
#	AceFemAssembly/AssembleFEGeneral.hpp
#	AceFemAssembly/tests/FE_Test_decl.hpp
#	AceFemAssembly/tests/FE_Test_def.hpp
#	AceFemAssembly/tests/FE_Test.hpp
CACHE INTERNAL "")

SET(AceFemAssembly_SOURCES
	AceFemAssembly/AssembleFE.cpp
#	AceFemAssembly/specific/AssembleFEAceLaplace.cpp
	AceFemAssembly/AssembleFEFactory.cpp
#	AceFemAssembly/tests/FE_Test.cpp
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
