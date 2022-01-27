SET(AceFemAssemblyTestFE_HEADERS
	AceFemAssembly/TestFE/FE_Test_decl.hpp
	AceFemAssembly/TestFE/FE_Test_def.hpp
	AceFemAssembly/TestFE/FE_Test.hpp
CACHE INTERNAL "")

SET(AceFemAssemblyTestFE_SOURCES
	AceFemAssembly/TestFE/FE_Test.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${AceFemAssemblyTestFE_HEADERS})
  STRING(REPLACE "TestFE/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/AceFemAssembly/TestFE"
  COMPONENT ${PACKAGE_NAME}
)
