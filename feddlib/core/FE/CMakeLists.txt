SET(FE_HEADERS
	FE/Domain_decl.hpp
	FE/Domain_def.hpp
	FE/Domain.hpp
	FE/EdgeElements.hpp
    FE/EntitiesOfElements.hpp
	FE/Elements.hpp
	FE/FE.hpp
	FE/FE_decl.hpp
	FE/FE_def.hpp
	FE/FiniteElement.hpp
    FE/sms.hpp
CACHE INTERNAL "")

SET(FE_SOURCES
	FE/Domain.cpp
    FE/EdgeElements.cpp
    FE/Elements.cpp
    FE/EntitiesOfElements.cpp
	FE/FE.cpp
    FE/FiniteElement.cpp    
    FE/SMSUtility.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${FE_HEADERS})
  STRING(REPLACE "FE/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/FE"
  COMPONENT ${PACKAGE_NAME}
)
