set(FE_HEADERS
	FE/Domain_decl.hpp
	FE/Domain_def.hpp
	FE/Domain.hpp
	FE/EdgeElements.hpp
    FE/EntitiesOfElements.hpp
	FE/Elements.hpp
	FE/FE.hpp
	FE/FE_decl.hpp
	FE/FE_def.hpp
	FE/FE_ElementAssembly.hpp
	FE/FE_ElementAssembly_decl.hpp
	FE/FE_ElementAssembly_def.hpp
	FE/FiniteElement.hpp
	FE/TriangleElements.hpp
    FE/sms.hpp
	FE/Helper.hpp
CACHE INTERNAL "")

set(FE_SOURCES
	FE/Domain.cpp
    FE/EdgeElements.cpp
    FE/Elements.cpp
    FE/EntitiesOfElements.cpp
	FE/FE.cpp
    FE/FE_ElementAssembly.cpp
    FE/FiniteElement.cpp  
	FE/TriangleElements.cpp  
    FE/SMSUtility.cpp
	FE/Helper.cpp
CACHE INTERNAL "")


set(LOCAL_HEADERS)
foreach(INC ${FE_HEADERS})
  string(REPLACE "FE/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
