set(AceFemAssembly_HEADERS
	AceFemAssembly/AssembleFE_decl.hpp
	AceFemAssembly/AssembleFE_def.hpp
	AceFemAssembly/AssembleFE.hpp
	AceFemAssembly/AssembleFEBlock_decl.hpp
	AceFemAssembly/AssembleFEBlock_def.hpp
	AceFemAssembly/AssembleFEBlock.hpp
	AceFemAssembly/AssembleFEFactory_decl.hpp
	AceFemAssembly/AssembleFEFactory_def.hpp
	AceFemAssembly/AssembleFEFactory.hpp
CACHE INTERNAL "")

set(AceFemAssembly_SOURCES
	AceFemAssembly/AssembleFE.cpp
	AceFemAssembly/AssembleFEBlock.cpp
	AceFemAssembly/AssembleFEFactory.cpp
CACHE INTERNAL "")


set(LOCAL_HEADERS)
foreach(INC ${AceFemAssembly_HEADERS})
  string(REPLACE "AceFemAssembly/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
