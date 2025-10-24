set(General_HEADERS
	General/BCBuilder.hpp
	General/BCBuilder_decl.hpp
	General/BCBuilder_def.hpp
	General/DefaultTypeDefs.hpp
	General/ExporterParaView.hpp
	General/ExporterParaView_decl.hpp
	General/ExporterParaView_def.hpp
	General/ExporterTxt.hpp
	General/SmallMatrix.hpp
	General/InputToOutputMappingClass.hpp
	General/InputToOutputMappingClass_decl.hpp
	General/InputToOutputMappingClass_def.hpp
	General/DifferentiableFuncClass.hpp
	General/DifferentiableFuncClass_decl.hpp
	General/DifferentiableFuncClass_def.hpp	
	General/HDF5Import.hpp
	General/HDF5Import_decl.hpp
	General/HDF5Import_def.hpp
	General/HDF5Export.hpp
	General/HDF5Export_decl.hpp
	General/HDF5Export_def.hpp
	General/HDF5Import.cpp
	General/HDF5Export.cpp
CACHE INTERNAL "")

set(General_SOURCES
	General/BCBuilder.cpp
    General/ExporterParaView.cpp
    General/ExporterTxt.cpp
	General/InputToOutputMappingClass.cpp
	General/DifferentiableFuncClass.cpp
CACHE INTERNAL "")


set(LOCAL_HEADERS)
foreach(INC ${General_HEADERS})
  string(REPLACE "General/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()
