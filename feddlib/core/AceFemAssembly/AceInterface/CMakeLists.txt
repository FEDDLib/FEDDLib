set(AceFemAssemblyAceInterface_HEADERS
	AceFemAssembly/AceInterface/NeoHookQuadraticTets.hpp
  AceFemAssembly/AceInterface/NeoHookQuadraticTets2.hpp
    AceFemAssembly/AceInterface/NeoHookQuadraticTets3.hpp
   AceFemAssembly/AceInterface/DeformationDiffusionNeoHook.hpp
  AceFemAssembly/AceInterface/sms.h
  AceFemAssembly/AceInterface/bdata.h
  AceFemAssembly/AceInterface/pdata6.h
  AceFemAssembly/AceInterface/cdata.h
  AceFemAssembly/AceInterface/eldata.h
  AceFemAssembly/AceInterface/iofile.h
  AceFemAssembly/AceInterface/rdata.h
  AceFemAssembly/AceInterface/strnum.h
  AceFemAssembly/AceInterface/hdata.h
  AceFemAssembly/AceInterface/prstrs.h
  AceFemAssembly/AceInterface/comblk.h
  AceFemAssembly/AceInterface/sdata.h
  AceFemAssembly/AceInterface/pointer.h
  AceFemAssembly/AceInterface/counts.h
  AceFemAssembly/AceInterface/tdata.h
  AceFemAssembly/AceInterface/sensitivity.h
CACHE INTERNAL "")

set(AceFemAssemblyAceInterface_SOURCES
	AceFemAssembly/AceInterface/NeoHookQuadraticTets.f
  AceFemAssembly/AceInterface/NeoHookQuadraticTets_cbinding.f
  AceFemAssembly/AceInterface/NeoHookQuadraticTets2.f
  AceFemAssembly/AceInterface/NeoHookQuadraticTets2_cbinding.f
    AceFemAssembly/AceInterface/NeoHookQuadraticTets3.f
  AceFemAssembly/AceInterface/NeoHookQuadraticTets3_cbinding.f
  AceFemAssembly/AceInterface/DeformationDiffusionNeoHook.f
  AceFemAssembly/AceInterface/DeformationDiffusionNeoHook_cbinding.f
	AceFemAssembly/AceInterface/SMSUtility.f
CACHE INTERNAL "")

set(LOCAL_HEADERS)
foreach(INC ${AceFemAssemblyAceInterface_HEADERS})
  string(REPLACE "AceInterface/" "" LOCAL_INC ${INC})
  append_set(LOCAL_HEADERS ${LOCAL_INC})
endforeach()

install(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/AceFemAssembly/AceInterface"
  COMPONENT ${PACKAGE_NAME}
)
