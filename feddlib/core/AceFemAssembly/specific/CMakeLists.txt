SET(AceFemAssemblySpecific_HEADERS
	AceFemAssembly/specific/AssembleFEAceLaplace_decl.hpp
	AceFemAssembly/specific/AssembleFEAceLaplace_def.hpp
	AceFemAssembly/specific/AssembleFEAceLaplace.hpp
	AceFemAssembly/specific/AssembleFENavierStokes_decl.hpp
	AceFemAssembly/specific/AssembleFENavierStokes_def.hpp
	AceFemAssembly/specific/AssembleFENavierStokes.hpp
	AceFemAssembly/specific/AssembleFENavierStokesNonNewtonian_decl.hpp
	AceFemAssembly/specific/AssembleFENavierStokesNonNewtonian_def.hpp
	AceFemAssembly/specific/AssembleFENavierStokesNonNewtonian.hpp
	AceFemAssembly/specific/AssembleFEAceLinElas_decl.hpp
	AceFemAssembly/specific/AssembleFEAceLinElas_def.hpp
	AceFemAssembly/specific/AssembleFEAceLinElas.hpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas_decl.hpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas_def.hpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas.hpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas2_decl.hpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas2_def.hpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas2.hpp
	AceFemAssembly/specific/AssembleFE_SCI_NH_decl.hpp
	AceFemAssembly/specific/AssembleFE_SCI_NH_def.hpp
	AceFemAssembly/specific/AssembleFE_SCI_NH.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK_decl.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK_def.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation_def.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation.hpp
CACHE INTERNAL "")

SET(AceFemAssemblySpecific_SOURCES
	AceFemAssembly/specific/AssembleFEAceLaplace.cpp
	AceFemAssembly/specific/AssembleFENavierStokes.cpp
	AceFemAssembly/specific/AssembleFENavierStokesNonNewtonian.cpp
	AceFemAssembly/specific/AssembleFEAceLinElas.cpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas.cpp
	AceFemAssembly/specific/AssembleFEAceNonLinElas2.cpp
	AceFemAssembly/specific/AssembleFE_SCI_NH.cpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK.cpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation.cpp
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
