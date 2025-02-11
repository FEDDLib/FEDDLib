SET(AceFemAssemblySpecific_HEADERS
	AceFemAssembly/specific/AssembleFE_Laplace_decl.hpp
	AceFemAssembly/specific/AssembleFE_Laplace_def.hpp
	AceFemAssembly/specific/AssembleFE_Laplace.hpp
  AceFemAssembly/specific/AssembleFENonLinLaplace_decl.hpp
  AceFemAssembly/specific/AssembleFENonLinLaplace_def.hpp
  AceFemAssembly/specific/AssembleFENonLinLaplace.hpp
	AceFemAssembly/specific/AssembleFENavierStokes_decl.hpp
	AceFemAssembly/specific/AssembleFENavierStokes_def.hpp
	AceFemAssembly/specific/AssembleFENavierStokes.hpp
	AceFemAssembly/specific/AssembleFEGeneralizedNewtonian_decl.hpp
	AceFemAssembly/specific/AssembleFEGeneralizedNewtonian_def.hpp
	AceFemAssembly/specific/AssembleFEGeneralizedNewtonian.hpp
	AceFemAssembly/specific/AssembleFE_LinElas_decl.hpp
	AceFemAssembly/specific/AssembleFE_LinElas_def.hpp
	AceFemAssembly/specific/AssembleFE_LinElas.hpp
	AceFemAssembly/specific/AssembleFE_NonLinElas_decl.hpp
	AceFemAssembly/specific/AssembleFE_NonLinElas_def.hpp
	AceFemAssembly/specific/AssembleFE_NonLinElas.hpp
	AceFemAssembly/specific/AssembleFE_NonLinElas2_decl.hpp
	AceFemAssembly/specific/AssembleFE_NonLinElas2_def.hpp
	AceFemAssembly/specific/AssembleFE_NonLinElas2.hpp
	AceFemAssembly/specific/AssembleFE_SCI_NH_decl.hpp
	AceFemAssembly/specific/AssembleFE_SCI_NH_def.hpp
	AceFemAssembly/specific/AssembleFE_SCI_NH.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK_decl.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK_def.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation_def.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/CarreauYasuda_decl.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/CarreauYasuda_def.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/CarreauYasuda.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/PowerLaw_decl.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/PowerLaw_def.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/PowerLaw.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/Dimless_Carreau_decl.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/Dimless_Carreau_def.hpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/Dimless_Carreau.hpp
CACHE INTERNAL "")

SET(AceFemAssemblySpecific_SOURCES
	AceFemAssembly/specific/AssembleFE_Laplace.cpp
  AceFemAssembly/specific/AssembleFENonLinLaplace.cpp
	AceFemAssembly/specific/AssembleFENavierStokes.cpp
	AceFemAssembly/specific/AssembleFEGeneralizedNewtonian.cpp
	AceFemAssembly/specific/AssembleFE_LinElas.cpp
	AceFemAssembly/specific/AssembleFE_NonLinElas.cpp
	AceFemAssembly/specific/AssembleFE_NonLinElas2.cpp
	AceFemAssembly/specific/AssembleFE_SCI_NH.cpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK.cpp
	AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation.cpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/CarreauYasuda.cpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/PowerLaw.cpp
	AceFemAssembly/specific/GeneralizedNewtonianModels/Dimless_Carreau.cpp
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
