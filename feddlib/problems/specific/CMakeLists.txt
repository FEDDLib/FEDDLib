SET(problems_specific_HEADERS
    #specific/Elasticity_decl.hpp
    #specific/Elasticity_def.hpp
    #specific/Elasticity.hpp
    specific/Geometry_decl.hpp
    specific/Geometry_def.hpp
    specific/Geometry.hpp
    specific/Laplace_decl.hpp
 	specific/Laplace_def.hpp
 	specific/Laplace.hpp
    specific/LaplaceBlocks_decl.hpp
 	specific/LaplaceBlocks_def.hpp
 	specific/LaplaceBlocks.hpp
 	 specific/DiffusionReaction_decl.hpp
 	specific/DiffusionReaction_def.hpp
 	specific/DiffusionReaction.hpp
    #specific/LaplacePseudoCoupled_decl.hpp
    #specific/LaplacePseudoCoupled_def.hpp
    #specific/LaplacePseudoCoupled.hpp
    specific/LinElas_decl.hpp
    specific/LinElas_def.hpp
    specific/LinElas.hpp
    specific/LinElasAssFE_decl.hpp
    specific/LinElasAssFE_def.hpp
    specific/LinElasAssFE.hpp
#    specific/LinElasFirstOrder_decl.hpp
#    specific/LinElasFirstOrder_def.hpp
#    specific/LinElasFirstOrder.hpp
    specific/NavierStokes_decl.hpp
	specific/NavierStokes_def.hpp
	specific/NavierStokes.hpp
	specific/NavierStokesAssFE_decl.hpp
	specific/NavierStokesAssFE_def.hpp
	specific/NavierStokesAssFE.hpp
    specific/NonLinElasticity_decl.hpp
    specific/NonLinElasticity_def.hpp
    specific/NonLinElasticity.hpp
    specific/NonLinLaplace_decl.hpp
    specific/NonLinLaplace_def.hpp
    specific/NonLinLaplace.hpp
    specific/NonLinElasAssFE_decl.hpp
    specific/NonLinElasAssFE_def.hpp
    specific/NonLinElasAssFE.hpp
    specific/NonLinTPM_decl.hpp
    specific/NonLinTPM_def.hpp
    specific/NonLinTPM.hpp
	specific/Stokes_decl.hpp
	specific/Stokes_def.hpp
	specific/Stokes.hpp
    specific/TPM_decl.hpp
    specific/TPM_def.hpp
    specific/TPM.hpp
    specific/FSI_decl.hpp
    specific/FSI_def.hpp
    specific/FSI.hpp
CACHE INTERNAL "")

SET(problems_specific_SOURCES
	#specific/Elasticity.cpp
    specific/Geometry.cpp
	specific/Laplace.cpp
	specific/LaplaceBlocks.cpp
    #specific/LaplacePseudoCoupled.cpp
    specific/LinElas.cpp
    specific/LinElasAssFE.cpp
#    specific/LinElasFirstOrder.cpp
	specific/DiffusionReaction.cpp
	specific/NavierStokes.cpp
	specific/NavierStokesAssFE.cpp
    specific/NonLinElasticity.cpp
    specific/NonLinElasAssFE.cpp
    specific/NonLinLaplace.cpp
    specific/NonLinTPM.cpp
    specific/Stokes.cpp
    specific/TPM.cpp
    specific/FSI.cpp
CACHE INTERNAL "")


SET(LOCAL_HEADERS)
FOREACH(INC ${Problems_HEADERS})
  STRING(REPLACE "specific/" "" LOCAL_INC ${INC})
  SET(LOCAL_HEADERS ${LOCAL_HEADERS} ${LOCAL_INC})
ENDFOREACH()

INSTALL(
  FILES ${LOCAL_HEADERS}
  DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/feddlib/core/problems/specific"
  COMPONENT ${PACKAGE_NAME}
)
