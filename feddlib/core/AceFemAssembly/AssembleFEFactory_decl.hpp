#ifndef ASSEMBLEFEFACTORY_DECL_hpp
#define ASSEMBLEFEFACTORY_DECL_hpp


#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFEAceLaplace.hpp"

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFEFactory {
  public:


	typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;
        typedef Teuchos::RCP<AssembleFE_Type> AssembleFEPtr_Type;

	AssembleFEFactory();

	AssembleFEPtr_Type build( string problemType, int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params);



};
}
#endif




