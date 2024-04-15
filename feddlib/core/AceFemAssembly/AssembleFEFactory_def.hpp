#ifndef ASSEMBLEFEFACTORY_DEF_hpp
#define ASSEMBLEFEFACTORY_DEF_hpp

#include "AssembleFEFactory_decl.hpp"

namespace FEDD {


template <class SC, class LO, class GO, class NO>
AssembleFEFactory<SC,LO,GO,NO>::AssembleFEFactory(){

}


template <class SC, class LO, class GO, class NO>
typename AssembleFEFactory<SC,LO,GO,NO>::AssembleFEPtr_Type AssembleFEFactory<SC,LO,GO,NO>::build(string problemType, int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple)
{
	AssembleFEPtr_Type assembleFE;

	AssembleFEPtr_Type assembleFESpecific;

	if(problemType == "Laplace"){
		//AssembleFE_Laplace<SC,LO,GO,NO> assembleFESpecific  = new AssembleFE_Laplace<SC,LO,GO,NO>(flag,nodesRefConfig, params);
		Teuchos::RCP<AssembleFE_Laplace<SC,LO,GO,NO>> assembleFESpecific(new AssembleFE_Laplace<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
  else if (problemType == "NonLinearLaplace") {
    	Teuchos::RCP<AssembleFENonLinLaplace<SC, LO, GO, NO>> assembleFESpecific(new AssembleFENonLinLaplace<SC, LO, GO, NO>(flag, nodesRefConfig,
                                                    params, tuple));
    	assembleFE = assembleFESpecific;
  }
	else if(problemType == "NavierStokes"){
		Teuchos::RCP<AssembleFENavierStokes<SC,LO,GO,NO>> assembleFESpecific(new AssembleFENavierStokes<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "GeneralizedNewtonian"){
		Teuchos::RCP<AssembleFEGeneralizedNewtonian<SC,LO,GO,NO>> assembleFESpecific(new AssembleFEGeneralizedNewtonian<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "LinearElasticity"){
		Teuchos::RCP<AssembleFE_LinElas<SC,LO,GO,NO>> assembleFESpecific(new AssembleFE_LinElas<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "NonLinearElasticity"){
		Teuchos::RCP<AssembleFE_NonLinElas<SC,LO,GO,NO>> assembleFESpecific(new AssembleFE_NonLinElas<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "NonLinearElasticity2"){
		Teuchos::RCP<AssembleFE_NonLinElas2<SC,LO,GO,NO>> assembleFESpecific(new AssembleFE_NonLinElas2<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "SCI_NH"){
		Teuchos::RCP<AssembleFE_SCI_NH<SC,LO,GO,NO>> assembleFESpecific(new AssembleFE_SCI_NH<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	// Structure interaction model established by Klemens Uhlmann
	else if(problemType == "SCI_SMC_MLCK"){
		Teuchos::RCP<AssembleFE_SCI_SMC_MLCK<SC,LO,GO,NO>> assembleFESpecific(new AssembleFE_SCI_SMC_MLCK<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "SCI_SMC_Active_Growth_Reorientation"){
		Teuchos::RCP<AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO>> assembleFESpecific(new AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else
    		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No specific implementation for your request.");


	return assembleFE;
};

}
#endif
