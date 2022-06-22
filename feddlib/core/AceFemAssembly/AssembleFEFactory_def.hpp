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
		//AssembleFEAceLaplace<SC,LO,GO,NO> assembleFESpecific  = new AssembleFEAceLaplace<SC,LO,GO,NO>(flag,nodesRefConfig, params);
		Teuchos::RCP<AssembleFEAceLaplace<SC,LO,GO,NO>> assembleFESpecific(new AssembleFEAceLaplace<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "NavierStokes"){
		Teuchos::RCP<AssembleFENavierStokes<SC,LO,GO,NO>> assembleFESpecific(new AssembleFENavierStokes<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "NavierStokesNonNewtonian"){
		Teuchos::RCP<AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>> assembleFESpecific(new AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "LinearElasticity"){
		Teuchos::RCP<AssembleFEAceLinElas<SC,LO,GO,NO>> assembleFESpecific(new AssembleFEAceLinElas<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}
	else if(problemType == "NonLinearElasticity"){
		Teuchos::RCP<AssembleFEAceNonLinElas<SC,LO,GO,NO>> assembleFESpecific(new AssembleFEAceNonLinElas<SC,LO,GO,NO>(flag,nodesRefConfig, params,tuple) );
		assembleFE = assembleFESpecific;
	}

	else
    		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No specific implementation for your request.");


	return assembleFE;
};

}
#endif
