#ifndef ASSEMBLEFEFACTORY_DEF_hpp
#define ASSEMBLEFEFACTORY_DEF_hpp

#include "AssembleFEFactory_decl.hpp"

namespace FEDD {

/*!
 \brief Constructor

*/
template <class SC, class LO, class GO, class NO>
AssembleFEFactory<SC,LO,GO,NO>::AssembleFEFactory(){

}

/*!
 \brief We only need on function to build assembleFE, where we define the problem type and the AssembleFE object we extend to a specific problem.

@param[in] problemType Type of specfic problem we focus on
@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem
*/
template <class SC, class LO, class GO, class NO>
typename AssembleFEFactory<SC,LO,GO,NO>::AssembleFEPtr_Type AssembleFEFactory<SC,LO,GO,NO>::build(string problemType, int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params)
{
	AssembleFEPtr_Type assembleFE;

	AssembleFEPtr_Type assembleFESpecific;

	if(problemType == "Laplace"){
		cout<< " builng Laplace elements " << endl;
		//AssembleFEAceLaplace<SC,LO,GO,NO> assembleFESpecific  = new AssembleFEAceLaplace<SC,LO,GO,NO>(flag,nodesRefConfig, params);
		Teuchos::RCP<AssembleFEAceLaplace<SC,LO,GO,NO>> assembleFESpecific(new AssembleFEAceLaplace<SC,LO,GO,NO>(flag,nodesRefConfig, params) );
		assembleFE = assembleFESpecific;

	}
	else if(problemType == "NavierStokes"){
		//AssembleFESpecificNavierStokes_Type assembeFESpecific  = new AssembleFeAceNavierStokes(flag,nodesRefConfig, parameters);
	}
	else
    		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No specific implementation for your request.");


	return assembleFE;
};

}
#endif

