#ifndef ASSEMBLEFE_DEF_hpp
#define ASSEMBLEFE_DEF_hpp

#include "AssembleFE_decl.hpp"

namespace FEDD {


template <class SC, class LO, class GO, class NO>
AssembleFE<SC,LO,GO,NO>::AssembleFE(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params)
{
	flag_=flag;
	nodesRefConfig_ = nodesRefConfig;

	timeStep_ =0. ;

	params_=params;

	// Reading through parameterlist
	dim_= params_->sublist("Parameter").get("Dimension",-1);
	FEType1_= params_->sublist("Parameter").get("Discretization1","none");
	//FEType2_= params_->sublist("Parameter").get("Discretization1","none");
	dofs1_= params_->sublist("Parameter").get("Dofs1",-1);
	//dofs2_= params_->sublist("Parameter").get("Dofs2",-1);
	//timeProblem_= params_->sublist("Timestepping Parameter").get("Timeproblem",false);

	checkParameter();

	solution_ = vec_dbl_Type(dofs1_);
/// Element Numbering for triangular elements:
/*!
    - Triangle numbering

                    2
                  * *
                *   *
              4	    5
            *       *
          *         *
        1 * * 3 * * 0
------------------------------------------------------------------------------------
*/
/*!
    - Tetrahedral numbering

                Face 1          Face2               Face 3            Face 4
                    2      2 * * 9 * * 3        3 * * 9 * * 2          	 3
                  * *      *          *          *          * 		 	* *
                *   *      *        *             *        *          *   *
              5	 	6      6      7                8      5         8	  7
            *       *      *    *                   *    *        *       *
          *         *      *  *                      *  *       *         *
        1 * * 4 * * 0       0                         1       1 * * 4 * * 0
------------------------------------------------------------------------------------
*/
}


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::checkParameter( ){
	TEUCHOS_TEST_FOR_EXCEPTION(dim_==-1, std::runtime_error, "Dimension not initialized");
	TEUCHOS_TEST_FOR_EXCEPTION(dofs1_==-1, std::runtime_error, "Dofs1 not initialized");
	//TEUCHOS_TEST_FOR_EXCEPTION(dofs2_==-1, std::runtime_error, "Dofs2 not initialized");
	TEUCHOS_TEST_FOR_EXCEPTION(FEType1_=="none", std::runtime_error, "FEType1 not initialized");
	//TEUCHOS_TEST_FOR_EXCEPTION(FEType2_=="none", std::runtime_error, "FEType2 not initialized");
	//TEUCHOS_TEST_FOR_EXCEPTION(timeProblem_==-1, std::runtime_error, "Dofs not initialized");

};


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::updateParams( ParameterListPtr_Type params){
	params_ = params;

};


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::advanceInTime( double dt){
	timeStep_ = timeStep_ + dt;

};


template <class SC, class LO, class GO, class NO>
double AssembleFE<SC,LO,GO,NO>::getTimestep(){
	return timeStep_ ;

};


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::updateSolution( vec_dbl_Type solution){

	TEUCHOS_TEST_FOR_EXCEPTION(solution_.size() != solution.size(), std::runtime_error, "Dofs of solutions is not the same");
	solution_ = solution;

};


template <class SC, class LO, class GO, class NO>
vec_dbl_Type AssembleFE<SC,LO,GO,NO>::getSolution( ){
	return solution_;

};


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::preProcessing( ){


};


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::postProcessing( ){


};


template <class SC, class LO, class GO, class NO>
int AssembleFE<SC,LO,GO,NO>::getDim( ){
	return dim_;

};


template <class SC, class LO, class GO, class NO>
vec2D_dbl_Type AssembleFE<SC,LO,GO,NO>::getNodesRefConfig( ){
	return nodesRefConfig_;

};


}
#endif
