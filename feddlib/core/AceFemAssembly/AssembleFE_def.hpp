#ifndef ASSEMBLEFE_DEF_hpp
#define ASSEMBLEFE_DEF_hpp

#include "AssembleFE_decl.hpp"

namespace FEDD {

/*!

 \brief Constructor 

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem

*/
template <class SC, class LO, class GO, class NO>
AssembleFE<SC,LO,GO,NO>::AssembleFE(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params)
{ 
	flag_=flag;
	nodesRefConfig_ = nodesRefConfig;

	timeStep_ =0. ;
	
	params_=params;

	// Reading through parameterlist
	dim_= params_->sublist("Parameter").get("dim",-1);
	FEType1_= params_->sublist("Parameter").get("Discretization1","none");
	FEType2_= params_->sublist("Parameter").get("Discretization1","none");
	dofs1_= params_->sublist("Parameter").get("Dofs1",-1);
	dofs2_= params_->sublist("Parameter").get("Dofs2",-1);
	timeProblem_= params_->sublist("Timestepping Parameter").get("Timeproblem",false);
	
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
/*!

 \brief Checking if all necessary parameters are set

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::checkParameter( ){
	TEUCHOS_TEST_FOR_EXCEPTION(dim_=-1, std::runtime_error, "Dimension not initialized");
	TEUCHOS_TEST_FOR_EXCEPTION(dofs1_=-1, std::runtime_error, "Dofs1 not initialized");
	TEUCHOS_TEST_FOR_EXCEPTION(dofs2_=-1, std::runtime_error, "Dofs2 not initialized");
	TEUCHOS_TEST_FOR_EXCEPTION(FEType1_=="none", std::runtime_error, "FEType1 not initialized");
	TEUCHOS_TEST_FOR_EXCEPTION(FEType2_=="none", std::runtime_error, "FEType2 not initialized");
	//TEUCHOS_TEST_FOR_EXCEPTION(timeProblem_==-1, std::runtime_error, "Dofs not initialized");

};



/*!

 \brief Setting external parameters 

@param[in] params Update parameters that can be set externaly 

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::updateParams( ParameterListPtr_Type params){
	params_ = params;

};

/*!

 \brief Advancing in time. Useful for timedependent problems

@param[in] dt Timestepping length

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::advanceInTime( double dt){
	timeStep_ = timeStep_ + dt;

};

/*!

 \brief Returning current timestep

@param[out] timeStep_ Current time step

*/
template <class SC, class LO, class GO, class NO>
double AssembleFE<SC,LO,GO,NO>::getTimestep(){
	return timeStep_ ;

};

/*!

 \brief Updating solution. Useful for timedependent problems

@param[in] solution New solution we want to insert

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::updateSolution( vec_dbl_Type solution){
        
	TEUCHOS_TEST_FOR_EXCEPTION(solution_.size() != solution.size(), std::runtime_error, "Dofs of solutions is not the same");
	solution_ = solution;

};

/*!

 \brief Returning solution

@param[out] solution_ 

*/
template <class SC, class LO, class GO, class NO>
vec_dbl_Type AssembleFE<SC,LO,GO,NO>::getSolution( ){
	return solution_;

};

/*!

 \brief Preprocessing. TBD.

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::preProcessing( ){
	

};

/*!

 \brief Postprocessing. TBD.


*/
template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::postProcessing( ){


};

/*!

 \brief Returning Dimension

@param[out] solution_ 

*/
template <class SC, class LO, class GO, class NO>
int AssembleFE<SC,LO,GO,NO>::getDim( ){
	return dim_;

};

/*!

 \brief Returning nodes of reference configuration

@param[out] solution_ 

*/
template <class SC, class LO, class GO, class NO>
vec2D_dbl_Type AssembleFE<SC,LO,GO,NO>::getNodesRefConfig( ){
	return nodesRefConfig_;

};


}
#endif

