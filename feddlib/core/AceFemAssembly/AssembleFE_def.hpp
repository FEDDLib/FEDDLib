#ifndef ASSEMBLEFE_DEF_hpp
#define ASSEMBLEFE_DEF_hpp

namespace FEDD {


template <class SC, class LO, class GO, class NO>
AssembleFE<SC,LO,GO,NO>::AssembleFE(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
rhsVec_(0),
jacobian_(0),
solution_(0)
{
	flag_=flag;
	nodesRefConfig_ = nodesRefConfig;

	timeStep_ =0. ;
	newtonStep_ =0 ;
	globalElementID_=-1; // First not set


	params_=params;

	// Reading through parameterlist
	dim_= params_->sublist("Parameter").get("Dimension",-1);

	timeIncrement_= params_->sublist("Timestepping Parameter").get("dt",0.1);

	diskTuple_= tuple;
	
	checkParameters();

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

                Face 1          Face2               Face 3          Face 4
                    2      2 * * 9 * * 3        3 * * 9 * * 2          	    3
                  * *      *          *          *          *             * *
                *   *      *        *             *        *            *   *
              5     6      6      7                8      5           8     7
            *       *      *    *                   *    *          *       *
          *         *      *  *                      *  *         *         *
        1 * * 4 * * 0       0                         1         1 * * 4 * * 0
------------------------------------------------------------------------------------
*/
}


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::checkParameters( ){
	TEUCHOS_TEST_FOR_EXCEPTION(dim_==-1, std::runtime_error, "Dimension not initialized");
};


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::updateParams( ParameterListPtr_Type params){
	params_ = params;

};


template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::advanceInTime( double dt){
	timeIncrement_ = dt;
	timeStep_ = timeStep_ + dt;
};

template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::advanceNewtonStep(){
	newtonStep_ = newtonStep_+1 ;

};


template <class SC, class LO, class GO, class NO>
double AssembleFE<SC,LO,GO,NO>::getTimeStep(){
	return timeStep_ ;

};

template <class SC, class LO, class GO, class NO>
int AssembleFE<SC,LO,GO,NO>::getNewtonStep(){
	return newtonStep_ ;

};

template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::updateSolution( vec_dbl_Type solution){

	//TEUCHOS_TEST_FOR_EXCEPTION(solution_.size() != solution.size(), std::runtime_error, "Dofs of solutions is not the same");
	this->solution_.reset( new vec_dbl_Type (solution.size(),0.) );

	for(int i=0; i< solution.size();i++)
		(*solution_)[i] = solution[i];

};




template <class SC, class LO, class GO, class NO>
vec_dbl_ptr_Type AssembleFE<SC,LO,GO,NO>::getSolution( ){
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
