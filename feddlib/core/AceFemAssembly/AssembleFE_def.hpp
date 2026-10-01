#ifndef ASSEMBLEFE_DEF_hpp
#define ASSEMBLEFE_DEF_hpp

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AssembleFE<SC, LO, GO, NO>::AssembleFE(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : rhsVec_(0), jacobian_(0), solution_(0) {
    flag_ = flag;
	nodesRefConfig_ = nodesRefConfig;

    timeStep_ = 0.;
    newtonStep_ = 0;
    globalElementID_ = -1; // First not set

    params_ = params;

	// Reading through parameterlist
    dim_ = params_->sublist("Parameter").get("Dimension", -1);

    timeIncrement_ = params_->sublist("Timestepping Parameter").get("dt", 0.1); // as develop: stationary problems set no dt

    diskTuple_ = tuple;
	
	checkParameters();

    // Checking for restart. In case of restart we need to adjust the current time step.
    bool restart = params_->sublist("Timestepping Parameter").get("Restart",false);
	  if(restart){
		  timeStep_ = params_->sublist("Timestepping Parameter").get("Time step", 0.0) ;//- timeIncrement_;
      // we need to look for the correct first time increment
      int numSegments = params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

      for(int i=1; i <= numSegments; i++){

            double startTime = params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
            
            if(startTime-1e-10 < timeStep_)
              timeIncrement_ = params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
      }
	
    }
    historyImported_=false; // If we restart from a previous solution, we also ne to import the history values.

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
void AssembleFE<SC,LO,GO,NO>::advanceInTime( Teuchos::RCP<TimeSteppingTools> timeSteppingTool){
  timeIncrement_ = timeSteppingTool->get_dt();
  timeStep_ = timeSteppingTool->currentTime();
};

template <class SC, class LO, class GO, class NO>
void AssembleFE<SC,LO,GO,NO>::synchronizeTime(Teuchos::RCP<TimeSteppingTools> timeSteppingTool){
  timeIncrement_ = timeSteppingTool->get_dt();
  timeStep_ = timeSteppingTool->currentTime();
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
	this->solution_.reset( new vec_dbl_Type (solution.size(),0.) ); // as develop: elements need not allocate it

	for(int i=0; i< solution.size();i++){
		(*this->solution_)[i] = solution[i];
	}

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

template <class SC, class LO, class GO, class NO>
void AssembleFE<SC, LO, GO, NO>::setLocalHistory(vec_dbl_Type history) {
	  TEUCHOS_TEST_FOR_EXCEPTION(history.size() != history_.size(), std::runtime_error, "Input history and current history have different length. History input " << history.size() << " history current " <<history_.size() );
    this->history_ = history;
    historyImported_=true;
};
template <class SC, class LO, class GO, class NO>
void AssembleFE<SC, LO, GO, NO>::setLocalHistoryUpdated(vec_dbl_Type history) {
	  TEUCHOS_TEST_FOR_EXCEPTION(history.size() != historyUpdated_.size(), std::runtime_error, "Input history and current history have different length. History input " << history.size() << " history current " <<historyUpdated_.size() );
    this->historyUpdated_ = history;
  
};


}
#endif
