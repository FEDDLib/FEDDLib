#ifndef LINEARPROBLEM_DEF_hpp
#define LINEARPROBLEM_DEF_hpp
#include "LinearProblem_decl.hpp"

/*!
 Definition of LinearProblem
 
 @brief  LinearProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
using namespace std;
namespace FEDD {
template<class SC,class LO,class GO,class NO>
LinearProblem<SC,LO,GO,NO>::LinearProblem(CommConstPtr_Type comm):
Problem<SC,LO,GO,NO>(comm,true)
{}

template<class SC,class LO,class GO,class NO>
LinearProblem<SC,LO,GO,NO>::LinearProblem(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm ):
Problem<SC,LO,GO,NO>(parameterList, comm)
{}
    
template<class SC,class LO,class GO,class NO>
LinearProblem<SC,LO,GO,NO>::~LinearProblem(){

}
}
#endif
