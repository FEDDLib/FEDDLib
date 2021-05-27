#ifndef MinPrecProblem_DEF_hpp
#define MinPrecProblem_DEF_hpp
#include "MinPrecProblem_decl.hpp"

/*!
 Definition of MinPrecProblem
 @brief  MinPrecProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


namespace FEDD {
template<class SC,class LO,class GO,class NO>
MinPrecProblem<SC,LO,GO,NO>::MinPrecProblem(ParameterListPtr_Type pl, CommConstPtr_Type comm):
Problem<SC,LO,GO,NO>(pl, comm)
{

}
    
template<class SC,class LO,class GO,class NO>
MinPrecProblem<SC,LO,GO,NO>::~MinPrecProblem( )
{
    
}

template<class SC,class LO,class GO,class NO>
void MinPrecProblem<SC,LO,GO,NO>::initializeSystem(BlockMatrixPtr_Type system){
    this->system_ = system;
}
    
template<class SC,class LO,class GO,class NO>
void MinPrecProblem<SC,LO,GO,NO>::initializeDomains(DomainConstPtr_vec_Type& domainPtr_vec){
    this->domainPtr_vec_ = domainPtr_vec;
}
template<class SC,class LO,class GO,class NO>
void MinPrecProblem<SC,LO,GO,NO>::initializeLinSolverBuilder(LinSolverBuilderPtr_Type linearSolverBuilder){
    this->linearSolverBuilder_ = linearSolverBuilder;
}
}
#endif
