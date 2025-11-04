#ifndef LAPLACE_def_hpp
#define LAPLACE_def_hpp

/*!
 Definition of Laplace
 
 @brief Laplace
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"

namespace FEDD {

template<class SC,class LO,class GO,class NO>
Laplace<SC,LO,GO,NO>::Laplace(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList, bool vectorLaplace):
Problem<SC,LO,GO,NO>(parameterList, domain->getComm()),
vectorLaplace_(vectorLaplace)
{
 
    this->addVariable( domain , FEType , "u" , 1);
    this->dim_ = this->getDomain(0)->getDimension();
}

template<class SC,class LO,class GO,class NO>
Laplace<SC,LO,GO,NO>::~Laplace(){

}
    
template<class SC,class LO,class GO,class NO>
void Laplace<SC,LO,GO,NO>::info(){
    this->infoProblem();
}

template<class SC,class LO,class GO,class NO>
void Laplace<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if (this->verbose_)
        std::cout << "-- Assembly Laplace ... " << std::flush;

    MatrixPtr_Type A;
    vec_dbl_Type funcParameter(1,0.);
    if (vectorLaplace_){
        A = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyLaplaceVecField(this->dim_, this->domain_FEType_vec_.at(0), 2, A );
    }
    else{
        A = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyLaplace(this->dim_, this->domain_FEType_vec_.at(0), 2, A );
    }
    
    this->system_->addBlock(A,0,0);
    
    this->assembleSourceTerm( 0. );
    
    this->addToRhs( this->sourceTerm_ );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

template<class SC,class LO,class GO,class NO>
typename Laplace<SC,LO,GO,NO>::MatrixPtr_Type Laplace<SC,LO,GO,NO>::getMassMatrix() const{
	
    MatrixPtr_Type A;
	A = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
	this->feFactory_->assemblyMass(this->dim_,this->domain_FEType_vec_.at(0),"Scalar", A);

	return A;

}
    
}
#endif
