#ifndef DIFFUSIONREACTION_def_hpp
#define DIFFUSIONREACTION_def_hpp
#include "DiffusionReaction_decl.hpp"
/*!
 Definition of Diffusion Reaction Equation
 
 @brief Diffusion Reaction Equation
 @author Lea Sassmannshausen
 @version 1.0
 @copyright LS
 */

namespace FEDD {

template<class SC,class LO,class GO,class NO>
DiffusionReaction<SC,LO,GO,NO>::DiffusionReaction(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList, vec2D_dbl_Type diffusionTensor,  RhsFunc_Type reactionFunc, bool vectorDiffusion):
Problem<SC,LO,GO,NO>(parameterList, domain->getComm()),
vectorDiffusion_(vectorDiffusion),
A_(),
u_rep_(),
reactionFunc_()
{
 
    this->addVariable( domain , FEType , "u" , 1);
    this->dim_ = this->getDomain(0)->getDimension();
	
	diffusionTensor_ = diffusionTensor;

    funcParameter_.push_back(this->parameterList_->sublist("Parameter").get("E0",1.0));
    funcParameter_.push_back(this->parameterList_->sublist("Parameter").get("E1",0.5));

	reactionFunc_ = reactionFunc;

    u_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapRepeated() ) );

	// Test for exception!!

}

template<class SC,class LO,class GO,class NO>
DiffusionReaction<SC,LO,GO,NO>::~DiffusionReaction(){

}
    
template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::info(){
    this->infoProblem();
}

/*!
    \brief assemble constant matrices, that remain the same

*/
template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::assembleConstantMatrices( std::string type ) const{
    
    if (this->verbose_)
        std::cout << "-- Assembly Laplace with Diffusion Tensor ... " << std::flush;

    A_.reset(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );

    this->feFactory_->assemblyLaplaceDiffusion(this->dim_, this->domain_FEType_vec_.at(0), 2, A_, this->diffusionTensor_ );

    if (this->system_.is_null())
        this->system_.reset(new BlockMatrix_Type(1));

    this->system_->addBlock(A_,0,0);
    
    //this->assembleSourceTerm( 0. );
    //this->addToRhs( this->sourceTerm_ );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

/*!
    In case of a nonlinear problem we call reassemble for each newton step. 
    For now diffusion-reaction is implemented for linear reaction term.
*/

template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if (type=="") {
        if (this->verbose_)
            std::cout << "-- Assembly Diffusion ... " << std::endl;

        assembleConstantMatrices();
        if (this->verbose_)
            std::cout << "done -- " << std::endl;
    }
    //else
    //    reAssemble( type );
    // Here we would enter the reaction component. Since it is dependent on u it would need an update. So probably an update RHS thing.

}

template<class SC,class LO,class GO,class NO>
typename DiffusionReaction<SC,LO,GO,NO>::MatrixPtr_Type DiffusionReaction<SC,LO,GO,NO>::getMassMatrix() const{
	
    MatrixPtr_Type A;
	A = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
	this->feFactory_->assemblyMass(this->dim_,this->domain_FEType_vec_.at(0),"Scalar", A);

	return A;

}
    
}
#endif
