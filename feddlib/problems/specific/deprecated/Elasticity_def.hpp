#ifndef ELASTICITY_def_hpp
#define ELASTICITY_def_hpp
#include "Elasticity_decl.hpp"
/*!
 Definition of Elasticity
 
 @brief Elasticity
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template<class SC,class LO,class GO,class NO>
Elasticity<SC,LO,GO,NO>::Elasticity(const DomainConstPtr_Type  &domain, std::string FEType, ParameterListPtr_Type parameterList)
{
    if(parameterList->sublist("Parameter").get("Material model","linear")=="linear")
        linearProblem_ = Teuchos::rcp( new LinElas_Type( domain, FEType, parameterList ) );
    else
        nonLinearProblem_ = Teuchos::rcp( new NonLinElas_Type( domain, FEType, parameterList ) );
}
template<class SC,class LO,class GO,class NO>
Elasticity<SC,LO,GO,NO>::~Elasticity(){
    linearProblem_.reset();
    nonLinearProblem_.reset();
}

template<class SC,class LO,class GO,class NO>
void Elasticity<SC,LO,GO,NO>::info(){
    this->infoProblem();
    if(!nonLinearProblem_.is_null())
        this->infoNonlinProblem();
}
    
template<class SC,class LO,class GO,class NO>
void Elasticity<SC,LO,GO,NO>::assemble(){
    if(!nonLinearProblem_.is_null())
        nonLinearProblem_->assemble();
    else
        linearProblem_->assemble();
}

template<class SC,class LO,class GO,class NO>
void Elasticity<SC,LO,GO,NO>::reAssemble(std::string type) const {
    
    if(!nonLinearProblem_.is_null())
        nonLinearProblem_->reAssemble( type );
    else{
        if(this->getVerbose())
            std::cout << "-- Nothing to reassemble for linear elasticity --" << std::endl;
    }
        
}

template<class SC,class LO,class GO,class NO>
void Elasticity<SC,LO,GO,NO>::reAssembleExtrapolation(){

    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only Newton/NOX implemented for nonlinear material models!");

}

template<class SC,class LO,class GO,class NO>
void Elasticity<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                            const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                            ) const
{
    if(!nonLinearProblem_.is_null())
        nonLinearProblem_->evalModelImpl( inArgs, outArgs );
    else{
        if(this->getVerbose())
            std::cout << "-- Nothing to evaluate for linear elasticity --" << std::endl;
    }
        
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > Elasticity<SC,LO,GO,NO>::create_W_op() const
{
    
    if(!nonLinearProblem_.is_null())
        return nonLinearProblem_->create_W_op( );
    else{
        if(this->getVerbose())
            std::cout << "-- No create_W_op for linear elasticity --" << std::endl;
        Teuchos::RCP<Thyra::LinearOpBase<SC> > W_opDummy;
        return W_opDummy;
    }
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > Elasticity<SC,LO,GO,NO>::create_W_prec() const
{
    
    if(!nonLinearProblem_.is_null())
        return nonLinearProblem_->create_W_prec( );
    else{
        if(this->getVerbose())
            std::cout << "-- No create_W_prec for linear elasticity --" << std::endl;
        Teuchos::RCP<Thyra::PreconditionerBase<SC> > W_precDummy;
        return W_precDummy;
    }
}

template<class SC,class LO,class GO,class NO>
void Elasticity<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type) const{
    
    if(!nonLinearProblem_.is_null())
        return nonLinearProblem_->calculateNonLinResidualVec( type );
    else{
        if(this->getVerbose())
            std::cout << "-- No computation of nonlinear residual for linear elasticity --" << std::endl;
    }
}
}
#endif
