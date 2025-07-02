#ifndef NonLinElasticity_def_hpp
#define NonLinElasticity_def_hpp
#include "NonLinElasticity_decl.hpp"
/*!
 Definition of NonLinElasticity
 
 @brief NonLinElasticity
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


namespace FEDD {
template<class SC,class LO,class GO,class NO>
NonLinElasticity<SC,LO,GO,NO>::NonLinElasticity(const DomainConstPtr_Type  &domain, std::string FEType, ParameterListPtr_Type parameterList):
NonLinearProblem<SC,LO,GO,NO>( parameterList, domain->getComm() ),
u_rep_()
{
    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    this->initNOXParameters();
    this->addVariable( domain , FEType , "u" , domain->getDimension());
    this->dim_ = this->getDomain(0)->getDimension();
    
    C_ = this->parameterList_->sublist("Parameter").get("C",1.);
    
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
//    poissonRatio_ = this->parameterList_->sublist("Parameter").get("Poisson Ratio",0.4);
//    mue_ = this->parameterList_->sublist("Parameter").get("Mu",2.0e+6);
//    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstante \lambda
//    E_ = mue_*2.*(1. + poissonRatio_);
//    lambda_ = (poissonRatio_*E_)/((1 + poissonRatio_)*(1 - 2*poissonRatio_));
    
}

template<class SC,class LO,class GO,class NO>
NonLinElasticity<SC,LO,GO,NO>::~NonLinElasticity(){

}

template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::info(){
    this->infoProblem();
    this->infoNonlinProblem();
}
    
template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::assemble(std::string type) const{
    
    if (type == ""){
        if (this->verbose_)
            std::cout << "-- Assembly nonlinear elasticity ... " << std::flush;

        MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        
        this->feFactory_->assemblyEmptyMatrix(A);

        this->system_.reset(new BlockMatrix_Type(1));
        this->system_->addBlock( A, 0, 0 );
                
        double density = this->parameterList_->sublist("Parameter").get("Density",1000.);
        std::string sourceType = 	this->parameterList_->sublist("Parameter").get("Source Type","volume");

        this->assembleSourceTerm( 0. );
        if(sourceType == "volume")
            this->sourceTerm_->scale(density);
        
        this->addToRhs( this->sourceTerm_ );

        this->setBoundariesRHS();
                
        
        this->solution_->putScalar(0.);
        
        u_rep_ = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ));
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true);
        
        if (this->verbose_)
            std::cout << "done -- " << std::endl;
        
        this->reAssemble("Newton-Residual");
    }
    else
        this->reAssemble(type);
}

template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::reAssemble(std::string type) const {
    std::string material_model = this->parameterList_->sublist("Parameter").get("Material model","Neo-Hooke");

    if (this->verbose_)
        std::cout << "-- Reassembly nonlinear elasticity with material model " << material_model <<" ("<<type <<") ... " << std::flush;
    
    if (type=="Newton-Residual") {
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);

        u_rep_->importFromVector(u, true);
        
        MultiVectorPtr_Type f = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
        MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyElasticityJacobianAndStressAceFEM(this->dim_, this->getDomain(0)->getFEType(), W, f, u_rep_, this->parameterList_, C_);
        
        MultiVectorPtr_Type fUnique = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
        fUnique->putScalar(0.);
        fUnique->exportFromVector( f, true, "Add" );

        this->residualVec_->addBlock( fUnique, 0 );
        this->system_->addBlock( W, 0, 0 );

        //        MultiVectorPtr_Type f = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
//        this->feFactory_->assemblyElasticityStressesAceFEM(this->dim_, this->getDomain(0)->getFEType(), f, u_rep_, material_model, E_, nu_, C_);
//        MultiVectorPtr_Type fUnique = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
//        fUnique->putScalar(0.);
//        fUnique->importFromVector( f, true, "Add" );
//        this->residualVec_->addBlock( fUnique, 0 );
    }
    else if(type=="Newton"){ //we already assemble the new tangent when we calculate the stresses above
        
    }
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}
    
template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const
{
    
}

template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions){


    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only Newton/NOX implemented for nonlinear material models!");

}


template<class SC,class LO,class GO,class NO>
void NonLinElasticity<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
    this->reAssemble("Newton-Residual");
    if (!type.compare("standard")){
        this->residualVec_->update(-1.,*this->rhs_,1.);
        //if ( !this->sourceTerm_.is_null() )
        //    this->residualVec_->update(-1.,*this->sourceTerm_,1.);
    }
    else if(!type.compare("reverse")){
        this->residualVec_->update(1.,*this->rhs_,-1.); // this = -1*this + 1*rhs
        //if ( !this->sourceTerm_.is_null() )
        //    this->residualVec_->update(1.,*this->sourceTerm_,1.);
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unknown type for residual computation.");
    }
    
    // this might be set again by the TimeProblem after adding of M*u
    this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    
}
}
#endif
