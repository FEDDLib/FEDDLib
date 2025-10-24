#ifndef NonLinTPM_def_hpp
#define NonLinTPM_def_hpp
#include "NonLinTPM_decl.hpp"

/*!
 Definition of NonLinTPM
 
 @brief NonLinTPM
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

    
template<class SC,class LO,class GO,class NO>
NonLinTPM<SC,LO,GO,NO>::NonLinTPM(const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList):
NonLinearProblem<SC,LO,GO,NO>(parameterList, domainVelocity->getComm())
{

    this->addVariable( domainVelocity , FETypeVelocity , "u" , domainVelocity->getDimension());
    this->addVariable( domainPressure , FETypePressure , "p" , 1);
    this->dim_ = this->getDomain(0)->getDimension();
}

template<class SC,class LO,class GO,class NO>
NonLinTPM<SC,LO,GO,NO>::~NonLinTPM(){

}

template<class SC,class LO,class GO,class NO>
void NonLinTPM<SC,LO,GO,NO>::info(){
    this->infoProblem();
}
    

//template<class SC,class LO,class GO,class NO>
//void NonLinTPM<SC,LO,GO,NO>::assemble( std::string type ){
//    if (type == "") {
//        if (this->verbose_)
//            std::cout << "-- Assembly ... " << std::flush;
//        
//        std::string tpmType = this->getParameterList()->sublist("Parameter").get("TPM Type","Biot");
//        
//        MapConstPtr_Type mapRepeatedConst1 = this->getDomain(0)->getMapRepeated();
//        MapConstPtr_Type mapRepeatedConst2 = this->getDomain(1)->getMapRepeated();
//        MapPtr_Type mapRepeated1 = Teuchos::rcp_const_cast<Map_Type>(mapRepeatedConst1);
//        MapPtr_Type mapRepeated2 = Teuchos::rcp_const_cast<Map_Type>(mapRepeatedConst2);
//
//        MatrixPtr_Type A00(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
//        MatrixPtr_Type A01(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
//        MatrixPtr_Type A10(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
//        MatrixPtr_Type A11(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
//        
//        MultiVectorPtr_Type F0(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
//        MultiVectorPtr_Type F1(new MultiVector_Type( this->getDomain(1)->getMapRepeated(), 1 ) );
//        
//        if (u_repNewton_.is_null()){
//            u_repNewton_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() )) ;
//            u_repNewton_->putScalar(0.);
//        }
//        if (p_repNewton_.is_null()){
//            p_repNewton_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() )) ;
//            p_repNewton_->putScalar(0.);
//        }
//        
//        if (u_repTime_.is_null()){
//            u_repTime_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() )) ;
//            u_repTime_->putScalar(0.);
//        }
//        if (p_repTime_.is_null()){
//            p_repTime_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() )) ;
//            p_repTime_->putScalar(0.);
//        }
//        if(tpmType == "Biot" || tpmType == "Biot-StVK"){
//            this->feFactory_->assemblyAceGenTPM(A00,
//                                                A01,
//                                                A10,
//                                                A11,
//                                                F0,
//                                                F1,
//                                                mapRepeated1,
//                                                mapRepeated2,
//                                                this->parameterList_,
//                                                u_repNewton_,
//                                                p_repNewton_,
//                                                u_repTime_,
//                                                p_repTime_,
//                                                true,
//                                                false);
//        }
//        
//        this->system_.reset(new BlockMatrix_Type(2));
//        this->system_->addBlock( A00, 0, 0 );
//        this->system_->addBlock( A01, 0, 1 );
//        this->system_->addBlock( A10, 1, 0 );
//        this->system_->addBlock( A11, 1, 1 );
//        
//        this->initializeVectors();
//        this->initializeVectorsNonLinear();
//        
//        MultiVectorPtr_Type F0Unique = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique() ) );
//        MultiVectorPtr_Type F1Unique = Teuchos::rcp(new MultiVector_Type( this->getDomain(1)->getMapUnique() ) );
//        F0Unique->exportFromVector( F0, false, "Add" );
//        F1Unique->exportFromVector( F1, false, "Add" );
//        
//        this->rhs_->addBlock( F0Unique, 0 );
//        this->rhs_->addBlock( F1Unique, 1 );
//            
//        this->assembleSourceTerm( 0. );
//        this->addToRhs( this->sourceTerm_ );
//        
//        if (this->verbose_)
//            std::cout << "done -- " << std::endl;
//    } else {
//        this->reAssmble(type);
//    }
//}
    

template<class SC,class LO,class GO,class NO>
void NonLinTPM<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
    this->assemble("Newton-Residual");
    if (type == "external" ){// we assume reverse computation: Ax-b but we need to scale the sourceTerm with -1
        this->residualVec_->update(1.,*this->rhs_,0.);        
        if ( !this->sourceTerm_.is_null() ){
            this->residualVec_->update(-1.,*this->sourceTerm_,1.);
        }
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unknown type for residual computation.");
    }
    this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    
}
    
template<class SC,class LO,class GO,class NO>
void NonLinTPM<SC,LO,GO,NO>::reAssemble( BlockMultiVectorPtr_Type previousSolution ) const {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "This should not be used.");
    MultiVectorConstPtr_Type u_prev = previousSolution->getBlock(0);
    MultiVectorConstPtr_Type p_prev = previousSolution->getBlock(1);
    u_repTime_->importFromVector(u_prev, true);
    p_repTime_->importFromVector(p_prev, true);
}


template<class SC,class LO,class GO,class NO>
void NonLinTPM<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if (this->verbose_)
        std::cout << "-- Assembly nonlinear TPM " <<" ("<<type <<") ... " << std::flush;
    
    std::string tpmType = this->parameterList_->sublist("Parameter").get("TPM Type","Biot");
    
    if (type=="Newton-Residual") {
        type = "Assemble";
    }
    
    if (type == "FirstAssemble") {
        if (u_repNewton_.is_null()){
            u_repNewton_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() )) ;
            u_repNewton_->putScalar(0.);
        }
        if (p_repNewton_.is_null()){
            p_repNewton_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() )) ;
            p_repNewton_->putScalar(0.);
        }
        
        if (u_repTime_.is_null()){
            u_repTime_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() )) ;
            u_repTime_->putScalar(0.);
        }
        if (p_repTime_.is_null()){
            p_repTime_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() )) ;
            p_repTime_->putScalar(0.);
        }
    }
    
    
    if (type == "SetSolutionInTime") {
            TEUCHOS_TEST_FOR_EXCEPTION( u_repTime_.is_null(), std::runtime_error, "u_repTime_ not initialized.");
            TEUCHOS_TEST_FOR_EXCEPTION( p_repTime_.is_null(), std::runtime_error, "p_repTime_ not initialized.");
            MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
            MultiVectorConstPtr_Type p = this->solution_->getBlock(1);
            u_repTime_->importFromVector(u, true);
            p_repTime_->importFromVector(p, true);
    }
    else if (type == "SetSolutionNewton") {
        TEUCHOS_TEST_FOR_EXCEPTION( u_repNewton_.is_null(), std::runtime_error, "u_repNewton_ not initialized.");
        TEUCHOS_TEST_FOR_EXCEPTION( p_repNewton_.is_null(), std::runtime_error, "p_repNewton_ not initialized.");

        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        MultiVectorConstPtr_Type p = this->solution_->getBlock(1);
        u_repNewton_->importFromVector(u, true);
        p_repNewton_->importFromVector(p, true);
    }
    else if( type == "FirstAssemble" || type == "Assemble" || type == "OnlyUpdate" || type == "AssembleAndUpdate"){
        if (this->verbose_)
            std::cout << "-- External assembly ... " << std::flush;
    
        MapConstPtr_Type mapRepeatedConst1 = this->getDomain(0)->getMapRepeated();
        MapConstPtr_Type mapRepeatedConst2 = this->getDomain(1)->getMapRepeated();
        MapPtr_Type mapRepeated1 = Teuchos::rcp_const_cast<Map_Type>(mapRepeatedConst1);
        MapPtr_Type mapRepeated2 = Teuchos::rcp_const_cast<Map_Type>(mapRepeatedConst2);
                
        MatrixPtr_Type A00(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type A01(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type A10(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type A11(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );

        
        MultiVectorPtr_Type F0(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
        MultiVectorPtr_Type F1(new MultiVector_Type( this->getDomain(1)->getMapRepeated(), 1 ) );
        
        bool update(type == "FirstAssemble" || type == "Assemble" || type == "AssembleAndUpdate");
        bool updateHistory(type == "OnlyUpdate" || type == "AssembleAndUpdate");
        
        if(tpmType == "Biot" || tpmType == "Biot-StVK"){
            this->feFactory_->assemblyAceGenTPM(A00,
                                                A01,
                                                A10,
                                                A11,
                                                F0,
                                                F1,
                                                mapRepeated1,
                                                mapRepeated2,
                                                this->parameterList_,
                                                u_repNewton_,
                                                p_repNewton_,
                                                u_repTime_,
                                                p_repTime_,
                                                update,
                                                updateHistory);
        }
        if (update) {
            this->system_->addBlock( A00, 0, 0 );
            this->system_->addBlock( A01, 0, 1 );
            this->system_->addBlock( A10, 1, 0 );
            this->system_->addBlock( A11, 1, 1 );

            
            MultiVectorPtr_Type F0Unique = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique() ) );
            MultiVectorPtr_Type F1Unique = Teuchos::rcp(new MultiVector_Type( this->getDomain(1)->getMapUnique() ) );
            F0Unique->exportFromVector( F0, false, "Add" );
            F1Unique->exportFromVector( F1, false, "Add" );
            
            this->rhs_->addBlock( F0Unique, 0 );
            this->rhs_->addBlock( F1Unique, 1 );
        }
    }
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
    
}
    
template<class SC,class LO,class GO,class NO>
void NonLinTPM<SC,LO,GO,NO>::reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const
{
    
}
    
template<class SC,class LO,class GO,class NO>
void NonLinTPM<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions){
    
    
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only Newton for NonLinTPM!");
    
}

template<class SC,class LO,class GO,class NO>
void NonLinTPM<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                  const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                                  ) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented! NonLinTPM evalModelImpl(...)");
   
}
 
template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > NonLinTPM<SC,LO,GO,NO>::create_W_op() const
{
    Teuchos::RCP<Thyra::LinearOpBase<SC> > dummy;
    return dummy;
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > NonLinTPM<SC,LO,GO,NO>::create_W_prec() const
{
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > dummy;
    return dummy;
}

}
#endif
