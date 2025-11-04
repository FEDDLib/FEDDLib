#ifndef LaplaceBlocks_def_hpp
#define LaplaceBlocks_def_hpp

/*!
 Definition of LaplaceBlocks
 
 @brief LaplaceBlocks
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"


namespace FEDD {

template<class SC,class LO,class GO,class NO>
LaplaceBlocks<SC,LO,GO,NO>::LaplaceBlocks( const DomainConstPtr_Type &domain1, const DomainConstPtr_Type &domain2, std::string FEType1, std::string FEType2, ParameterListPtr_Type parameterList ):
Problem<SC,LO,GO,NO>(parameterList, domain1->getComm())
{
 
    this->addVariable( domain1 , FEType1 , "u1" , 1);
    this->addVariable( domain2 , FEType2 , "u2" , 1);
    this->dim_ = this->getDomain(0)->getDimension();
}

template<class SC,class LO,class GO,class NO>
LaplaceBlocks<SC,LO,GO,NO>::~LaplaceBlocks(){

}
    
template<class SC,class LO,class GO,class NO>
void LaplaceBlocks<SC,LO,GO,NO>::info(){
    this->infoProblem();
}

template<class SC,class LO,class GO,class NO>
void LaplaceBlocks<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if (this->verbose_)
        std::cout << "-- Assembly LaplaceBlocks ... " << std::flush;

    MatrixPtr_Type A1;
    MatrixPtr_Type A2;

   
    A1 = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyLaplace(this->dim_, this->domain_FEType_vec_.at(0), 2, A1, true, 0);

    A2 = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(1)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyLaplace(this->dim_, this->domain_FEType_vec_.at(1), 2, A2, true , 1);

    this->system_->addBlock( A1, 0, 0 );
    this->system_->addBlock( A2, 1, 1 );
        
    this->assembleSourceTerm( 0. );
    
    this->addToRhs( this->sourceTerm_ );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

//template<class SC,class LO,class GO,class NO>
//void LaplaceBlocks<SC,LO,GO,NO>::assembleSourceTerm(double time){
//
//    MultiVectorPtr_Type FERhs1;
//    MultiVectorPtr_Type FERhs2;
//
//    vec_dbl_Type funcParameter(1,0.);
//    FERhs1 = Teuchos::rcp(new MultiVector_Type( this->domainPtr_vec_.at(0)->getMapRepeated() ) );
//    this->feFactory_->assemblyRHS(this->dim_,
//                                  this->domain_FEType_vec_.at(0),
//                                  FERhs1,
//                                  "Scalar",
//                                  this->rhsFuncVec_[0],
//                                  funcParameter,
//                                  0);
//
//    FERhs2 = Teuchos::rcp(new MultiVector_Type( this->domainPtr_vec_.at(1)->getMapRepeated() ) );
//    this->feFactory_->assemblyRHS(this->dim_,
//                                  this->domain_FEType_vec_.at(1),
//                                  FERhs2,
//                                  "Scalar",
//                                  this->rhsFuncVec_[0],
//                                  funcParameter,
//                                  1);
//
//
//    this->sourceTerm_->getBlockNonConst(0)->putScalar(0.);
//
//    this->sourceTerm_->getBlockNonConst(0)->exportFromVector( FERhs1, false, "Add" );
//
//    this->sourceTerm_->getBlockNonConst(1)->putScalar(0.);
//
//    this->sourceTerm_->getBlockNonConst(1)->exportFromVector( FERhs2, false, "Add" );
//
//}
    
//template<class SC,class LO,class GO,class NO>
//void LaplaceBlocks<SC,LO,GO,NO>::assembleRhs(double time){
//
//    MultiVectorPtr_Type FERhs1;
//    MultiVectorPtr_Type FERhs2;
//
//    vec_dbl_Type funcParameter(1,0.);
//    FERhs1 = Teuchos::rcp(new MultiVector_Type( this->domainPtr_vec_.at(0)->getMapRepeated() ) );
//    this->feFactory_->assemblyRHS(this->dim_,
//                                  this->domain_FEType_vec_.at(0),
//                                  FERhs1,
//                                  "Scalar",
//                                  this->rhsFuncVec_[0],
//                                  funcParameter,
//                                  0);
//    
//    FERhs2 = Teuchos::rcp(new MultiVector_Type( this->domainPtr_vec_.at(1)->getMapRepeated() ) );
//    this->feFactory_->assemblyRHS(this->dim_,
//                                  this->domain_FEType_vec_.at(1),
//                                  FERhs2,
//                                  "Scalar",
//                                  this->rhsFuncVec_[0],
//                                  funcParameter,
//                                  1);
//    
//    
//    this->sourceTerm_->getBlockNonConst(0)->putScalar(0.);
//    
//    this->sourceTerm_->getBlockNonConst(0)->exportFromVector( FERhs1, false, "Add" );
//
//    this->sourceTerm_->getBlockNonConst(1)->putScalar(0.);
//    
//    this->sourceTerm_->getBlockNonConst(1)->exportFromVector( FERhs2, false, "Add" );
//    
//}

//template<class SC,class LO,class GO,class NO>
//int LaplaceBlocks<SC,LO,GO,NO>::SetupPreconditioner(BMat_ptr_Type systemPrec){
//
//#ifdef ASSERTS_WARNINGS
//    MYASSERT(this->boolProblemAssembled_,"Problem not assembled!");
//#endif
//
//    Teuchos::RCP<PrecondManagerFROSch<SC,LO,GO,NO> > precondManager(new PrecondManagerFROSch<SC,LO,GO,NO>());
//
//    precondManager->BuildPreconditioner( this->XpetraPrec_, this->XpetraMatrix_, systemPrec, this->domainPtr_vec_.at(0), this->bcFactory_ , this->parameterList_ , this->comm_);
//    
//    
//    return 0;
//};

//template<class SC,class LO,class GO,class NO>
//int LaplaceBlocks<SC,LO,GO,NO>::SetupPreconditioner( BMat_ptr_Type systemPrec, ThyraConstLinOpPtr_Type thyraMatrix, ThyraPrecPtr_Type thyraPreconditioner, LinSolverBuilderPtr_Type linearSolverBuilder ) const{
//    
//#ifdef ASSERTS_WARNINGS
//    MYASSERT(this->boolProblemAssembled_,"Problem not assembled!");
//#endif
//    
//    Teuchos::RCP<PrecondManagerFROSch<SC,LO,GO,NO> > precondManager(new PrecondManagerFROSch<SC,LO,GO,NO>());
//    
//    if (!thyraPreconditioner.is_null() && !thyraMatrix.is_null() && !linearSolverBuilder.is_null()) {
//        precondManager->BuildPreconditioner( thyraPreconditioner, thyraMatrix, linearSolverBuilder ,systemPrec, this->domainPtr_vec_, this->bcFactory_ , this->parameterList_ , this->PListThyraSolver_, this->comm_, this->PrecondtionerIsBuilt_);
//    }
//    else{
//        precondManager->BuildPreconditioner( this->ThyraPreconditioner_, this->ThyraMatrix_, this->LinearSolverBuilder_ ,systemPrec, this->domainPtr_vec_, this->bcFactory_ , this->parameterList_ , this->PListThyraSolver_, this->comm_, this->PrecondtionerIsBuilt_);
//    }
//
//    
//    
//    return 0;
//};
//
//template<class SC,class LO,class GO,class NO>
//void LaplaceBlocks<SC,LO,GO,NO>::initSolutionWithBC(){
//    
////    this->bcFactory_->SetRHS(this->system_, this->solution_, 0.);
//    
//}
}
#endif
