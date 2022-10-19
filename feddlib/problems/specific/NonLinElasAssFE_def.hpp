#ifndef NonLinElasAssFE_def_hpp
#define NonLinElasAssFE_def_hpp
#include "NonLinElasAssFE_decl.hpp"
/*!
 Definition of NonLinElasticity
 
 @brief NonLinElasticity
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


namespace FEDD {
template<class SC,class LO,class GO,class NO>
NonLinElasAssFE<SC,LO,GO,NO>::NonLinElasAssFE(const DomainConstPtr_Type  &domain, std::string FEType, ParameterListPtr_Type parameterList):
NonLinearProblem<SC,LO,GO,NO>( parameterList, domain->getComm() ),
u_rep_()
{
    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    this->initNOXParameters();
    this->addVariable( domain , FEType , "u" , domain->getDimension());
    this->dim_ = this->getDomain(0)->getDimension();
    
    C_ = this->parameterList_->sublist("Parameter").get("C",1.);
    
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    eModVec_.reset(new MultiVector_Type(domain->getElementMap(),1));
}

template<class SC,class LO,class GO,class NO>
NonLinElasAssFE<SC,LO,GO,NO>::~NonLinElasAssFE(){

}

template<class SC,class LO,class GO,class NO>
void NonLinElasAssFE<SC,LO,GO,NO>::info(){
    this->infoProblem();
    this->infoNonlinProblem();
}
    
template<class SC,class LO,class GO,class NO>
void NonLinElasAssFE<SC,LO,GO,NO>::assemble(std::string type) const{
    
    if (type == ""){
        if (this->verbose_)
            std::cout << "-- Assembly nonlinear elasticity ... " << std::flush;

        MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        
        this->feFactory_->assemblyEmptyMatrix(A);

        this->system_.reset(new BlockMatrix_Type(1));
        this->system_->addBlock( A, 0, 0 );
                
        double density = this->parameterList_->sublist("Parameter").get("Density",1000.);
        string sourceType = 	this->parameterList_->sublist("Parameter").get("Source Type","volume");

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
void NonLinElasAssFE<SC,LO,GO,NO>::reAssemble(std::string type) const {
    std::string material_model = this->parameterList_->sublist("Parameter").get("Material model","Neo-Hooke");

    if (this->verbose_)
        std::cout << "-- Reassembly nonlinear elasticity with material model " << material_model <<" ("<<type <<") ... " << std::flush;
    
    if (type=="Newton-Residual") {
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);

        u_rep_->importFromVector(u, true);
        
        MultiVectorPtr_Type f = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
        f->putScalar(0.);
        MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->system_->addBlock( W, 0, 0 );  
        MultiVectorPtr_Type fRep = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) ); 
        fRep->putScalar(0.);   
        this->residualVec_->addBlock( fRep, 0 ); 
        if(this->parameterList_->sublist("Parameter").get("SCI",false) == true)
       		 this->feFactory_->assemblyNonLinearElasticity(this->dim_, this->getDomain(0)->getFEType(),2, this->dim_, u_rep_, this->system_, this->residualVec_, this->parameterList_,this->getDomain(0), this->eModVec_,true);
       	else 
 			this->feFactory_->assemblyNonLinearElasticity(this->dim_, this->getDomain(0)->getFEType(),2, this->dim_, u_rep_, this->system_, this->residualVec_, this->parameterList_,true);
        fRep->scale(-1.);
        
        MultiVectorPtr_Type fUnique = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
        fUnique->putScalar(0.);
        fUnique->exportFromVector( fRep, true, "Add" );

        this->residualVec_->addBlock( fUnique, 0 );

 //        MultiVectorPtr_Type f = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
//        this->feFactory_->assemblyElasticityStressesAceFEM(this->dim_, this->getDomain(0)->getFEType(), f, u_rep_, material_model, E_, nu_, C_);
//        MultiVectorPtr_Type fUnique = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
//        fUnique->putScalar(0.);
//        fUnique->importFromVector( f, true, "Add" );
//        this->residualVec_->addBlock( fUnique, 0 );
    }
    else if(type=="Newton"){ //we already assemble the new tangent when we calculate the stresses above
        
//        MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), 10 ) );
//        this->feFactory_->assemblyElasticityJacobianAceFEM(this->dim_, this->getDomain(0)->getFEType(), W, u_rep_, material_model, E_, nu_, C_);
//        this->system_->addBlock( W, 0, 0 );
    }
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}
    
template<class SC,class LO,class GO,class NO>
void NonLinElasAssFE<SC,LO,GO,NO>::reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const
{
    
}

template<class SC,class LO,class GO,class NO>
void NonLinElasAssFE<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions){


    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only Newton/NOX implemented for nonlinear material models!");

}

template<class SC,class LO,class GO,class NO>
void NonLinElasAssFE<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                            const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                            ) const
{
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;
    using Teuchos::ArrayView;
    using Teuchos::Array;
    
    TEUCHOS_TEST_FOR_EXCEPTION( this->solution_->getBlock(0)->getMap()->getUnderlyingLib() != "Tpetra", std::runtime_error, "Use of NOX only supports Tpetra. Epetra support must be implemented.");
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    TEUCHOS_TEST_FOR_EXCEPTION( inArgs.get_x().is_null(), std::logic_error, "inArgs.get_x() is null.");
    
    RCP< const Thyra::VectorBase< SC > > vecThyra = inArgs.get_x();
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    RCP< Thyra::VectorBase< SC > > vecThyraNonConst = rcp_const_cast<Thyra::VectorBase< SC > >(vecThyra);
    
    this->solution_->fromThyraMultiVector(vecThyraNonConst);
    
    const RCP<Thyra::MultiVectorBase<SC> > f_out = outArgs.get_f();
    const RCP<Thyra::LinearOpBase<SC> > W_out = outArgs.get_W_op();
    const RCP<Thyra::PreconditionerBase<SC> > W_prec_out = outArgs.get_W_prec();
    
    typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;
    typedef Xpetra::Matrix<SC,LO,GO,NO> XpetraMatrix_Type;
    typedef RCP<XpetraMatrix_Type> XpetraMatrixPtr_Type;
    typedef RCP<const XpetraMatrix_Type> XpetraMatrixConstPtr_Type;
    
    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);
    
    if ( fill_f || fill_W || fill_W_prec ) {
        
        // ****************
        // Get the underlying xpetra objects
        // ****************
        if (fill_f) {
            
            this->calculateNonLinResidualVec("standard");
            
            RCP<Thyra::MultiVectorBase<SC> > f_thyra = this->getResidualVector()->getThyraMultiVector();
            f_out->assign(*f_thyra);
        }
        
        XpetraMatrixPtr_Type W;
        if (fill_W) {
            
            this->reAssemble("Newton");
            
            this->setBoundariesSystem();
            
            RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
            RCP<TpetraMatrix_Type> W_tpetraMat = rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);
            
            XpetraMatrixConstPtr_Type W_systemXpetra = this->getSystem()->getBlock( 0, 0 )->getXpetraMatrix();
            
            XpetraMatrixPtr_Type W_systemXpetraNonConst = rcp_const_cast<XpetraMatrix_Type>(W_systemXpetra);
            Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_systemXpetraNonConst);
            Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            Teuchos::RCP<TpetraMatrix_Type> tpetraMatXpetra = xTpetraMat.getTpetra_CrsMatrixNonConst();
            
            W_tpetraMat->resumeFill();
            
            for (auto i=0; i<tpetraMatXpetra->getMap()->getLocalNumElements(); i++) {
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;
                tpetraMatXpetra->getLocalRowView( i, indices, values);
                W_tpetraMat->replaceLocalValues( i, indices, values);
            }
            W_tpetraMat->fillComplete();
            
        }
        
        if (fill_W_prec) {
            this->setupPreconditioner( "Monolithic" );
            
            // ch 26.04.19: After each setup of the preconditioner we check if we use a two-level precondtioner with multiplicative combination between the levels.
            // If this is the case, we need to pre apply the coarse level to the residual(f_out).
            
            std::string levelCombination = this->parameterList_->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive");
            if (!levelCombination.compare("Multiplicative")) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Multiplicative Level Combination is not supported for NOX. In general we need to pre-apply the coarse problem. This must be implemented here.");
            }
            
        }
    }
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > NonLinElasAssFE<SC,LO,GO,NO>::create_W_op() const
{
    
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = this->system_->getThyraLinOp();
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
    return W_op;
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > NonLinElasAssFE<SC,LO,GO,NO>::create_W_prec() const
{
    this->initializeSolverBuilder();
    this->initializePreconditioner();
    
    Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  this->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst= Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);
    
    return thyraPrecNonConst;
}

template<class SC,class LO,class GO,class NO>
void NonLinElasAssFE<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
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
