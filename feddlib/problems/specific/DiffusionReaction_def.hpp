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

    MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
    // Here we insert the assembly of the reaction part.
    
    vec_dbl_Type param(0);
    for(int i=0; i< funcParameter_.size(); i++){
        param.push_back(funcParameter_[i]);
    }

    this->feFactory_->assemblyLinearReactionTerm( this->dim_, this->domain_FEType_vec_.at(0), N,  true, param, reactionFunc_);
    N->scale(-1.0);
    MatrixPtr_Type AN = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );

    N->addMatrix(1.,AN,0.);
    A_->addMatrix(1.,AN,0.);

    AN->fillComplete( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getMapUnique() );

    if (this->system_.is_null())
        this->system_.reset(new BlockMatrix_Type(1));

    this->system_->addBlock(AN,0,0);
    
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

}

/*template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const
{

}

// This needs to be adjusted according to the nonlinear reaction part

template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::reAssemble(std::string type) const {

    
    if (this->verbose_)
        std::cout << "-- Reassembly Reaction-Diffusion ("<< type <<") ... " << std::flush;
    
	// Sensible input is the reaction function. Might distinguish between linear and nonlinear Reaction term.
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    
    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );

    vec_dbl_Type param(0);
		for(int i=0; i< funcParameter_.size(); i++){
			param.push_back(funcParameter_[i]);
		}
    if (type=="FixedPoint") {
        
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);

        u_rep_->importFromVector(u, true);

        MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
		// Here we insert the assembly of the reaction part.
		
        this->feFactory_->assemblyReactionTerm( this->dim_, this->domain_FEType_vec_.at(0), N, u_rep_, true, param, reactionFunc_);
        
		// Adding A to ANW
        A_->addMatrix(1.,ANW,0.);
		// Addind N to ANW
        N->addMatrix(1.,ANW,1.);
    }
    else if(type=="Newton"){ // We assume that reAssmble("FixedPoint") was already called for the current iterate
        MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyDReactionTerm( this->dim_, this->domain_FEType_vec_.at(0), W, u_rep_, true, param, reactionFunc_);
        W->resumeFill();

        W->fillComplete( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getMapUnique());
        
        A_->addMatrix(1.,ANW,0.);
        W->addMatrix(1.,ANW,1.);

    }
    ANW->fillComplete( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getMapUnique() );
    
    this->system_->addBlock( ANW, 0, 0 );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions){

    if (this->verbose_)
        std::cout << "-- Reassembly Reaction Diffusion (Extrapolation) ... " << std::flush;

    
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);

    if (previousSolutions.size()>=2) {

        MultiVectorPtr_Type extrapolatedVector = Teuchos::rcp( new MultiVector_Type( previousSolutions[0]->getBlock(0) ) );

        extrapolatedVector->update( -1., *previousSolutions[1]->getBlock(0), 2. );

        u_rep_->importFromVector(extrapolatedVector, true);
    }
    else if(previousSolutions.size()==1){
        MultiVectorConstPtr_Type u = previousSolutions[0]->getBlock(0);
        u_rep_->importFromVector(u, true);
    }
    else if (previousSolutions.size()==0){
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true);
    }

    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

    MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

	vec_dbl_Type param(0);
	for(int i=0; i< funcParameter_.size(); i++){
		param.push_back(funcParameter_[i]);
	}
    this->feFactory_->assemblyReactionTerm( this->dim_, this->domain_FEType_vec_.at(0), N, u_rep_, true, param, reactionFunc_);


    N->resumeFill();
    N->scale(density);
    N->fillComplete( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getMapUnique());

    A_->addMatrix(1.,ANW,0.);
    N->addMatrix(-1.,ANW,1.);

    ANW->fillComplete( this->getDomain(0)->getMapUnique(), this->getDomain(0)->getMapUnique());

    this->system_->addBlock( ANW, 0, 0 );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
    this->reAssemble("FixedPoint");

    // We need to account for different parameters of time discretizations here
    // This is ok for bdf with 1.0 scaling of the system. Would be wrong for Crank-Nicolson - might be ok now for CN
    if (this->coeff_.size() == 0)
        this->system_->apply( *this->solution_, *this->residualVec_ );
    else
        this->system_->apply( *this->solution_, *this->residualVec_, this->coeff_ );
    

    if (!type.compare("standard")){
        this->residualVec_->update(-1.,*this->rhs_,1.);
//        if ( !this->sourceTerm_.is_null() )
//            this->residualVec_->update(-1.,*this->sourceTerm_,1.);
        // this might be set again by the TimeProblem after addition of M*u
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
        
    }
    else if(!type.compare("reverse")){
        this->residualVec_->update(1.,*this->rhs_,-1.); // this = -1*this + 1*rhs
//        if ( !this->sourceTerm_.is_null() )
//            this->residualVec_->update(1.,*this->sourceTerm_,1.);
        // this might be set again by the TimeProblem after addition of M*u
        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
        
    }
}

template<class SC,class LO,class GO,class NO>
void DiffusionReaction<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                        const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs ) const
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
        bool system_updated = false;
        if (fill_f) {

            this->calculateNonLinResidualVec("standard");

            system_updated = true;
            Teuchos::RCP<Thyra::MultiVectorBase<SC> > f_thyra = this->getResidualVector()->getThyraMultiVector();
            f_out->assign(*f_thyra);
        }

        XpetraMatrixPtr_Type W;
        if (fill_W) {

            this->reAssemble("Newton");

            this->setBoundariesSystem();

            Teuchos::RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
            Teuchos::RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);

            XpetraMatrixConstPtr_Type W_systemXpetra = this->getSystem()->getBlock( 0, 0 )->getXpetraMatrix();

            XpetraMatrixPtr_Type W_systemXpetraNonConst = rcp_const_cast<XpetraMatrix_Type>(W_systemXpetra);
            Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_systemXpetraNonConst);
            Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            Teuchos::RCP<TpetraMatrix_Type> tpetraMatXpetra = xTpetraMat.getTpetra_CrsMatrixNonConst();

            W_tpetraMat->resumeFill();

            for (auto i=0; i<tpetraMatXpetra->getMap()->getNodeNumElements(); i++) {
                ArrayView< const LO > indices;
                ArrayView< const SC > values;
                tpetraMatXpetra->getLocalRowView( i, indices, values);
                W_tpetraMat->replaceLocalValues( i, indices.size(), values.getRawPtr(), indices.getRawPtr() );
            }
            W_tpetraMat->fillComplete();

        }

        if (fill_W_prec) {
            this->setupPreconditioner( "Monolithic" );

            // ch 26.04.19: After each setup of the preconditioner we check if we use a two-level precondtioner with multiplicative combination between the levels.
            // If this is the case, we need to pre apply the coarse level to the residual(f_out).

            std::string levelCombination = this->parameterList_->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive");
            if (!levelCombination.compare("Multiplicative")) {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Multiplicative Level Combination is not supported for NOX.");
//                ParameterListPtr_Type solverPList = this->getLinearSolverBuilder()->getNonconstParameterList();
//
//                solverPList->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);
//
//                Teuchos::RCP<const Thyra::LinearOpBase<SC> > thyra_linOp = this->getPreconditionerConst()->getThyraPrecConst()->getUnspecifiedPrecOp();
//
//                f_out->describe(*out,Teuchos::VERB_EXTREME);
//                vecThyraNonConst->describe(*out,Teuchos::VERB_EXTREME);
//                Thyra::apply( *thyra_linOp, Thyra::NOTRANS, *f_out, vecThyraNonConst.ptr() );
//                solverPList->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",false);
            }

        }
    }
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > DiffusionReaction<SC,LO,GO,NO>::create_W_op() const
{   
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = this->system_->getThyraLinOp();
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
    return W_op;
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > DiffusionReaction<SC,LO,GO,NO>::create_W_prec() const
{
    this->initializeSolverBuilder();
    this->initializePreconditioner();
    
    Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  this->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst= Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);
    
    return thyraPrecNonConst;
}
*/


template<class SC,class LO,class GO,class NO>
typename DiffusionReaction<SC,LO,GO,NO>::MatrixPtr_Type DiffusionReaction<SC,LO,GO,NO>::getMassMatrix() const{
	
    MatrixPtr_Type A;
	A = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
	this->feFactory_->assemblyMass(this->dim_,this->domain_FEType_vec_.at(0),"Scalar", A);

	return A;

}
    
}
#endif
