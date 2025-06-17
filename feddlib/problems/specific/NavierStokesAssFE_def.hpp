#ifndef NAVIERSTOKESASSFE_def_hpp
#define NAVIERSTOKESASSFE_def_hpp
#include "NavierStokesAssFE_decl.hpp"

/*!
 Definition of Navier-Stokes

 @brief Navier-Stokes
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

/*void sxOne2D(double* x, double* res, double t, double* parameter){

    res[0] = 1.;
    res[1] = 0.;
    return;
}

void syOne2D(double* x, double* res, double t, double* parameter){

    res[0] = 0.;
    res[1] = 1.;

    return;
}
void sDummyFunc(double* x, double* res, double t, double* parameter){

    return;
}

double OneFunction(double* x, int* parameter)
{
    return 1.0;
}*/

namespace FEDD {



template<class SC,class LO,class GO,class NO>
NavierStokesAssFE<SC,LO,GO,NO>::NavierStokesAssFE( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList ):
NonLinearProblem<SC,LO,GO,NO>( parameterList, domainVelocity->getComm() ),
A_(),
pressureIDsLoc(new vec_int_Type(2)),
u_rep_(),
p_rep_(),
viscosity_element_() // Added here also viscosity field 
{

    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    this->initNOXParameters();

    this->addVariable( domainVelocity , FETypeVelocity , "u" , domainVelocity->getDimension());
    this->addVariable( domainPressure , FETypePressure , "p" , 1);
    this->dim_ = this->getDomain(0)->getDimension();

    u_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );

    p_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );

    if (parameterList->sublist("Parameter").get("Calculate Coefficients",false)) {
        vec2D_dbl_ptr_Type vectmpPointsPressure = domainPressure->getPointsUnique();
        vec2D_dbl_Type::iterator it;
        int front = -1;
        int back = -1;
        if (domainPressure->getDimension() == 2) {
            it = find_if(vectmpPointsPressure->begin(), vectmpPointsPressure->end(),
                    [&] (const std::vector<double>& a){
                        if (a.at(0) >= 0.15-1.e-12 && a.at(0) <= 0.15+1.e-12
                            && a.at(1) >= 0.2-1.e-12 && a.at(1) <= 0.2+1.e-12) {
                            return true;
                        }
                        else {
                            return false;
                        }
                    });

            if (it != vectmpPointsPressure->end()) {
                front = distance(vectmpPointsPressure->begin(),it);
            }
            it = find_if(vectmpPointsPressure->begin(), vectmpPointsPressure->end(),
                    [&] (const std::vector<double>& a){
                        if (a.at(0) >= 0.25-1.e-12 && a.at(0) <= 0.25+1.e-12
                            && a.at(1) >= 0.2-1.e-12 && a.at(1) <= 0.2+1.e-12) {
                            return true;
                        }
                        else {
                            return false;
                        }
                    });

            if (it != vectmpPointsPressure->end()) {
                back = distance(vectmpPointsPressure->begin(),it);
            }
            pressureIDsLoc->at(0) = front;
            pressureIDsLoc->at(1) = back;
        }
        else if(domainPressure->getDimension() == 3){
#ifdef ASSERTS_WARNINGS
            MYASSERT(false,"Not implemented to calc coefficients in 3D!");
#endif
        }
    }

}

template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::info(){
    this->infoProblem();
    this->infoNonlinProblem();
}

template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if (type=="") {
        if (this->verbose_)
            std::cout << "-- Assembly Navier-Stokes ... " << std::endl;

        assembleConstantMatrices();
        
        if (this->verbose_)
            std::cout << "done -- " << std::endl;
    }
    else
        reAssemble( type );

};

template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::assembleConstantMatrices() const{
    
    if (this->verbose_)
        std::cout << "-- Assembly constant matrices Navier-Stokes ... " << std::flush;
    
    double viscosity = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    
    // Egal welcher Wert, da OneFunction nicht von parameter abhaengt
    int* dummy;
    
    A_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    
 	MapConstPtr_Type pressureMap;
    if ( this->getDomain(1)->getFEType() == "P0" )
        pressureMap = this->getDomain(1)->getElementMap();
    else
        pressureMap = this->getDomain(1)->getMapUnique();
    
	if (this->system_.is_null())
        this->system_.reset(new BlockMatrix_Type(2));

	if (this->residualVec_.is_null())
        this->residualVec_.reset(new BlockMultiVector_Type(2));
    
    MatrixPtr_Type B(new Matrix_Type( pressureMap, this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    MatrixPtr_Type BT(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
    MatrixPtr_Type C(new Matrix_Type( pressureMap,1));


	this->system_->addBlock(A_,0,0);
	this->system_->addBlock(BT,0,1);
	this->system_->addBlock(B,1,0);
	this->system_->addBlock(C,1,1);

	this->feFactory_->assemblyNavierStokes(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(1)->getFEType(), 2, this->dim_,1,u_rep_,p_rep_,this->system_,this->residualVec_,this->coeff_, this->parameterList_,false, "Jacobian", true/*call fillComplete*/);

    if ( !this->getFEType(0).compare("P1") ) {
        C.reset(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyBDStabilization( this->dim_, "P1", C, true);
        C->resumeFill();
        C->scale( -1. / ( viscosity * density ) );
        C->fillComplete( pressureMap, pressureMap );
        
        this->system_->addBlock( C, 1, 1 );
    }



    
#ifdef FEDD_HAVE_TEKO
    if ( !this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic").compare("Teko") ) {
        if (!this->parameterList_->sublist("General").get("Assemble Velocity Mass",false)) {
            MatrixPtr_Type Mvelocity(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            //
            this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity, true );
            //
            this->getPreconditionerConst()->setVelocityMassMatrix( Mvelocity );
            if (this->verbose_)
                std::cout << "\nVelocity mass matrix for LSC block preconditioner is assembled." << std::endl;
        } else {
            if (this->verbose_)
                std::cout << "\nVelocity mass matrix for LSC block preconditioner not assembled." << std::endl;
        }
    }
#endif
    std::string precType = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( precType == "Diagonal" || precType == "Triangular" ) {
        MatrixPtr_Type Mpressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        
        this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(1), "Scalar", Mpressure, true );
        SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
        Mpressure->scale(-1./kinVisco);
        this->getPreconditionerConst()->setPressureMassMatrix( Mpressure );
    }
    if (this->verbose_)
        std::cout << " Call Reassemble FixedPoint and Newton to allocate the Matrix pattern " << std::endl;

    // this->reAssemble("FixedPoint");
    this->reAssemble("Newton");
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
    
};
    

template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const
{

}
    
   

template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::reAssemble(std::string type) const {

    
    if (this->verbose_)
        std::cout << "-- Reassembly Navier-Stokes ("<< type <<") ... " << std::flush;
    
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    
    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

    MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
    u_rep_->importFromVector(u, true);

    MultiVectorConstPtr_Type p = this->solution_->getBlock(1);
    p_rep_->importFromVector(p, true); 
   

   if (type=="Rhs") {

   		this->system_->addBlock(ANW,0,0);

        this->feFactory_->assemblyNavierStokes(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(1)->getFEType(), 2, this->dim_,1,u_rep_,p_rep_,this->system_, this->residualVec_,this->coeff_,this->parameterList_, true, "FixedPoint",  true);  // We can also change that to "Jacobian"     
 		this->feFactory_->assemblyNavierStokes(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(1)->getFEType(), 2, this->dim_,1,u_rep_,p_rep_,this->system_, this->residualVec_,this->coeff_,this->parameterList_, true, "Rhs",  true);

    }
    else if (type=="FixedPoint" ) {

   		this->system_->addBlock(ANW,0,0);
		this->feFactory_->assemblyNavierStokes(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(1)->getFEType(), 2, this->dim_,1,u_rep_,p_rep_,this->system_, this->residualVec_,this->coeff_,this->parameterList_, true, "FixedPoint",  true);
    }
	else if(type=="Newton"){ 
        
        this->system_->addBlock(ANW,0,0);
		this->feFactory_->assemblyNavierStokes(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(1)->getFEType(), 2, this->dim_,1,u_rep_,p_rep_,this->system_,this->residualVec_, this->coeff_,this->parameterList_, true,"Jacobian", true);

    }
	
    this->system_->addBlock(ANW,0,0);

    if (this->verbose_)
        std::cout << "done -- " << std::endl;
 
    
}


template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
	//this->reAssemble("FixedPoint");
    this->reAssemble("Rhs");

    // We need to additionally add the residual component for the stabilization block, as it is not part of the AssembleFE routines and calculated externally here.
    if(this->getDomain(0)->getFEType() == "P1"){

        MultiVectorPtr_Type residualPressureTmp = Teuchos::rcp(new MultiVector_Type( this->getDomain(1)->getMapUnique() ));

         this->system_->getBlock(1,1)->apply( *(this->solution_->getBlock(1)), *residualPressureTmp);
           
        this->residualVec_->getBlockNonConst(1)->update(1.,*residualPressureTmp,1.);


    }
    // We need to account for different parameters of time discretizations here
    // This is ok for bdf with 1.0 scaling of the system. Would be wrong for Crank-Nicolson - might be ok now for CN

    /*if (this->coeff_.size() == 0)
        this->system_->apply( *this->solution_, *this->residualVec_ );
    else
        this->system_->apply( *this->solution_, *this->residualVec_, this->coeff_ );*/
    
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
void NavierStokesAssFE<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                              const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                              ) const
{
    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( !type.compare("Monolithic"))
        evalModelImplMonolithic( inArgs, outArgs );
    else if ( !type.compare("Teko")){
#ifdef FEDD_HAVE_TEKO
        evalModelImplBlock( inArgs, outArgs );
#else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Teko not found! Build Trilinos with Teko.");
#endif
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner/solver type.");
}

/*!
	\brief Monolithic Approach for Nonlinear Solver NOX. Input. Includes calculation of the residual vector and update (reAssembly) of non constant matrices with new solution.
		   ResidualVec and SystemMatrix of this class are then converted into the corresponding Thyra/Tpetra objects for Solver.



*/
template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::evalModelImplMonolithic(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                        const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs ) const
{


    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;
    using Teuchos::ArrayView;
    using Teuchos::Array;
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
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
    typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;
    typedef RCP<const TpetraMatrix_Type> TpetraMatrixConstPtr_Type;
 
    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);


    if ( fill_f || fill_W || fill_W_prec ) {

        // ****************
        // Get the underlying xpetra objects
        // ****************
        if (fill_f) {

            this->calculateNonLinResidualVec("standard"); // Calculating residual Vector

			// Changing the residualVector into a ThyraMultivector

            Teuchos::RCP<Thyra::MultiVectorBase<SC> > f_thyra = this->getResidualVector()->getThyraMultiVector();
            f_out->assign(*f_thyra);
        }
        TpetraMatrixPtr_Type W;
        if (fill_W) {
            this->reAssemble("Newton"); // ReAssembling matrices with updated u  in this class

            this->setBoundariesSystem(); // setting boundaries to the system

			Teuchos::RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
            Teuchos::RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);
            
            TpetraMatrixConstPtr_Type W_systemTpetra = this->getSystem()->getMergedMatrix()->getTpetraMatrix();           
            TpetraMatrixPtr_Type W_systemTpetraNonConst = rcp_const_cast<TpetraMatrix_Type>(W_systemTpetra);
            
            //Tpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_systemXpetraNonConst);
            //Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            
            Teuchos::RCP<TpetraMatrix_Type> tpetraMatTpetra = W_systemTpetraNonConst; //xTpetraMat.getTpetra_CrsMatrixNonConst();
            W_tpetraMat->resumeFill();

            for (auto i=0; i<tpetraMatTpetra->getMap()->getLocalNumElements(); i++) {
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;
                tpetraMatTpetra->getLocalRowView( i, indices, values);
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
/*!
	\brief Block Approach for Nonlinear Solver NOX. Input. Includes calculation of the residual vector and update (reAssembly) of non constant matrices with new solution.
		   ResidualVec and SystemMatrix of this class are then converted into the corresponding Thyra/Tpetra objects for Solver.



*/
#ifdef FEDD_HAVE_TEKO
template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::evalModelImplBlock(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                   const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs ) const
{


    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;
    using Teuchos::ArrayView;
    using Teuchos::Array;

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    TEUCHOS_TEST_FOR_EXCEPTION( inArgs.get_x().is_null(), std::logic_error, "inArgs.get_x() is null.");

    RCP< const Thyra::VectorBase< SC > > vecThyra = inArgs.get_x();
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    RCP< Thyra::VectorBase< SC > > vecThyraNonConst = rcp_const_cast<Thyra::VectorBase< SC > >(vecThyra);

    RCP< Thyra::ProductVectorBase< SC > > vecThyraBlock = rcp_dynamic_cast<Thyra::ProductVectorBase< SC > > (vecThyraNonConst);

    this->solution_->getBlockNonConst(0)->fromThyraMultiVector( vecThyraBlock->getNonconstVectorBlock(0) );
    this->solution_->getBlockNonConst(1)->fromThyraMultiVector( vecThyraBlock->getNonconstVectorBlock(1) );

    const RCP<Thyra::MultiVectorBase<SC> > f_out = outArgs.get_f();
    const RCP<Thyra::LinearOpBase<SC> > W_out = outArgs.get_W_op();
    const RCP<Thyra::PreconditionerBase<SC> > W_prec_out = outArgs.get_W_prec();

    typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
    typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;
    typedef RCP<const TpetraMatrix_Type> TpetraMatrixConstPtr_Type;

    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);

    if ( fill_f || fill_W || fill_W_prec ) {

        // ****************
        // Get the underlying xpetra objects
        // ****************
        if (fill_f) {

            this->calculateNonLinResidualVec("standard");

            Teko::MultiVector f0;
            Teko::MultiVector f1;
            f0 = this->getResidualVector()->getBlockNonConst(0)->getThyraMultiVector();
            f1 = this->getResidualVector()->getBlockNonConst(1)->getThyraMultiVector();

            std::vector<Teko::MultiVector> f_vec; f_vec.push_back(f0); f_vec.push_back(f1);

            Teko::MultiVector f = Teko::buildBlockedMultiVector(f_vec);

            f_out->assign(*f);
        }

        TpetraMatrixPtr_Type W;
        if (fill_W) {

            typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;

            this->reAssemble("Newton");

            this->setBoundariesSystem();

            RCP<ThyraBlockOp_Type> W_blocks = rcp_dynamic_cast<ThyraBlockOp_Type>(W_out);
            RCP<const ThyraOp_Type> W_block00 = W_blocks->getBlock(0,0);
            RCP<ThyraOp_Type> W_block00NonConst = rcp_const_cast<ThyraOp_Type>( W_block00 );
            RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator( W_block00NonConst );

            RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);

            TpetraMatrixConstPtr_Type W_matrixTpetra = this->getSystem()->getBlock(0,0)->getTpetraMatrix();
            TpetraMatrixPtr_Type W_matrixTpetraNonConst = rcp_const_cast<TpetraMatrix_Type>(W_matrixTpetra);
            RCP<TpetraMatrix_Type> tpetraMatTpetra = W_matrixTpetraNonConst;

            W_tpetraMat->resumeFill();

            for (auto i=0; i<tpetraMatTpetra->getMap()->getLocalNumElements(); i++) {
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;
                tpetraMatTpetra->getLocalRowView( i, indices, values);
                W_tpetraMat->replaceLocalValues( i, indices, values);
            }
            W_tpetraMat->fillComplete();

        }

        if (fill_W_prec) {
            if (stokesTekoPrecUsed_){
                this->setupPreconditioner( "Teko" );
            }
            else
                stokesTekoPrecUsed_ = true;

            // ch 26.04.19: After each setup of the preconditioner we check if we use a two-level precondtioner with multiplicative combination between the levels.
            // If this is the case, we need to pre apply the coarse level to the residual(f_out).

            ParameterListPtr_Type tmpSubList = sublist( sublist( sublist( sublist( this->parameterList_, "Teko Parameters" ) , "Preconditioner Types" ) , "Teko" ) , "Inverse Factory Library" );

            std::string levelCombination1 = tmpSubList->sublist( "FROSch-Velocity" ).get("Level Combination","Additive");
            std::string levelCombination2 = tmpSubList->sublist( "FROSch-Pressure" ).get("Level Combination","Additive");

            if ( !levelCombination1.compare("Multiplicative") || !levelCombination2.compare("Multiplicative") ) {

                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Multiplicative Level Combination is not supported for NOX.");
                ParameterListPtr_Type solverPList = this->getLinearSolverBuilder()->getNonconstParameterList();

//                    pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);
//
//                    Teuchos::RCP<const Thyra::LinearOpBase<SC> > thyra_linOp = this->getPreconditioner()->getThyraPrec()->getUnspecifiedPrecOp();
//                    Thyra::apply( *thyra_linOp, Thyra::NOTRANS, *thyraB, thyraX.ptr() );
//                    pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",false);


            }
        }
    }
}
#endif

template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::calculateNonLinResidualVecWithMeshVelo(std::string type, double time, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const{


   // this->reAssembleFSI( "FixedPoint", u_minus_w, P );
    
    // We need to account for different parameters of time discretizations here
    // This is ok for bdf with 1.0 scaling of the system. Would be wrong for Crank-Nicolson
    
    this->system_->apply( *this->solution_, *this->residualVec_ );
//    this->residualVec_->getBlock(0)->writeMM("Ax.mm");
//    this->rhs_->getBlock(0)->writeMM("nsRHS.mm");
    if (!type.compare("standard")){
        this->residualVec_->update(-1.,*this->rhs_,1.);
        if ( !this->sourceTerm_.is_null() )
            this->residualVec_->update(-1.,*this->sourceTerm_,1.);
    }
    else if(!type.compare("reverse")){
        this->residualVec_->update(1.,*this->rhs_,-1.); // this = -1*this + 1*rhs
        if ( !this->sourceTerm_.is_null() )
            this->residualVec_->update(1.,*this->sourceTerm_,1.);
    }
    
    // this might be set again by the TimeProblem after addition of M*u
    this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    
//    this->residualVec_->getBlock(0)->writeMM("b_Ax.mm");
    
}

    
template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > NavierStokesAssFE<SC,LO,GO,NO>::create_W_op() const
{
    //this->reAssemble("FixedPoint");
    //this->reAssemble("Newton");

    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( !type.compare("Monolithic"))
        return create_W_op_Monolithic( );
    else if ( !type.compare("Teko")){
#ifdef FEDD_HAVE_TEKO
        return create_W_op_Block( );
#else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Teko not found! Build Trilinos with Teko.");
#endif
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner/solver type.");
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > NavierStokesAssFE<SC,LO,GO,NO>::create_W_op_Monolithic() const
{
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = this->system_->getThyraLinOp();
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
    return W_op;
}

#ifdef FEDD_HAVE_TEKO
template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > NavierStokesAssFE<SC,LO,GO,NO>::create_W_op_Block() const
{

    Teko::LinearOp thyraF = this->system_->getBlock(0,0)->getThyraLinOp();
    Teko::LinearOp thyraBT = this->system_->getBlock(0,1)->getThyraLinOp();
    Teko::LinearOp thyraB = this->system_->getBlock(1,0)->getThyraLinOp();

    if (!this->system_->blockExists(1,1)){
        MatrixPtr_Type dummy = Teuchos::rcp( new Matrix_Type( this->system_->getBlock(1,0)->getMap(), 1 ) );
        dummy->fillComplete();
        this->system_->addBlock( dummy, 1, 1 );
    }

    Teko::LinearOp thyraC = this->system_->getBlock(1,1)->getThyraLinOp();

    Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = Thyra::block2x2(thyraF,thyraBT,thyraB,thyraC);
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
    return W_op;
}
#endif

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > NavierStokesAssFE<SC,LO,GO,NO>::create_W_prec() const
{

    this->initializeSolverBuilder();

    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    this->setBoundariesSystem();

    if (!type.compare("Teko")) { //
        this->setupPreconditioner( type );
        stokesTekoPrecUsed_ = false;
    }
    else{
        this->setupPreconditioner( type );
    }

    Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  this->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst = Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);

    return thyraPrecNonConst;

}
/*template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::reAssembleFSI(std::string type, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const {
    
    if (this->verbose_)
        std::cout << "-- Reassembly Navier-Stokes ("<< type <<") for FSI ... " << std::flush;
    
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);

    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    if (type=="FixedPoint") {
        
        MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyAdvectionVecField( this->dim_, this->domain_FEType_vec_.at(0), N, u_minus_w, true );
        
        N->resumeFill();
        N->scale(density);
        N->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
        A_->addMatrix(1.,ANW,0.);

        N->addMatrix(1.,ANW,1.);
        // P must be scaled correctly in FSI
        P->addMatrix(1.,ANW,1.);


    }
    else if(type=="Newton"){
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "reAssembleFSI should only be called for FPI-System.");
    }
    ANW->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );

    this->system_->addBlock( ANW, 0, 0 );

    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}*/


template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions){

    if (this->verbose_)
        std::cout << "-- Reassembly Navier-Stokes (Extrapolation) ... " << std::flush;

    
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

    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

    MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyAdvectionVecField( this->dim_, this->domain_FEType_vec_.at(0), N, u_rep_, true );

    N->resumeFill();
    N->scale(density);
    N->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());

    A_->addMatrix(1.,ANW,0.);
    N->addMatrix(1.,ANW,1.);

    ANW->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );

    this->system_->addBlock( ANW, 0, 0 );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}






// Use converged solution to compute viscosity field
template<class SC,class LO,class GO,class NO>
void NavierStokesAssFE<SC,LO,GO,NO>::computeSteadyPostprocessingViscosity_Solution() {
    

    MultiVectorConstPtr_Type u = this->solution_->getBlock(0); // solution_ is initialized in problem_def.hpp so the most general class
    u_rep_->importFromVector(u, true); // this is the current velocity solution at the nodes - distributed at the processors with repeated values
   
    MultiVectorConstPtr_Type p = this->solution_->getBlock(1);
    p_rep_->importFromVector(p, true);  // this is the current pressure solution at the nodes - distributed at the processors with repeated values
    
    // Reset here the viscosity so at this moment this makes only sense to call at the end of a simulation
    // to visualize viscosity field
    // @ToDo Add possibility for transient problem to save viscosity solution in each time step
    viscosity_element_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getElementMap() ) );
    this->feFactory_->computeSteadyViscosityFE_CM(this->dim_, this->getDomain(0)->getFEType(), this->getDomain(1)->getFEType(), this->dim_,1,u_rep_,p_rep_,this->parameterList_);        
  
    Teuchos::RCP<const MultiVector<SC,LO,GO,NO>> exportSolutionViscosityAssFE = this->feFactory_->const_output_fields->getBlock(0); // For now we assume that viscosity is always saved in first block
    viscosity_element_->importFromVector(exportSolutionViscosityAssFE, true);  
  
  

}







}

#endif
