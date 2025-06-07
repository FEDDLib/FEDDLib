#ifndef NAVIERSTOKES_def_hpp
#define NAVIERSTOKES_def_hpp
#include "NavierStokes_decl.hpp"

/*!
 Definition of Navier-Stokes

 @brief Navier-Stokes
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

void sxOne2D(double* x, double* res, double t, double* parameter){

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
}

namespace FEDD {



template<class SC,class LO,class GO,class NO>
NavierStokes<SC,LO,GO,NO>::NavierStokes( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList ):
NonLinearProblem<SC,LO,GO,NO>( parameterList, domainVelocity->getComm() ),
A_(),
pressureIDsLoc(new vec_int_Type(2)),
u_rep_()
{

    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    this->initNOXParameters();

    this->addVariable( domainVelocity , FETypeVelocity , "u" , domainVelocity->getDimension());
    this->addVariable( domainPressure , FETypePressure , "p" , 1);
    this->dim_ = this->getDomain(0)->getDimension();

    u_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );

    if (parameterList->sublist("Parameter").get("Calculate Coefficients",false)) {
        vec2D_dbl_ptr_Type vectmpPointsPressure = domainPressure->getPointsUnique();
        vec2D_dbl_Type::iterator it;
        int front = -1;
        int back = -1;
        if (domainPressure->getDimension() == 2) {
            it = find_if(vectmpPointsPressure->begin(), vectmpPointsPressure->end(),
                    [&] (const vector<double>& a){
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
                    [&] (const vector<double>& a){
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
void NavierStokes<SC,LO,GO,NO>::info(){
    this->infoProblem();
    this->infoNonlinProblem();
}

template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::assemble( std::string type ) const{
    
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
void NavierStokes<SC,LO,GO,NO>::assembleConstantMatrices() const{
    
    if (this->verbose_)
        std::cout << "-- Assembly constant matrices Navier-Stokes ... " << std::flush;
    
    double viscosity = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    
    // Egal welcher Wert, da OneFunction nicht von parameter abhaengt
    int* dummy;
    
    A_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    
    if ( this->parameterList_->sublist("Parameter").get("Symmetric gradient",false) )
        this->feFactory_->assemblyStress(this->dim_, this->domain_FEType_vec_.at(0), A_, OneFunction, dummy, true);
    else
        this->feFactory_->assemblyLaplaceVecField(this->dim_, this->domain_FEType_vec_.at(0), 2, A_, true);
    
    A_->resumeFill();
    
    A_->scale(viscosity);
    A_->scale(density);
    
    A_->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());

    if (this->system_.is_null())
        this->system_.reset(new BlockMatrix_Type(2));
    
    this->system_->addBlock( A_, 0, 0 );
    assembleDivAndStab();
    
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
    string precType = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( precType == "Diagonal" || precType == "Triangular" ) {
        MatrixPtr_Type Mpressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        
        this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(1), "Scalar", Mpressure, true );
        SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
        Mpressure->scale(-1./kinVisco);
        this->getPreconditionerConst()->setPressureMassMatrix( Mpressure );
    }

    if (this->verbose_)
        std::cout << " Call Reassemble FixedPoint and Newton to allocate the Matrix pattern " << std::endl;
    
    // This was moved here from 'create_W_op'.
    // Here it will definetly be called before create_W_op and create_W_prec is called.
    this->reAssemble("FixedPoint");
    this->reAssemble("Newton");

    if (this->verbose_)
        std::cout << "done -- " << std::endl;
    
};
    
template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::assembleDivAndStab() const{
    
    double viscosity = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    
    // Egal welcher Wert, da OneFunction nicht von parameter abhaengt

    MatrixPtr_Type BT(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
    
    MapConstPtr_Type pressureMap;
    if ( this->getDomain(1)->getFEType() == "P0" )
        pressureMap = this->getDomain(1)->getElementMap();
    else
        pressureMap = this->getDomain(1)->getMapUnique();
    
    MatrixPtr_Type B(new Matrix_Type( pressureMap, this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    
    MatrixPtr_Type C;
    
    this->feFactory_->assemblyDivAndDivTFast(this->dim_, this->getFEType(0), this->getFEType(1), 2, B, BT, this->getDomain(0)->getMapVecFieldUnique(), pressureMap, true );
    
    B->resumeFill();
    BT->resumeFill();
    
    B->scale(-1.);
    BT->scale(-1.);
    
    B->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), pressureMap );
    BT->fillComplete( pressureMap, this->getDomain(0)->getMapVecFieldUnique() );
    
    this->system_->addBlock( BT, 0, 1 );
    this->system_->addBlock( B, 1, 0 );
    
    if ( !this->getFEType(0).compare("P1") ) {
        C.reset(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyBDStabilization( this->dim_, "P1", C, true);
        C->resumeFill();
        C->scale( -1. / ( viscosity * density ) );
        C->fillComplete( pressureMap, pressureMap );
        
        this->system_->addBlock( C, 1, 1 );
    }
    //else 
    //    C.reset(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
    //    this->feFactory_->assemblyEmptyMatrix(C);
               
        //this->system_->addBlock( C, 1, 1 );

};

template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const
{

}
    
template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::reAssembleFSI(std::string type, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const {
    
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
}
    

template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::reAssemble(std::string type) const {

    
    if (this->verbose_)
        std::cout << "-- Reassembly Navier-Stokes ("<< type <<") ... " << std::flush;
    
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    
    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    if (type=="FixedPoint") {
        
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true);

        MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyAdvectionVecField( this->dim_, this->domain_FEType_vec_.at(0), N, u_rep_, true );
        
        N->resumeFill();
        N->scale(density);
        N->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
        
        A_->addMatrix(1.,ANW,0.);
        N->addMatrix(1.,ANW,1.);
    }
    else if(type=="Newton"){ // We assume that reAssmble("FixedPoint") was already called for the current iterate
        MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyAdvectionInUVecField( this->dim_, this->domain_FEType_vec_.at(0), W, u_rep_, true );
        W->resumeFill();
        W->scale(density);
        W->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
        this->system_->getBlock( 0, 0 )->addMatrix(1.,ANW,0.);
        W->addMatrix(1.,ANW,1.);
    }
    ANW->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );
    
    this->system_->addBlock( ANW, 0, 0 );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions){

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

//template<class SC,class LO,class GO,class NO>
//int NavierStokes<SC,LO,GO,NO>::ComputeDragLift(vec_dbl_ptr_Type &values){
//
//    int dimension = this->domainPtr_vec_.at(0)->GetDimension();
//    MultiVector_ptr_vec_ptr_Type sol_unique_vec(new std::vector<MultiVector_ptr_Type>(0));
//    this->system_->FillSplitVector64(*this->solution_, sol_unique_vec);
//    this->system_->BuildRepeatedVectorBlocks(sol_unique_vec);
//    Vector_ptr_Type u_rep(new Epetra_Vector(*(this->system_->GetRepeatedVec(0))));
//    Teuchos::RCP<Epetra_FECrsMatrix> 	N(new Epetra_FECrsMatrix(Epetra_DataAccess::Copy,*(this->domainPtr_vec_.at(0)->GetMapXDimUnique()),10));
//    u_rep.reset(new Epetra_Vector(*(this->system_->GetRepeatedVec(0))));
//    N.reset(new Epetra_FECrsMatrix(Epetra_DataAccess::Copy,*(this->domainPtr_vec_.at(0)->GetMapXDimUnique()),10));
//    this->feFactory_->AssemblyAdvectionXDim(dimension, this->domain_FEType_vec_.at(0), 7, N, u_rep /* u */, setZeros);
//
//    Teuchos::RCP<Epetra_CrsMatrix>   	AN(new Epetra_CrsMatrix(Epetra_DataAccess::Copy,*(this->domainPtr_vec_.at(0)->GetMapXDimUnique()),10));
//
//    EpetraExt::MatrixMatrix::Add(*A_,false,1.,*AN,1.);
//    EpetraExt::MatrixMatrix::Add(*N,false,1.,*AN,1.);
//
//    AN.reset(new Epetra_CrsMatrix(Epetra_DataAccess::Copy,*(this->domainPtr_vec_.at(0)->GetMapXDimUnique()),10));
//    EpetraExt::MatrixMatrix::Add(*A_,false,1.,*AN,1.);
//    EpetraExt::MatrixMatrix::Add(*N,false,1.,*AN,1.);
//    AN->FillComplete();
//
//
//    Teuchos::RCP<Epetra_FECrsMatrix> 	B_T (new Epetra_FECrsMatrix(Epetra_DataAccess::Copy,*(this->domainPtr_vec_.at(0)->GetMapXDimUnique()),5));
//    Teuchos::RCP<Epetra_FECrsMatrix> 	B (new Epetra_FECrsMatrix(Epetra_DataAccess::Copy,*(this->domainPtr_vec_.at(1)->GetMapUnique()),5));
//    this->feFactory_->AssemblyDivergence(dimension, this->domain_FEType_vec_.at(0), this->domain_FEType_vec_.at(1), 2, B, B_T, setZeros, this->domainPtr_vec_.at(0)->GetMapXDimUnique(), this->domainPtr_vec_.at(1)->GetMapUnique());
//    B_T->Scale(-1.);
//    Teuchos::RCP<BlockElement> 	BE_AN(new BlockElement(AN));
//    Teuchos::RCP<BlockElement> 	BE_B_T(new BlockElement(B_T));
//
//    BE_AN.reset(new BlockElement(AN));
//    BE_B_T.reset(new BlockElement(B_T));
//
//    this->system_->ReplaceBlock(BE_AN,0,0);
//    this->system_->ReplaceBlock(BE_B_T,0,1);
//
//    Teuchos::RCP<BCBuilder> bCFactoryDrag(new BCBuilder(sublist(this->parameterList_,"Parameter")));
//    Teuchos::RCP<BCBuilder> bCFactoryLift(new BCBuilder(sublist(this->parameterList_,"Parameter")));
//
//    bCFactoryDrag->AddBC(sxOne2D, 4, 0, this->domainPtr_vec_.at(0), "Dirichlet", dimension);
//    bCFactoryDrag->AddBC(sDummyFunc, 666, 1, this->domainPtr_vec_.at(1), "Neumann", 1);
//
//    bCFactoryLift->AddBC(syOne2D, 4, 0, this->domainPtr_vec_.at(0), "Dirichlet", dimension);
//    bCFactoryLift->AddBC(sDummyFunc, 666, 1, this->domainPtr_vec_.at(1), "Neumann", 1);
//
//    Teuchos::RCP<Epetra_Vector>  dragVec(new Epetra_Vector(*(*this->solution_)(0)));
//    Teuchos::RCP<Epetra_Vector>	liftVec(new Epetra_Vector(*(*this->solution_)(0)));
//    dragVec->PutScalar(0.);
//    liftVec->PutScalar(0.);
//    bCFactoryDrag->SetRHS(this->system_, dragVec);
//    bCFactoryLift->SetRHS(this->system_, liftVec);
//
//    Teuchos::RCP<Epetra_Vector>	mat_sol(new Epetra_Vector(*(*this->solution_)(0)));
//    this->system_->Apply(*this->solution_,*mat_sol);
//    mat_sol->Scale(-1.);
//    double dragCoeff;
//    double liftCoeff;
//    mat_sol->Dot(*dragVec,&dragCoeff);
//    mat_sol->Dot(*liftVec,&liftCoeff);
//    values->at(0) = dragCoeff;
//    values->at(1) = liftCoeff;
//    if (this->verbose_) {
//        cout<< "Not scaled drag coefficient: " << dragCoeff<< endl;
//        cout<< "Not scaled lift coefficient: " << liftCoeff<< endl;
//    }
//    Teuchos::RCP<Epetra_Vector>     pressureSolutuion( new Epetra_Vector(*((*(sol_unique_vec->at(1)))(0))));
//    double p1 = numeric_limits<double>::min();
//    double p2 = numeric_limits<double>::min();
//    if (pressureIDsLoc->at(0)>-1) {
//        p1 = (*pressureSolutuion)[pressureIDsLoc->at(0)];
//
//    }
//    if (pressureIDsLoc->at(1)>-1) {
//        p2 = (*pressureSolutuion)[pressureIDsLoc->at(1)];
//    }
//    this->comm_->Barrier();
//    double res;
//    this->comm_->MaxAll(&p1,&res,1);
//    values->at(2) = res;
//    this->comm_->MaxAll(&p2,&res,1);
//    values->at(3) = res;
//    return 0;
//}

//template<class SC,class LO,class GO,class NO>
//typename NavierStokes<SC,LO,GO,NO>::MultiVector_Type NavierStokes<SC,LO,GO,NO>::GetExactSolution(double time){
//#ifdef ASSERTS_WARNINGS
//    MYASSERT(false,"no analytic solution.");
//#endif
//    return *this->solution_;
//}


//template<class SC,class LO,class GO,class NO>
//void NavierStokes<SC,LO,GO,NO>::set_x0(const Teuchos::ArrayView<const SC> &x0_in){
//#ifdef TEUCHOS_DEBUG
//    TEUCHOS_ASSERT_EQUALITY(xSpace_->dim(), x0_in.size());
//#endif
//    Thyra::DetachedVectorView<SC> x0(x0_);
//    x0.sv().values()().assign(x0_in);
//}

template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
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
void NavierStokes<SC,LO,GO,NO>::evalModelImplMonolithic(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
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
void NavierStokes<SC,LO,GO,NO>::evalModelImplBlock(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
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
void NavierStokes<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
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
void NavierStokes<SC,LO,GO,NO>::calculateNonLinResidualVecWithMeshVelo(std::string type, double time, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const{


    this->reAssembleFSI( "FixedPoint", u_minus_w, P );
    
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
Teuchos::RCP<Thyra::LinearOpBase<SC> > NavierStokes<SC,LO,GO,NO>::create_W_op() const
{
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
Teuchos::RCP<Thyra::LinearOpBase<SC> > NavierStokes<SC,LO,GO,NO>::create_W_op_Monolithic() const
{
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = this->system_->getThyraLinOp();
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
    return W_op;
}

#ifdef FEDD_HAVE_TEKO
template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > NavierStokes<SC,LO,GO,NO>::create_W_op_Block() const
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
Teuchos::RCP<Thyra::PreconditionerBase<SC> > NavierStokes<SC,LO,GO,NO>::create_W_prec() const
{
    this->initializeSolverBuilder();

    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    this->setBoundariesSystem();

    if (!type.compare("Teko")) { //
        this->setupPreconditioner( type );
        stokesTekoPrecUsed_ = false;
    }
    else{
        this->setupPreconditioner( type ); // initializePreconditioner( type );
    }
    

    Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  this->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst = Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);

    return thyraPrecNonConst;

}
}

#endif
