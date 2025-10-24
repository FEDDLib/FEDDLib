#ifndef NAVIERSTOKES_def_hpp
#define NAVIERSTOKES_def_hpp
#include "NavierStokes_decl.hpp"

#ifndef NAVIER_STOKES_START
#define NAVIER_STOKES_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Assemble Navier-Stokes:") + std::string(S))));
#endif

#ifndef NAVIER_STOKES_STOP
#define NAVIER_STOKES_STOP(A) A.reset();
#endif

/*!
 Definition of Navier-Stokes

 @brief Navier-Stokes
 @authors Christian Hochmuth, Lea Saßmannshausen
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

void zeroDirichletBC(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;

    return;
}

double OneFunction(double* x, int* parameter)
{
    return 1.0;
}

void dummyFuncRhs(double* x, double* res, double* parameters){

    // parameters[1] contains the surface flag, parameter[0] contains the inlet flag
    if(parameters[1]==parameters[0])
        res[0]=1;
    else
        res[0] = 0.;

    return;
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

    if(!this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("PCD")
        || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("LSC")
        || !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("PCD") 
        || !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("LSC-Pressure-Laplace") )
    { 

        int flagOutlet = this->parameterList_->sublist("General").get("Flag Outlet Fluid", 3);
        int flagInterface = this->parameterList_->sublist("General").get("Flag Interface", 6);

        this->bcFactoryPCD_.reset(new BCBuilder<SC,LO,GO,NO>( ));
        this->bcFactoryPCD_->addBC(zeroDirichletBC, flagOutlet, 0, Teuchos::rcp_const_cast<Domain_Type>( domainPressure ), "Dirichlet", 1);
        // this->bcFactoryPCD_->addBC(zeroDirichletBC, flagInterface, 0, Teuchos::rcp_const_cast<Domain_Type>( domainPressure ), "Dirichlet", 1);
        // this->bcFactoryPCD_->addBC(zeroDirichletBC, 9, 0, Teuchos::rcp_const_cast<Domain_Type>( domainPressure ), "Dirichlet", 1);
        // this->bcFactoryPCD_->addBC(zeroDirichletBC, 10, 0, Teuchos::rcp_const_cast<Domain_Type>( domainPressure ), "Dirichlet", 1);

    } 

    if(this->parameterList_->sublist("General").get("Augmented Lagrange",false))  
        augmentedLagrange_ = true;

    // Establish the non zero pattern of the system matrix in (0,0) block
    NNZ_A_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    establishNNZPattern();

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
    else if(type=="UpdateTime"){
        this->newtonStep_ = 0;
        // timeSteppingTool_->t_ = timeSteppingTool_->t_ + timeSteppingTool_->dt_prev_;
    }
    else if(type =="UpdateConvectionDiffusionOperator")
        this->updateConvectionDiffusionOperator();
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
    if ( !this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic").compare("Teko") 
    || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("PCD")
    || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("LSC")) {

        // ###############################################
        // LSC Preconditioner
        // Constructing velocity mass matrix
        // If the Velocity Mass Matrix is the identity matrix, 
        // it results in the BFBt preconditioner
        if (!this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","None").compare("LSC")
         || !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","None").compare("LSC-Pressure-Laplace")
         || !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","None").compare("SIMPLE")
         || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("LSC")) {
                        
            MatrixPtr_Type Mvelocity(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // Constructing velocity mass matrix
            if(this->parameterList_->sublist("Parameter").get("BFBT",false)){
                if(this->verbose_)
                    std::cout << "\n Setting M_u to be the identity Matrix to use BFBT preconditioner " << std::endl;

                this->feFactory_->assemblyIdentity( Mvelocity );
            }
            else{ 
                this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity, true );
            }
            // Adding the velocity mass matrix Mu to the preconditioner
            this->getPreconditionerConst()->setVelocityMassMatrix( Mvelocity );

           if (this->verbose_)
                std::cout << "\n Velocity mass matrix for LSC block preconditioner is assembled and used for the preconditioner." << std::endl;

            // For LSC-Pressure-Laplace approach we add also the Laplaian on the pressure space to the preconditioner
            if(!this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("LSC-Pressure-Laplace")
                || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("LSC"))
            {
                MatrixPtr_Type Lp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
                this->feFactory_->assemblyLaplace( this->dim_, this->domain_FEType_vec_.at(1), 2, Lp, true );
                
                BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));
                bcBlockMatrix->addBlock(Lp,0,0);
                bcFactoryPCD_->setSystemScaled(bcBlockMatrix); // Setting boundary information where the Diagonal entry is kept 

                this->getPreconditionerConst()->setPressureLaplaceMatrix( Lp);
            }
        } 
        // ###############################################
        // PCD Preconditioner
        // For the PCD preconditioner we additionally need 
        // assemble the pressure mass matrix, the pressure
        // Laplacian and the pressure convection diffusion
        // opertor.
        else if(!this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("PCD") 
        || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("PCD") ){
            
            // ###############################################
            // Velocity mass matrix: Currently this is set to not have an error in preconditioner. PLEASE FIX
            MatrixPtr_Type Mvelocity(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity, true );
            this->getPreconditionerConst()->setVelocityMassMatrix( Mvelocity );
            // ###############################################


            // Pressure mass matrix
            MatrixPtr_Type Mpressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(1), "Scalar", Mpressure, true ); 
            Mp_= Mpressure;
            this->getPreconditionerConst()->setPressureMassMatrix( Mpressure );
            // --------------------------------------------------------------------------------------------

            // --------------------------------------------------------------------------------------------
            // Pressure Laplace matrix
            SC density = this->parameterList_->sublist("Parameter").get("Density",1.); // Ap need to be scaled with viscosity

            MatrixPtr_Type Lp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyLaplace( this->dim_, this->domain_FEType_vec_.at(1), 2, Lp, true );
            // Lp->resumeFill();
            // Lp->scale(density);
            // Lp->fillComplete();
            Ap_.reset(new Matrix_Type(Lp)); // Setting Ap_ as Lp without any BC
            
            // Adding boundary information to pressure Laplace operator
            BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));
            bcBlockMatrix->addBlock(Lp,0,0);

            bcFactoryPCD_->setSystemScaled(bcBlockMatrix); // Setting boundary information where the Diagonal entry is kept 
            this->getPreconditionerConst()->setPressureLaplaceMatrix( Lp);   // Adding pressure laplacian to preconditioner
            // --------------------------------------------------------------------------------------------

            // --------------------------------------------------------------------------------------------
            // Pressure convection diffusion operator 
            MatrixPtr_Type Kp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            // --------------------------------------------------------------------------------------------
            // Advection component
            MatrixPtr_Type AdvPressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyAdvectionVecFieldScalar( this->dim_, this->domain_FEType_vec_.at(1), this->domain_FEType_vec_.at(0),AdvPressure, u_rep_, true ); 
           
            // Diffusion component: \nu * \Delta
            MatrixPtr_Type Ap2(new Matrix_Type( Ap_) );

            SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.); // Ap need to be scaled with viscosity

            Ap2->resumeFill();
            Ap2->scale(kinVisco);
            // Ap2->scale(density);
            Ap2->fillComplete(); 
            
            
            MatrixPtr_Type K_robin(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow()*2 ) );          
            int flagInlet =this->parameterList_->sublist("General").get("Flag Inlet Fluid", 2);

            vec_dbl_Type funcParameter(1,flagInlet);
            funcParameter.push_back(0.0); // Dummy for flag
            this->feFactory_->assemblySurfaceRobinBC(this->dim_, this->getDomain(1)->getFEType(),this->getDomain(0)->getFEType(),u_rep_,K_robin, funcParameter, dummyFuncRhs,this->parameterList_);
            K_robin->addMatrix(-1.,Kp,1.); // adding robin boundary condition to to Kp
            
            // Adding laplace and convetion-diffusion operator to Kp
            Ap2->addMatrix(1.,Kp,1.); // adding advection to diffusion
            AdvPressure->addMatrix(1.,Kp,1.); // adding advection to diffusion
            
            // Kp->scale(density);
            Kp->fillComplete();

            bcBlockMatrix->addBlock(Kp,0,0);
            bcFactoryPCD_->setSystemScaled(bcBlockMatrix);

            // Adding Kp to the preconditioner
            this->getPreconditionerConst()->setPCDOperator( Kp );  
            // --------------------------------------------------------------------------------------------

        }
    }
#endif
    std::string precType = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( precType == "Diagonal" || precType == "Triangular" ) {
        MatrixPtr_Type Mpressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        
        this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(1), "Scalar", Mpressure, true );
        SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.);

        if(augmentedLagrange_){
            double gamma = this->parameterList_->sublist("General").get("Gamma",1.0);
            Mpressure->scale(-1./(kinVisco+gamma));
        }
        else{
            Mpressure->scale(-1./kinVisco);
        }
        this->getPreconditionerConst()->setPressureMassMatrix( Mpressure );
    }

    if (this->verbose_)
        std::cout << " Allocate the Matrix pattern " << std::endl;
    
    // After the constant matrices are assembled, we establish the matrix pattern for the (0,0) block for the advection term
    MatrixPtr_Type A_withNNZ = Teuchos::rcp( new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    A_->addMatrix(1.,A_withNNZ,0.); // We add the previously computed laplacian
    NNZ_A_->addMatrix(1.,A_withNNZ,1.); // We add the zero matrix containing the nnz pattern

    A_withNNZ->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
    this->system_->addBlock( A_withNNZ, 0, 0 ); // We replace the (0,0) block with the new matrix containing the laplacian and the nnz pattern for the advection term

    if (this->verbose_)
        std::cout << "done -- " << std::endl;
    
};
    

template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::updateConvectionDiffusionOperator() const{

    if ( !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("PCD") 
                || !this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic").compare("PCD")) 
    {
        std::cout << "NavierStokes<SC,LO,GO,NO>::updateConvectionDiffusionOperator() set to TRUE" << std::endl;
        NAVIER_STOKES_START(ReassemblePCD," Reassembling Matrix for PCD ");
      
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true);

        // PCD Operator  
        MatrixPtr_Type Fp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        // --------------------------------------------------------------------------------------------
        // Advection component
        MatrixPtr_Type AdvPressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyAdvectionVecFieldScalar( this->dim_, this->domain_FEType_vec_.at(1), this->domain_FEType_vec_.at(0),AdvPressure, u_rep_, true ); 
        
        // Diffusion component: \nu * \Delta
        MatrixPtr_Type Ap2(new Matrix_Type( Ap_ ) ); // We use A_p which we already stored
        SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
        SC density = this->parameterList_->sublist("Parameter").get("Density",1.);

        Ap2->resumeFill();
        Ap2->scale(kinVisco);
        Ap2->fillComplete();

        // ---------------------
        // Robin boundary
        MatrixPtr_Type Kext(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow()*2 ) );          
        int flagInlet =this->parameterList_->sublist("General").get("Flag Inlet Fluid", 2);
        vec_dbl_Type funcParameter(1,flagInlet);
        funcParameter.push_back(0.0); // Dummy for flag

        this->feFactory_->assemblySurfaceRobinBC(this->dim_, this->getDomain(1)->getFEType(),this->getDomain(0)->getFEType(),u_rep_,Kext, funcParameter, dummyFuncRhs,this->parameterList_);
        Kext->addMatrix(-1.,Fp,1.); // adding advection to diffusion
        
        // Adding laplace an convection together
        Ap2->addMatrix(1.,Fp,1.); // adding advection to diffusion
        AdvPressure->addMatrix(1.,Fp,1.); // adding advection to diffusion

        // Finally if we deal with a transient problem we additionally add the Mass term 1/delta t M_p
        ///TODO: Extract parameters from timestepping tool.
        if(this->parameterList_->sublist("Timestepping Parameter").get("dt",-1.)> -1 ){ // In case we have a timeproblem
            MatrixPtr_Type Mp2(new Matrix_Type( Mp_ ) );
            double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",-1.);
            Mp2->resumeFill();
            if(this->parameterList_->sublist("Timestepping Parameter").get("BDF",1) == 1) // BDF 1
                Mp2->scale(1./dt);
            else if(this->parameterList_->sublist("Timestepping Parameter").get("BDF",1) == 2) // BDF 1
                Mp2->scale(3./(2.*dt));
            else
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "PCD operator for transient problems only defined for BDF-1 and BDF-2.");

            // Mp2->scale(density);
            Mp2->fillComplete();
            Mp2->addMatrix(1.,Fp,1.);
        }
        // Fp->scale(1./density);
        Fp->fillComplete();

        BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));

        bcBlockMatrix->addBlock(Fp,0,0);   
        bcFactoryPCD_->setSystemScaled(bcBlockMatrix); // We set the boundary conditions into Fp. Both Ap and Fp have Dirichlet zero bc on the outlet
                                                       // Addionally, a robin bc is applied to the inlet bc for Fp. 
                                                       // Note, if no surfaces are available, the matrix Kext, containg the robin bc is zero,
                                                       // so no robin bc is applied.

        this->getPreconditionerConst()->setPCDOperator( Fp );      

        NAVIER_STOKES_STOP(ReassemblePCD);       
    }
}
template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::assembleDivAndStab() const{
    
    double viscosity = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    
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
    
    if ( !this->getFEType(0).compare("P1") ||  !this->getFEType(0).compare("Q1") ) {
        C.reset(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyBDStabilization( this->dim_, this->getFEType(0), C, true);
        C->resumeFill();
        C->scale( -1. / ( viscosity * density ) ); // scaled with dynamic viscosity with mu = nu*rho
        C->fillComplete( pressureMap, pressureMap );
        
        this->system_->addBlock( C, 1, 1 );
    }

    // Implementation of augmented lagrange. This works in theory, but the resulting matrix changes the nnz pattern and has a wider FE Stenciln than before. This can present an issue for the algebraic overlap.
    // Compare for example to 'ANALYSIS OF AUGMENTED LAGRANGIAN-BASED PRECONDITIONERS FOR THE STEADY INCOMPRESSIBLE NAVIER–STOKES EQUATIONS, MICHELE BENZI AND ZHEN WANG' for the theory behind this.
    // In augmented lagrange, the term \gamma B^T M_p^{-1} B is added to the velocity-velocity block, where M_p is the pressure mass matrix.
    if(augmentedLagrange_){
        NAVIER_STOKES_START(AssembleAugmentedLagrangianComponent,"AssembleDivAndStab: AL - Assemble BT Mp B");

        double gamma = this->parameterList_->sublist("General").get("Gamma",1.0);

        MatrixPtr_Type Mp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(1), "Scalar", Mp, true ); 

        MatrixPtr_Type MpInv = Mp->buildDiagonalInverse("Diagonal");

        MpInv->resumeFill();
        MpInv->scale(gamma);
        MpInv->fillComplete();


        MatrixPtr_Type BT_M(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension()*this->getDomain(0)->getApproxEntriesPerRow() ) );
        BT_M->Multiply(BT,false,MpInv,false);

        BT_Mp_ = BT_M;

        MatrixPtr_Type BT_M_B(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension()*this->getDomain(0)->getApproxEntriesPerRow() ) );
        BT_M_B->Multiply(BT_M,false,B,false);

        BT_Mp_B_ = BT_M_B;

        NAVIER_STOKES_STOP(AssembleAugmentedLagrangianComponent);
    }
   
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

    // If we use augmented lagrange, the matrix fe-stencil increases and we need to allocate more nnz entries
    int allocationFactor = 1;
    if(augmentedLagrange_)
        allocationFactor = 3;
    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), allocationFactor*this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
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
    
    if(augmentedLagrange_)
        BT_Mp_B_->addMatrix(1.,ANW,1.);

    ANW->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );

    this->system_->addBlock( ANW, 0, 0 );
 
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

// This function assembles only the non zero pattern of the block (0,0) of the system matrix
// This is done by using a zero solution vector zeroVec
// Then, the advection matrices are assembled and added to an empty matrix ANW
// The resulting matrix ANW contains then the non zero pattern 
// This is then stored as NNZ_A_
template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::establishNNZPattern() const {

   
    if (this->verbose_)
        std::cout << "-- Establish NNZ Pattern Navier-Stokes ... " << std::flush;
    
    MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        
    MultiVectorPtr_Type zeroVec = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated(), 1 ) );
    zeroVec->putScalar(0.0);

    MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyAdvectionVecField( this->dim_, this->domain_FEType_vec_.at(0), N, zeroVec, true );
    
    N->addMatrix(1.,ANW,0.);
     
    MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyAdvectionInUVecField( this->dim_, this->domain_FEType_vec_.at(0), W, zeroVec, true );
 
    W->addMatrix(1.,ANW,1.);
    
    ANW->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );

    NNZ_A_= ANW;
 
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
void NavierStokes<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    
    if (this->verbose_)
        std::cout << "-- NavierStokes::calculateNonLinResidualVec ("<< type <<") ... " << std::flush;
    
    this->reAssemble("FixedPoint");
    // We need to account for different parameters of time discretizations here
    // This is ok for bdf with 1.0 scaling of the system. Would be wrong for Crank-Nicolson - might be ok now for CN
    if (this->coeff_.size() == 0)
        this->system_->apply( *this->solution_, *this->residualVec_ );
    else
        this->system_->apply( *this->solution_, *this->residualVec_, this->coeff_ );
    
    // The additional component needs to be acconted for in residual as well due to augmented lagrange
    if(augmentedLagrange_){
        MultiVectorPtr_Type rhsAL = Teuchos::rcp( new MultiVector_Type( this->residualVec_->getBlock(0) ) );
        BT_Mp_->apply( *this->residualVec_->getBlock(1), *rhsAL );
        this->residualVec_->getBlockNonConst(0)->update(1.,*rhsAL,1.);

    }
    if (!type.compare("standard")){
        this->residualVec_->update(-1.,*this->rhs_,1.);
        // this might be set again by the TimeProblem after addition of M*u
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
        
    }
    else if(!type.compare("reverse")){
        this->residualVec_->update(1.,*this->rhs_,-1.); // this = -1*this + 1*rhs
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
        
}

}

#endif
