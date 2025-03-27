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

void drag2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 1.;
    res[1] = 0.;
    
    return;
}

void drag3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void lift2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 1.;
    
    return;
}

void lift3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 1.;
    res[2] = 0.;
    
    return;
}
void sDummyFunc(double* x, double* res, double t, double* parameter){

    return;
}

double OneFunction(double* x, int* parameter)
{
    return 1.0;
}

using namespace std;
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::outArg;

namespace FEDD {



template<class SC,class LO,class GO,class NO>
NavierStokes<SC,LO,GO,NO>::NavierStokes( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList ):
NonLinearProblem<SC,LO,GO,NO>( parameterList, domainVelocity->getComm() ),
A_(),
pressureIDsLoc(new vec_int_Type(2)),
u_rep_(),
p_rep_()
{

    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    this->initNOXParameters();

    this->addVariable( domainVelocity , FETypeVelocity , "u" , domainVelocity->getDimension());
    this->addVariable( domainPressure , FETypePressure , "p" , 1);
    this->dim_ = this->getDomain(0)->getDimension();

    u_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    p_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );

    this->newtonStep_=0;

    
    if(this->parameterList_->sublist("Timestepping Parameter").get("dt",-1.)> 0)
        timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));

    if(this->dim_ ==3){
        // Values we need to estimate RE and CFL in 3D
        domainVelocity->getMesh()->calcDiamTetraeder();
        domainVelocity->getMesh()->calcRhoTetraeder();
        domainVelocity->getMesh()->determineLongestEdge();
        // Reynolds number and CFL number estimations
        exporterTxtCFLMax_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtCFLMax_->setup( "CFL_max", this->comm_ );
        exporterTxtReMax_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtReMax_->setup( "Re_max", this->comm_ );
        
        exporterTxtCFLMin_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtCFLMin_->setup( "CFL_min", this->comm_ );
        // exporterTxtReMin_ = Teuchos::rcp(new ExporterTxt () );
        // exporterTxtReMin_->setup( "Re_min", this->comm_ );
        // ---------------------------------------------
    }

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
    if ( parameterList->sublist("General").get("Export drag and lift",false) ){
        exporterTxtDrag_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtDrag_->setup( "drag_force", this->comm_ );
        exporterTxtLift_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtLift_->setup( "lift_force", this->comm_ );
    }
    if ( parameterList->sublist("Parameter").get("Set Zeros",false) ){
        double eps = parameterList->sublist("Parameter").get("Zeros Tolerance",1.e-13);
        this->feFactory_->doSetZeros(eps);
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
    else if(type=="UpdateTime"){
        this->newtonStep_ = 0;
        timeSteppingTool_->t_ = timeSteppingTool_->t_ + timeSteppingTool_->dt_prev_;

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

   // In case of a monolithic preconditioner and a P2-P1 discretization we have the option to correct the pressure to have mean value = 0. This way, generally, we can improve scalabilty and results. 
    // The real correction is then done via projection in the Overlapping Operator of FROSch,here we only assemble a as \int p dx . a is assembled as a column vector but in the Dissertation of C. Hochmuth defined as row.
    if(this->parameterList_->sublist("Parameter").get("Use Pressure Correction",false) && (!this->getFEType(0).compare("P2") || (!this->getFEType(0).compare("Q2") && !this->getFEType(1).compare("Q1"))) && !this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic").compare("Monolithic")){ 
        // Projection vector a: \int p dx, for pressure component and 0 for velocity.
        BlockMultiVectorPtr_Type projection(new BlockMultiVector_Type (2));

        MultiVectorPtr_Type P(new MultiVector_Type( this->getDomain(1)->getMapUnique(), 1 ) );

        this->feFactory_->assemblyPressureMeanValue( this->dim_,this->getFEType(1),P) ;

        // Velocity component is set to zero, such that the projection vector only influences the pressure part
        MultiVectorPtr_Type vel0(new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
        vel0->putScalar(0.);

        // Adding components to projection vector 
        projection->addBlock(vel0,0);
        projection->addBlock(P,1);

        // Setting projection vector in preconditioner to later pass to paramterlist in FROSch
        this->getPreconditionerConst()->setPressureProjection( projection );    

        if (this->verbose_)
            std::cout << "\n 'Use pressure correction' was set to 'true'. This requieres a version of Trilinos of that includes pressure correction in the FROSch_OverlappingOperator!!" << std::endl;  

    }
    else if(this->parameterList_->sublist("Parameter").get("Use Pressure Correction",false) && (!this->getFEType(0).compare("P2") || (!this->getFEType(0).compare("Q2") && !this->getFEType(1).compare("Q1"))) && !this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic").compare("Teko")){ 
        // Projection vector a: \int p dx, for pressure component and 0 for velocity.
        BlockMultiVectorPtr_Type projection(new BlockMultiVector_Type (1));

        MultiVectorPtr_Type P(new MultiVector_Type( this->getDomain(1)->getMapUnique(), 1 ) );

        this->feFactory_->assemblyPressureMeanValue( this->dim_,this->getFEType(1),P) ;

        // Velocity component is set to zero, such that the projection vector only influences the pressure part
        // MultiVectorPtr_Type vel0(new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
        // vel0->putScalar(0.);

        // Adding components to projection vector 
        projection->addBlock(P,0);

        // Setting projection vector in preconditioner to later pass to paramterlist in FROSch
        this->getPreconditionerConst()->setPressureProjection( projection );    

        if (this->verbose_)
            std::cout << "\n 'Use pressure correction' was set to 'true'. This requieres a version of Trilinos of that includes pressure correction in the FROSch_OverlappingOperator!!" << std::endl;  

    }
    
#ifdef FEDD_HAVE_TEKO
    if ( !this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic").compare("Teko") 
    || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("PCD")
    || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("LSC")) {

        if (!this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("LSC")
         || !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("LSC-Pressure-Laplace")
         || !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("SIMPLE")
         || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("LSC")) {
            MatrixPtr_Type Mvelocity(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            //
            int extraDeg = this->parameterList_->sublist("Parameter").get("Extra Deg",0);
            if(this->parameterList_->sublist("Parameter").get("BFBT",false)){
                if(this->verbose_)
                    std::cout << "\n Setting M_u to be the identity Matrix to use BFBT preconditioner " << std::endl;

                this->feFactory_->assemblyIdentity( Mvelocity, true );
                Mvelocity->resumeFill();
                Mvelocity->fillComplete();
            }
            else{ // For whatever reason, when we have a stationary problem a higher degree for the quadrature improves results
                if(this->parameterList_->sublist("Timestepping Parameter").get("dt",-1.)> 0 ) // In case we have a timeproblem
                    this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity, true,extraDeg );
                else
                    this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity, true,extraDeg );
            }
            // this->getComm()->barrier();

            // MatrixPtr_Type Mvelocity0(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity0, true,0 );
            // Mvelocity0->writeMM("M_0");
                
            // this->getComm()->barrier();

            //  MatrixPtr_Type Mvelocity1(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity1, true,1 );
            // Mvelocity1->writeMM("M_1");
            // //
            // this->getComm()->barrier();

            // MatrixPtr_Type Mvelocity2(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity2, true,2 );
            // Mvelocity2->writeMM("M_2");

        
            //
            BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));
            if(this->parameterList_->sublist("Parameter").get("BC in LSC Mu",false)){
                double epsilon = this->parameterList_->sublist("Parameter").get("Scaling Mu Matrix",0.0);
                bcBlockMatrix->addBlock(Mvelocity,0,0);
                this->bcFactory_->setSystemScaled(bcBlockMatrix,1.0+epsilon); // setSystemScaled(bcBlockMatrix); 
            }
            //
            this->getPreconditionerConst()->setVelocityMassMatrix( Mvelocity );

           if (this->verbose_)
                std::cout << "\n Velocity mass matrix for LSC block preconditioner is assembled and used for the preconditioner." << std::endl;

            MatrixPtr_Type Lp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyLaplace( this->dim_, this->domain_FEType_vec_.at(1), 2, Lp, true );//assemblyIdentity(Lp); //
                       
            bcBlockMatrix->addBlock(Lp,0,0);
            double eps = this->parameterList_->sublist("Parameter").get("Scaling Ap Matrix",0.0);

            this->bcFactoryPressureLaplace_->setSystemScaled(bcBlockMatrix, 1.0+eps ); 
            this->getPreconditionerConst()->setPressureLaplaceMatrix( Lp);

            // Weighting Vector for Scaling Matrix H
            if(this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("LSC").sublist("Strategy Settings").get("Use W-Scaling",false)||
            this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("LSC").get("Use W-Scaling",false))
            {
                MultiVectorPtr_Type W(new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) );
                double epsilon = this->parameterList_->sublist("Parameter").get("Scaling W Matrix",0.1);
                this->feFactory_->assemblyWeightedMatrix( this->dim_,this->getFEType(0), epsilon,0 ,W, this->parameterList_) ;
                this->getPreconditionerConst()->setWScaling( W );

                // ExporterPtr_Type Exporter = Teuchos::rcp(new Exporter_Type());
                
                // DomainConstPtr_Type dom = this->getDomain(0);
                // std::string varName = "W";
                
                // MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( dom->getMesh() );

                // Exporter->setup(varName, meshNonConst, this->getFEType(0));

                // MultiVectorConstPtr_Type exportVector = W;
                
                // Exporter->addVariable( exportVector, "W_Scaling", "Vector", this->dim_, dom->getMapUnique() );

                // Exporter->save(0.);

                if (this->verbose_)
                    std::cout << "\n Computed W-Scaling Vector for LSC and added to preconditioner." << std::endl;

            }

        } 
        
        if(!this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("PCD") 
        || !this->parameterList_->sublist("General").get("Preconditioner Method","Diagonal").compare("PCD") ){
            
            // Velocity mass matrix
            MatrixPtr_Type Mvelocity(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(0), "Vector", Mvelocity, true,0 );
            this->getPreconditionerConst()->setVelocityMassMatrix( Mvelocity );


              // Pressure mass matrix
            MatrixPtr_Type Mpressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(1), "Scalar", Mpressure, true,2 ); //assemblyIdentity(Mpressure);//
            Mp_= Mpressure;
            this->getPreconditionerConst()->setPressureMass( Mpressure );
            // --------------------------------------------------------------------------------------------

            // Pressure Laplace matrix
            MatrixPtr_Type Lp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyLaplace( this->dim_, this->domain_FEType_vec_.at(1), 2, Lp, true );//assemblyIdentity(Lp); //
            Ap_.reset(new Matrix_Type(Lp)); // Setting Ap_ as Lp without any BC
        
            // Adding Boundary Conditions
            BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));
            bcBlockMatrix->addBlock(Lp,0,0);
            double epsilon = this->parameterList_->sublist("Parameter").get("Scaling Ap Matrix",0.0);

            this->bcFactoryPressureLaplace_->setSystemScaled(bcBlockMatrix,1.+ epsilon ); 
            this->getPreconditionerConst()->setPressureLaplaceMatrix( Lp);
            //gitLp->writeMM("A_p");
            // --------------------------------------------------------------------------------------------

            // PCD Operator  
            MatrixPtr_Type Kp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            // --------------------------------------------------------------------------------------------
            // Advection component
            MatrixPtr_Type AdvPressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyAdvectionVecFieldScalar( this->dim_, this->domain_FEType_vec_.at(1), this->domain_FEType_vec_.at(0),AdvPressure, u_rep_, true ); 
           
            // Diffusion component: \nu * \Delta
            MatrixPtr_Type Ap2(new Matrix_Type( Ap_) );
            //this->feFactory_->assemblyLaplace( this->dim_, this->domain_FEType_vec_.at(1), 2, Ap2, true );//assemblyIdentity(Lp);
            SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
            Ap2->resumeFill();
            Ap2->scale(kinVisco);
            Ap2->fillComplete(); 
            // ---------------------
            // if(this->parameterList_->sublist("Parameter").get("Fp-Ap Option 1",false)){ // Setting in Ap2 the boundaries of Lp
            //     BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));
            //     bcBlockMatrix->addBlock(Ap2,0,0);
            //     this->bcFactoryPressureLaplace_->setSystemScaled(bcBlockMatrix); 
            // }
            
            // if(this->parameterList_->sublist("Parameter").get("Fp-Ap Option 2",false)){  // Setting in Ap2 the boundaries of Fp
            //     BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));
            //     bcBlockMatrix->addBlock(Ap2,0,0);
            //     this->bcFactoryPressureFp_->setSystemScaled(bcBlockMatrix);     
            // }
            
            if(this->parameterList_->sublist("Parameter").get("Robin BC",false)){
                MatrixPtr_Type Kext(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow()*2 ) );          
                vec_dbl_Type funcParameter(1,kinVisco);
                this->feFactory_->assemblySurfaceRobinBC(this->dim_, this->getDomain(1)->getFEType(),this->getDomain(0)->getFEType(),u_rep_,Kext, funcParameter, this->rhsFuncVec_[0],this->parameterList_);
                Kext->addMatrix(-1.,Kp,1.); // adding advection to diffusion
            }

            // Adding laplace an convection together
            Ap2->addMatrix(1.,Kp,1.); // adding advection to diffusion
            AdvPressure->addMatrix(1.,Kp,1.); // adding advection to diffusion
            
            Kp->fillComplete();

            bcBlockMatrix->addBlock(Kp,0,0);
            epsilon = this->parameterList_->sublist("Parameter").get("Scaling Fp Matrix",0.0);
            this->bcFactoryPressureFp_->setSystemScaled(bcBlockMatrix,1.+epsilon);

            this->getPreconditionerConst()->setPCDOperator( Kp );  

        }
    }
#endif
    string precType = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( precType == "Diagonal" || precType == "Triangular" || precType == "PCD" || precType == "LSC"  ) {
        MatrixPtr_Type Mpressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        
        this->feFactory_->assemblyMass( this->dim_, this->domain_FEType_vec_.at(1), "Scalar", Mpressure, true );
        SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
        Mpressure->scale(-1./kinVisco);
        this->getPreconditionerConst()->setPressureMassMatrix( Mpressure ); // FOR PCD THIS IS DUMMY
    }
    
    
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

    MatrixPtr_Type B1(new Matrix_Type( pressureMap, this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    MatrixPtr_Type BT1(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyDivAndDivTFast(this->dim_, this->getFEType(0), this->getFEType(1), 2, B1, BT1, this->getDomain(0)->getMapVecFieldUnique(), pressureMap, true );

    B_ = B1; 
    BT_ = BT1;
    
    if ( !this->getFEType(0).compare("P1") || !this->getFEType(0).compare("Q1") ) {
        C.reset(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyBDStabilization( this->dim_, this->getFEType(1), C, true);
        C->resumeFill();
        C->scale( -1. / ( viscosity * density ) );
        C->fillComplete( pressureMap, pressureMap );
        //C->print();
        this->system_->addBlock( C, 1, 1 );
    }
    // else 
    // {
    //     C.reset(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
    //     this->feFactory_->assemblyEmptyMatrix(C);
               
    //     this->system_->addBlock( C, 1, 1 );
    // }
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

        if ( !this->parameterList_->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE").compare("PCD") 
                || !this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic").compare("PCD")) {
        
            NAVIER_STOKES_START(ReassemblePCD," Reassembling Matrix for PCD ")
          
            // // --------------------------------------------------------------------------------------------
            BlockMatrixPtr_Type bcBlockMatrix(new BlockMatrix_Type (1));

            // PCD Operator  
            MatrixPtr_Type Fp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            // --------------------------------------------------------------------------------------------
            // Advection component
            MatrixPtr_Type AdvPressure(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
            this->feFactory_->assemblyAdvectionVecFieldScalar( this->dim_, this->domain_FEType_vec_.at(1), this->domain_FEType_vec_.at(0),AdvPressure, u_rep_, true ); 
           
            // Diffusion component: \nu * \Delta
            MatrixPtr_Type Ap2(new Matrix_Type( Ap_ ) ); // We use A_p which we already stored
            //this->feFactory_->assemblyLaplace( this->dim_, this->domain_FEType_vec_.at(1), 2, Ap2, true );//assemblyIdentity(Lp);
            SC kinVisco = this->parameterList_->sublist("Parameter").get("Viscosity",1.);
            Ap2->resumeFill();
            Ap2->scale(kinVisco);
            Ap2->fillComplete(); 
            // ---------------------
          
            
            if(this->parameterList_->sublist("Parameter").get("Robin BC",false)){
                MatrixPtr_Type Kext(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow()*2 ) );          
                vec_dbl_Type funcParameter(1,kinVisco);
                this->feFactory_->assemblySurfaceRobinBC(this->dim_, this->getDomain(1)->getFEType(),this->getDomain(0)->getFEType(),u_rep_,Kext, funcParameter, this->rhsFuncVec_[0],this->parameterList_);
                Kext->addMatrix(-1.,Fp,1.); // adding advection to diffusion
            }
            // Setting boundary conditions in Fp
           

            // Adding laplace an convection together
            Ap2->addMatrix(1.,Fp,1.); // adding advection to diffusion
            AdvPressure->addMatrix(1.,Fp,1.); // adding advection to diffusion
  
            // Finally if we deal with a transient problem we additionally add the Mass term 1/delta t M_p
            if(this->parameterList_->sublist("Timestepping Parameter").get("dt",-1.)> -1 ){ // In case we have a timeproblem
                MatrixPtr_Type Mp2(new Matrix_Type( Mp_ ) );
                double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",-1.);
                Mp2->resumeFill();
                if(this->parameterList_->sublist("Timestepping Parameter").get("BDF",1) < 2) // BDF 1
                    Mp2->scale(1./dt);
                else // BDF 2
                    Mp2->scale(3./(2.*dt));
                Mp2->fillComplete();

                bcBlockMatrix->addBlock(Mp2,0,0);
                this->bcFactoryPressureLaplace_->setSystemScaled(bcBlockMatrix);

                Mp2->addMatrix(1.,Fp,1.);

            }
            Fp->fillComplete();

            bcBlockMatrix->addBlock(Fp,0,0);   
            double epsilon = this->parameterList_->sublist("Parameter").get("Scaling Fp Matrix",0.0);
            this->bcFactoryPressureFp_->setSystemScaled(bcBlockMatrix,1.+epsilon); 


            this->getPreconditionerConst()->setPCDOperator( Fp );       
            NAVIER_STOKES_STOP(ReassemblePCD);       
        }
    }
    else if(type=="Newton"){ 
        // if(this->parameterList_->sublist("Parameter").get("Symmetric BC",true))
        //     this->reAssemble("FixedPoint");

        // We assume that reAssmble("FixedPoint") was already called for the current iterate
        MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyAdvectionInUVecField( this->dim_, this->domain_FEType_vec_.at(0), W, u_rep_, true );
        W->resumeFill();
        W->scale(density);
        W->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
        this->system_->getBlock( 0, 0 )->addMatrix(1.,ANW,0.);
        W->addMatrix(1.,ANW,1.);
        W_ = W;

        
    }
    if(!ANW->isFillComplete())
        ANW->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );
    
    this->system_->addBlock( ANW, 0, 0 );
    
    MultiVectorConstPtr_Type p = this->solution_->getBlock(1);
    p_rep_->importFromVector(p, true);

    Teuchos::Array<SC> norm(1);
    p_rep_->norm2(norm());

    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}


template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const{
    if (this->verbose_)
        std::cout << "--- calculateNonLinResidualVec ("<< type <<") ... " << std::flush;
    

    // We need to account for different parameters of time discretizations here
    // This is ok for bdf with 1.0 scaling of the system. Would be wrong for Crank-Nicolson - might be ok now for CN
    if(!type.compare("standard") || !type.compare("reverse")){

        this->reAssemble("FixedPoint");

        this->bcFactory_->setSystem(this->system_);
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
    else if(!type.compare("rhs_W")){
        // this->reAssemble("Newton");
        BlockMatrixPtr_Type sysW(new BlockMatrix_Type (2));
        sysW->addBlock(W_,0,0);  
        
        // W_->print();
        this->bcFactory_->setSystem(sysW);
        sysW->apply( *this->solution_, *this->residualVec_ );
        
        // this->residualVec_->update(0.,*this->residualVec_,-1.); // this = -1*this + 1*rhs

        // this->residualVec_->scale(-1.);
        this->bcFactory_->setRHS(this->residualVec_);
        // this->residualVec_->scale(-1.);

        // this->solution_->print();
        // this->residualVec_->print();

        // this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );    
    
    }
    // We need to account for different parameters of time discretizations here
    // This is ok for bdf with 1.0 scaling of the system. Would be wrong for Crank-Nicolson - might be ok now for CN
     else if(!type.compare("residual_W")){
        this->reAssemble("FixedPoint");

        this->bcFactory_->setSystem(this->system_);
        this->system_->apply( *this->solution_, *this->residualVec_ );
       
        
        this->residualVec_->update(1.,*this->rhs_,-1.); // this = -1*this + 1*rhs

        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );    
        
    }


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

template<class SC,class LO,class GO,class NO>
void NavierStokes<SC,LO,GO,NO>::computeValuesOfInterestAndExport(){

   if ( this->parameterList_->sublist("General").get("Export drag and lift",false) ) {
        int dim = this->dim_;
        TEUCHOS_TEST_FOR_EXCEPTION( this->parameterList_->sublist("Parameter").get("Criterion","Residual") == "Update",  std::runtime_error, "Wrong nonlinear criterion to calculate the drag coefficient. The last system is the Newton system but we need the fixed point system. Either use Criterion=Residual or implement for Criterion=Update." );
        
        //TEUCHOS_TEST_FOR_EXCEPTION( this->hasSourceTerm(),  std::runtime_error, "We need to substract the additional source term: drag = < F*u + B_T*p + C1_T*lamba - f, v >" );
        
        Teuchos::Array<SC> drag(1);
        Teuchos::Array<SC> lift(1);
        
        BlockMultiVectorPtr_Type uDrag = Teuchos::rcp( new BlockMultiVector_Type( this->getSolution() ) );
        BlockMultiVectorPtr_Type uLift = Teuchos::rcp( new BlockMultiVector_Type( this->getSolution() ) );
        // should be the last fixed point system without boundary conditions or the last extrapolation system without boundary values.
        // We need to reassemble B and BT, because we might have set Dirichlet boundary conditions in BT (less likely in B)
        this->assembleDivAndStab();
        this->reAssemble("FixedPoint");

        this->getSystem()->apply( *this->getSolution(), *uDrag );
        this->getSystem()->apply( *this->getSolution(), *uLift );
        
        //MultiVectorPtr_Type C1T_lambda = Teuchos::rcp( new MultiVector_Type( this->getSolution()->getBlock(0) ) );
        //this->system_->getBlock(0,3)->apply( *this->getSolution()->getBlock(3), *C1T_lambda );
        
        //uDrag->getBlockNonConst(0)->update( 1., *C1T_lambda, 1. ); // velocity + C1_T * lambda
        //uLift->getBlockNonConst(0)->update( 1., *C1T_lambda, 1. ); // velocity + C1_T * lambda
        
        BCPtr_Type bcFactoryDrag = Teuchos::rcp( new BC_Type( ) );
        BCPtr_Type bcFactoryLift = Teuchos::rcp( new BC_Type( ) );
        
        DomainConstPtr_Type domainVelocityConst = this->getDomain(0);
        DomainPtr_Type domainVelocity = Teuchos::rcp_const_cast<Domain_Type>(domainVelocityConst);
        if( dim == 2 ){
            bcFactoryDrag->addBC(drag2D, 4, 0, domainVelocity, "Dirichlet", dim); // obstacle
            bcFactoryDrag->addBC(drag2D, 5, 0, domainVelocity, "Dirichlet", dim); // interface; check main fsi for matching flags at the obstacle and interface
            bcFactoryLift->addBC(lift2D, 4, 0, domainVelocity, "Dirichlet", dim);
            bcFactoryLift->addBC(lift2D, 5, 0, domainVelocity, "Dirichlet", dim);
        }
        else if( dim == 3 ){
            bcFactoryDrag->addBC(drag3D, 4, 0, domainVelocity, "Dirichlet", dim); // check main fsi for matching
            bcFactoryDrag->addBC(drag3D, 6, 0, domainVelocity, "Dirichlet", dim); // check main fsi for matching flags at the obstacle and interface
            bcFactoryLift->addBC(lift3D, 4, 0, domainVelocity, "Dirichlet", dim);
            bcFactoryLift->addBC(lift3D, 6, 0, domainVelocity, "Dirichlet", dim);
        }
        
        BlockMultiVectorPtr_Type vD = Teuchos::rcp( new BlockMultiVector_Type( this->getSolution() ) );
        BlockMultiVectorPtr_Type vL = Teuchos::rcp( new BlockMultiVector_Type( this->getSolution() ) );
        
        vD->putScalar(0.);
        vL->putScalar(0.);
        
        bcFactoryDrag->setRHS( vD );
        bcFactoryLift->setRHS( vL );
        
        uDrag->dot( vD, drag() );
        uLift->dot( vL, lift() );
        
       double density = this->getParameterList()->sublist("Parameter").get("Density",1.);
       double D = 0.1;
       double H = 0.41;
       double uMean = this->getParameterList()->sublist("Parameter").get("MeanVelocity",1.0);
       
       drag[0] *= -(2./(density*uMean*uMean*D*H));
       lift[0] *= -(2./(density*uMean*uMean*D*H));

        // drag[0] *= -1.;
        // lift[0] *= -1.;
        
        exporterTxtDrag_->exportData( drag[0] );
        exporterTxtLift_->exportData( lift[0] );
    }

    if ( this->parameterList_->sublist("General").get("Export RE and CFL",true) ) {

        MultiVectorPtr_Type Re(new MultiVector_Type( this->getDomain(0)->getElementMap(), 1 ) );
        MultiVectorPtr_Type CFL(new MultiVector_Type( this->getDomain(0)->getElementMap(), 1 ) );

        this->feFactory_->assemblyCFLandRe(this->dim_,this->domain_FEType_vec_.at(0), this->parameterList_,CFL,Re,u_rep_);

        //Re->print();
        Teuchos::Array<SC> norm(1);
        CFL->normInf(norm());               
        exporterTxtCFLMax_->exportData( norm[0] );


        Teuchos::ArrayRCP<  SC > CFLArray = CFL->getDataNonConst(0);
        double minimum = 1000.;
        for (int i=0 ; i< CFLArray.size(); i++)
            if(minimum > CFLArray[i])
                minimum = CFLArray[i];

        reduceAll<int, double> (*this->getComm(), REDUCE_MIN, minimum, outArg (minimum));

        exporterTxtCFLMin_->exportData( minimum );

        Re->normInf(norm());               
        exporterTxtReMax_->exportData( norm[0] );
        // exporterTxtReMin_->exportData( reMin[0] );
        if ( exporterRe_.is_null() && this->parameterList_->sublist("General").get("Plot RE and CFL",false)){
            exporterRe_ = Teuchos::rcp(new Exporter_Type());
            
            DomainConstPtr_Type dom = this->getDomain(0);
            std::string varName = "RE_CFL_Element";
            
            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( dom->getMesh() );
            exporterRe_->setup(varName, meshNonConst,"P0");

            MultiVectorConstPtr_Type exportVector = Re;
            
            exporterRe_->addVariable( exportVector, "RE", "Scalar", 1, dom->getElementMap() );

            MultiVectorConstPtr_Type exportVector2 = CFL;
            exporterRe_->addVariable( exportVector2, "CFL", "Scalar", 1, dom->getElementMap() );

        }

        if (!exporterRe_.is_null()){
            MultiVectorConstPtr_Type exportVector = Re;
            MultiVectorConstPtr_Type exportVector2 = CFL;

            this->exporterRe_->updateVariables(exportVector, "RE");
            this->exporterRe_->updateVariables(exportVector2, "CFL");
            
            double exportTime = 0.0;
            if(this->parameterList_->sublist("Timestepping Parameter").get("dt",-1.)> 0)
                exportTime=this->timeSteppingTool_->currentTime() ;

            this->exporterRe_->save(exportTime );
        
        }

    }
}

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
    else if ( !type.compare("Teko") || !type.compare("Diagonal") || !type.compare("PCD") || !type.compare("LSC")){
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

            this->calculateNonLinResidualVec("standard"); // Calculating residual Vector

			// Changing the residualVector into a ThyraMultivector

            Teuchos::RCP<Thyra::MultiVectorBase<SC> > f_thyra = this->getResidualVector()->getThyraMultiVector();
            f_out->assign(*f_thyra);
        }

        XpetraMatrixPtr_Type W;
        if (fill_W) {

            this->reAssemble("Newton"); // ReAssembling matrices with updated u  in this class

            this->setBoundariesSystem(); // setting boundaries to the system
            // this->bcFactory_->setDirichletColumn(this->getSystem()->getBlock(1,0),false);
            // this->bcFactory_->setDirichletColumn(this->getSystem()->getBlock(0,0),true);
			// Changing the system Matrix into a tpetra Matrix (block matrices have 'getXpetraMatrix' feature)
            Teuchos::RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
            Teuchos::RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);

            XpetraMatrixConstPtr_Type W_systemXpetra = this->getSystem()->getMergedMatrix()->getXpetraMatrix(); // The current system matrix of this class

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

        if (fill_W_prec ) {
        
            if (stokesMonoPrecUsed_){
                int newtonLimit = this->parameterList_->sublist("Parameter").get("newtonLimit",2);
                if(this->newtonStep_ < newtonLimit || this->parameterList_->sublist("Parameter").get("Rebuild Preconditioner every Newton Iteration",true) )
                {
                    this->setupPreconditioner( "Monolithic" );
                }
                else{
                    if (this->verbose_)
                        cout << " ############ Skipping preconditioner reconstruction #############" << endl;
                }
            }
            else
                stokesMonoPrecUsed_ = true;
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
            this->newtonStep_ ++; 


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

            Teko::MultiVector f0;
            Teko::MultiVector f1;
            f0 = this->getResidualVector()->getBlockNonConst(0)->getThyraMultiVector();
            f1 = this->getResidualVector()->getBlockNonConst(1)->getThyraMultiVector();

            std::vector<Teko::MultiVector> f_vec; f_vec.push_back(f0); f_vec.push_back(f1);

            Teko::MultiVector f = Teko::buildBlockedMultiVector(f_vec);

            f_out->assign(*f);
        }

        XpetraMatrixPtr_Type W;
        if (fill_W) {

            typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;

            this->reAssemble("Newton");

            this->setBoundariesSystem();
            // this->bcFactory_->setDirichletColumn(this->getSystem()->getBlock(1,0),false);
            // this->bcFactory_->setDirichletColumn(this->getSystem()->getBlock(0,0),true);
            RCP<ThyraBlockOp_Type> W_blocks = rcp_dynamic_cast<ThyraBlockOp_Type>(W_out);
            RCP<const ThyraOp_Type> W_block00 = W_blocks->getBlock(0,0);
            RCP<ThyraOp_Type> W_block00NonConst = rcp_const_cast<ThyraOp_Type>( W_block00 );
            RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator( W_block00NonConst );

            RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);

            XpetraMatrixConstPtr_Type W_matrixXpetra = this->getSystem()->getBlock(0,0)->getXpetraMatrix();
            XpetraMatrixPtr_Type W_matrixXpetraNonConst = rcp_const_cast<XpetraMatrix_Type>(W_matrixXpetra);
            Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_matrixXpetraNonConst);
            Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            RCP<TpetraMatrix_Type> tpetraMatXpetra = xTpetraMat.getTpetra_CrsMatrixNonConst();

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
            std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");

            if (stokesTekoPrecUsed_){
                int newtonLimit = this->parameterList_->sublist("Parameter").get("newtonLimit",2);
                if(this->newtonStep_ < newtonLimit || this->parameterList_->sublist("Parameter").get("Rebuild Preconditioner every Newton Iteration",true) )
                {
                    this->setupPreconditioner( type );
                }
                else{
                    if (this->verbose_)
                        cout << " ############ Skipping preconditioner reconstruction #############" << endl;
                }
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
            this->newtonStep_ ++; 

        }
    }
}
#endif

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
    this->reAssemble("FixedPoint");
    this->reAssemble("Newton");

    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( !type.compare("Monolithic"))
        return create_W_op_Monolithic( );
    else if ( !type.compare("Teko") || !type.compare("Diagonal") || !type.compare("PCD") || !type.compare("LSC") ){
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
    
    // BlockMatrixPtr_Type system = this->getSystem();
    
    // Teuchos::RCP<const ThyraBlockOp_Type> W_opBlocksConst = system->getThyraLinBlockOp();
    // Teuchos::RCP<ThyraBlockOp_Type> W_opBlocks = Teuchos::rcp_const_cast<ThyraBlockOp_Type >(W_opBlocksConst);
    // Teuchos::RCP<ThyraOp_Type> W_op = Teuchos::rcp_dynamic_cast<ThyraOp_Type >(W_opBlocks);



    return W_op;
}
#endif

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > NavierStokes<SC,LO,GO,NO>::create_W_prec() const
{

    this->initializeSolverBuilder();

    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    this->setBoundariesSystem();
    // this->bcFactory_->setDirichletColumn(this->getSystem()->getBlock(1,0),false);
    // this->bcFactory_->setDirichletColumn(this->getSystem()->getBlock(0,0),true);

    if (!type.compare("Teko") || !type.compare("Diagonal") || !type.compare("Triangular") || !type.compare("PCD") || !type.compare("LSC")) { //
        this->setupPreconditioner( type );
        stokesTekoPrecUsed_ = false;
    }
    else{
        this->setupPreconditioner( type ); // initializePreconditioner( type );
        stokesMonoPrecUsed_ = false;
   }
    

    Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  this->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst = Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);

    return thyraPrecNonConst;

}
}

#endif
