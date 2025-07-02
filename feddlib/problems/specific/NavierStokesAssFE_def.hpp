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
