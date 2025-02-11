#ifndef STOKES_def_hpp
#define STOKES_def_hpp
#include "Stokes_decl.hpp"
/*!
 Definition of Stokes
 
 @brief Stokes
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
void ZeroDirichlet(double* x, double* res, double t, double* parameters){
    
    res[0] = 0.;
    
    return;
}
    
double OneFunction(double* x, int* parameter)
{
    return 1.0;
}
    
template<class SC,class LO,class GO,class NO>
Stokes<SC,LO,GO,NO>::Stokes(const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList):
Problem<SC,LO,GO,NO>(parameterList, domainVelocity->getComm())
{

    this->addVariable( domainVelocity , FETypeVelocity , "u" , domainVelocity->getDimension());
    this->addVariable( domainPressure , FETypePressure , "p" , 1);
    
    this->dim_ = this->getDomain(0)->getDimension();
}

template<class SC,class LO,class GO,class NO>
Stokes<SC,LO,GO,NO>::~Stokes(){

}

template<class SC,class LO,class GO,class NO>
void Stokes<SC,LO,GO,NO>::info(){
    this->infoProblem();
}
    
template<class SC,class LO,class GO,class NO>
void Stokes<SC,LO,GO,NO>::assemble( std::string type ) const{
    
    if (this->verbose_)
        std::cout << "-- Assembly ... " << std::flush;

    double viscosity = this->parameterList_->sublist("Parameter").get("Viscosity",1.);

    MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
    MatrixPtr_Type BT(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
    
    MapConstPtr_Type pressureMap;
    if ( this->getDomain(1)->getFEType() == "P0" )
        pressureMap = this->getDomain(1)->getElementMap();
    else
        pressureMap = this->getDomain(1)->getMapUnique();

    MatrixPtr_Type B(new Matrix_Type( pressureMap, this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    
    MatrixPtr_Type C;
    if (this->verbose_)
        std::cout << " A ... " << std::flush;
    int* dummy;
    if ( this->parameterList_->sublist("Parameter").get("Symmetric gradient",false) )
        this->feFactory_->assemblyStress(this->dim_, this->domain_FEType_vec_.at(0), A, OneFunction, dummy, true);
    else
        this->feFactory_->assemblyLaplaceVecField(this->dim_, this->domain_FEType_vec_.at(0), 2, A, true);

    if (this->verbose_)
        std::cout << "B and B^T ... " << std::flush;

    this->feFactory_->assemblyDivAndDivT(this->dim_, this->getFEType(0), this->getFEType(1), 2, B, BT, this->getDomain(0)->getMapVecFieldUnique(), pressureMap, true );
    
    A->resumeFill();
    B->resumeFill();
    BT->resumeFill();
    
    A->scale(viscosity);
    B->scale(-1.);
    BT->scale(-1.);

    A->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
    B->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), pressureMap );
    BT->fillComplete( pressureMap, this->getDomain(0)->getMapVecFieldUnique() );
    
    this->system_.reset(new BlockMatrix_Type(2));
  
    this->system_->addBlock( A, 0, 0 );
    this->system_->addBlock( BT, 0, 1 );
    this->system_->addBlock( B, 1, 0 );
    
//    this->initializeVectors();
    
    if ( !this->getFEType(0).compare("P1") ) {
        C.reset(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        this->feFactory_->assemblyBDStabilization( this->dim_, "P1", C, true);
        C->resumeFill();
        C->scale( -1./viscosity );
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
        
        Mpressure->resumeFill();
        Mpressure->scale(-1./viscosity);
        Mpressure->fillComplete( pressureMap, pressureMap );
        this->getPreconditionerConst()->setPressureMassMatrix( Mpressure );
    }

    this->assembleSourceTerm( 0. );
    this->addToRhs( this->sourceTerm_ );
    //this->rhs_->print();
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
    
}

}
#endif
