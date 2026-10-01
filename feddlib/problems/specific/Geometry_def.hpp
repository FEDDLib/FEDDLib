#ifndef GEOMETRY_def_hpp
#define GEOMETRY_def_hpp

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"

namespace FEDD {
// Funktion fuer die rechte Seite der DGL in 2D
void ZeroFErhsFunc2D(double* x, double* result, double* parameters)
{
    // Wir setzen die rechte Seite f_vec = (0, 0)
    result[0] = 0.0;
    result[1] = 0.0;
    return;
}

// Funktion fuer die rechte Seite der DGL in 3D
void ZeroFErhsFunc3D(double* x, double* result, double* parameters)
{
    // Wir setzen die rechte Seite f_vec = (0, 0, 0)
    result[0] = 0.0;
    result[1] = 0.0;
    result[2] = 0.0;
    return;
}

double HeuristicScaling(double* x, double* parameter)
{
    if(x[0] < parameter[0]) // Wenn die Distanz des Schwerpunkts des Elements zum Interface kleiner als 0.03 ist
    {
        return parameter[1];
    }
    else
    {
        return 1.0;
    }
}


template<class SC,class LO,class GO,class NO>
Geometry<SC,LO,GO,NO>::Geometry(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList):
Problem_Type(parameterList,domain->getComm())
{
    this->addVariable( domain , FEType , "d_f" , domain->getDimension());
    this->dim_ = this->getDomain(0)->getDimension();
}


template<class SC,class LO,class GO,class NO>
Geometry<SC,LO,GO,NO>::~Geometry()
{

}

template<class SC,class LO,class GO,class NO>
void Geometry<SC,LO,GO,NO>::info(){
    this->infoProblem();
}
    
template<class SC,class LO,class GO,class NO>
void Geometry<SC,LO,GO,NO>::assemble( std::string type ) const
{
    
    MatrixPtr_Type H = Teuchos::rcp( new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

    if (this->parameterList_->sublist("Parameter").get("Model","Laplace")=="Laplace"){
        if (this->verbose_)
            std::cout << "-- Assembly Geometry (scaled Laplace)... " << std::flush;
        double distanceLaplace = this->parameterList_->sublist("Parameter").get("Distance Laplace",0.03);
        double coefficientLaplace = this->parameterList_->sublist("Parameter").get("Coefficient Laplace",1000.);
        if (this->verbose_){
                std::cout << "\n    Distance Laplace = " << distanceLaplace << "  coefficient Laplace = " << coefficientLaplace << " ... " <<  std::flush;
        }
        vec_dbl_Type parameter(2);
        parameter[0] = distanceLaplace;
        parameter[1] = coefficientLaplace;
        
        this->feFactory_->assemblyLaplaceXDim(this->dim_, this->getDomain(0)->getFEType(), H, HeuristicScaling, &(parameter.at(0)) );
    }
    else if (this->parameterList_->sublist("Parameter").get("Model","Laplace")=="Elasticity"){
        if (this->verbose_)
            std::cout << "-- Assembly Geometry (linear elasticity)... " << std::flush;
        double poissonRatio = this->parameterList_->sublist("Parameter").get("Poisson Ratio",0.3);
        double mu = this->parameterList_->sublist("Parameter").get("Mu",2.0e+6);
        double youngModulus = mu*2.*(1 + poissonRatio);
        double lambda = (poissonRatio*youngModulus)/((1 + poissonRatio)*(1 - 2*poissonRatio));
        if (this->verbose_){
            std::cout << "\n    Poisson ratio = " << poissonRatio << "  mu = "<<mu << "  lambda = "<< lambda << "  E = " << youngModulus << "... " << std::flush;
        }
        
        this->feFactory_->assemblyLinElasXDim( this->dim_, this->getDomain(0)->getFEType(), H, lambda, mu );
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown model for geometry.");
    
    this->system_->addBlock( H, 0, 0 );
    
    this->assembleSourceTerm( 0./*time*/ );

    this->addToRhs( this->sourceTerm_ );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

}
#endif
