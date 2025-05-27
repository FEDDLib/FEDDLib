#ifndef LinElasFirstOrder_def_hpp
#define LinElasFirstOrder_def_hpp
#include "LinElasFirstOrder_decl.hpp"
namespace FEDD {

template<class SC,class LO,class GO,class NO>
LinElasFirstOrder<SC,LO,GO,NO>::LinElasFirstOrder(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList):
Problem_Type(parameterList, domain->getComm())
{
    // Bem.: Hier benutzen wir auch direkt den Konstruktor der Klasse Problem, wo z.B. parameterList auf parameterList_ gesetzt wird

    // d_s steht displacement (d) der Struktur (_s). Im Allgemeinen vektorwertig, daher domainDisplacement->GetDimension().
    // Andernfalls benutzen wir stattdessen 1

    // Siehe Problem_def.hpp fuer die Funktion AddVariable()
    this->addVariable( domain , FEType , "d_s" ,domain->getDimension() );
    this->addVariable( domain , FEType , "v_s" ,domain->getDimension() );
    this->dim_ = this->getDomain(0)->getDimension();
}

template<class SC,class LO,class GO,class NO>
LinElasFirstOrder<SC,LO,GO,NO>::~LinElasFirstOrder()
{

}

template<class SC,class LO,class GO,class NO>
void LinElasFirstOrder<SC,LO,GO,NO>::info(){
    this->infoProblem();
}

template<class SC,class LO,class GO,class NO>
void LinElasFirstOrder<SC,LO,GO,NO>::assemble( std::string type ) const
{
    if(this->verbose_)
        std::cout << "-- Assembly linear elasticity first order ... " << std::flush;

    // Hole die Dichte \rho (density) und die Paramter \nu (Poisson-ratio) und \mu (zweite Lamé-Konstante)
    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
    double poissonRatio = this->parameterList_->sublist("Parameter").get("Poisson Ratio",0.4);
    double mu = this->parameterList_->sublist("Parameter").get("Mu",2.0e+6);

    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstanten \lambda
    double youngModulus = mu*2.*(1 + poissonRatio);
    double lambda = (poissonRatio*youngModulus)/((1 + poissonRatio)*(1 - 2*poissonRatio));

    // Initialisiere die Steifigkeitsmatrix. Das letzte Argument gibt die (ungefaehre) Anzahl an Eintraege pro Zeile an
    MatrixPtr_Type K = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
    // MatrixPtr_Type K = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapVecFieldUnique(), 10 ) );

    // Assembliere die Steifigkeitsmatrix. Die 2 gibt degree an, d.h. die Ordnung der Quadraturformel, die benutzt werden soll.
    this->feFactory_->assemblyLinElasXDim( this->dim_, this->getDomain(0)->getFEType(), K, lambda, mu );
    K->resumeFill();
    K->scale(-1.);
    K->fillComplete();
    
    MatrixPtr_Type I = Teuchos::rcp(new Matrix_Type( this->getDomain(1)->getMapVecFieldUnique(), 1 ) );
    this->feFactory_->assemblyIdentity( I );
    I->resumeFill();
    I->scale(-density); // use density aswell, since every diagonal mass matrix in TimeProblem will be scaled with density aswell
    I->fillComplete();

    // Fuege die Steifikeitsmatrix als Blockeintrag an der Stelle (1,1) (in C dann (0,0)) in die Blockmatrix hinein.
    this->system_->addBlock( K, 0, 0 );
    this->system_->addBlock( I, 1, 1 );
    
    // Initialisiere den Loesungsvektor mit der DomainMap und die rechte Seite RHS mit der RangeMap.
    // int nmbVectors = 1 in Problem_decl.hpp
//    this->initializeVectors();
    
    this->assembleSourceTerm( 0. );

    this->rhs_->addBlock( this->sourceTerm_->getBlockNonConst(0), 0 );

    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

//template<class SC,class LO,class GO,class NO>
//void LinElasFirstOrder<SC,LO,GO,NO>::assembleSourceTerm(double time)
//{
//    
//    double force = this->parameterList_->sublist("Parameter").get("Volume force",0.);
//    
//    MultiVectorPtr_Type sourceTerm = Teuchos::rcp(new MultiVector_Type( this->domainPtr_vec_.at(0)->getMapVecFieldRepeated() ) );
//    vec_dbl_Type funcParameter(3,0.);
//    
//    funcParameter[0] = 0.; // degree of function
//    funcParameter[1] = time;
//    funcParameter[2] = force;
//    
//    this->feFactory_->assemblyRHS(this->dim_, this->domain_FEType_vec_.at(0), sourceTerm, "Vector", this->rhsFuncVec_[0], funcParameter);
//    
//    this->sourceTerm_->getBlockNonConst(0)->putScalar(0.);
//
//    this->sourceTerm_->getBlockNonConst(0)->exportFromVector( sourceTerm, true, "Add" );
//
//    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
//
//    this->sourceTerm_->getBlockNonConst(0)->scale(density);
// }

}
#endif
