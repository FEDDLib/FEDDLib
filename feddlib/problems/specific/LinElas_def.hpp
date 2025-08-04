#ifndef LINELAS_def_hpp
#define LINELAS_def_hpp
#include "LinElas_decl.hpp"
namespace FEDD {
// Funktion fuer den homogenen Dirichletrand
// void ZeroDirichlet(double* x, double* result, double t, double* parameters)
// {
//     result[0] = 0.;
//     return;
// }

// TODO: Fuer FSI auf 0 und 0 abaendern!!!!!!!!!!
// Funktion fuer die rechte Seite der DGL in 2D
void rhsFunc2D(double* x, double* result, double* parameters)
{
    // Wir setzen die rechte Seite g_vec als g_vec = (0, g), mit g = -2.
    result[0] = 0.;
    result[1] = 0.;
//    if (parameters[0]<0.1) {
//        result[1] = -.2;
//    }

    return;
}

// Funktion fuer die rechte Seite der DGL in 3D
void rhsFunc3D(double* x, double* result, double* parameters)
{
    // Wir setzen die rechte Seite g_vec als g_vec = (0, g, 0), mit g = -2.
    result[0] = 0.;
    // result[1] = -2.;
    result[1] = 0.;
    result[2] = 0.;
    return;
}


template<class SC,class LO,class GO,class NO>
LinElas<SC,LO,GO,NO>::LinElas(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList):
Problem_Type(parameterList, domain->getComm() )
{
    // Bem.: Hier benutzen wir auch direkt den Konstruktor der Klasse Problem, wo z.B. parameterList auf parameterList_ gesetzt wird

    // d_s steht displacement (d) der Struktur (_s). Im Allgemeinen vektorwertig, daher domainDisplacement->GetDimension().
    // Andernfalls benutzen wir stattdessen 1

    // Siehe Problem_def.hpp fuer die Funktion AddVariable()
    this->addVariable( domain , FEType , "d_s" ,domain->getDimension() );
    this->dim_ = this->getDomain(0)->getDimension();

    d_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );

}

template<class SC,class LO,class GO,class NO>
LinElas<SC,LO,GO,NO>::~LinElas()
{

}

template<class SC,class LO,class GO,class NO>
void LinElas<SC,LO,GO,NO>::info(){
    this->infoProblem();
}

template<class SC,class LO,class GO,class NO>
void LinElas<SC,LO,GO,NO>::assemble( std::string type ) const
{
    if(this->verbose_)
        std::cout << "-- Assembly linear elasticity ... " << std::flush;

    // Hole die Dichte \rho (density) und die Paramter \nu (Poisson-ratio) und \mu (zweite Lamé-Konstante)
    double density = this->parameterList_->sublist("Parameter").get("Density",1000.);
    
    double poissonRatio = this->parameterList_->sublist("Parameter").get("Poisson Ratio",0.4);
    double mu = this->parameterList_->sublist("Parameter").get("Mu",2.0e+6);

    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstanten \lambda
    double youngModulus = mu*2.*(1 + poissonRatio);

    double lambda = (poissonRatio*youngModulus)/((1 + poissonRatio)*(1 - 2*poissonRatio));

    // Initialisiere die Steifigkeitsmatrix. Das letzte Argument gibt die (ungefaehre) Anzahl an Eintraege pro Zeile an
    MatrixPtr_Type K = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    // MatrixPtr_Type K = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapVecFieldUnique(), 10 ) );

    // Assembliere die Steifigkeitsmatrix. Die 2 gibt degree an, d.h. die Ordnung der Quadraturformel, die benutzt werden soll.
#ifdef FEDD_HAVE_ACEGENINTERFACE
    bool useInterface = this->parameterList_->sublist("Parameter").get("Use AceGen Interface", true);
    if(this->getFEType(0) =="P2" && useInterface){
        MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
        d_rep_->importFromVector(d, true);

        this->system_->addBlock( K, 0, 0 );
        // Assembliere die Steifigkeitsmatrix. Die 2 gibt degree an, d.h. die Ordnung der Quadraturformel, die benutzt werden soll.
        this->feFactory_->assemblyLinearElasticity(this->dim_, this->getDomain(0)->getFEType(),2, this->dim_, d_rep_, this->system_, this->rhs_, this->parameterList_,false, "Jacobian", true);
    }
    else 
        this->feFactory_->assemblyLinElasXDim( this->dim_, this->getDomain(0)->getFEType(), K, lambda, mu );
#else
    this->feFactory_->assemblyLinElasXDim( this->dim_, this->getDomain(0)->getFEType(), K, lambda, mu );
#endif
    // Setup fuer die linke Seite des zu loesdenen GLS. Beachte, dass system_ via system_() im Standardkonstruktor (von der Klasse Problem)
    // initialisiert worden ist, hier wird dieser also erst richtig initialisiert.
    // Ein Objekt der Klasse Bmat ist eine Blockmatrix; also ist system_ eine Blockmatrix (Objekt von BMat)
    
    // Fuege die Steifikeitsmatrix als Blockeintrag an der Stelle (1,1) (in C dann (0,0)) in die Blockmatrix hinein.
    this->system_->addBlock( K, 0, 0 );
    
    this->assembleSourceTerm( 0. );
    //this->sourceTerm_->scale(density);
    this->addToRhs( this->sourceTerm_ );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

//template<class SC,class LO,class GO,class NO>
//void LinElas<SC,LO,GO,NO>::assembleSourceTerm(double time)
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
//    this->sourceTerm_->getBlockNonConst(0)->exportFromVector( sourceTerm, false, "Add" );
//
//    double density = this->parameterList_->sublist("Parameter").get("Density",1.);
//
//    this->sourceTerm_->getBlockNonConst(0)->scale(density);
//    // Noch mit \rho = density skalieren
//    // feRhs->scale(density);
//    // BlockMultiVectorPtr_Type blockMVSourceT = Teuchos::rcp( new BlockMultiVector_Type( 1 ) );
//    // blockMVSourceT->addBlock( feRhs, 0 );
//
//    // Setze die rechte Seite auf das Attribut des SourceTerms aus dem Problem
//    // Dieses Attribut wird in der Zeitintegration benoetigt, vgl. DAESolverInTime
//    // this->sourceTerm_ = blockMVSourceT;
//
////    MultiVectorPtr_Type feRhsUnique = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldUnique() ) );
////
////    feRhsUnique->importFromVector( feRhs, false, "Add" );
////
////    // Noch mit \rho = density skalieren
////    feRhsUnique->scale(density);
////
////    this->sourceTerm_->addBlock( feRhsUnique, 0 );
//}

}
#endif
