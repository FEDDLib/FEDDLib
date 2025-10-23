#ifndef LINELASASSFE_def_hpp
#define LINELASASSFE_def_hpp
#include "LinElas_decl.hpp"
namespace FEDD {


template<class SC,class LO,class GO,class NO>
LinElasAssFE<SC,LO,GO,NO>::LinElasAssFE(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList):
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
LinElasAssFE<SC,LO,GO,NO>::~LinElasAssFE()
{

}

template<class SC,class LO,class GO,class NO>
void LinElasAssFE<SC,LO,GO,NO>::info(){
    this->infoProblem();
}

template<class SC,class LO,class GO,class NO>
void LinElasAssFE<SC,LO,GO,NO>::assemble( std::string type ) const
{
    if(this->verbose_)
        std::cout << "-- Assembly linear elasticity ... " << std::flush;
    // Initialisiere die Steifigkeitsmatrix. Das letzte Argument gibt die (ungefaehre) Anzahl an Eintraege pro Zeile an
    MatrixPtr_Type K = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    // MatrixPtr_Type K = Teuchos::rcp(new Matrix_Type( this->domainPtr_vec_.at(0)->getMapVecFieldUnique(), 10 ) );

	MultiVectorConstPtr_Type d = this->solution_->getBlock(0);

    d_rep_->importFromVector(d, true);

    this->system_->addBlock( K, 0, 0 );
    // Assembliere die Steifigkeitsmatrix. Die 2 gibt degree an, d.h. die Ordnung der Quadraturformel, die benutzt werden soll.
    this->feFactory_->assemblyLinearElasticity(this->dim_, this->getDomain(0)->getFEType(),2, this->dim_, d_rep_, this->system_, this->rhs_, this->parameterList_,false, "Jacobian", true);

    // Setup fuer die linke Seite des zu loesdenen GLS. Beachte, dass system_ via system_() im Standardkonstruktor (von der Klasse Problem)
    // initialisiert worden ist, hier wird dieser also erst richtig initialisiert.
    // Ein Objekt der Klasse Bmat ist eine Blockmatrix; also ist system_ eine Blockmatrix (Objekt von BMat)
     
    double density = this->parameterList_->sublist("Parameter").get("Density",1000.);
    std::string sourceType = this->parameterList_->sublist("Parameter").get("Source Type","volume");
    this->assembleSourceTerm( 0. );
    if(sourceType == "volume")
        this->sourceTerm_->scale(density);
    this->addToRhs( this->sourceTerm_ );
    
    if (this->verbose_)
        std::cout << "done -- " << std::endl;
}

}
#endif
