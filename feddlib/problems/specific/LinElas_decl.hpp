#ifndef LINELAS_decl_hpp
#define LINELAS_decl_hpp
#include "feddlib/problems/abstract/Problem.hpp"

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class LinElas : public Problem<SC,LO,GO,NO>  {

public:
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;

    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;
    
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;


    // Konstruktor
    LinElas(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList);
    
    // Destruktor
    ~LinElas();

    virtual void info();
    
    // Assemblierung der Matrix, hier: \sigma = 2\mu \epsilon + \lambda \div(u) I, mit \epsilon = 0.5 ( \grad u + (\grad u)^T )
    virtual void assemble( std::string type = "" ) const;

    // Assemblierung der rechten Seite der DGL (RHS), welche in AdvanceinTime() genutzt wird;
    // also fuer zeitabhaengige Probleme ist
    // Beachte: In AdvanceTime() wird in jedem Zeitschritt der SourceTerm neu berechnet, auch
    // wenn dieser konstant ist!
    // TODO: Baue in Problem_def.hpp sowas wie boolHasTimeDependentSourceTerm ein, mit default auf false.
//    void assembleSourceTerm(double time);

    virtual void getValuesOfInterest( vec_dbl_Type& values ){}
    
    virtual void computeValuesOfInterestAndExport() {}
    // Steifigkeitsmatrix des Problems der linearen Elastizitaet gegeben wie in assemble().
    // Moeglicherweise nicht noetig (vgl. Laplace.hpp)
    // Falls es doch irgendwann benutzt wird, denke daran den Konstruktor zu aendern (vgl. Stokes.hpp)
    // Teuchos::RCP<Matrix_Type> 	K_;
//    virtual void assembleExternal( std::string type ){}
private:
    MultiVectorPtr_Type d_rep_;

};
}

#endif
