#ifndef GEOMETRY_decl_hpp
#define GEOMETRY_decl_hpp
#include "feddlib/problems/abstract/Problem.hpp"

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Geometry : public Problem<SC,LO,GO,NO>  {


public:
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;


    // Konstruktor
    Geometry(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList);

    // Destruktor
    ~Geometry();

    virtual void info();
    
    // Assemblierung des Fortsetzungsoperators; zunaechst erstmal diskret harmonisch mit heuristischer Skalierung
    virtual void assemble( std::string type = "" ) const;

    virtual void getValuesOfInterest( vec_dbl_Type& values ){}
    
    virtual void computeValuesOfInterestAndExport() {}

//    virtual void assembleExternal( std::string type ){}
private:

};
}
#endif
