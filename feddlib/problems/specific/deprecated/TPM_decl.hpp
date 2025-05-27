#ifndef TPM_decl_hpp
#define TPM_decl_hpp
#include "feddlib/problems/abstract/Problem.hpp"
/*!
 Declaration of TPM
 
 @brief TPM
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class TPM : public Problem<SC,LO,GO,NO> {
    
public:
    
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    typedef typename Problem_Type::Map_Type Map_Type;    
    typedef typename Problem_Type::MapPtr_Type MapPtr_Type;
    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;

    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;
    
    TPM( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList );
    
    ~TPM();
    
    virtual void info();
    
    virtual void assemble( std::string type ) const;
    
//    virtual void assembleExternal( std::string type );

    virtual void getValuesOfInterest( vec_dbl_Type& values ){};
    
    virtual void computeValuesOfInterestAndExport() {};
    
  private:
    mutable MultiVectorPtr_Type u_repNewton_;
    mutable MultiVectorPtr_Type p_repNewton_;
    mutable MultiVectorPtr_Type u_repTime_;
    mutable MultiVectorPtr_Type p_repTime_;
    /*####################*/

};
}
#endif
