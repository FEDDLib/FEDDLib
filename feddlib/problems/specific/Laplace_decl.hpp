#ifndef LAPLACE_decl_hpp
#define LAPLACE_decl_hpp
#include "feddlib/problems/abstract/Problem.hpp"
/*!
 Declaration of Laplace
 
 @brief Laplace
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Laplace : public Problem<SC,LO,GO,NO>  {
    
public:
    
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    Laplace(const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList, bool vectorLaplace = false );
    
    ~Laplace();
    
    virtual void info();
    
    virtual void assemble( std::string type = "" ) const;
    
    virtual void getValuesOfInterest( vec_dbl_Type& values ){}

	MatrixPtr_Type getMassMatrix() const; // new for calculating L2-Error
    
    virtual void computeValuesOfInterestAndExport() {}
//    virtual int SetupPreconditioner(BMat_ptr_Type systemPrec, ThyraConstLinOpPtr_Type thyraMatrix=Teuchos::null, ThyraPrecPtr_Type thyraPreconditioner = Teuchos::null, LinSolverBuilderPtr_Type linearSolverBuilder = Teuchos::null) const;

private:
    /*####################*/
    bool vectorLaplace_;
};
}
#endif
