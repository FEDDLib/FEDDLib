#ifndef LaplaceBlocks_decl_hpp
#define LaplaceBlocks_decl_hpp
#include "feddlib/problems/abstract/Problem.hpp"
/*!
 Declaration of LaplaceBlocks
 
 @brief LaplaceBlocks
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class LaplaceBlocks : public Problem<SC,LO,GO,NO>  {
    
public:
    
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    LaplaceBlocks( const DomainConstPtr_Type &domain1, const DomainConstPtr_Type &domain2, std::string FEType1, std::string FEType2, ParameterListPtr_Type parameterList );
    
    ~LaplaceBlocks();
    
    virtual void info();
    
    virtual void assemble( std::string type = "" ) const;
    
    virtual void getValuesOfInterest( vec_dbl_Type& values ){}
    
    virtual void computeValuesOfInterestAndExport() {}
//    virtual int SetupPreconditioner(BMat_ptr_Type systemPrec, ThyraConstLinOpPtr_Type thyraMatrix=Teuchos::null, ThyraPrecPtr_Type thyraPreconditioner = Teuchos::null, LinSolverBuilderPtr_Type linearSolverBuilder = Teuchos::null) const;
//    
//    virtual void assembleExternal( std::string type ){}
    
    void setPrecLists( ParameterListPtr_Type p1list, ParameterListPtr_Type p2list ){ precList1_ = p1list; precList2_ = p2list;}
    
    ParameterListPtr_Type getPrecList1(){return precList1_;}
    
    ParameterListPtr_Type getPrecList2(){return precList2_;}
    
private:
    ParameterListPtr_Type precList1_;
    ParameterListPtr_Type precList2_;
    /*####################*/
};
}
#endif
