#ifndef NonLinElasticity_decl_hpp
#define NonLinElasticity_decl_hpp
#include "feddlib/problems/abstract/NonLinearProblem.hpp"
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_ModelEvaluatorBase_decl.hpp>
/*!
 Declaration of NonLinElasticity
 
 @brief NonLinElasticity
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD{
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class NonLinElasticity : public NonLinearProblem<SC,LO,GO,NO>  {
    
public:
    //! @name Public Types
    //@{
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;
    
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;
    
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;
    
    typedef NonLinearProblem<SC,LO,GO,NO> NonLinearProblem_Type;
    typedef typename NonLinearProblem_Type::BlockMultiVectorPtrArray_Type BlockMultiVectorPtrArray_Type;
    
    typedef typename NonLinearProblem_Type::TpetraMatrix_Type TpetraMatrix_Type;
    
    typedef typename NonLinearProblem_Type::ThyraVecSpace_Type ThyraVecSpace_Type;
    typedef typename NonLinearProblem_Type::ThyraVec_Type ThyraVec_Type;
    typedef typename NonLinearProblem_Type::ThyraOp_Type ThyraOp_Type;
    
    typedef typename NonLinearProblem_Type::TpetraOp_Type TpetraOp_Type;

    typedef Tpetra::CrsMatrix<SC, LO, GO, NO> tpetra_matrix;

    
    //@}
    
    //! @name Constructor/Destructor
    //@{
    NonLinElasticity( const DomainConstPtr_Type  &domain, std::string FEType, ParameterListPtr_Type parameterList );
    //@}
    ~NonLinElasticity();

    virtual void info();
    
    virtual void assemble( std::string type = "" ) const;

    void reAssemble(std::string type) const;

    void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const override{}
    
    virtual void reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const;
    
    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions);
    
    virtual void calculateNonLinResidualVec(std::string type, double time=0.) const;
    
    void getValuesOfInterest( vec_dbl_Type& values ) override {}
    
    void computeValuesOfInterestAndExport() override {}
    
//    virtual void assembleExternal( std::string type ){}
        
private:
    
    mutable MultiVectorPtr_Type u_rep_;
    double E_;
    double mue_;
    double C_;
    double poissonRatio_;
    double lambda_;
    /*####################*/

};
}
#endif
