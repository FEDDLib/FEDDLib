#ifndef NonLinTPM_decl_hpp
#define NonLinTPM_decl_hpp
#include "feddlib/problems/abstract/Problem.hpp"
/*!
 Declaration of NonLinTPM
 
 @brief NonLinTPM
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class NonLinTPM : public NonLinearProblem<SC,LO,GO,NO>  {
    
public:
    
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    typedef typename Problem_Type::Map_Type Map_Type;    
    typedef typename Problem_Type::MapPtr_Type MapPtr_Type;
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

    NonLinTPM( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList );
    
    ~NonLinTPM();
    
    virtual void info();
    
    virtual void assemble(std::string type = "") const;    
        
    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const;
    
    virtual void reAssemble( MatrixPtr_Type& massmatrix, std::string type ) const;// not needed
        
    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions);// not needed
    
    virtual void calculateNonLinResidualVec(std::string type, double time=0.) const;
    
    virtual void getValuesOfInterest( vec_dbl_Type& values ){}
    
    virtual void computeValuesOfInterestAndExport() {}
    
//    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() const;// not needed
    
//    Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() const;// not needed
    
//    virtual void assembleExternal( std::string type ) {}
    
    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() const;
    
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() const;
    
private:
    
    virtual void evalModelImpl(
                               const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                               const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                               ) const;
    
    mutable MultiVectorPtr_Type u_repNewton_;
    mutable MultiVectorPtr_Type p_repNewton_;
    mutable MultiVectorPtr_Type u_repTime_;
    mutable MultiVectorPtr_Type p_repTime_;
    /*####################*/

};
}
#endif
