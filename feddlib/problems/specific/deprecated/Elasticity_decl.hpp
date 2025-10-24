#ifndef ELASTICITY_decl_hpp
#define ELASTICITY_decl_hpp

#include "feddlib/problems/abstract/LinearProblem.hpp"
#include "feddlib/problems/abstract/NonLinearProblem.hpp"

namespace FEDD{
/*!
 Declaration of Elasticity
 
 @brief Elasticity
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Elasticity : public NonLinearProblem<SC,LO,GO,NO>  {  
// We want Elasticity to be a general nonlinear problem, however if the underlying elasticity problem is linear we still use the reassmble functions but they dont do anything in this case.
public:
    typedef LinearProblem<SC,LO,GO,NO> LinearProblem_Type;
    typedef NonLinearProblem<SC,LO,GO,NO> NonLinearProblem_Type;
    typedef Teuchos::RCP<LinearProblem_Type> LinearProblemPtr_Type;
    typedef Teuchos::RCP<NonLinearProblem_Type> NonLinearProblemPtr_Type;
    typedef LinElas<SC,LO,GO,NO> LinElas_Type;
    typedef NonLinElasticity<SC,LO,GO,NO> NonLinElas_Type;
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;
    
    Elasticity( const DomainConstPtr_Type  &domain, std::string FEType, ParameterListPtr_Type parameterList );

    ~Elasticity();

    virtual void info();
    
    virtual void assemble();

    virtual void reAssemble(std::string type="Newton") const;
    
    virtual void reAssembleExtrapolation( );
    
    virtual void calculateNonLinResidualVec(std::string type)const;
        
    virtual void getValuesOfInterest( vec_dbl_Type& values ){}
    
    virtual void computeValuesOfInterestAndExport() {}
    
//    virtual void assembleExternal( std::string type ){}

    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() const;
    
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() const;
    
private:
    
    virtual void evalModelImpl(
                               const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                               const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                               ) const;

    /*####################*/
    LinearProblemPtr_Type linearProblem_;
    NonLinearProblemPtr_Type nonLinearProblem_;
};
}
#endif
