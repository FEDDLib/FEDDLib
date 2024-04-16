#ifndef NAVIERSTOKESASSFE_decl_hpp
#define NAVIERSTOKESASSFE_decl_hpp
#include "feddlib/problems/abstract/NonLinearProblem.hpp"
#include "Xpetra_ThyraUtils.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include <Thyra_ProductVectorBase.hpp>
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_ModelEvaluatorBase_decl.hpp>
/*!
 Declaration of Navier-Stokes

 @brief Navier-Stokes
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD{

    /*!
     Declaration of Navier-Stokes

     @brief Navier-Stokes
     @author Christian Hochmuth
     @version 1.0
     @copyright CH
     */

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class NavierStokesAssFE : public NonLinearProblem<SC,LO,GO,NO>  {

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
    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef NonLinearProblem<SC,LO,GO,NO> NonLinearProblem_Type;
    typedef typename NonLinearProblem_Type::BlockMultiVectorPtrArray_Type BlockMultiVectorPtrArray_Type;

    typedef typename NonLinearProblem_Type::TpetraMatrix_Type TpetraMatrix_Type;

    typedef typename NonLinearProblem_Type::ThyraVecSpace_Type ThyraVecSpace_Type;
    typedef typename NonLinearProblem_Type::ThyraVec_Type ThyraVec_Type;
    typedef typename NonLinearProblem_Type::ThyraOp_Type ThyraOp_Type;
    typedef Thyra::BlockedLinearOpBase<SC> ThyraBlockOp_Type;

    typedef typename NonLinearProblem_Type::TpetraOp_Type TpetraOp_Type;
    //@}

    //! @name Constructor/Destructor
    //@{
    NavierStokesAssFE( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList );
    //@}

    virtual void info();

    virtual void assemble( std::string type = "" ) const;
    
    void assembleConstantMatrices() const;
    
    void assembleDivAndStab() const;
    
    void reAssemble( std::string type ) const;
    
    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const{};
    
    //void reAssembleFSI(std::string type, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const;
    
    virtual void reAssemble(MatrixPtr_Type& massmatrix, std::string type ) const;

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions);

    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const; //standard or reverse
    
    void calculateNonLinResidualVecWithMeshVelo(std::string type, double time, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const;
//    virtual int ComputeDragLift(vec_dbl_ptr_Type &values);

    virtual void getValuesOfInterest( vec_dbl_Type& values ){};
    
    virtual void computeValuesOfInterestAndExport() {};

//    virtual void assembleExternal( std::string type ){};
    /*####################*/

    void computeSteadyPostprocessingViscosity_Solution(); // Compute the viscosity based on the current velocity solution and save it inside viscosity_element_


    mutable MatrixPtr_Type 	A_;
    vec_int_ptr_Type pressureIDsLoc;
    MultiVectorPtr_Type u_rep_;
    MultiVectorPtr_Type p_rep_;

    MultiVectorPtr_Type viscosity_element_; //In case of a generalized-Newtonian fluid one can compute viscosity for visualization

private:
    mutable bool stokesTekoPrecUsed_; //Help variable to signal that we constructed the initial preconditioner for NOX with the Stokes system and we do not need to compute it if fill_W_prec is called for the first time. However, the preconditioner is only correct if a Stokes system is solved in the first nonlinear iteration. This only affects the block preconditioners of Teko
    /*####################*/

public:

    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() const;
    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Monolithic() const;
#ifdef FEDD_HAVE_TEKO
    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Block() const;
#endif
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() const;

private:

    virtual void evalModelImpl(
                       const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                       const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                       ) const;

    void evalModelImplMonolithic(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                 const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;

#ifdef FEDD_HAVE_TEKO
    void evalModelImplBlock(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                            const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;
#endif

};
}
#endif
