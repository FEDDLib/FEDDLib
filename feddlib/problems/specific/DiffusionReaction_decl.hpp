#ifndef DIFFUSION_decl_hpp
#define DIFFUSION_decl_hpp

#include "feddlib/problems/abstract/Problem.hpp"
#include "feddlib/problems/abstract/NonLinearProblem.hpp"
#include <Thyra_ProductVectorBase.hpp>
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_ModelEvaluatorBase_decl.hpp>
#include "feddlib/core/FEDDCore.hpp"

/*!
 Declaration of DiffusionReaction
 
 @brief Diffusion
 @author Lea Sassmannshausen
 @version 1.0
 @copyright LS
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class DiffusionReaction : public Problem<SC,LO,GO,NO> {
    
public:
    
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    
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
    typedef Thyra::BlockedLinearOpBase<SC> ThyraBlockOp_Type;
    
    typedef typename NonLinearProblem_Type::TpetraOp_Type TpetraOp_Type;

   // typedef boost::function<void(double* x,double* res, const double* parameters)> Reac_func_Type;

    DiffusionReaction (const DomainConstPtr_Type &domain, std::string FEType, ParameterListPtr_Type parameterList, vec2D_dbl_Type diffusionTensor,  RhsFunc_Type reactionFunc, bool vectorDiffusion=false);
    
    ~DiffusionReaction();
    
    virtual void info();
    
    void assembleConstantMatrices( std::string type = "" ) const;

    virtual void assemble( std::string type = "" ) const;

    virtual void getValuesOfInterest( vec_dbl_Type& values ){}

	MatrixPtr_Type getMassMatrix() const; // new for calculating L2-Error
    
    virtual void computeValuesOfInterestAndExport() {}

    /*void reAssemble( std::string type = "" ) const;

    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const{}
    
    //    virtual int SetupPreconditioner(BMat_ptr_Type systemPrec, ThyraConstLinOpPtr_Type thyraMatrix=Teuchos::null, ThyraPrecPtr_Type thyraPreconditioner = Teuchos::null, LinSolverBuilderPtr_Type linearSolverBuilder = Teuchos::null) const;

	 virtual void reAssemble(MatrixPtr_Type& massmatrix, std::string type ) const;

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions);

    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const; //standard or reverse*/


    mutable MatrixPtr_Type 	A_;
    MultiVectorPtr_Type u_rep_;
	vec_dbl_Type funcParameter_;


    //Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() const;
    
    //Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() const;

private:
    /*####################*/
    bool vectorDiffusion_;
	vec2D_dbl_Type diffusionTensor_;
	RhsFunc_Type reactionFunc_;

	/*virtual void evalModelImpl(
                       const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                       const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                       ) const;*/

    
};
}
#endif
