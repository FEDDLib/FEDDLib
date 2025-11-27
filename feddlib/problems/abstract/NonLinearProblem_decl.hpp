#ifndef NONLINEARPROBLEM_DECL_hpp
#define NONLINEARPROBLEM_DECL_hpp

#include <Thyra_StateFuncModelEvaluatorBase.hpp>
#include <Thyra_ProductVectorBase.hpp>
#include <Teko_Utilities.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/problems/abstract/Problem.hpp"

/*!
 Declaration of NonLinearProblem

 @brief  NonLinearProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD{
template<class LO_, class GO_, class NO_>
class BlockMap;

template <class SC = default_sc,
          class LO = default_lo,
          class GO = default_go,
          class NO = default_no>
class NonLinearProblem : public Problem<SC,LO,GO,NO> , public Thyra::StateFuncModelEvaluatorBase<SC> {

public:

    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;
    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;
    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;
    

    typedef Tpetra::Map<LO,GO,NO> TpetraMap_Type;
    typedef Teuchos::RCP<TpetraMap_Type> TpetraMapPtr_Type;
    typedef Teuchos::RCP<const TpetraMap_Type> TpetraMapConstPtr_Type;
    typedef const TpetraMapConstPtr_Type TpetraMapConstPtrConst_Type;


    typedef BlockMap<LO,GO,NO> BlockMap_Type;
    typedef Teuchos::RCP<BlockMap_Type> BlockMapPtr_Type;
    typedef Teuchos::RCP<const BlockMap_Type> BlockMapConstPtr_Type;
    typedef Teuchos::Array<BlockMultiVectorPtr_Type> BlockMultiVectorPtrArray_Type;

    using ThyraTypes = ThyraTypedefs<SC>;
    using TpetraTypes = TpetraTypedefs<SC,LO,GO,NO>;
    using ThyraVecSpace_Type = typename ThyraTypes::ThyraVecSpace_Type;
    typedef Teuchos::RCP<const ThyraVecSpace_Type> ThyraVecSpaceConstPtr_Type;
    using ThyraVec_Type = typename ThyraTypes::ThyraVec_Type;
    using TpetraMatrix_Type = typename TpetraTypes::TpetraMatrix_Type;
    using ThyraOp_Type = typename ThyraTypes::ThyraOp_Type;
    using TpetraOp_Type = typename TpetraTypes::TpetraOp_Type;
    typedef Thyra::BlockedLinearOpBase<SC> ThyraBlockOp_Type;

    /// @brief Constructor
    /// @param comm 
    NonLinearProblem(CommConstPtr_Type comm);

    /// @brief Constructor with parameterlist
    /// @param parameterList 
    /// @param comm 
    NonLinearProblem(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm);

    ~NonLinearProblem();

    virtual void info() = 0;

    /// @brief Information about the non-linear problem
    void infoNonlinProblem();

    /// @brief Initialisation of the non-linear problem with system, vectors, and Thyra vector spaces for NOX
    /// @param nmbVectors 
    void initializeProblem(int nmbVectors=1);
    
    /// @brief assemble of type exectuted by the derived specific non-linear problem classes
    /// @param type for example Newton
    virtual void assemble( std::string type = "" ) const = 0;

    /// @brief  Virtual class to extract values of interest that are computed during the solve
    /// @param values 
    virtual void getValuesOfInterest( vec_dbl_Type& values ) = 0;
    
    /// @brief Solving the non-linear problem and updating the solution
    /// @param criterion Update or Residual
    /// @param criterionValue The actual value of the criterion
    /// @return number of linear iterations
    int solveAndUpdate( const std::string& criterion , double& criterionValue );

    /// @brief This is where the linear solve specifically happens
    /// @return Number of linear iterations
    int solveUpdate( );


    /// @brief Reassemble with previous solution. I think it is not used anymore. @TODO: Look into this.
    /// @param previousSolution 
    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const = 0;
    
    void reAssembleAndFill( BlockMatrixPtr_Type bMat, std::string type="FixedPoint" );

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions) = 0;

    /// @brief Initialisation of the non-linear vectors like residual and previous solution
    /// @param nmbVectors 
    void initializeVectorsNonLinear(int nmbVectors=1);

    /// @brief Calculate the 2-norm of the residual vector
    /// @return Value of the norm of the residual
    double calculateResidualNorm() const;

    /// @brief Virtual function which is implemented in the specific non-linear problem classes to calculate the non-linear residual vector
    /// @param type standard or reverse depending on Newton formulation, e.g. as in NOX or FEDDLib-Newton
    /// @param time current timestep
    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const = 0; //type=standard or reverse

    /// @brief Calculate the non-linear residual vector with given coefficients for time-dependent problems (if used for timeproblem)
    virtual void calculateNonLinResidualVec(SmallMatrix<double>& coeff, std::string type="standard", double time=0., BlockMatrixPtr_Type systemMass = Teuchos::null); //type=standard or reverse
    
    /// @brief Get the residual vector
    /// @return residual vector
    BlockMultiVectorPtr_Type getResidualVector() const;
    
    /// @brief Get previous solution. Needed for time-dependent problems
    /// @return 
    BlockMultiVectorPtr_Type getPreviousSolution() const{ return previousSolution_; }

    virtual Thyra::ModelEvaluatorBase::InArgs<SC> getNominalValues() const;

    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > get_x_space() const;

    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > get_f_space() const;

    virtual ::Thyra::ModelEvaluatorBase::InArgs<SC> createInArgs() const;

    void initNOXParameters( );

    void initVectorSpaces( );

    void initVectorSpacesMonolithic( );

    void initVectorSpacesBlock( );

    virtual ::Thyra::ModelEvaluatorBase::OutArgs<SC> createOutArgsImpl() const;

    
    void setNonlinearIterationStep(int newtonStep)  { this->newtonStep_ = newtonStep;}   // For SolveFixedPoint etc. we need to set the Newton Step manually
    int  getNonlinearIterationStep() const override { return newtonStep_; }              // This overrides the default in Problem to provide the actual Newton step

    
    double nonLinearTolerance_;
    BlockMultiVectorPtr_Type    previousSolution_;
    mutable BlockMultiVectorPtr_Type    residualVec_;
    SmallMatrix<double> coeff_;// coefficients for a time-dependent problem

    mutable int newtonStep_;

    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() const;
    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Monolithic() const;
#ifdef FEDD_HAVE_TEKO
    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Block() const;
#endif
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() const;


private:
    mutable bool precInitOnly_; //Help variable to signal that we constructed the initial preconditioner 
                                // for NOX already and we do not need to compute it if fill_W_prec is 
                                // called for the first time. However, the preconditioner is only correct 
                                // if a linear system is solved in the first nonlinear iteration. 


    Thyra::ModelEvaluatorBase::InArgs<SC> nominalValues_;

    Thyra::ModelEvaluatorBase::InArgs<SC> prototypeInArgs_;
    Thyra::ModelEvaluatorBase::OutArgs<SC> prototypeOutArgs_;

    Teuchos::RCP<const ThyraVecSpace_Type> xSpace_;
    Teuchos::RCP<const ThyraVecSpace_Type> fSpace_;

    Teuchos::RCP<ThyraVec_Type> x0_;

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
