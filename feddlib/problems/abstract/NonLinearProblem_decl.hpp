#ifndef NONLINEARPROBLEM_DECL_hpp
#define NONLINEARPROBLEM_DECL_hpp

#include "Problem.hpp"
#include <Thyra_StateFuncModelEvaluatorBase.hpp>


/*!
 Declaration of NonLinearProblem

 @brief  NonLinearProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD{
template<class SC_, class LO_, class GO_, class NO_>
class Problem;
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

    typedef Thyra::VectorSpaceBase<SC> ThyraVecSpace_Type;
    typedef Teuchos::RCP<const ThyraVecSpace_Type> ThyraVecSpaceConstPtr_Type;
    typedef Thyra::VectorBase<SC> ThyraVec_Type;
    typedef Tpetra::CrsMatrix<SC, LO, GO, NO> TpetraMatrix_Type;
    typedef Thyra::LinearOpBase<SC> ThyraOp_Type;
    typedef Tpetra::Operator<SC,LO,GO,NO> TpetraOp_Type;
    typedef Thyra::BlockedLinearOpBase<SC> ThyraBlockOp_Type;

    NonLinearProblem(CommConstPtr_Type comm);

    NonLinearProblem(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm);

    ~NonLinearProblem();

    virtual void info() = 0;

    void infoNonlinProblem();

    void initializeProblem(int nmbVectors=1);
    
    virtual void assemble( std::string type = "" ) const = 0;

    virtual void getValuesOfInterest( vec_dbl_Type& values ) = 0;
    
    int solveAndUpdate( const std::string& criterion , double& criterionValue );

    int solveUpdate( );

//    virtual void reAssemble(std::string type="FixedPoint") const = 0;

    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const = 0;
    
    void reAssembleAndFill( BlockMatrixPtr_Type bMat, std::string type="FixedPoint" );

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions) = 0;
    //MultiVector_ptr_vec_ptr_Type allPreviousSolutions
//    virtual int ComputeDragLift(vec_dbl_ptr_Type &values) = 0;    

    void initializeVectorsNonLinear(int nmbVectors=1);

    double calculateResidualNorm() const;

    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const = 0; //type=standard or reverse

    // if used for timeproblem
    virtual void calculateNonLinResidualVec(SmallMatrix<double>& coeff, std::string type="standard", double time=0.); //type=standard or reverse
    
    BlockMultiVectorPtr_Type getResidualVector() const;
    
    BlockMultiVectorPtr_Type getPreviousSolution() const{ return previousSolution_; };

    virtual Thyra::ModelEvaluatorBase::InArgs<SC> getNominalValues() const;

    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > get_x_space() const;

    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > get_f_space() const;

    virtual ::Thyra::ModelEvaluatorBase::InArgs<SC> createInArgs() const;

    void initNOXParameters( );

    void initVectorSpaces( );

    void initVectorSpacesMonolithic( );

    void initVectorSpacesBlock( );

    // virtual void evalModelImpl(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                            //    const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const = 0;

    virtual ::Thyra::ModelEvaluatorBase::OutArgs<SC> createOutArgsImpl() const;
    
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

    ::Thyra::ModelEvaluatorBase::InArgs<SC> prototypeInArgs_;
    ::Thyra::ModelEvaluatorBase::OutArgs<SC> prototypeOutArgs_;

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
