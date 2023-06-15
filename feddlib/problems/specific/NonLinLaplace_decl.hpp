#ifndef NonLinLaplace_decl_hpp
#define NonLinLaplace_decl_hpp
#include "feddlib/problems/abstract/NonLinearProblem.hpp"
#include <Thyra_ModelEvaluatorBase_decl.hpp>
#include <Thyra_PreconditionerBase.hpp>
/*!
 Declaration of NonLinLaplace

 @brief NonLinLaplace
 @author Kyrill Ho
 @version 1.0
 @copyright KH
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go,
          class NO = default_no>
class NonLinLaplace : public NonLinearProblem<SC, LO, GO, NO> {

public:
  //! @name Public Types
  //@{
  typedef Problem<SC, LO, GO, NO> Problem_Type;
  typedef typename Problem_Type::Matrix_Type Matrix_Type;
  typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;

  typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;

  typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
  typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;

  typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
  typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
  typedef
      typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;

  typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
  typedef
      typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

  typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
  typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

  typedef NonLinearProblem<SC, LO, GO, NO> NonLinearProblem_Type;
  typedef typename NonLinearProblem_Type::BlockMultiVectorPtrArray_Type
      BlockMultiVectorPtrArray_Type;

  typedef typename NonLinearProblem_Type::TpetraMatrix_Type TpetraMatrix_Type;

  typedef typename NonLinearProblem_Type::ThyraVecSpace_Type ThyraVecSpace_Type;
  typedef typename NonLinearProblem_Type::ThyraVec_Type ThyraVec_Type;
  typedef typename NonLinearProblem_Type::ThyraOp_Type ThyraOp_Type;

  typedef typename NonLinearProblem_Type::TpetraOp_Type TpetraOp_Type;

  typedef Tpetra::CrsMatrix<SC, LO, GO, NO> tpetra_matrix;

  //@}

  //! @name Constructor/Destructor
  //@{
  NonLinLaplace(const DomainConstPtr_Type &domain, std::string FEType,
                ParameterListPtr_Type parameterList);
  //@}
  ~NonLinLaplace();

  virtual void info();

  virtual void assemble(std::string type = "") const;

  void initAssemble() const;

  void reAssemble(std::string type) const;

  virtual void reAssemble(BlockMultiVectorPtr_Type previousSolution) const {};

  virtual void reAssemble(MatrixPtr_Type &massmatrix,
                          std::string type) const {};

  virtual void
  reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions);

  virtual void calculateNonLinResidualVec(std::string type,
                                          double time = 0.) const;

  virtual void getValuesOfInterest(vec_dbl_Type &values){};

  virtual void computeValuesOfInterestAndExport(){};

  //    virtual void assembleExternal( std::string type ){};

  Teuchos::RCP<Thyra::LinearOpBase<SC>> create_W_op() const;

  Teuchos::RCP<Thyra::PreconditionerBase<SC>> create_W_prec() const;

private:
  virtual void
  evalModelImpl(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;
  mutable MultiVectorPtr_Type u_rep_;
};
} // namespace FEDD
#endif
