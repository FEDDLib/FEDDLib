#ifndef ASSEMBLEFENONLINLAPLACE_DECL_hpp
#define ASSEMBLEFENONLINLAPLACE_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/FE/Helper.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go,
          class NO = default_no>
class AssembleFENonLinLaplace : public AssembleFE<SC, LO, GO, NO> {
  public:
    typedef Matrix<SC, LO, GO, NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

    typedef MultiVector<SC, LO, GO, NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

    typedef AssembleFE<SC, LO, GO, NO> AssembleFE_Type;

    /*!
     \brief Assemble the element Jacobian matrix.
     \return the element Jacobian matrix
    */
    void assembleJacobian() override;

    /*!
     \brief Assemble the element right hand side vector.
     \return the element right hand side vector
    */
    void assembleRHS() override;

    /*!
    \brief Assemble the element Jacobian matrix.
    @param[in] block ID i
    */
    void assembleJacobianBlock(LO i) override {};

  protected:
    AssembleFENonLinLaplace(
        int flag, vec2D_dbl_Type nodesRefConfig,
        ParameterListPtr_Type parameters,

        tuple_disk_vec_ptr_Type tuple);
  private:
    void assemblyNonLinLaplacian(SmallMatrixPtr_Type &elementMatrix);

    void buildTransformation(SmallMatrix<SC> &B);

    void applyBTinv(vec3D_dbl_ptr_Type &dPhiIn, vec3D_dbl_Type &dPhiOut,
                    SmallMatrix<SC> &Binv);

    friend class AssembleFEFactory<SC, LO, GO, NO>;
    string FEType_;
    int dofs_;
    int numNodes_;
    int dofsElement_;
};

} // namespace FEDD
#endif
