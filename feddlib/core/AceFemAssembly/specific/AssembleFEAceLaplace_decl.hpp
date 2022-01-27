#ifndef ASSEMBLEFEACELAPLACE_DECL_hpp
#define ASSEMBLEFEACELAPLACE_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFEAceLaplace : public AssembleFE<SC,LO,GO,NO> {
  public:
   
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

	typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

	typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;

	AssembleFEAceLaplace(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type parameters);

	virtual void assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix);
	virtual void assemblyRHS(MultiVectorPtr_Type &elementVector);

	virtual void assemblyJacobian(MatrixPtr_Type &A) {};
	virtual void assemblyMass(MatrixPtr_Type &A) {};

	
 };

}
#endif

