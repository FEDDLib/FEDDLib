#ifndef ASSEMBLEFEACELAPLACE_DECL_hpp
#define ASSEMBLEFEACELAPLACE_DECL_hpp

#include "AssembleFE_decl.hpp"

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFEAceLaplace : public AssembleFE<SC,LO,GO,NO> {
  public:
   
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

	typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;

	AssembleFEAceLaplace(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type parameters);

	virtual void assemblyLaplacian(MatrixPtr_Type &elementMatrix) const;
	virtual void assemblyRHS(VectorPtr_Type &elementVector) const;

   private:
 	virtual void assemblyLaplacian(MatrixPtr_Type &A) =0;

	virtual void assemblyAdvectionVecField(MatrixPtr_Type &A) =0;
	virtual void assemblyDivAndDivT(MatrixPtr_Type &A) =0;*/
 }

}


