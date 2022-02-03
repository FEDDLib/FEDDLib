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

	/*!
	 \brief Assemble the element Jacobian matrix.
	 \return the element Jacobian matrix
	*/
	virtual SmallMatrixPtr_Type assembleJacobian();

	/*!
	 \brief Assemble the element right hand side vector.
	 \return the element right hand side vector
	*/
	virtual vec_dbl_Type assembleRHS();	


   private:
	void assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix);

	void gradPhi(	int Dimension,
                    int intFE,
                    int i,
                    vec_dbl_Type &QuadPts,
                    vec_dbl_ptr_Type &value);
    
    /*! Most of the quadrature formulas can be found in http://code-aster.org/doc/v11/en/man_r/r3/r3.01.01.pdf 01/2021  */
    void getQuadratureValues(int Dimension,
                            int Degree,
                            vec2D_dbl_ptr_Type &QuadPts,
                            vec_dbl_ptr_Type &QuadW,
                            std::string FEType);
    
    int getDPhi(	vec3D_dbl_ptr_Type &DPhi,
                	vec_dbl_ptr_Type &weightsDPhi,
                    int Dimension,
                    std::string FEType,
                    int Degree);

	void buildTransformation(SmallMatrix<SC>& B);

    UN determineDegree(UN dim,
                       std::string FEType,
                       UN degFunc);

    int getPhi(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree,
               			    std::string FETypeQuadPoints="");

	void phi(int dim,
			  int intFE,
			  int i,
			  vec_dbl_Type &p,
			  double* value);

    void applyBTinv(vec3D_dbl_ptr_Type& dPhiIn,
                    vec3D_dbl_Type& dPhiOut,
                    SmallMatrix<SC>& Binv);
	
 };

}
#endif

