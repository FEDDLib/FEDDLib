#ifndef  AssembleFE_NonLinElas2_DECL_hpp
#define  AssembleFE_NonLinElas2_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class  AssembleFE_NonLinElas2 : public AssembleFE<SC,LO,GO,NO> {
  public:
   
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

	typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

	typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;


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
	/*!
		\brief Update the parameter read from the ParameterList.
		@param[in] Parameter as read from the xml file
	*/
    void updateParameter(std::string type, double value) override;

   protected:
	 AssembleFE_NonLinElas2(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type parameters,   tuple_disk_vec_ptr_Type tuple); 
   private:
	void assemblyNonLinElas(SmallMatrixPtr_Type &elementMatrix);

    friend class AssembleFEFactory<SC,LO,GO,NO>; // Must have for specfic classes

	
	double E_ ; 
   	double lambda_;
	double poissonRatio_;
	std::string FEType_ ; // FEType of Disk

	int dofs_ ; // Degrees of freedom per node

	int numNodes_ ; // Number of nodes of element

	int dofsElement_; // "Dimension of return matrix"

 };

}
#endif

