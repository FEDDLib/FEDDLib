#ifndef ASSEMBLEFEACENAVIERSTOKES_DECL_hpp
#define ASSEMBLEFEACENAVIERSTOKES_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/Helper.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFEAceNavierStokes : public AssembleFE<SC,LO,GO,NO> {
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
	*/
	virtual void assembleJacobian();

	/*!
	 \brief Assemble the element right hand side vector.
	*/
	virtual void assembleRHS();	

	SmallMatrixPtr_Type getFixedPointMatrix(){return ANB_;};

	void setCoeff(SmallMatrix_Type coeff);
   protected:


   private:

	/*!

	 \brief Constructor for AssembleFEAceNavierStokes

	@param[in] flag Flag of element
	@param[in] nodesRefConfig Nodes of element in reference configuration
	@param[in] params Parameterlist for current problem
	@param[in] tuple vector of element information tuples. 
	*/
	AssembleFEAceNavierStokes(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type parameters,tuple_disk_vec_ptr_Type tuple); 

	/*!

	 \brief Assembly function for vector values laplacian \f$ \int_T \nabla v \cdot \nabla u ~dx\f$ 
	@param[in] &elementMatrix

	*/
	void assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix);

	/*!

	 \brief Assembly advection vector field \f$ \int_T \nabla v \cdot u(\nabla u) ~dx\f$ 
	@param[in] &elementMatrix

	*/
	void assemblyAdvection(SmallMatrixPtr_Type &elementMatrix);
	
	/*!
	 \brief Assembly advection vector field in u  
	@param[in] &elementMatrix

	*/
	void assemblyAdvectionInU(SmallMatrixPtr_Type &elementMatrix); 

	/*!

	 \brief Assembly \f$ \int_T  div(v) p ~dx\f$ / \f$ \int_T  div(u) q ~dx\f$
	@param[in] &elementMatrix

	*/
	void assemblyDivAndDivT(SmallMatrixPtr_Type &elementMatrix);



    friend class AssembleFEFactory<SC,LO,GO,NO>; // Must have for specfic classes

	void buildTransformation(SmallMatrix<SC>& B);

	void applyBTinv(vec3D_dbl_ptr_Type& dPhiIn,
		            vec3D_dbl_Type& dPhiOut,
		            SmallMatrix<SC>& Binv);

	//tuple_disk_vec_ptr_Type returnTuple(); /// @todo return tuple in case or check tuple

    /// @todo Why do we need dofs1_ and dofs2_ in the abstract class? I think, we should think about a general framework for this
	/// \todo Put into Parameterlist.
    int dofsVelocity_;
    int dofsPressure_;

	string FETypeVelocity_;
	string FETypePressure_;

	int numNodesVelocity_;
	int numNodesPressure_;

	int dofsElementVelocity_ ;
	int dofsElementPressure_ ;

	int dofsElement_;

	vec_dbl_Type solutionVelocity_;
	vec_dbl_Type solutionPressure_;

	SmallMatrixPtr_Type constantMatrix_;
	SmallMatrixPtr_Type ANB_;

	SmallMatrix_Type coeff_;

	double viscosity_ ;
   	 double density_ ;

	
 };

}
#endif

