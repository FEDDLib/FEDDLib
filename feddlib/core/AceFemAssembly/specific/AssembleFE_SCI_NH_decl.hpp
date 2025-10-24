#ifndef ASSEMBLEFE_SCI_NH_DECL_hpp
#define ASSEMBLEFE_SCI_NH_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFE_SCI_NH : public AssembleFE<SC,LO,GO,NO> {
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
        void assembleJacobianBlock(LO i) override {}

		void advanceInTime( double dt) override;

    protected:
        AssembleFE_SCI_NH(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple);

    private:

        void assembleDeformationDiffusionNeoHook(SmallMatrixPtr_Type &elementMatrix);

        friend class AssembleFEFactory<SC,LO,GO,NO>; // Must have for specfic classes

        double E0_;
		double E1_;
		double poissonRatio_;
		double c1_;
		double D0_;
		double m_;
	    
	    std::string FEType_ ; // FEType of Disk

	    int dofsSolid_ ; // Degrees of freedom per node
		int dofsChem_;
	    int numNodesSolid_ ; // Number of nodes of element
		int numNodesChem_ ; // Number of nodes of element


	    int dofsElement_; // "Dimension of return matrix"
		
		int dofOrdering_; // Order of DOFs: 
						  // dofOrdering = 1 -> 'u1 v1 w1 c1 u2 v2 w2 c2 ... un vn wn cn'
						  // dofOrdering = 2 -> 'u1 v1 w1 u2 v2 w2 ... un vn wn c1 c2 c3 ... cn'
		vec_dbl_Type solution_n_;
		vec_dbl_Type solution_n1_; 	

		vec2D_dbl_Type timeParametersVec_;
    	double numSegments_ ;			  
};

}
#endif //ASSEMBLEFE_SCI_NH_DECL_hpp
