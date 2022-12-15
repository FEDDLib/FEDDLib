#ifndef AssembleFEAceDeformDiffu22_DECL_hpp
#define AssembleFEAceDeformDiffu22_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/Helper.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

//#ifdef FEDD_HAVE_ACEGENINTERFACE
//#include "aceinterface.h"
//#include "ace2.h"
//#endif

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFEAceDeformDiffu2 : public AssembleFE<SC,LO,GO,NO> {
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
	    virtual void assembleJacobian();

        /*!
	    \brief Assemble the element right hand side vector.
	    \return the element right hand side vector
	    */
	    virtual void assembleRHS();

    protected:
        AssembleFEAceDeformDiffu2(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple);

    private:

        void assembleDeformationDiffusionNeoHook(SmallMatrixPtr_Type &elementMatrix);

        friend class AssembleFEFactory<SC,LO,GO,NO>; // Must have for specfic classes
	    
	    string FEType_ ; // FEType of Disk

	    int dofsSolid_ ; // Degrees of freedom per node
		int dofsChem_;
	    int numNodesSolid_ ; // Number of nodes of element
		int numNodesChem_ ; // Number of nodes of element


	    int dofsElement_; // "Dimension of return matrix"
		
		int dofOrdering_; // Order of DOFs: 
						  // dofOrdering = 1 -> 'u1 v1 w1 c1 u2 v2 w2 c2 ... un vn wn cn'
						  // dofOrdering = 2 -> 'u1 v1 w1 u2 v2 w2 ... un vn wn c1 c2 c3 ... cn'

		/*
		fA -Fibre angle_1
		$[Lambda]$C50 -LambdaC50_2
		$[Gamma]$3 -Gamma3_3
		$[Lambda]$BarCDotMax -LambdaBarCDotMax_4
		$[Lambda]$BarCDotMin -LambdaBarCDotMin_5
		$[Gamma]$2 -Gamma2_6
		$[Gamma]$1 -Gamma1_7
		$[Eta]$1 -Eta1_8
		Ca50 -Ca50_9
		k2 -K2_10
		k5 -K5_11
		k3 -K3_12
		k4 -K4_13
		k7 -K7_14
		$[Kappa]$C -KappaC_15
		$[Beta]$1 -Beta1_16
		$[Mu]$a -MuA_17
		$[Alpha]$ -Alpha_18
		$[Epsilon]$1 -Epsilon1_19
		$[Epsilon]$2 -Epsilon2_20
		c1 -C1_21
		$[Alpha]$1 -Alpha1_22
		$[Alpha]$2 -Alpha2_23
		p1 -P1_24
		p3 -P3_25
		c50 -C50_26
		d0 -D0_27
		m -M_28
		startTime -StartTime_29 // wann das Material aktiv wird
		$[Rho]$0 -Density_30
		*/
		double fA_;
		double lambdaC50_;
		double gamma3_;
		double lambdaBarCDotMax_;
		double lambdaBarCDotMin_;
		double gamma2_;
		double gamma1_;
		double eta1_;
		double ca50_;
		double k2_;
		double k5_;
		double k3_;
		double k4_;
		double k7_;
		double kappaC_;
		double beta1_;
		double muA_;
		double alpha_;
		double epsilon1_;
		double epsilon2_;
		double c1_;
		double alpha1_;
		double alpha2_;
		double p1_;
		double p3_;
		double c50_;
		double d0_;
		double m_;
		double startTime_;
		double rho_;

		int iCode_; // Integration Code
		vec_dbl_Type historyUpdated_;

		vec_dbl_Type solutionC_n_;
		vec_dbl_Type solutionC_n1_; 
 


};

}
#endif //AssembleFEAceDeformDiffu2_DECL_hpp
