#ifndef AssembleFE_SCI_SMC_Active_Growth_Reorientation_DECL_hpp
#define AssembleFE_SCI_SMC_Active_Growth_Reorientation_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/AssembleFEBlock.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#ifdef FEDD_HAVE_ACEGENINTERFACE
//#include <aceinterface.h>
//#include <ace2.h>
//#endif

/*!
\class AssembleFE_SCI_SMC_Active_Growth_Reorientation
	Coupled deformation diffusion problem with smooth-muscle model with active response, growth and reorientation
	Derived from AssembleFE base class
	Active response with MLCK and MLCP
*/


namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFE_SCI_SMC_Active_Growth_Reorientation : public AssembleFE<SC,LO,GO,NO> {
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
	    void assembleJacobian() override;

        /*!
	    \brief Assemble the element right hand side vector.
	    */
	    void assembleRHS() override;

 		/*!
	    \brief Assemble block parts of the element Jacobian matrix.
	    \return the element Jacobian matrix of block i 
	    */
		void assembleJacobianBlock(LO i) override{};

		void advanceInTime(double dt) override;

        void getMassMatrix(SmallMatrixPtr_Type &massMatrix ){massMatrix=massMatrix_;};

    protected:
        AssembleFE_SCI_SMC_Active_Growth_Reorientation(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple);

    private:

        void assemble_SCI_SMC_Active_Growth_Reorientation(SmallMatrixPtr_Type &elementMatrix);


        friend class AssembleFEFactory<SC,LO,GO,NO>; // Must have for specfic classes
	    
		std::string FEType_ ; // FEType of Disk

		SmallMatrixPtr_Type massMatrix_;

	    int dofsSolid_ ; // Degrees of freedom per node
		int dofsChem_;
	    int numNodesSolid_ ; // Number of nodes of element
		int numNodesChem_ ; // Number of nodes of element


	    int dofsElement_; // "Dimension of return matrix"
		
		int dofOrdering_; // Order of DOFs: 
						  // dofOrdering = 1 -> 'u1 v1 w1 c1 u2 v2 w2 c2 ... un vn wn cn'
						  // dofOrdering = 2 -> 'u1 v1 w1 u2 v2 w2 ... un vn wn c1 c2 c3 ... cn'

		
		
		int iCode_; // Integration Code
		int historyLength_; // Length of history vector
		
		vec_dbl_Type historyUpdated_;
		vec_dbl_Type history_;

		vec_dbl_Type solutionC_n_;
		vec_dbl_Type solutionC_n1_; 

		vec2D_dbl_Type timeParametersVec_;
    	double numSegments_ ;
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
		$[Gamma]$6 -Gamma6_24"
		"$[Lambda]$P50 -LambdaP50_25"
		"kDotMin -KDotMin_26"
		"$[Zeta]$1 -Zeta1_27"
		"kDotMax -KDotMax_28"
		"$[Gamma]$4 -Gamma4_29"
		"$[Lambda]$BarDotPMin -LambdaBarDotPMin_30"
		"$[Lambda]$BarDotPMax -LambdaBarDotPMax_31"
		"$[Gamma]$5 -Gamma5_32"
		"$[Zeta]$2 -Zeta2_33"
		"$[CapitalDelta]$$[Lambda]$BarPMin -DeltaLambdaBarPMin_34"
		"p1 -P1_35"
		"p3 -P3_36"
		"c50 -C50_37"
		"d0 -D0_38"
		"m -M_39"
		"activeStartTime -ActiveStartTime_40"
		"k$[Eta]$Plus -kEtaPlus_41"
		"m$[Eta]$Plus -mEtaPlus_42"
		"growthStartTime -growthStartTime_43"
		"reorientationStartTime -reorientationStartTime_44"
		"growthEndTime -growthEndTime_45"
		"reorientationEndTime -reorientationEndTime_46"
		"k$[Theta]$Plus -KThetaPlus_47"
		"k$[Theta]$Minus -KThetaMinus_48"
		"m$[Theta]$Plus -MThetaPlus_49"
		"m$[Theta]$Minus -MThetaMinus_50"
		"$[Theta]$Plus1 -ThetaPlus1_51"
		"$[Theta]$Plus2 -ThetaPlus2_52"
		"$[Theta]$Plus3 -ThetaPlus3_53"
		"$[Theta]$Minus1 -ThetaMinus1_54"
		"$[Theta]$Minus2 -ThetaMinus2_55"
		"$[Theta]$Minus3 -ThetaMinus3_56"
		"$[Rho]$ -Density_57"	
	*/

		double fA_;
		double lambdaC50_;
		double gamma3_;
		double lambdaBarCDotMax_;
		double lambdaBarCDotMin_;
		double gamma2_;
		double gamma1_;
		double eta_;
		double ca50_;
		double k2_;
		double k5_;
		double k3_;
		double k4_;
		double k7_;
		double kappa_;
		double beta1_;
		double muA_;
		double beta2_;
		double alpha2_;
		double alpha3_;
		double alpha1_;
		double alpha4_;
		double alpha5_;
		double gamma6_ ;
		double lambdaP50_;
		double kDotMin_ ;
	    double zeta1_ ;
		double kDotMax_ ;
		double gamma4_ ;
		double lambdaBarDotPMin_ ;
		double lambdaBarDotPMax_ ;
		double gamma5_;
		double zeta2_;
		double DeltaLambdaBarPMin_ ;
		double p1_;
		double p3_;
		double c50_;
		double d0_;
		double m_;
		double activeStartTime_;
		double kEtaPlus_ ;
		double mEtaPlus_ ;
		double growthStartTime_ ;
		double reorientationStartTime_ ;
		double growthEndTime_ ;
		double reorientationEndTime_ ;
		double kThetaPlus_ ;
		double kThetaMinus_ ;
		double mThetaPlus_ ;
		double mThetaMinus_ ;
		double thetaPlus1_ ;
		double thetaPlus2_ ;
		double thetaPlus3_ ;
		double thetaMinus1_ ;
		double thetaMinus2_ ;
		double thetaMinus3_ ;
		double kMin_;
		double rho_;
};

}
#endif //AssembleFE_SCI_SMC_Active_Growth_Reorientation_DECL_hpp
