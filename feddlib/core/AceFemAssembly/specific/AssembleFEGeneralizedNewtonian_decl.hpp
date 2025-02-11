#ifndef AssembleFEGeneralizedNewtonian_DECL_hpp
#define AssembleFEGeneralizedNewtonian_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFENavierStokes.hpp"
#include "feddlib/core/FE/Helper.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/General/DifferentiableFuncClass.hpp"
// Add the generalized Newtonian Fluid models for the viscosity
#include "feddlib/core/AceFemAssembly/specific/GeneralizedNewtonianModels/CarreauYasuda.hpp"
#include "feddlib/core/AceFemAssembly/specific/GeneralizedNewtonianModels/PowerLaw.hpp"
#include "feddlib/core/AceFemAssembly/specific/GeneralizedNewtonianModels/Dimless_Carreau.hpp"

namespace FEDD
{

	template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
	class AssembleFEGeneralizedNewtonian : public AssembleFENavierStokes<SC, LO, GO, NO>
	{
	public:
		typedef Matrix<SC, LO, GO, NO> Matrix_Type;
		typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

		typedef SmallMatrix<SC> SmallMatrix_Type;
		typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

		typedef MultiVector<SC, LO, GO, NO> MultiVector_Type;
		typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

		typedef AssembleFE<SC, LO, GO, NO> AssembleFE_Type;

		typedef DifferentiableFuncClass<SC, LO, GO, NO> DifferentiableFuncClass_Type;
		typedef Teuchos::RCP<DifferentiableFuncClass_Type> DifferentiableFuncClassPtr_Type;

		typedef InputToOutputMappingClass<SC, LO, GO, NO> InputToOutputMappingClass_Type;
		typedef Teuchos::RCP<InputToOutputMappingClass_Type> InputToOutputMappingClassPtr_Type;
		// smart pointer inside we need a type

		/*!
		 \brief Assemble the element Jacobian matrix.
		*/
		void assembleJacobian() override;

		/*!
		 \brief Assemble the element right hand side vector.
		*/
		void assembleRHS() override;

		/*!
		  \brief Assemble the element Jacobian matrix.
	 	  @param[in] block ID i
	    */
		void assembleJacobianBlock(LO i){TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No implementation");};

		/*!
	    \brief Compute the viscosity for an element depending on the knwon velocity solution.
		*/
		void computeLocalconstOutputField() override;

		/*
			\brief Assembly of FixedPoint- Matrix (System Matrix K with current u) 
	     */
	    void assembleFixedPoint();

	   SmallMatrixPtr_Type getFixedPointMatrix(){return this->ANB_;};

	protected:
		std::string shearThinningModel;
		int dofsElementViscosity_;

		/*!

		 \brief Constructor for AssembleFEAceNavierStokes

		@param[in] flag Flag of element
		@param[in] nodesRefConfig Nodes of element in reference configuration
		@param[in] params Parameterlist for current problem
		@param[in] tuple vector of element information tuples.
		*/
		AssembleFEGeneralizedNewtonian(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type parameters, tuple_disk_vec_ptr_Type tuple);

		/*!

		 \brief Assembly function for shear stress tensor which includes viscosity function \f$ \int_T \nabla v : (2\eta(\Dot{\gamma}(u))) D(u) ~dx\f$, which is a highly nonlinear Term
		@param[in] &elementMatrix
		*/
		void assemblyStress(SmallMatrixPtr_Type &elementMatrix);

		/*!

		 \brief Assembly function for directional derivative contribution of shear stress tensor term (see function for detailed formula)
		@param[in] &elementMatrix
		*/
		void assemblyStressDev(SmallMatrixPtr_Type &elementMatrix);
		/*!

		/*!

		 \brief Assembly function for neumann boundary term - If we want to have same outflow boundary condition as for the Newtonian case (Laplace operator) we have to subtract boundary integral term see
		 @article{Pacheco2021,
		abstract = {The matter of appropriate boundary conditions for open or truncated outflow regions in internal flow is still focus of discussion and research. In most practical applications, one can at best estimate mean pressure values or flow rates at such outlets. In the context of finite element methods, it is known that enforcing mean pressures through the pseudo-tractions arising from the Laplacian Navier-Stokes formulation yields accurate, physically consistent solutions. Nevertheless, when generalised Newtonian fluid models are considered, the resulting non-uniform viscosity fields render the classical Laplacian formulation inadequate. Thus, it is common practice to use the socalled stress-divergence formulation with natural boundary conditions known for causing nonphysical outflow behaviour. In order to overcome such a limitation, this work presents a novel mixed variational formulation that can be seen as a generalisation of the Laplacian Navier-Stokes form to fluids with shear-rate-dependent viscosity, as appearing in hemodynamic and polymeric flows. By appropriately manipulating the viscous terms in the variational formulation and employing a simple projection of the constitutive law, it is possible to devise a formulation with the desired natural boundary conditions and low computational complexity. Several numerical examples are presented to showcase the potential of our method, revealing improved accuracy and robustness in comparison with the state of the art.},
		author = {Pacheco, Douglas R Q and Steinbach, Olaf},
		doi = {10.51375/ijcvse.2021.1.6},
		file = {:home/user/.local/share/data/Mendeley Ltd./Mendeley Desktop/Downloaded/Pacheco, Steinbach - 2021 - On outflow boundary conditions in finite element simulations of non-Newtonian internal flows.pdf:pdf},
		journal = {International Journal of Computing and Visualization in Science and Engineering},
		mendeley-groups = {DoctoralThesis/Numerical Fluid Dynamics/NonNewtonian/Boundary_Conditions},
		number = {April},
		title = {{On outflow boundary conditions in finite element simulations of non-Newtonian internal flows}},
		year = {2021}
		}
		@param[in] &elementMatrix
		*/
		void assemblyOutflowNeumannBoundaryTerm(SmallMatrixPtr_Type &elementMatrix){TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Will be added when normal evaluation on the finite element is finished");} ;

		/*!

		/*!
		 \brief Assembly function for directional derivative contribution  neumann boundary term
		@param[in] &elementMatrix
		*/
		void assemblyOutflowNeumannBoundaryTermDev(SmallMatrixPtr_Type &elementMatrix){TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Will be added when normal evaluation on the finite element is finished");};

		/*!
		\brief Computation of shear rate using the current velocity solution at the nodes and the derivative of the shape function
		@param[in]  dPhiTrans             - Derivative of shape function
		@param[in]  dim 		          - dimension of the problem
		@param[in, out] gammaDot 		  - shear rate value depending on inputs

		*/
		void computeShearRate(vec3D_dbl_Type dPhiTrans, vec_dbl_ptr_Type &gammaDot, int dim);

		friend class AssembleFEFactory<SC, LO, GO, NO>; // Must have for specfic classes

		InputToOutputMappingClassPtr_Type viscosityModel; // viscosity Model can be in theory any Input to output mapping

	private:
	};

}
#endif
