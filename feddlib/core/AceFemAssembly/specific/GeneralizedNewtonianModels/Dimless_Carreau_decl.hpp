#ifndef DIMLESS_CARREAU_DECL_hpp
#define DIMLESS_CARREAU_DECL_hpp

#include "feddlib/core/General/DifferentiableFuncClass.hpp"
// #include "feddlib/core/AceFemAssembly/Helper.hpp"
#include "feddlib/core/FEDDCore.hpp"

namespace FEDD
{

   /*!
   \class Dimless_Carreau
   \brief This class is derived from the abstract class DifferentiableFuncClass and should provide functionality to evaluate the viscosity function specified by Carreau model in dimensionless form that why we need to multiply with a reference viscosity to obtain the actual viscosity value(see [1])
     \tparam SC The scalar type. So far, this is always double, but having it as a template parameter would allow flexibily, e.g., for using complex instead
     \tparam LO The local ordinal type. The is the index type for local indices
     \tparam GO The global ordinal type. The is the index type for global indices
     \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

 In general, there are several equations used to describe the material behavior of blood.
 Here (in this folder), we implement generalized Newtonian constitutive equations that capture the shear-thinning behavior of blood.
 The chosen shear-thinning model, such as the Power-Law or Carreau-Yasuda model, provides a function to update viscosity,
 which is no longer constant but depends on the shear rate (and other fixed constant parameters).

 In our FEM code, we need an update function for viscosity, depending on the chosen model.
 If we apply Newton's method or NOX, we also require the directional derivative of the viscosity function,
 which also depends on the chosen model. Additionally, we need functions to set the required parameters
 of the chosen model and a function to print the parameter values.

 The material parameters can be provided through a Teuchos::ParameterList object,
 which contains all the parameters specified in the input file 'parametersProblem.xml'.
 The structure of the input file and, consequently,
 the resulting parameter list can be chosen freely. The FEDDLib will handle reading the parameters from the
 file and making them available.


IMPORTANT: Regard e.g. following paper for computation of derivative
@article{article_he,
author = {He, Xin and Neytcheva, Maya and Vuik, C.},
year = {2015},
month = {06},
pages = {33-58},
title = {On Preconditioning of Incompressible Non-Newtonian Flow Problems},
volume = {33},
journal = {Journal of Computational Mathematics},
doi = {10.4208/jcm.1407-m4486}
}

In our Finite Element Method (FEM) assembly, we require the directional derivative of the shear stress term with respect to the velocity vector for Newton's method. Since our viscosity function depends nonlinearly on the shear rate, which in turn depends on velocity,
we must account for this contribution. Thus, when assembling the directional derivative term (in AssembleFEGeneralizedNewtonian), we call the function "this->viscosityModel->computeDerivative."
It's crucial to note that through the chain rule, we obtain different contributions when computing the directional derivative of eta(gamma_dot) with respect to velocity.
The individual terms are, for instance, detailed in the literature paper above.  In general one obtain:

d (eta( gamma_Dot(v + eps*delta_v) ) )/ d(eps) |eps=0  =  d(eta)/ d(gamma_Dot) * d(gamma_Dot)*d(Pi_||) * d(Pi_||)/ d(D) : d (D(v+eps*delta_v))/d(eps) |eps=0
Here Pi_|| is the second invariant of the strain-rate tensor D(v).
In total we get four contributions:
IMPORTANT: Consider that the contraction operator : is used in the formula above
 1. d(eta)/ d(gamma_Dot)               - Compute the derivative of eta with respect to gamma_Dot -- For e.g. Carreau-Yasuda model just compute it analytically
 2. d(gamma_Dot)*d(Pi_||)              - The shear rate can be defined in terms of the second invariant of D by means of gamma_Dot = 2 \sqrt( -  Pi_||) -- So the derivative is just d(gamma_Dot)*d(Pi_||) = - 1/ sqrt( -  Pi_||) = -2/ gamma_Dot
 3. d(Pi_||)/ d(D)                     - Derivative of second invariant of D with respect to D -- For incompressible fluids, the first invariant is zero, and the result can be found in the literature. However, as it depends on velocity, this contribution will be implemented in the assemblyFESpecific
 4. d (D(v+eps*delta_v))/d(eps) |eps=0 - Directional derivative of D with respect to v - also implemented in assemblyFESpecific


For our viscosity models we define our derivative function "computeDerivative" to return the result of d(eta)/ d(gamma_Dot) * d(gamma_Dot)*d(Pi_||), which encompasses the first two contributions dependent on the chosen viscosity model.
*/

   template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
   class Dimless_Carreau : public DifferentiableFuncClass<SC, LO, GO, NO>
   {
   public:
      typedef MultiVector<SC, LO, GO, NO> MultiVector_Type;
      typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
      typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

      typedef DifferentiableFuncClass<SC, LO, GO, NO> DifferentiableFuncClass_Type;

      // Inherited Function from base abstract class

      /*!
       \brief Each constitutive model includes different material parameters which will be specified in parametersProblem.xml
              This function should set the specififc needed parameters for each model to the defined values.
      @param[in] ParameterList as read from the xml file (maybe redundant)
      */
      void setParams(ParameterListPtr_Type params) override;
      /*!
       \brief Update the viscosity according to a chosen shear thinning generalized newtonian constitutive equation. Viscosity depends on spatial coordinates due to its dependency on velocity gradients
       @param[in] params as read from the xml file (maybe redundant)
       @param[in] shearRate scalar value of computed shear rate
       @param[in,out] viscosity value of viscosity
      */
      void evaluateMapping(ParameterListPtr_Type params, double shearRate, double &viscosity) override;
      /*!
       \brief For Newton method and NOX we need additional term in Jacobian considering directional derivative of our functional formulation.
               IMPORTANT: Here we implement the contribution of this: d(eta)/ d(gamma_Dot) * d(gamma_Dot)*d(Pi_||).
              For each constitutive model the function looks different and will be defined inside this function
       @param[in] params as read from the xml file (maybe redundant)
       @param[in] shearRate scalar value of computed shear rate
       @param[in,out] res scalar value of   d(eta)/ d(gamma_Dot) * d(gamma_Dot)*d(Pi_||)
      */
      void evaluateDerivative(ParameterListPtr_Type params, double shearRate, double &res) override;

      /*!
      \brief Print parameter values used in model at runtime
      */
      void echoInformationMapping() override;

      // New Added Functions
      /*!
       \brief Get the current viscosity value
       \return the scalar value of viscosity
       */
      double getViscosity() { return viscosity_; };

      /*!

      \brief Constructor for Dimless_Carreau
      @param[in] parameters Parameterlist for current problem	*/
      Dimless_Carreau(ParameterListPtr_Type parameters);

   private:
      double viscosity_;
      std::string shearThinningModel_; // for printing out which model is actually used
      //!
      double characteristicTime; // corresponds to \lambda in the formulas in the literature here dimensionless
      double fluid_index_n;      // corresponds to n in the formulas being the power-law index
      double nu_0;               // is the zero shear-rate viscosity here dimensionless
      double nu_infty;           // is the infnite shear-rate viscosity here dimensionless
      double shear_rate_limitZero;
      double reference_viscosity; // to obtain actual viscosity we have to multiply with reference viscosity
   };

}
#endif
