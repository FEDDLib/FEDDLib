#ifndef DIMLESS_CARREAU_DEF_hpp
#define DIMLESS_CARREAU_DEF_hpp

#include "Dimless_Carreau_decl.hpp"

namespace FEDD {

template <class SC, class LO, class GO, class NO>
Dimless_Carreau<SC,LO,GO,NO>::Dimless_Carreau(ParameterListPtr_Type params):
DifferentiableFuncClass<SC,LO,GO,NO>(params)
{
    this->params_=params;
    // Reading through parameterlist
    shearThinningModel_ = this->params_->sublist("Material").get("ShearThinningModel","");
	characteristicTime  = this->params_->sublist("Material").get("Dimless_CharacteristicTime_Lambda",0.);            // corresponds to \lambda in the formulas in the literature
    fluid_index_n       = this->params_->sublist("Material").get("Carreau_FluidIndex_n",1.);                         // corresponds to n in the formulas being the power-law index
    nu_0                = this->params_->sublist("Material").get("Dimless_ZeroShearRateViscosity_eta0",0.);          // is the zero shear-rate viscosity
    nu_infty            = this->params_->sublist("Material").get("Dimless_InftyShearRateViscosity_etaInfty",0.);     // is the infnite shear-rate viscosity
    reference_viscosity = this->params_->sublist("Material").get("Reference_Viscosity",0.); 
    shear_rate_limitZero= this->params_->sublist("Material").get("Numerical_ZeroValue_ShearRate",1e-8);  
    
    viscosity_ = 0.;
    //	TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "No discretisation Information for Velocity in Navier Stokes Element." );

	
}

template <class SC, class LO, class GO, class NO>
void Dimless_Carreau<SC,LO,GO,NO>::evaluateMapping(ParameterListPtr_Type params, double shearRate, double &viscosity) {
	
    viscosity = this->reference_viscosity*( this->nu_infty +(this->nu_0-this->nu_infty)*(pow(1.0+pow(this->characteristicTime*shearRate,2.0)    , (this->fluid_index_n-1.0)/2.0 )));
    this->viscosity_ = viscosity; // Multiply with reference viscosity to obtain actual viscosity with dimension of reference viscosity [Pa s]
}


                
template <class SC, class LO, class GO, class NO>
void Dimless_Carreau<SC,LO,GO,NO>::evaluateDerivative(ParameterListPtr_Type params, double shearRate, double &res) {
	
// The function is composed of d_eta/ d_GammaDot * d_GammaDot/ D_Tau with d_GammaDot * d_GammaDot/ D_Tau= - 2/GammaDot
// So for a Carreau-like Fluid we do not get the problem that the shear rate is in denominator
res = (-2.0)*(this->nu_0-this->nu_infty)*(this->fluid_index_n-1.0)*pow(this->characteristicTime, 2)*pow(1.0+pow(this->characteristicTime*shearRate,2.0)    , ((this->fluid_index_n-3.0)/2.0) );
res = this->reference_viscosity*res;

}


template <class SC, class LO, class GO, class NO>
void Dimless_Carreau<SC,LO,GO,NO>::setParams(ParameterListPtr_Type params){
    this->params_=params;
    // Reading through parameterlist
    shearThinningModel_ = this->params_->sublist("Material").get("ShearThinningModel","");
	characteristicTime  = this->params_->sublist("Material").get("Dimless_CharacteristicTime_Lambda",0.);            // corresponds to \lambda in the formulas in the literature
    fluid_index_n       = this->params_->sublist("Material").get("Carreau_FluidIndex_n",1.);                         // corresponds to n in the formulas being the power-law index
    nu_0                = this->params_->sublist("Material").get("Dimless_ZeroShearRateViscosity_eta0",0.);          // is the zero shear-rate viscosity
    nu_infty            = this->params_->sublist("Material").get("Dimless_InftyShearRateViscosity_etaInfty",0.);     // is the infnite shear-rate viscosity    shear_rate_limitZero= this->params_->sublist("Material").get("Numerical_ZeroValue_ShearRate",1e-8);  
    reference_viscosity = this->params_->sublist("Material").get("Reference_Viscosity",0.); 

 }




template <class SC, class LO, class GO, class NO>
void Dimless_Carreau<SC,LO,GO,NO>::echoInformationMapping(){
            std::cout << "************************************************************ "  <<std::endl;
            std::cout << "-- Chosen material model ..." << this->shearThinningModel_ << " --- "  <<std::endl;
            std::cout << "-- Specified material parameters:" <<std::endl;
            std::cout << "-- eta_0:"     <<  this->nu_0 <<std::endl;
            std::cout << "-- eta_Infty:" <<  this->nu_infty << std::endl;
            std::cout << "-- Fluid index n:" << this->fluid_index_n << std::endl;
            std::cout << "-- Reference viscosity:" << this->reference_viscosity <<std::endl;
            std::cout << "-- Characteristic time lambda:" << this->characteristicTime << std::endl;
            std::cout << "************************************************************ "  <<std::endl;
  }




}
#endif

