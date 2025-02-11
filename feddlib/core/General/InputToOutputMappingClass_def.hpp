#ifndef InputToOutputMappingClass_DEF_hpp
#define InputToOutputMappingClass_DEF_hpp

#include "InputToOutputMappingClass_decl.hpp"

namespace FEDD {


template <class SC, class LO, class GO, class NO>
InputToOutputMappingClass<SC,LO,GO,NO>::InputToOutputMappingClass(ParameterListPtr_Type params)
{

	params_=params;

}

template <class SC, class LO, class GO, class NO>
void InputToOutputMappingClass<SC,LO,GO,NO>::updateParams( ParameterListPtr_Type params){
	params_ = params;
};

}
#endif
