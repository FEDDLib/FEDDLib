#ifndef DIFFERENTIABLEFUNCCLASS_DEF_hpp
#define DIFFERENTIABLEFUNCCLASS_DEF_hpp

#include "DifferentiableFuncClass_decl.hpp"

namespace FEDD {


template <class SC, class LO, class GO, class NO>
DifferentiableFuncClass<SC,LO,GO,NO>::DifferentiableFuncClass(ParameterListPtr_Type params):InputToOutputMappingClass<SC,LO,GO,NO>(params)
{

	params_=params;

}


}
#endif
