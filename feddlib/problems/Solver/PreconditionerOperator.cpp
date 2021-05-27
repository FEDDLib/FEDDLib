#include "PreconditionerOperator_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "PreconditionerOperator_def.hpp"
namespace FEDD {
template class PreconditionerOperator<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
