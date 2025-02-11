#include "AssembleFENonLinLaplace_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFENonLinLaplace_def.hpp"
namespace FEDD {
template class AssembleFENonLinLaplace<default_sc, default_lo, default_go,
                                       default_no>;
}
#endif // HAVE_EXPLICIT_INSTANTIATION
