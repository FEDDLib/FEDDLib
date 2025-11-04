#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFE_Laplace_decl.hpp"
#include "AssembleFE_Laplace_def.hpp"
namespace FEDD {
    template class AssembleFE_Laplace<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

