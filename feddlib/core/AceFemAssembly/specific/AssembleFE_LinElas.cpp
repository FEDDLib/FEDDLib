#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFE_LinElas_decl.hpp"
#include "AssembleFE_LinElas_def.hpp"
namespace FEDD {
    template class AssembleFE_LinElas<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

