#include "AssembleFE_NonLinElas_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFE_NonLinElas_def.hpp"
namespace FEDD {
    template class AssembleFE_NonLinElas<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

