#include "AssembleFEAceLinElas_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFEAceLinElas_def.hpp"
namespace FEDD {
    template class AssembleFEAceLinElas<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

