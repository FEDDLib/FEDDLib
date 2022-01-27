#include "AssembleFEAceLaplace_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFEAceLaplace_def.hpp"
namespace FEDD {
    template class AssembleFEAceLaplace<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

