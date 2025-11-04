#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFE_decl.hpp"
#include "AssembleFE_def.hpp"
namespace FEDD {
    template class AssembleFE<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

