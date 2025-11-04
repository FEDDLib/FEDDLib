#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFE_SCI_NH_decl.hpp"
#include "AssembleFE_SCI_NH_def.hpp"
namespace FEDD {
    template class AssembleFE_SCI_NH<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
