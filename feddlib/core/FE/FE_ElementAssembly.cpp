#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "FE_ElementAssembly_decl.hpp"
#include "FE_ElementAssembly_def.hpp"
namespace FEDD {
    template class FE_ElementAssembly<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

