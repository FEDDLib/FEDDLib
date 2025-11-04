#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "FE_decl.hpp"
#include "FE_def.hpp"
namespace FEDD {
    template class FE<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

