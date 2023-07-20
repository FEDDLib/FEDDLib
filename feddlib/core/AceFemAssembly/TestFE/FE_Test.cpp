#include "FE_Test_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "FE_Test_def.hpp"
namespace FEDD {
    template class FE_Test<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

