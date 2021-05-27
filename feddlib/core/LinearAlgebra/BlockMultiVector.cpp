#include "BlockMultiVector_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "BlockMultiVector_def.hpp"
namespace FEDD {
    template class BlockMultiVector<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

