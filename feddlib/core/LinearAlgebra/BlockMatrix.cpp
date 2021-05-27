
#include "BlockMatrix_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "BlockMatrix_def.hpp"
namespace FEDD {
    template class BlockMatrix<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

