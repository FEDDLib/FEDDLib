#include "LaplaceBlocks_decl.hpp"
#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "LaplaceBlocks_def.hpp"
//template class LaplaceBlocks<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
namespace FEDD {
template class LaplaceBlocks<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION