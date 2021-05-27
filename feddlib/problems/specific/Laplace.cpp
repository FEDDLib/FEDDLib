#include "Laplace_decl.hpp"
#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "Laplace_def.hpp"
//template class Laplace<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
namespace FEDD {
template class Laplace<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION