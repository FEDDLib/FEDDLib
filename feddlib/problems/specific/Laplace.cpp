#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "Laplace_decl.hpp"
#include "Laplace_def.hpp"
namespace FEDD {
template class Laplace<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
