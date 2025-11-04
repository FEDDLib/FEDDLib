#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "Preconditioner_decl.hpp"
#include "Preconditioner_def.hpp"
namespace FEDD {
template class Preconditioner<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
