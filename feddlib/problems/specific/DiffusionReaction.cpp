#include "DiffusionReaction_decl.hpp"
#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "DiffusionReaction_def.hpp"
//template class Diffusion<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
namespace FEDD {
template class DiffusionReaction<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
