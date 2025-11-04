#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AdaptiveMeshRefinement_decl.hpp"
#include "AdaptiveMeshRefinement_def.hpp"
namespace FEDD {
    template class AdaptiveMeshRefinement<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

