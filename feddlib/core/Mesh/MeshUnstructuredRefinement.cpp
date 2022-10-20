#include "MeshUnstructuredRefinement_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "MeshUnstructuredRefinement_def.hpp"
namespace FEDD {
    template class MeshUnstructuredRefinement<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

