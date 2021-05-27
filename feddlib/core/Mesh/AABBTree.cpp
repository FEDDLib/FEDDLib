#include "AABBTree_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AABBTree_def.hpp"
namespace FEDD {
    template class AABBTree<default_sc, default_lo, default_go, default_no>; // What to do here?
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
