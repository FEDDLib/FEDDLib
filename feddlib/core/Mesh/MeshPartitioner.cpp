#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "MeshPartitioner_decl.hpp"
#include "MeshPartitioner_def.hpp"
namespace FEDD {
    template class MeshPartitioner<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

