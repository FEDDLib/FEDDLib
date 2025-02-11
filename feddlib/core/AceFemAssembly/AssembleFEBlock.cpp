#include "AssembleFEBlock_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFEBlock_def.hpp"
namespace FEDD {
    template class AssembleFEBlock<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

