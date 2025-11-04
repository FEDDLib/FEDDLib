#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "AssembleFEFactory_decl.hpp"
#include "AssembleFEFactory_def.hpp"
namespace FEDD {
    template class AssembleFEFactory<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

