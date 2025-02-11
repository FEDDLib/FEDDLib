#include "NavierStokesAssFE_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "NavierStokesAssFE_def.hpp"
namespace FEDD{
    template class NavierStokesAssFE<default_sc, default_lo, default_go, default_no>;
}
#endif
