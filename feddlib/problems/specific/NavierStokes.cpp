#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "NavierStokes_decl.hpp"
#include "NavierStokes_def.hpp"
namespace FEDD{
    template class NavierStokes<default_sc, default_lo, default_go, default_no>;
}
#endif
