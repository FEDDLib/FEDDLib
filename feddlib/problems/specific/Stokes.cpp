#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "Stokes_decl.hpp"
#include "Stokes_def.hpp"
namespace FEDD {
    template class Stokes<default_sc, default_lo, default_go, default_no>;
}
#endif
