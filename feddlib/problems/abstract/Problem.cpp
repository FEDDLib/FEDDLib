#include "Problem_decl.hpp"
#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "Problem_def.hpp"
namespace FEDD{
    template class Problem<default_sc, default_lo, default_go, default_no>;
}
#endif