#include "Elasticity_decl.hpp"
#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "Elasticity_def.hpp"
namespace FEDD{
template class Elasticity<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION