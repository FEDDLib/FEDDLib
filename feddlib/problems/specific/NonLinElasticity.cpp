#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "NonLinElasticity_decl.hpp"
#include "NonLinElasticity_def.hpp"
namespace FEDD{
template class NonLinElasticity<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
