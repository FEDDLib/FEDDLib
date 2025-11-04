#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "HDF5Export_decl.hpp"
#include "HDF5Export_def.hpp"
namespace FEDD {
template class HDF5Export<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
