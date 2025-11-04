#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "HDF5Import_decl.hpp"
#include "HDF5Import_def.hpp"
namespace FEDD {
template class HDF5Import<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
