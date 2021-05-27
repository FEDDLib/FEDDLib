#include "Matrix_decl.hpp"

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "Matrix_def.hpp"
namespace FEDD {
    template class Matrix<default_sc, default_lo, default_go, default_no>;
}
#endif  // HAVE_EXPLICIT_INSTANTIATION

