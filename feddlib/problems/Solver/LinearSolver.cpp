#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "LinearSolver_decl.hpp"
#include "LinearSolver_def.hpp"
namespace FEDD{
template class LinearSolver<default_sc, default_lo, default_go, default_no>;
}
#endif
