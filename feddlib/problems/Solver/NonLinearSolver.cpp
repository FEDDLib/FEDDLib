#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "NonLinearSolver_decl.hpp"
#include "NonLinearSolver_def.hpp"
namespace FEDD{
    template class NonLinearSolver<default_sc, default_lo, default_go, default_no>;
}
#endif
