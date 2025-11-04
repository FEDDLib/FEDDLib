#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "NonLinearProblem_decl.hpp"
#include "NonLinearProblem_def.hpp"
namespace FEDD{
    template class NonLinearProblem<default_sc, default_lo, default_go, default_no>;
}
#endif
