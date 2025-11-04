#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "LinearProblem_decl.hpp"
#include "LinearProblem_def.hpp"
namespace FEDD {
    template class LinearProblem<default_sc, default_lo, default_go, default_no>;
}
#endif
