#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "TimeProblem_decl.hpp"
#include "TimeProblem_def.hpp"
namespace FEDD {
    template class TimeProblem<default_sc, default_lo, default_go, default_no>;
}

#endif
