#include "TimeProblem_decl.hpp"
#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "TimeProblem_def.hpp"
//    template class Problem<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
namespace FEDD {
    template class TimeProblem<default_sc, default_lo, default_go, default_no>;
}

#endif