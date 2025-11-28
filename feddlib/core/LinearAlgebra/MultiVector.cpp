#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifdef HAVE_EXPLICIT_INSTANTIATION
#include "MultiVector_decl.hpp"
#include "MultiVector_def.hpp"
namespace FEDD {
    template class MultiVector<default_sc, default_lo, default_go, default_no>;
    template class MultiVector<default_lo, default_lo, default_go, default_no>;
#ifdef DEFAULT_GO_IS_LONG_LONG
    // TODO: this instantiation of MultiVector is used in amr. Is this necessary or would MultiVector<LO, LO, GO, NO> suffice?                                                                                                                           
    template class MultiVector<default_go, default_lo, default_go, default_no>;
#endif
}
#endif  // HAVE_EXPLICIT_INSTANTIATION
#endif
