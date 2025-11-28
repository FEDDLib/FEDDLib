#ifndef DEFAULTTYPEDEFS_hpp
#define DEFAULTTYPEDEFS_hpp
#include <TpetraCore_config.h>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>

typedef double default_sc;
typedef int default_lo;
#if defined HAVE_TPETRA_INT_LONG_LONG
typedef long long default_go;
#define DEFAULT_GO_IS_LONG_LONG
#elif !defined HAVE_TPETRA_INT_LONG_LONG && defined HAVE_TPETRA_INT_INT // TODO: [JK] Can this be removed?
typedef int default_go;
#else
typedef long default_go;
#endif
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType default_no;

#endif
