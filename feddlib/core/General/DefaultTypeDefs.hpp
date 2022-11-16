#ifndef DEFAULTTYPEDEFS_hpp
#define DEFAULTTYPEDEFS_hpp
#include <TpetraCore_config.h>
#include <Kokkos_DefaultNode.hpp>

typedef double default_sc;
typedef int default_lo;
#if defined HAVE_TPETRA_INT_LONG_LONG
typedef long long default_go;
#elif !defined HAVE_TPETRA_INT_LONG_LONG && defined HAVE_TPETRA_INT_INT && defined HAVE_XPETRA_EPETRA
typedef int default_go;
#else
typedef long default_go;
#endif
typedef KokkosClassic::DefaultNode::DefaultNodeType default_no;
#endif
