#ifndef ACEGENINTERFACECHECK_hpp
#define ACEGENINTERFACECHECK_hpp

#include "feddlib/core/core_config.h"

#include <Teuchos_GlobalMPISession.hpp>

#include <cstdlib>
#include <iostream>

namespace FEDD {

/*!
 \brief For cases that need the AceGen elements (Interface2).

 Returns true if FEDDLib was built with the AceGen interface. Otherwise rank 0 explains why the case
 cannot run and false is returned, so that the case stops with an error right away instead of failing
 later in the element assembly. Call it after the MPI session has been created:

     if (!FEDD::aceGenInterfaceAvailable())
         return EXIT_FAILURE;
*/
inline bool aceGenInterfaceAvailable()
{
#ifdef FEDD_HAVE_ACEGENINTERFACE
    return true;
#else
    if (Teuchos::GlobalMPISession::getRank() == 0)
        std::cerr << "This case needs FEDDLib built with the AceGen interface (Interface2); "
                     "configure FEDDLib with -D TPL_ENABLE_AceGENInterface=ON." << std::endl;
    return false;
#endif
}

}

#endif
