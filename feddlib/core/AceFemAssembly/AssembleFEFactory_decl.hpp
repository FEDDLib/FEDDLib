#ifndef ASSEMBLEFEFACTORY_DECL_hpp
#define ASSEMBLEFEFACTORY_DECL_hpp


#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFEAceLaplace.hpp"

namespace FEDD {

    /*!
    \class AssembleFEFactory
    \brief This class allows for constructing AssembleFE objects.

    \tparam SC The scalar type. So far, this is always double, but having it as a template parameter would allow flexibily, e.g., for using complex instead
    \tparam LO The local ordinal type. The is the index type for local indices
    \tparam GO The global ordinal type. The is the index type for global indices (this should actually by removed since the class should operate only on element level)
    \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

    This class provides a function for constructing objects of classes derived from AssembleFE. In particular, an AssembleFE object is constructed as follows:
    \code
    AssembleFE assmeblyFe = AssembleFEFactory<>::build("problemType",flag,nodesRefConfig,params);
    \endcode
    If a new element assembly is implemented, it has to be added to the build() function. Only then, this element assembly can be executed within the FEDDLib.
    */
    template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
    class AssembleFEFactory {
    public:

        typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;
        typedef Teuchos::RCP<AssembleFE_Type> AssembleFEPtr_Type;

        AssembleFEFactory();

        AssembleFEPtr_Type build( string problemType, int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params);

    };

}
#endif
