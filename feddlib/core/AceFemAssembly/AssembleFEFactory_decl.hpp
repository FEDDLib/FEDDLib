#ifndef ASSEMBLEFEFACTORY_DECL_hpp
#define ASSEMBLEFEFACTORY_DECL_hpp


#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/AssembleFEBlock.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_Laplace.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFENonLinLaplace.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_LinElas.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_NonLinElas.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_NonLinElas2.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFENavierStokes.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFEGeneralizedNewtonian.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_SCI_NH_decl.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_SCI_SMC_MLCK_decl.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp"

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
    template <class SC , class LO , class GO , class NO>
    class AssembleFEFactory {
    public:

        typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;
        typedef Teuchos::RCP<AssembleFE_Type> AssembleFEPtr_Type;
        /*!
         \brief Constructor

        */
        AssembleFEFactory();

        /*!
         \brief We only need this one function to build assembleFE, where we define the problem type and the AssembleFE object we extend to a specific problem.

         @param[in] problemType Type of specfic problem we focus on
         @param[in] flag Flag of element
         @param[in] nodesRefConfig Nodes of element in reference configuration
         @param[in] params Parameterlist for current problem
        */
        AssembleFEPtr_Type build( string problemType, int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple);

    };

}
#endif
