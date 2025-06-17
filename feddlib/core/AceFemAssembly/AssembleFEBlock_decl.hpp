
#ifndef AssembleFEBlock_DECL_hpp
#define AssembleFEBlock_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/FE/Helper.hpp"
#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"


namespace FEDD {
	
    /*!
    \class AssembleFEBlock
    \brief This abstract class defining the interface for any type of element assembly rountines in the FEDDLib.

    \tparam SC The scalar type. So far, this is always double, but having it as a template parameter would allow flexibily, e.g., for using complex instead
    \tparam LO The local ordinal type. The is the index type for local indices
    \tparam GO The global ordinal type. The is the index type for global indices
    @todo This should actually be removed since the class should operate only on element level)
    \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

    Any new assembly routine on element level should implemented following the interface provided in this class. During the setup of a specific boundary value problem one AssembleFEBlock object will be constructed using the AssembleFEBlockFactory for each finite element. This is can be understood roughly as follows:
    \code
    for (int i=1; i<numElements; i++) {
        AssembleFEBlock assmeblyFe[i] = AssembleFEBlockFactory<>::build("problemType",flag,nodesRefConfig,params,tuple);
    }
    \endcode
    It is not possible to construct an AssembleFEBlock object without using the AssembleFEBlockFactory since the constructor is protected and hence not directly accessible.

    Similar to constructing the AssembleFEBlock, all other member functions will be called automatically by the FEDDLib during the program flow. For instance, the assembly of the element Jacobian matrices will be performed:
    \code
    for (int i=1; i<numElements; i++) {
        assmeblyFe[i].assembleJacobian();
        Matrix_Type elementJacobian[i] = assmeblyFe[i].getJacobian();
    }
    \endcode
    A specific implementation of a class derived from AssembleFEBlock can only interact with the FEDDLib by implementing the public member functions in AssembleFEBlock for
    - Construction
    - Assmebly of the Jacobian and right hand side
    - Getting the Jacobian and right hand side
    - Upating the solution
    - ...

    They will be automatically executed as the construction and assembly of the Jacobian; see above.

    If additional public member functions are added, they will not be executed from the FEDDLib. Therefore, we only allow for adding additional protected or private functions.

    Upon construction, the FEDDLib will provide some information, such as
    - The element flag
    - The coordinates of the finite element nodes
    - ...

    Additional parameters, such as material parameters, can provided through a Teuchos::ParameterList object which will contain all the parameters specified in the input file `ABC.xml`. The structure of the input file and, hence, of the resulting parameter list can be chosen freely depending on the specific implementation of an element assembly. The FEDDLib will take care of reading the parameters from the file and making them available to every AssembleFEBlock object.
    */
    template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
    class AssembleFEBlock: public AssembleFE<SC,LO,GO,NO>  {
    public:


        typedef SmallMatrix<SC> SmallMatrix_Type;
        typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

	    // Teuchos:: Array anstatt vec_dbl_Type

        typedef AssembleFEBlock<SC,LO,GO,NO> AssembleFEBlock_Type;

        /*!
         \brief Assemble the element Jacobian matrix.
         \return the element Jacobian matrix
        */
        virtual void assembleJacobian() = 0;

        /*!
         \brief Assemble the element right hand side vector.
         \return the element right hand side vector
        */
        virtual void assembleRHS() = 0;

        virtual void assembleJacobianBlock(LO i) =0;
        //virtual void assembleMass(MatrixPtr_Type &A) =0;

    protected:

        /*!
         \brief Constructor
         @param[in] flag Flag of element
         @param[in] nodesRefConfig Nodes of element in reference configuration
         @param[in] params Parameterlist for current problem
		 @param[in] tuple vector of element information tuples. 
        */
        AssembleFEBlock(int flag,
                   vec2D_dbl_Type nodesRefConfig,
                   ParameterListPtr_Type parameters,
		   		tuple_disk_vec_ptr_Type tuple);

		//void readTuple(); /// @todo To have tuple information in basis class as well?
		//tuple_disk_vec_ptr_Type getTuple();  
    private:

        void assembleMonolithicSystem(SmallMatrixPtr_Type elementMatrix);
		
        int dimSystem_;

	    std::string FEType_ ; // FEType of Disk

	    int dofsSolid_ ; // Degrees of freedom per node
		int dofsChem_;
	    int numNodesSolid_ ; // Number of nodes of element
		int numNodesChem_ ; // Number of nodes of element


	    int dofsElement_; // "Dimension of return matrix"
		

        friend class AssembleFEFactory<SC,LO,GO,NO>;
    };
}
#endif
