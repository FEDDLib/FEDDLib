#ifndef ASSEMBLEFE_DECL_hpp
#define ASSEMBLEFE_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/FE/Helper.hpp"

namespace FEDD {

   template <class SC = default_sc,
              class LO = default_lo,
              class GO = default_go,
              class NO = default_no>
    class AssembleFEFactory;
	
    /*!
    \class AssembleFE
    \brief This abstract class defining the interface for any type of element assembly rountines in the FEDDLib.

    \tparam SC The scalar type. So far, this is always double, but having it as a template parameter would allow flexibily, e.g., for using complex instead
    \tparam LO The local ordinal type. The is the index type for local indices
    \tparam GO The global ordinal type. The is the index type for global indices
    @todo This should actually be removed since the class should operate only on element level)
    \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

    Any new assembly routine on element level should implemented following the interface provided in this class. During the setup of a specific boundary value problem one AssembleFE object will be constructed using the AssembleFEFactory for each finite element. This is can be understood roughly as follows:
    \code
    for (int i=1; i<numElements; i++) {
        AssembleFE assmeblyFe[i] = AssembleFEFactory<>::build("problemType",flag,nodesRefConfig,params,tuple);
    }
    \endcode
    It is not possible to construct an AssembleFE object without using the AssembleFEFactory since the constructor is protected and hence not directly accessible.

    Similar to constructing the AssembleFE, all other member functions will be called automatically by the FEDDLib during the program flow. For instance, the assembly of the element Jacobian matrices will be performed:
    \code
    for (int i=1; i<numElements; i++) {
        assmeblyFe[i].assembleJacobian();
        Matrix_Type elementJacobian[i] = assmeblyFe[i].getJacobian();
    }
    \endcode
    A specific implementation of a class derived from AssembleFE can only interact with the FEDDLib by implementing the public member functions in AssembleFE for
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

    Additional parameters, such as material parameters, can provided through a Teuchos::ParameterList object which will contain all the parameters specified in the input file `ABC.xml`. The structure of the input file and, hence, of the resulting parameter list can be chosen freely depending on the specific implementation of an element assembly. The FEDDLib will take care of reading the parameters from the file and making them available to every AssembleFE object.
    */
    template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
    class AssembleFE {
    public:


        typedef SmallMatrix<SC> SmallMatrix_Type;
        typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

	    // Teuchos:: Array anstatt vec_dbl_Type

        typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;

        /*!
         \brief Compute everything.
        */
        //virtual void compute() = 0;

        /*!
         \brief Assemble the element Jacobian matrix.
        */
        virtual void assembleJacobian() = 0;


        /*!
         \brief Assemble the element Jacobian matrix.
         @param[in] block ID i
        */
        virtual void assembleJacobianBlock(LO i) = 0;

        /*!
         \brief Assemble the element right hand side vector.
        */
        virtual void assembleRHS() = 0;

        /*!
         \brief Get the currently assembled element Jacobian matrix
         \return the element Jacobian matrix
        */
        SmallMatrixPtr_Type getJacobian() {return jacobian_;}; 

        /*!
         \brief Get the currently assembled element Jacobian matrix
         \return the element Jacobian matrix
        */
        SmallMatrixPtr_Type getJacobianBlock(LO i) {return jacobianBlock_;}; 

        /*!
         \brief Get the currently assembled right hand side vector.
         \return the element right hand side vector
        */
        vec_dbl_ptr_Type getRHS(){return rhsVec_;};

        //virtual void assembleMass(MatrixPtr_Type &A) =0;

        /*!
         \brief Check the input parameters from the constructor and the ParameterList for completeness and consistency.
        */
        virtual void checkParameters();

        /*!
         \brief Set or update the parameters read from the ParameterList.
         @param[in] ParameterList as read from the xml file
        */
        virtual void updateParams(ParameterListPtr_Type params);

        /*!
         \brief Update the parameter read from the ParameterList.
         @param[in] Parameter as read from the xml file
        */
        virtual void updateParameter(string type, double value) {};
        /*!
         \brief This function is called every time the FEDDLib proceeds from one to the next time step. The size of the time step will always be provided as input.
         @param[in] dt Timestepping length
        */
        virtual void advanceInTime(double dt);
        /*!
         \brief Get the time state of the object.
         \return the timestep
        */
        double getTimeStep();

        /*!
         \brief This function is called every time the FEDDLib proceeds from one to the next newton step. The size of the time step will always be provided as input. 
        */
        void advanceNewtonStep();

        /*!
         \brief Get the time state of the object.
         \return newtonStep.
        */
        int getNewtonStep();

        /*!
         \brief Update the solution vector.
         @todo We still have to fix the ordering of the dofs.
         @param[in] solution
        */
        void updateSolution(vec_dbl_Type solution);

        /*!
         \brief Get the current local solution vector.
         \return the solution vector.
        */
        vec_dbl_ptr_Type getSolution();

        /*!
         \brief This function is called in the beginning of each Newton step before actually assmblying anything.
        */
        void preProcessing();

        /*!
         \brief This function is called at the end of each Newton step after updating the solution vector.
        */
        void postProcessing();
		/// @todo PostProcessing: Teuchos::Array with values and one global Array with Strings and names

        /*!
         \brief Get the spatial dimension. (Typically 2 or 3)
         \return dimension.
        */
        int getDim();

        /*!
         \brief Return the coordnates of the finite element nodes.
         \return a 2D array with the coordnates for each nodes
         @todo How is the ordering?
        */
        vec2D_dbl_Type getNodesRefConfig();

        /*!
         @todo Still work in Progress with RHS and Mass Matrix
        */
        void addRHSFunc(RhsFunc_Type rhsFunc){ rhsFunc_ = rhsFunc;};

        /*!
         \brief Return vector of tupled with element based values. First column per tuple string with description, second column with corresponding value
			\return elementInformation
        */
		tuple_sd_vec_ptr_Type getTupleElement(){return elementIntormation_;};

        /*!
         \brief Returns the time increment. Required by AceGen implementation.
            \return timeIncrement
        */
        double getTimeIncrement(){return timeIncrement_;};

        void setGlobalElementID(GO goID){globalElementID_ = goID;};

        GO getGlobalElementID(){return globalElementID_;};

        /*!
        \brief E.g. In case of non-newtonian fluids the viscosity is not constant - Compute the viscosity for an element depending on the known velocity solution
        */
	    virtual void computeLocalconstOutputField() {TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "computeLocalconstOutputField not yet implemented"); };
          /*!

        /*!
         \brief Obtain value of resulting postprocessing field at nodes/ inside an element
         \return values
        */
        vec_dbl_Type getLocalconstOutputField() {return constOutputField_;};


    protected:

        /*!
         \brief Constructor
         @param[in] flag Flag of element
         @param[in] nodesRefConfig Nodes of element in reference configuration
         @param[in] params Parameterlist for current problem
		 @param[in] tuple vector of element information tuples. 
        */
        AssembleFE(int flag,
                   vec2D_dbl_Type nodesRefConfig,
                   ParameterListPtr_Type parameters,
		   		tuple_disk_vec_ptr_Type tuple);

		//void readTuple(); /// @todo To have tuple information in basis class as well?
		//tuple_disk_vec_ptr_Type getTuple();  

		SmallMatrixPtr_Type jacobian_;
   		SmallMatrixPtr_Type jacobianBlock_;

		vec_dbl_ptr_Type rhsVec_;

        RhsFunc_Type rhsFunc_;

        int dim_;

		tuple_disk_vec_ptr_Type diskTuple_;
		tuple_sd_vec_ptr_Type elementIntormation_;
        /// @todo Why "Reference Configuration"? 
        vec2D_dbl_Type nodesRefConfig_;
        bool timeProblem_;
        int flag_;
        double timeStep_ ; 
        int newtonStep_ ;
        ParameterListPtr_Type paramsMaterial_;
        ParameterListPtr_Type params_;
        vec_dbl_ptr_Type solution_ ;
        double timeIncrement_;
        GO globalElementID_;


        // This can be any postprocessing output field ddefined inside an element using converged solution
        vec_dbl_Type constOutputField_ ; // can be a vector with values on P1/ P2 nodes or just averaged element value

        friend class AssembleFEFactory<SC,LO,GO,NO>;
    };
}
#endif
