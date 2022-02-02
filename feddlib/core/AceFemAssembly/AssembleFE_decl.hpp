#ifndef ASSEMBLEFE_DECL_hpp
#define ASSEMBLEFE_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"

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
    @todo This should actually by removed since the class should operate only on element level)
    \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

    Any new assembly routine on element level should implemented following the interface provided in this class. During the setup of a specific boundary value problem one AssembleFE object will be constructed using the AssembleFEFactory for each finite element. This is can be understood roughly as follows:
    \code
    for (int i=1; i<numElements; i++) {
        AssembleFE assmeblyFe[i] = AssembleFEFactory<>::build("problemType",flag,nodesRefConfig,params);
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

        /// @todo This matrix type should not be needed
        typedef Matrix<SC,LO,GO,NO> Matrix_Type;
        typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

        /// @todo Can we use a vector type which does not need a map?
        typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
        typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

        typedef SmallMatrix<SC> SmallMatrix_Type;
        typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

        typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;

        /*!
         \brief Assemble the element Jacobian matrix.
         \return the element Jacobian matrix
        */
        virtual SmallMatrixPtr_Type &elementMatrix assembleJacobian() = 0;

        /*!
         \brief Assemble the element right hand side vector.
         \return the element right hand side vector
        */
        virtual MultiVectorPtr_Type &elementVector assembleRHS() = 0;

        /*!
         \brief Get the currently assembled element Jacobian matrix
         \return the element Jacobian matrix
        */
        virtual SmallMatrixPtr_Type &elementMatrix getJacobian() = 0;

        /*!
         \brief Get the currently assembled right hand side vector.
         \return the element right hand side vector
        */
        virtual MultiVectorPtr_Type &elementVector getRHS() = 0;

        //virtual void assembleMass(MatrixPtr_Type &A) =0;

        /*!
         \brief Check the input parameters from the constructor and the ParameterList for consistency.
        */
        void checkParameters();

        /*!
         \brief Update the ParameterList initially specified at construction time.
        */
        void updateParams( ParameterListPtr_Type params);

        /*!
         \brief This function is called every time the FEDDLib proceeds from one to the next time step. The size of the time step will always be provided as input.
         @param[in] dt Timestepping length
        */
        void advanceInTime(double dt);

        /*!
         \brief Get the time state of the object.
         \return the time.
        */
        double getTime();

        /*!
         \brief Update the solution vector.
         @todo We still have to fix the ordering of the dofs.
        */
        void updateSolution(vec_dbl_Type solution);

        /*!
         \brief Get the current local solution vector.
         \return the solution vector.
        */
        vec_dbl_Type getSolution();

        /*!
         \brief This function is called in the beginning of each Newton step before actually assmblying anything.
        */
        void preProcessing();

        /*!
         \brief This function is called at the end of each Newton step after updating the solution vector.
        */
        void postProcessing(); // gespeicherte Daten rein und rausziehen. Manche Sachen nachträglich nicht mehr ändern


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
         @todo Should this be "changeRHSFunction()" instead? Probably, we want to specificy a RHS per defaiult at construction. This function is only for changing it, right?
        */
        void addRHSFunc(RhsFunc_Type rhsFunc){ rhsFunc_ = rhsFunc;};

    protected:

        AssembleFE(int flag,
                   vec2D_dbl_Type nodesRefConfig,
                   ParameterListPtr_Type parameters);

        // For example if we consider an element with multiple discretizations i.e. P2-P1 Velocity-Pressure.
        // We might need more than one FEType or Degree of freedom information
        // Vectoren mit Informationen besser so abepsicher sieher z.B. bcBuilder.
        int numFEType_;
        string FEType1_;
        string FEType2_;
        int dofs1_;
        int dofs2_;

        RhsFunc_Type rhsFunc_;

        int dim_;

        vec2D_dbl_Type nodesRefConfig_;
        bool timeProblem_;
        int flag_;
        double timeStep_ ;
        ParameterListPtr_Type paramsMaterial_;// Including flags
        ParameterListPtr_Type params_;// Including flags
        vec_dbl_Type solution_ ;

        friend class AssembleFEFactory<SC,LO,GO,NO>;
    };
}
#endif
