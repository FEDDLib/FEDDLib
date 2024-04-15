#ifndef DIFFERENTIABLEFUNCCLASS_DECL_hpp
#define DIFFERENTIABLEFUNCCLASS_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/General/InputToOutputMappingClass.hpp"
//#include "feddlib/core/LinearAlgebra/Matrix.hpp" maybe we will need something like this for multidimensional?


namespace FEDD {


    /*!
    \class DifferentiableFuncClass
    \brief This abstract class is derived from the abstract class of general input to output mapping.
           It is defining the general concepts of a function and computing function evaluation and evaluation of its derivative.

    \tparam SC The scalar type. So far, this is always double, but having it as a template parameter would allow flexibily, e.g., for using complex instead
    \tparam LO The local ordinal type. The is the index type for local indices
    \tparam GO The global ordinal type. The is the index type for global indices
    @todo This should actually be removed since the class should operate only on element level)
    \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

    The material parameters can be provided through a Teuchos::ParameterList object which will contain 
    all the parameters specified in the input file `ABC.xml`.
    The structure of the input file and, hence, of the resulting parameter list can be chosen freely. The FEDDLib will take care of reading the parameters 
    from the file and making them available.
    */
    template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
    class DifferentiableFuncClass : public InputToOutputMappingClass<SC,LO,GO,NO> {
    public:

        typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
        typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
        typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
      
        // Inherited Function from base abstract class 
        /*!
         \brief Implement a mapping for evaluating output in dependence of given input and specified parameters for MultiVectorPtr_Type inputs/ outputs
         @param[in] params Parameterlist as read from the xml file (maybe redundant)
         @param[in] x input
         @param[in,out] res output
        */
        virtual void evaluateMapping(ParameterListPtr_Type params, MultiVectorConstPtr_Type input, MultiVectorPtr_Type &output) override { TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not yet implemented for MultiVectorConstPtr_Type"); };

        /*!
         \brief Computes value of derivative of defined function in evaluateMapping
         @param[in] params Parameterlist as read from the xml file (maybe redundant)
         @param[in] x Independent variable
         @param[in,out] res Dependent variable
        */
        virtual void evaluateDerivative(ParameterListPtr_Type params, MultiVectorConstPtr_Type x, MultiVectorPtr_Type &res) override { TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not yet implemented for MultiVectorConstPtr_Type"); };

       
        /*!
         \brief Function could include different parameters which will be specified in *.xml
                This function should set the needed parameters for defined function to the specified values.
        @param[in] ParameterList as read from the xml file (maybe redundant)
        */
        virtual void setParams(ParameterListPtr_Type params) override = 0;

        /*!
        \brief Print parameter values used in model at runtime
        */
        virtual void echoInformationMapping() override = 0;

        /*!
         \brief Implements a functional description for evaluating res in dependence of given dependent variables and specified parameters
         @param[in] params Parameterlist as read from the xml file (maybe redundant)
         @param[in] x Independent variable
         @param[in,out] res Dependent variable
        */
        virtual void evaluateMapping(ParameterListPtr_Type params, double x, double &res) override = 0;

        /*!
         \brief Computes value of derivative of defined function in evaluateMapping
         @param[in] params Parameterlist as read from the xml file (maybe redundant)
         @param[in] x Independent variable
         @param[in,out] res Dependent variable
        */
        virtual void evaluateDerivative(ParameterListPtr_Type params, double x, double &res) override = 0;



    protected:

        /*!
         \brief Constructor
         @param[in] params Parameterlist for current problem
        */
        DifferentiableFuncClass(ParameterListPtr_Type parameters);


        ParameterListPtr_Type params_;

    };
}
#endif
