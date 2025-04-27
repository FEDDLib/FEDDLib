#ifndef MinPrecProblem_DECL_hpp
#define MinPrecProblem_DECL_hpp

#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/problems/abstract/Problem.hpp"

/*!
 Declaration of MinPrecProblem

 @brief  MinPrecProblem, this class is a minimal interface to a Problem.
 It provides all the information needed to build a preconditioner:
 the underlying system in BlockMatrix form and the corresponding domains.
 It is a concrete problem and not abstract.
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


namespace FEDD {
template <class SC , class LO , class GO , class NO >
class Problem;
template <class SC , class LO , class GO , class NO >
class Preconditioner;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MinPrecProblem : public Problem<SC,LO,GO,NO> {

public:
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Domain_Type Domain_Type;
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::DomainConstPtr_vec_Type DomainConstPtr_vec_Type;

    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;

    typedef typename Problem_Type::Comm_Type Comm_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef typename Problem_Type::LinSolverBuilderPtr_Type  LinSolverBuilderPtr_Type;
    
    // hasSourceTerm and boolLinearProblem should be irrelevant, as this class should be only used when constructing a precondtioner
    MinPrecProblem(ParameterListPtr_Type pl, CommConstPtr_Type comm);

    ~MinPrecProblem();
    
    void initializeSystem(BlockMatrixPtr_Type system);

    void initializeDomains(DomainConstPtr_vec_Type& domainPtr_vec);

    void initializeLinSolverBuilder(LinSolverBuilderPtr_Type linearSolverBuilder);
    
    virtual void info(){
        if ( this->comm_->getRank() )
            std::cout<< "Minimal implementation of Problem. This object has sufficient information for the setup of FROSch." << std::endl;
    }
    
    virtual void assemble( std::string type="" )const {}
    
    virtual void getValuesOfInterest( vec_dbl_Type& values ){}

    virtual void computeValuesOfInterestAndExport() {}

//    virtual void assembleExternal( std::string type ){}

protected:

private:

};
}
#endif
