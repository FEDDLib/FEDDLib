#ifndef LINEARPROBLEM_DECL_hpp
#define LINEARPROBLEM_DECL_hpp
#include "Problem.hpp"

/*!
 Declaration of LinearProblem
 
 @brief  LinearProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template<class SC_, class LO_, class GO_, class NO_>
class Problem;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class LinearProblem : public Problem<SC,LO,GO,NO> {
//class LinearProblem : public Problem {
public:
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;
    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;
    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

    LinearProblem(CommConstPtr_Type comm);

    LinearProblem(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm);

    ~LinearProblem();

    virtual void assemble() = 0;

    virtual void getValuesOfInterest( vec_dbl_Type& values ) = 0;
    
private:


};
}


#endif
