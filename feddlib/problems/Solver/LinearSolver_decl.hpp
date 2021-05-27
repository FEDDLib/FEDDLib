#ifndef LinearSolver_DECL_hpp
#define LinearSolver_DECL_hpp

#include "feddlib/problems/problems_config.h"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/problems/abstract/TimeProblem.hpp"

#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_DefaultZeroLinearOp_decl.hpp>
/*!
 Declaration of LinearSolver
 
 @brief  LinearSolver
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template<class SC_, class LO_, class GO_, class NO_>
class Problem;
template<class SC_, class LO_, class GO_, class NO_>
class TimeProblem;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class LinearSolver {

public:
    
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;
    typedef TimeProblem<SC,LO,GO,NO> TimeProblem_Type;
    typedef Teuchos::RCP<const Thyra::LinearOpBase<SC> > ThyraLinOpConstPtr_Type;
    typedef Thyra::BlockedLinearOpBase<SC> ThyraLinOpBlock_Type;
    typedef Teuchos::RCP< ThyraLinOpBlock_Type > ThyraLinOpBlockPtr_Type;
    typedef Teuchos::RCP< const ThyraLinOpBlock_Type > ThyraLinOpBlockConstPtr_Type;
    typedef Teuchos::RCP<Thyra::PreconditionerBase<SC> > ThyraPrecPtr_Type;

    typedef typename Problem_Type::Matrix_Type Matrix_Type;    
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;    
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

//    typedef typename Problem_Type::Preconditioner_Type Preconditioner_Type;
//    typedef typename Problem_Type::PreconditionerPtr_Type PreconditionerPtr_Type;
    
    LinearSolver();
    
    ~LinearSolver();
    
    int solve(Problem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string type="Monolithic" );
    
    int solve(TimeProblem_Type* timeProblem, BlockMultiVectorPtr_Type rhs, std::string type="Monolithic" );
    
    int solveMonolithic(Problem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string type );
    
    int solveMonolithic(TimeProblem_Type* problem, BlockMultiVectorPtr_Type rhs );
    
#ifdef FEDD_HAVE_TEKO
    int solveTeko(Problem_Type* problem, BlockMultiVectorPtr_Type rhs );

    int solveTeko(TimeProblem_Type* problem, BlockMultiVectorPtr_Type rhs );
#endif
    int solveBlock(Problem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string precType );
    
    int solveBlock(TimeProblem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string precType );
    
private:

};
}
#endif
