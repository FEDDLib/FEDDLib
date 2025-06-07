#ifndef NONLINEARSOLVER_DECL_hpp
#define NONLINEARSOLVER_DECL_hpp
#include "feddlib/problems/problems_config.h"
#include "feddlib/problems/abstract/NonLinearProblem.hpp"
#include "feddlib/problems/abstract/TimeProblem.hpp"
#ifdef FEDD_HAVE_NOX
#include <NOX.H>
#include <NOX_Thyra.H>
#include <NOX_SolverStats.hpp>
#endif
/*!
 Declaration of NonLinearSolver
 
 @brief  NonLinearSolver
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD{
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class NonLinearSolver {
    
public:
    
    typedef Problem<SC,LO,NO,GO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;
    
    typedef TimeProblem<SC,LO,GO,NO> TimeProblem_Type;
    typedef NonLinearProblem<SC,LO,GO,NO> NonLinearProblem_Type;
    
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef Teuchos::RCP<ExporterTxt> ExporterTxtPtr_Type;
    
    NonLinearSolver();
    
    NonLinearSolver(std::string type);
    
    ~NonLinearSolver();
    
    void solve(NonLinearProblem_Type& problem);

    void solve(TimeProblem_Type& problem, double time=0., vec_dbl_ptr_Type valuesForExport = Teuchos::null );
    
	int getNonLinIts() {return nonLinearIts_;};
private:
#ifdef FEDD_HAVE_NOX
    void solveNOX(NonLinearProblem_Type& problem);
    
    void solveNOX( TimeProblem_Type& problem, vec_dbl_ptr_Type valuesForExport = Teuchos::null );
    
#endif
    void solveFixedPoint(NonLinearProblem_Type& problem);

    void solveNewton(NonLinearProblem_Type& problem);
    
    void solveFixedPoint(TimeProblem_Type& problem, double time);
    
    void solveNewton(TimeProblem_Type& problem, double time, vec_dbl_ptr_Type valuesForExport = Teuchos::null );

    void solveExtrapolation(TimeProblem_Type& problem, double time);
    
    std::string 	type_;

	int nonLinearIts_ =0;


};
}
#endif
