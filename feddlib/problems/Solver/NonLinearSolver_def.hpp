#ifndef NONLINEARSOLVER_DEF_hpp
#define NONLINEARSOLVER_DEF_hpp
#include "NonLinearSolver_decl.hpp"
/*!
 Definition of NonLinearSolver

 @brief  NonLinearSolver
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {


template<class SC,class LO,class GO,class NO>
NonLinearSolver<SC,LO,GO,NO>::NonLinearSolver():
type_("")
{}


template<class SC,class LO,class GO,class NO>
NonLinearSolver<SC,LO,GO,NO>::NonLinearSolver(string type):
type_(type)
{}

template<class SC,class LO,class GO,class NO>
NonLinearSolver<SC,LO,GO,NO>::~NonLinearSolver(){

}

template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solve(NonLinearProblem_Type &problem){

    if (!type_.compare("FixedPoint")) {
        solveFixedPoint(problem);
    }
    else if(!type_.compare("Newton")){
        solveNewton(problem);
    }
    else if(!type_.compare("NOX")){
#ifdef FEDD_HAVE_NOX
        solveNOX(problem);
#endif
    }

}

template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solve(TimeProblem_Type &problem, double time, vec_dbl_ptr_Type valuesForExport){

    if (!type_.compare("FixedPoint")) {
        solveFixedPoint(problem,time);
    }
    else if(!type_.compare("Newton")){
        solveNewton(problem,time, valuesForExport);
    }
    else if(!type_.compare("NOX")){
#ifdef FEDD_HAVE_NOX
        solveNOX(problem, valuesForExport);
#endif
    }
    else if(!type_.compare("Extrapolation")){
        solveExtrapolation(problem, time);
    }
}

#ifdef FEDD_HAVE_NOX
template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solveNOX(NonLinearProblem_Type &problem){

    bool verbose = problem.getVerbose();
    Teuchos::RCP<NonLinearProblem<SC,LO,GO,NO> > problemPtr = Teuchos::rcpFromRef(problem);
    Teuchos::RCP<Teuchos::ParameterList> p = sublist(problemPtr->getParameterList(),"ThyraSolver");
//    p->set("Preconditioner Type", "None"); // CH 16.04.19: preconditioner will be built seperately
//    sublist( sublist(p, "Linear Solver Types") , "Belos")->set("Left Preconditioner If Unspecified",true);
    problemPtr->getLinearSolverBuilder()->setParameterList(p);

    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory = problemPtr->getLinearSolverBuilder()->createLinearSolveStrategy("");

	//problemPtr->set_W_factory(lowsFactory);

    // Create the initial guess
    Teuchos::RCP<Thyra::VectorBase<SC> > initial_guess = problemPtr->getNominalValues().get_x()->clone_v();
    Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<SC>::zero());

    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = problemPtr->create_W_op();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > W_prec = problemPtr->create_W_prec();
    Teuchos::RCP<NOX::Thyra::Group> nox_group(new NOX::Thyra::Group(initial_guess,
                                                                    problemPtr.getConst(),
                                                                    W_op,
                                                                    lowsFactory.getConst(),
                                                                    W_prec,
                                                                    Teuchos::null,
                                                                    Teuchos::null,
                                                                    Teuchos::null));

    nox_group->computeF();

    // Create the NOX status tests and the solver
    // Create the convergence tests

    Teuchos::RCP<NOX::StatusTest::NormUpdate> updateTol =
        Teuchos::rcp(new NOX::StatusTest::NormUpdate( problemPtr->getParameterList()->sublist("Parameter").get("updateTol",1.0e-6) ) );
    
    Teuchos::RCP<NOX::StatusTest::RelativeNormF> relresid =
        Teuchos::rcp(new NOX::StatusTest::RelativeNormF( problemPtr->getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-4) ) );
    
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
        Teuchos::rcp(new NOX::StatusTest::NormWRMS(problemPtr->getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-4), problemPtr->getParameterList()->sublist("Parameter").get("absNonLinTol",1.0e-6)));
    
    Teuchos::RCP<NOX::StatusTest::NormF> absRes =
        Teuchos::rcp(new NOX::StatusTest::NormF( problemPtr->getParameterList()->sublist("Parameter").get("absNonLinTol",1.0e-6) ) );
    
    Teuchos::RCP<NOX::StatusTest::Combo> converged;
    
    if ( !problemPtr->getParameterList()->sublist("Parameter").get("Combo","AND").compare("AND") )
        converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    else if (!problemPtr->getParameterList()->sublist("Parameter").get("Combo","AND").compare("OR") )
        converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    
    if ( problemPtr->getParameterList()->sublist("Parameter").get("Use rel tol",true) )
        converged->addStatusTest(relresid);
    if ( problemPtr->getParameterList()->sublist("Parameter").get("Use update tol",false) )
        converged->addStatusTest(updateTol);
    if (problemPtr->getParameterList()->sublist("Parameter").get("Use WRMS",false))
        converged->addStatusTest(wrms);
    if (problemPtr->getParameterList()->sublist("Parameter").get("Use abs tol",false))
        converged->addStatusTest(absRes);

    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(problemPtr->getParameterList()->sublist("Parameter").get("MaxNonLinIts",10)));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    // Create nox parameter list
    Teuchos::RCP<Teuchos::ParameterList> nl_params = sublist(problemPtr->getParameterList(),"NOXSolver");

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(nox_group, combo, nl_params);
    NOX::StatusTest::StatusType solveStatus = solver->solve();
    double nonLinearIts = solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearSolves;
    double linearIts = solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearIterations;

    linearIts/=nonLinearIts;
    if (verbose){
        cout << "############################################################" << endl;
        cout << "### Total nonlinear iterations : " << nonLinearIts << "  with an average of " << linearIts << " linear iterations ###" << endl;
        cout << "############################################################" << endl;
    }
    
    if ( problemPtr->getParameterList()->sublist("Parameter").get("Cancel MaxNonLinIts",false) ) {
        TEUCHOS_TEST_FOR_EXCEPTION((int)nonLinearIts == problemPtr->getParameterList()->sublist("Parameter").get("MaxNonLinIts",10) ,std::runtime_error,"Maximum nonlinear Iterations reached. Problem might have converged in the last step. Still we cancel here.");
    }
	nonLinearIts_ = nonLinearIts;

}
    
template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solveNOX(TimeProblem_Type &problem, vec_dbl_ptr_Type valuesForExport){
    
    bool verbose = problem.getVerbose();
    Teuchos::RCP<TimeProblem_Type> problemPtr = Teuchos::rcpFromRef(problem);
    Teuchos::RCP<Teuchos::ParameterList> p = sublist(problemPtr->getParameterList(),"ThyraSolver");

    problemPtr->getLinearSolverBuilder()->setParameterList(p);
    
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> >
    lowsFactory = problemPtr->getLinearSolverBuilder()->createLinearSolveStrategy("");
    
    TEUCHOS_TEST_FOR_EXCEPTION(problemPtr->getSolution()->getNumVectors()>1, std::runtime_error, "With the current implementation NOX can only be used with 1 MultiVector column.");
    // Create the initial guess and fill with last solution
    Teuchos::RCP<Thyra::VectorBase<SC> > initialGuess = problemPtr->getNominalValues().get_x()->clone_v();
    // Try to convert to a ProductVB. If resulting pointer is not null we need to use the ProductMV below, otherwise it is a monolithic vector.
    Teuchos::RCP<Thyra::ProductVectorBase<SC> > initialGuessProd = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<SC> >(initialGuess);
    Teuchos::RCP<Thyra::MultiVectorBase<SC> > solMV;
    if (!initialGuessProd.is_null())
        solMV = problemPtr->getSolution()->getProdThyraMultiVector();
    else
        solMV = problemPtr->getSolution()->getThyraMultiVector();

    Thyra::assign(initialGuess.ptr(), *solMV->col(0));
    //Thyra::V_S(initialGuess.ptr(),Teuchos::ScalarTraits<SC>::zero());
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = problemPtr->create_W_op();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > W_prec = problemPtr->create_W_prec();
    Teuchos::RCP<NOX::Thyra::Group> nox_group(new NOX::Thyra::Group(initialGuess,
                                                                    problemPtr.getConst(),
                                                                    W_op,
                                                                    lowsFactory.getConst(),
                                                                    W_prec,
                                                                    Teuchos::null,
                                                                    Teuchos::null,
                                                                    Teuchos::null));
    
    nox_group->computeF();
    
   // Create the NOX status tests and the solver
    // Create the convergence tests
    Teuchos::RCP<NOX::StatusTest::NormUpdate> updateTol =
        Teuchos::rcp(new NOX::StatusTest::NormUpdate( problemPtr->getParameterList()->sublist("Parameter").get("updateTol",1.0e-6) ) );
    
    Teuchos::RCP<NOX::StatusTest::RelativeNormF> relresid =
        Teuchos::rcp(new NOX::StatusTest::RelativeNormF( problemPtr->getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-4) ) );
    
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms =
        Teuchos::rcp(new NOX::StatusTest::NormWRMS(problemPtr->getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-4), problemPtr->getParameterList()->sublist("Parameter").get("absNonLinTol",1.0e-6)));
    
    Teuchos::RCP<NOX::StatusTest::NormF> absRes =
        Teuchos::rcp(new NOX::StatusTest::NormF( problemPtr->getParameterList()->sublist("Parameter").get("absNonLinTol",1.0e-6) ) );
    
    Teuchos::RCP<NOX::StatusTest::Combo> converged;
    
    if ( !problemPtr->getParameterList()->sublist("Parameter").get("Combo","AND").compare("AND") )
        converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    else if (!problemPtr->getParameterList()->sublist("Parameter").get("Combo","AND").compare("OR") )
        converged = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    
    if ( problemPtr->getParameterList()->sublist("Parameter").get("Use rel tol",true) )
        converged->addStatusTest(relresid);
    if ( problemPtr->getParameterList()->sublist("Parameter").get("Use update tol",false) )
        converged->addStatusTest(updateTol);
    if (problemPtr->getParameterList()->sublist("Parameter").get("Use WRMS",false))
        converged->addStatusTest(wrms);
    if (problemPtr->getParameterList()->sublist("Parameter").get("Use abs tol",false))
        converged->addStatusTest(absRes);

    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
    Teuchos::rcp(new NOX::StatusTest::MaxIters(problemPtr->getParameterList()->sublist("Parameter").get("MaxNonLinIts",10)));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv =
    Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RCP<NOX::StatusTest::Combo> combo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);
    
    // Create nox parameter list
    Teuchos::RCP<Teuchos::ParameterList> nl_params = sublist(problemPtr->getParameterList(),"NOXSolver");
    
    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver =
    NOX::Solver::buildSolver(nox_group, combo, nl_params);
    NOX::StatusTest::StatusType solveStatus = solver->solve();
    
    double nonLinearIts = solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearSolves;
    double linearIts = solver->getSolverStatistics()->linearSolve.allNonlinearSolves_NumLinearIterations;
    
    linearIts/=nonLinearIts;
    if (verbose){
        cout << "############################################################" << endl;
        cout << "### Total nonlinear iterations : " << nonLinearIts << "  with an average of " << linearIts << " linear iterations ###" << endl;
        cout << "############################################################" << endl;
    }
    
    if ( problemPtr->getParameterList()->sublist("Parameter").get("Cancel MaxNonLinIts",false) ) {
        TEUCHOS_TEST_FOR_EXCEPTION((int)nonLinearIts == problemPtr->getParameterList()->sublist("Parameter").get("MaxNonLinIts",10) ,std::runtime_error,"Maximum nonlinear Iterations reached. Problem might have converged in the last step. Still we cancel here.");
    }
    
    if (!valuesForExport.is_null()) {
        if (valuesForExport->size() == 2){
            (*valuesForExport)[0] = linearIts;
            (*valuesForExport)[1] = nonLinearIts;
        }
    }
}
#endif

template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solveFixedPoint(NonLinearProblem_Type &problem){

    bool verbose = problem.getVerbose();
    TEUCHOS_TEST_FOR_EXCEPTION(problem.getRhs()->getNumVectors()!=1,std::logic_error,"We need to change the code for numVectors>1.");
    // -------
    // fix point iteration
    // -------
    double	gmresIts = 0.;
    double residual0 = 1.;
    double residual = 1.;
    
    double tol = problem.getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-6);
    int maxNonLinIts = problem.getParameterList()->sublist("Parameter").get("MaxNonLinIts",10);
    int nlIts=0;

    double criterionValue = 1.;
    std::string criterion = problem.getParameterList()->sublist("Parameter").get("Criterion","Residual");

    while ( nlIts < maxNonLinIts ) {

        problem.calculateNonLinResidualVec("reverse");

        problem.setBoundariesSystem();
        
        if (criterion=="Residual")
            residual = problem.calculateResidualNorm();
        
        if (nlIts==0)
            residual0 = residual;
    
        if (criterion=="Residual"){
            criterionValue = residual/residual0;
            if (verbose)
                cout << "### Fixed Point iteration : " << nlIts << "  relative nonlinear residual : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }


        gmresIts += problem.solveAndUpdate( criterion, criterionValue );
        nlIts++;
        if(criterion=="Update"){
            if (verbose)
                cout << "### Fixed Point iteration : " << nlIts << "  residual of update : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }
        // ####### end FPI #######
    }

    gmresIts/=nlIts;
    if (verbose)
        cout << "### Total FPI : " << nlIts << "  with average gmres its : " << gmresIts << endl;
    if ( problem.getParameterList()->sublist("Parameter").get("Cancel MaxNonLinIts",false) ) {
        TEUCHOS_TEST_FOR_EXCEPTION( nlIts == maxNonLinIts ,std::runtime_error,"Maximum nonlinear Iterations reached. Problem might have converged in the last step. Still we cancel here.");
    }
    
}

template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solveNewton( NonLinearProblem_Type &problem ){

    bool verbose = problem.getVerbose();

    TEUCHOS_TEST_FOR_EXCEPTION(problem.getRhs()->getNumVectors()!=1,std::logic_error,"We need to change the code for numVectors>1.")
    // -------
    // Newton
    // -------
    double	gmresIts = 0.;
    double residual0 = 1.;
    double residual = 1.;
    double tol = problem.getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-6);
    int nlIts=0;
    int maxNonLinIts = problem.getParameterList()->sublist("Parameter").get("MaxNonLinIts",10);
    double criterionValue = 1.;
    std::string criterion = problem.getParameterList()->sublist("Parameter").get("Criterion","Residual");

    while ( nlIts < maxNonLinIts ) {
        //this makes only sense for Navier-Stokes/Stokes, for other problems, e.g., non linear elasticity, it should do nothing.

        problem.calculateNonLinResidualVec("reverse");

        if (criterion=="Residual")
            residual = problem.calculateResidualNorm();

        problem.assemble("Newton");

        problem.setBoundariesSystem();

        if (nlIts==0)
            residual0 = residual;
        
        if (criterion=="Residual"){
            criterionValue = residual/residual0;
            if (verbose)
                cout << "### Newton iteration : " << nlIts << "  relative nonlinear residual : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }
        // PRINT INFOS
        gmresIts += problem.solveAndUpdate( criterion, criterionValue );
        nlIts++;
        if(criterion=="Update"){
            if (verbose)
                cout << "### Newton iteration : " << nlIts << "  residual of update : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }

        // ####### end FPI #######
    }

    gmresIts/=nlIts;
    if (verbose)
        cout << "### Total Newton iterations : " << nlIts << "  with average gmres its : " << gmresIts << endl;
    if ( problem.getParameterList()->sublist("Parameter").get("Cancel MaxNonLinIts",false) ) {
        TEUCHOS_TEST_FOR_EXCEPTION(nlIts == maxNonLinIts ,std::runtime_error,"Maximum nonlinear Iterations reached. Problem might have converged in the last step. Still we cancel here.");
    }
}

template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solveFixedPoint(TimeProblem_Type &problem, double time){

    bool verbose = problem.getVerbose();
    problem.setBoundariesRHS(time);
    TEUCHOS_TEST_FOR_EXCEPTION(problem.getRhs()->getNumVectors()!=1,std::logic_error,"We need to change the code for numVectors>1.")

    // -------
    // fix point iteration
    // -------
    double	gmresIts = 0.;
    double residual0 = 1.;
    double residual = 1.;
    double tol = problem.getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-6);
    int nlIts=0;
    int maxNonLinIts = problem.getParameterList()->sublist("Parameter").get("MaxNonLinIts",10);
    double criterionValue = 1.;
    std::string criterion = problem.getParameterList()->sublist("Parameter").get("Criterion","Residual");

    while ( nlIts < maxNonLinIts ) {
        
        problem.calculateNonLinResidualVec("reverse", time);

        if (criterion=="Residual")
            residual = problem.calculateResidualNorm();

        if (nlIts==0)
            residual0 = residual;
                    
        // Linearization of system matrix is done in calculateNonLinResidualVec
        // Now we need to combine it with the mass matrix
        problem.combineSystems();
        
        problem.setBoundariesSystem();
        
        if (criterion=="Residual"){
            criterionValue = residual/residual0;
            if (verbose)
                cout << "### Fixed Point iteration : " << nlIts << "  relative nonlinear residual : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }

        gmresIts += problem.solveAndUpdate( criterion, criterionValue );
        
        nlIts++;
        if(criterion=="Update"){
            if (verbose)
                cout << "### Fixed Point iteration : " << nlIts << "  residual of update : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }
        // ####### end FPI #######
    }
    
    gmresIts/=nlIts;
    if (verbose)
        cout << "### Total FPI : " << nlIts << "  with average gmres its : " << gmresIts << endl;
    if ( problem.getParameterList()->sublist("Parameter").get("Cancel MaxNonLinIts",false) ) {
        TEUCHOS_TEST_FOR_EXCEPTION( nlIts == maxNonLinIts ,std::runtime_error,"Maximum nonlinear Iterations reached. Problem might have converged in the last step. Still we cancel here.");
    }
}



template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solveNewton(TimeProblem_Type &problem, double time, vec_dbl_ptr_Type valuesForExport ){

    bool verbose = problem.getVerbose();
    problem.setBoundariesRHS(time);


    TEUCHOS_TEST_FOR_EXCEPTION(problem.getRhs()->getNumVectors()!=1,std::logic_error,"We need to change the code for numVectors>1.")
    
    // -------
    // Newton iteration
    // -------
    double	gmresIts = 0.;
    double residual0 = 1.;
    double residual = 1.;
    double tol = problem.getParameterList()->sublist("Parameter").get("relNonLinTol",1.0e-6);
    int nlIts=0;
    int maxNonLinIts = problem.getParameterList()->sublist("Parameter").get("MaxNonLinIts",10);
    double criterionValue = 1.;
    std::string criterion = problem.getParameterList()->sublist("Parameter").get("Criterion","Residual");
    std::string timestepping = problem.getParameterList()->sublist("Timestepping Parameter").get("Class","Singlestep");

    while ( nlIts < maxNonLinIts ) {
        if (timestepping == "External")
            problem.calculateNonLinResidualVec("external", time);
        else
            problem.calculateNonLinResidualVec("reverse", time);
        if (criterion=="Residual")
            residual = problem.calculateResidualNorm();
        
        if (nlIts==0)
            residual0 = residual;
        
        if (criterion=="Residual"){
            criterionValue = residual/residual0;
//            exporterTxt->exportData( criterionValue );
            if (verbose)
                cout << "### Newton iteration : " << nlIts << "  relative nonlinear residual : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }

        // Systems are combined in timeProblem.assemble("Newton") and then combined
        problem.assemble("Newton"); 

        problem.setBoundariesSystem();

        problem.getSystem()->writeMM("Assembled");


        if (timestepping == "External"){//AceGen
            gmresIts += problem.solveAndUpdate( "ResidualAceGen", criterionValue );
        //    exporterTxt->exportData( criterionValue );

            //problem.assembleExternal( "OnlyUpdate" );// update AceGEN internal variables
        }
        else
            gmresIts += problem.solveAndUpdate( criterion, criterionValue );
        
        nlIts++;

        //problem.getSolution()->getBlock(0)->print();
        if(criterion=="Update"){
            if (verbose)
                cout << "### Newton iteration : " << nlIts << "  residual of update : " << criterionValue << endl;
            if ( criterionValue < tol )
                break;
        }

        // ####### end FPI #######
    }

    gmresIts/=nlIts;
    if (verbose)
        cout << "### Total Newton iteration : " << nlIts << "  with average gmres its : " << gmresIts << endl;
    if ( problem.getParameterList()->sublist("Parameter").get("Cancel MaxNonLinIts",false) ) {
        TEUCHOS_TEST_FOR_EXCEPTION(nlIts == maxNonLinIts ,std::runtime_error,"Maximum nonlinear Iterations reached. Problem might have converged in the last step. Still we cancel here.");
    }
    if (!valuesForExport.is_null()) {
        if (valuesForExport->size() == 2){
            (*valuesForExport)[0] = gmresIts;
            (*valuesForExport)[1] = nlIts;
        }
       

    }
    
}


template<class SC,class LO,class GO,class NO>
void NonLinearSolver<SC,LO,GO,NO>::solveExtrapolation(TimeProblem<SC,LO,GO,NO> &problem, double time){

    bool verbose = problem.getVerbose();

    problem.assemble("Extrapolation");

    problem.setBoundaries(time); // Setting boundaries to system rhs. The rest of the rhs (e.g. M*u_t) must/should be implemented in DAESolver

    int	gmresIts = problem.solve( );

    if (verbose) {
        cout << "### GMRES Its : " << gmresIts << endl;
    }
}
}
#endif
