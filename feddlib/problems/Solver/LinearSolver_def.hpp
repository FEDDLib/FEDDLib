#ifndef LinearSolver_DEF_hpp
#define LinearSolver_DEF_hpp
#include "LinearSolver_decl.hpp"
/*!
 Definition of LinearSolver

 @brief  LinearSolver
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {


template<class SC,class LO,class GO,class NO>
LinearSolver<SC,LO,GO,NO>::LinearSolver(){}


template<class SC,class LO,class GO,class NO>
LinearSolver<SC,LO,GO,NO>::~LinearSolver(){}

template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solve(Problem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string type ){

    int its=0;
    if (!type.compare("Monolithic") || !type.compare("MonolithicConstPrec"))
        its = solveMonolithic( problem, rhs, type );
    else if (!type.compare("Teko")){
#ifdef FEDD_HAVE_TEKO
//        its = solveTeko( problem, rhs );
        its = solveBlock( problem, rhs, "Teko" );
#else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Teko not found! Build Trilinos with Teko.");
#endif
    }
    else if (type=="Diagonal" || type=="Triangular"){
        its = solveBlock( problem, rhs, type );
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown solver type.");

    return its;
}

template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solve(TimeProblem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string type ){

    int its=0;

    if (!type.compare("Monolithic"))
        its = solveMonolithic( problem, rhs );
    else if (!type.compare("Teko")){
#ifdef FEDD_HAVE_TEKO
//        its = solveTeko( problem, rhs );
        its = solveBlock( problem, rhs, "Teko" );
#else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Teko not found! Build Trilinos with Teko.");
#endif
    }
    else if(!type.compare("FaCSI") || type == "FaCSI-Teko" )
        its = solveBlock( problem, rhs, type );
    else if (type=="Diagonal" || type=="Triangular")
        its = solveBlock( problem, rhs, type );
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown solver type.");
    
    return its;
}

template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solveMonolithic(Problem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string type ){

    bool verbose(problem->getVerbose());
    int its=0;
    if (problem->getParameterList()->get("Zero Initial Guess",true)) {
        problem->getSolution()->putScalar(0.);
    }
    Teuchos::RCP<Thyra::MultiVectorBase<SC> >thyraX = problem->getSolution()->getThyraMultiVector();

    Teuchos::RCP<const Thyra::MultiVectorBase<SC> > thyraB;
    if (rhs.is_null())
        thyraB = problem->getRhs()->getThyraMultiVectorConst();
    else
        thyraB = rhs->getThyraMultiVectorConst();

    ParameterListPtr_Type pListThyraSolver = sublist( problem->getParameterList(), "ThyraSolver" );

    pListThyraSolver->setParameters( problem->getParameterList()->sublist("ThyraPreconditioner") );

    problem->getLinearSolverBuilder()->setParameterList(pListThyraSolver);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory = problem->getLinearSolverBuilder()->createLinearSolveStrategy("");

    bool iterativeSolve = !pListThyraSolver->get("Linear Solver Type", "Belos").compare("Belos");
    if (iterativeSolve && (type != "MonolithicConstPrec" || problem->getPreconditioner()->getThyraPrec().is_null()))
        problem->setupPreconditioner("Monolithic");

    if (!pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive").compare("Multiplicative")) {
        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);

        Teuchos::RCP<const Thyra::LinearOpBase<SC> > thyra_linOp = problem->getPreconditioner()->getThyraPrec()->getUnspecifiedPrecOp();
        Thyra::apply( *thyra_linOp, Thyra::NOTRANS, *thyraB, thyraX.ptr() );
        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",false);
    }

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > solver = lowsFactory->createOp();
//    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > solver = linearOpWithSolve(*lowsFactory, problem->getSystem()->getThyraLinOp());
    ThyraLinOpConstPtr_Type thyraMatrix = problem->getSystem()->getThyraLinOp();
    if ( iterativeSolve ) {
        ThyraPrecPtr_Type thyraPrec = problem->getPreconditioner()->getThyraPrec();
        Thyra::initializePreconditionedOp<SC>(*lowsFactory, thyraMatrix, thyraPrec.getConst(), solver.ptr());
    }
    else{
        Thyra::initializeOp<SC>(*lowsFactory, thyraMatrix, solver.ptr());
    }

    {
        Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*solver, Thyra::NOTRANS, *thyraB, thyraX.ptr());
        if (verbose)
            std::cout << status << std::endl;
        if ( !pListThyraSolver->get("Linear Solver Type","Belos").compare("Belos") )
            its = status.extraParameters->get("Belos/Iteration Count",0);
        else
            its = 0;

        problem->getSolution()->fromThyraMultiVector(thyraX);
    }

    return its;
}

template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solveMonolithic(TimeProblem_Type* timeProblem, BlockMultiVectorPtr_Type rhs){


    bool verbose(timeProblem->getVerbose());
    int its=0;
    // timeProblem->getSystem()->getBlock(0,0)->writeMM("System");
    ProblemPtr_Type problem = timeProblem->getUnderlyingProblem();

    if (problem->getParameterList()->get("Zero Initial Guess",true)) {
        problem->getSolution()->putScalar(0.);
    }

    ParameterListPtr_Type pListThyraSolver = sublist( problem->getParameterList(), "ThyraSolver" );

    Teuchos::RCP<Thyra::MultiVectorBase<SC> >thyraX = problem->getSolution()->getThyraMultiVector();
    Teuchos::RCP<const Thyra::MultiVectorBase<SC> > thyraB;

    if (rhs.is_null())
        thyraB = problem->getRhs()->getThyraMultiVectorConst();
    else
        thyraB = rhs->getThyraMultiVectorConst();
    
    pListThyraSolver->setParameters( problem->getParameterList()->sublist("ThyraPreconditioner") );

    problem->getLinearSolverBuilder()->setParameterList(pListThyraSolver);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory = problem->getLinearSolverBuilder()->createLinearSolveStrategy("");


    problem->setupPreconditioner( "Monolithic" );

    if (!pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive").compare("Multiplicative")) {
        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);

        Teuchos::RCP<const Thyra::LinearOpBase<SC> > thyra_linOp = problem->getPreconditioner()->getThyraPrec()->getUnspecifiedPrecOp();
        Thyra::apply( *thyra_linOp, Thyra::NOTRANS, *thyraB, thyraX.ptr() );
        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",false);
    }

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > solver = lowsFactory->createOp();
    //    solver = linearOpWithSolve(*lowsFactory, problem->getSystem()->getThyraLinOp());

    //timeProblem->combineSystems();
    // timeProblem->getSystemCombined()->getBlock(0,0)->writeMM("SystemCombined");

    ThyraLinOpConstPtr_Type thyraMatrix = timeProblem->getSystemCombined()->getThyraLinOp();

    // Printing the stiffness matrix for the first newton iteration
    // timeProblem->getSystemCombined()->writeMM("stiffnessMatrixWihtDirichlet");
    // timeProblem->getSystem()->writeMM("stiffnessMatrixFull");
    // rhs->writeMM("rhs");

    if ( !pListThyraSolver->get("Linear Solver Type","Belos").compare("Belos") ) {
        ThyraPrecPtr_Type thyraPrec = problem->getPreconditioner()->getThyraPrec();
        Thyra::initializePreconditionedOp<SC>(*lowsFactory, thyraMatrix, thyraPrec.getConst(), solver.ptr());
    }
    else{
        Thyra::initializeOp<SC>(*lowsFactory, thyraMatrix, solver.ptr());
    }
    {
        Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*solver, Thyra::NOTRANS, *thyraB, thyraX.ptr());
        if (verbose)
            std::cout << status << std::endl;
        problem->getSolution()->fromThyraMultiVector(thyraX);
        // problem->getSolution()->writeMM("solution");
        if ( !pListThyraSolver->get("Linear Solver Type","Belos").compare("Belos") ){
            its = status.extraParameters->get("Belos/Iteration Count",0);
            double achievedTol = status.extraParameters->get("Belos/Achieved Tolerance",-1.);
        }
        else
            its = 0;
    }
    return its;
}

#ifdef FEDD_HAVE_TEKO
template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solveTeko(Problem_Type* problem, BlockMultiVectorPtr_Type rhs ){

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    bool verbose(problem->getVerbose());
    int its=0;
    if (problem->getParameterList()->get("Zero Initial Guess",true)) {
        problem->getSolution()->putScalar(0.);
    }
    // typedef Teuchos::RCP< Thyra::ProductMultiVectorBase< double > > 	Teko::BlockedMultiVector
    // convert them to teko compatible sub vectors
    Teko::MultiVector x0_th = problem->getSolution()->getBlockNonConst(0)->getThyraMultiVector();
    Teko::MultiVector x1_th = problem->getSolution()->getBlockNonConst(1)->getThyraMultiVector();

    Teko::MultiVector b0_th;
    Teko::MultiVector b1_th;
    if (rhs.is_null()){
        b0_th = problem->getRhs()->getBlockNonConst(0)->getThyraMultiVector();
        b1_th = problem->getRhs()->getBlockNonConst(1)->getThyraMultiVector();
    }
    else{
        b0_th = rhs->getBlockNonConst(0)->getThyraMultiVector();
        b1_th = rhs->getBlockNonConst(1)->getThyraMultiVector();
    }

    std::vector<Teko::MultiVector> x_vec; x_vec.push_back(x0_th); x_vec.push_back(x1_th);
    std::vector<Teko::MultiVector> b_vec; b_vec.push_back(b0_th); b_vec.push_back(b1_th);

    Teko::MultiVector x = Teko::buildBlockedMultiVector(x_vec); // these will be used in the Teko solve
    Teko::MultiVector b = Teko::buildBlockedMultiVector(b_vec);

    ParameterListPtr_Type pListThyraSolver = sublist( problem->getParameterList(), "ThyraSolver" );

//    pListThyraSolver->setParameters( problem->getParameterList()->sublist("ThyraPreconditioner") );

    problem->getLinearSolverBuilder()->setParameterList(pListThyraSolver);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory = problem->getLinearSolverBuilder()->createLinearSolveStrategy("");

    problem->setupPreconditioner( "Teko" );

    //    if (!pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive").compare("Multiplicative")) {
    //        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);
    //
    //        Teuchos::RCP<const Thyra::LinearOpBase<SC> > thyra_linOp = problem->getPreconditioner()->getThyraPrec()->getUnspecifiedPrecOp();
    //        Thyra::apply( *thyra_linOp, Thyra::NOTRANS, *thyraB, thyraX.ptr() );
    //        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",false);
    //    }


    lowsFactory->setOStream(out);
//    lowsFactory->setVerbLevel(Teuchos::VERB_EXTREME);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > solver = lowsFactory->createOp();

    ThyraPrecPtr_Type thyraPrec = problem->getPreconditioner()->getThyraPrec();

    ThyraLinOpConstPtr_Type thyraMatrix = problem->getPreconditioner()->getTekoOp();

    Thyra::initializePreconditionedOp<SC>(*lowsFactory, thyraMatrix, thyraPrec.getConst(), solver.ptr());

    {
        Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*solver, Thyra::NOTRANS, *b, x.ptr());
        if (verbose)
            std::cout << status << std::endl;

        its = status.extraParameters->get("Belos/Iteration Count",0);
        Teuchos::RCP< Thyra::ProductMultiVectorBase< double > > xBlocked = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase< double > >( x );
        problem->getSolution()->getBlockNonConst(0)->fromThyraMultiVector( xBlocked->getNonconstMultiVectorBlock( 0 ) );
        problem->getSolution()->getBlockNonConst(1)->fromThyraMultiVector( xBlocked->getNonconstMultiVectorBlock( 1 ) );
    }

    return its;
}

template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solveTeko(TimeProblem_Type* timeProblem, BlockMultiVectorPtr_Type rhs ){

    bool verbose(timeProblem->getVerbose());
    int its=0;

    ProblemPtr_Type problem = timeProblem->getUnderlyingProblem();
    if (problem->getParameterList()->get("Zero Initial Guess",true)) {
        problem->getSolution()->putScalar(0.);
    }
    // typedef Teuchos::RCP< Thyra::ProductMultiVectorBase< double > > 	Teko::BlockedMultiVector
    // convert them to teko compatible sub vectors
    Teko::MultiVector x0_th = problem->getSolution()->getBlockNonConst(0)->getThyraMultiVector();
    Teko::MultiVector x1_th = problem->getSolution()->getBlockNonConst(1)->getThyraMultiVector();

    Teko::MultiVector b0_th;
    Teko::MultiVector b1_th;
    if (rhs.is_null()){
        b0_th = problem->getRhs()->getBlockNonConst(0)->getThyraMultiVector();
        b1_th = problem->getRhs()->getBlockNonConst(1)->getThyraMultiVector();
    }
    else{
        b0_th = rhs->getBlockNonConst(0)->getThyraMultiVector();
        b1_th = rhs->getBlockNonConst(1)->getThyraMultiVector();
    }

    std::vector<Teko::MultiVector> x_vec; x_vec.push_back(x0_th); x_vec.push_back(x1_th);
    std::vector<Teko::MultiVector> b_vec; b_vec.push_back(b0_th); b_vec.push_back(b1_th);

    Teko::MultiVector x = Teko::buildBlockedMultiVector(x_vec); // these will be used in the Teko solve
    Teko::MultiVector b = Teko::buildBlockedMultiVector(b_vec);

    ParameterListPtr_Type pListThyraSolver = sublist( problem->getParameterList(), "ThyraSolver" );

//    pListThyraSolver->setParameters( problem->getParameterList()->sublist("ThyraPreconditioner") );

    problem->getLinearSolverBuilder()->setParameterList(pListThyraSolver);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory = problem->getLinearSolverBuilder()->createLinearSolveStrategy("");

    problem->setupPreconditioner( "Teko" );

//    if (!pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive").compare("Multiplicative")) {
//        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);
//
//        Teuchos::RCP<const Thyra::LinearOpBase<SC> > thyra_linOp = problem->getPreconditioner()->getThyraPrec()->getUnspecifiedPrecOp();
//        Thyra::apply( *thyra_linOp, Thyra::NOTRANS, *thyraB, thyraX.ptr() );
//        pListThyraSolver->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",false);
//    }

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);

    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > solver = lowsFactory->createOp();

    ThyraPrecPtr_Type thyraPrec = problem->getPreconditioner()->getThyraPrec();

    ThyraLinOpConstPtr_Type thyraMatrix = problem->getPreconditioner()->getTekoOp();

    Thyra::initializePreconditionedOp<SC>(*lowsFactory, thyraMatrix, thyraPrec.getConst(), solver.ptr());

    {
        Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*solver, Thyra::NOTRANS, *b, x.ptr());
        if (verbose)
            std::cout << status << std::endl;

        its = status.extraParameters->get("Belos/Iteration Count",0);
        Teuchos::RCP< Thyra::ProductMultiVectorBase< double > > xBlocked = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase< double > >( x );
        problem->getSolution()->getBlockNonConst(0)->fromThyraMultiVector( xBlocked->getNonconstMultiVectorBlock( 0 ) );
        problem->getSolution()->getBlockNonConst(1)->fromThyraMultiVector( xBlocked->getNonconstMultiVectorBlock( 1 ) );
    }

    return its;
}
#endif
 
template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solveBlock(Problem_Type* problem, BlockMultiVectorPtr_Type rhs, std::string type ){
    typedef Thyra::DefaultZeroLinearOp<SC> ZeroOp_Type;
    typedef Teuchos::RCP<ZeroOp_Type> ZeroOpPtr_Type;

    int rank = problem->getComm()->getRank();
    bool verbose(problem->getVerbose());
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    int its=0;
    
    if (problem->getParameterList()->get("Zero Initial Guess",true)) {
        problem->getSolution()->putScalar(0.);
    }
    
    
    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraX = problem->getSolution()->getProdThyraMultiVector();
    
    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraRHS;
    if ( rhs.is_null() )
        thyraRHS = problem->getRhs()->getProdThyraMultiVector();
    else
        thyraRHS = rhs->getProdThyraMultiVector();
    
    ParameterListPtr_Type pListThyraSolver = sublist( problem->getParameterList(), "ThyraSolver" );
    
    problem->getLinearSolverBuilder()->setParameterList(pListThyraSolver);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory = problem->getLinearSolverBuilder()->createLinearSolveStrategy("");
    
    problem->setupPreconditioner( type );
    
    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);
    
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > solver = lowsFactory->createOp();
    
    ThyraPrecPtr_Type thyraPrec = problem->getPreconditioner()->getThyraPrec();

//    BlockMatrixPtr_Type system = problem->getSystem();
//    ThyraLinOpBlockConstPtr_Type thyraMatrixBlocks = system->getThyraLinBlockOp();
//    ThyraLinOpBlockPtr_Type thyraMatrixBlocksNonConst = Teuchos::rcp_const_cast<ThyraLinOpBlock_Type>( thyraMatrixBlocks );
//    
//    for (int i=0; i<system->size(); i++) {
//        for (int j=0; j<system->size(); j++) {
//            if ( !thyraMatrixBlocks->blockExists(i,j) ) {
//                ZeroOpPtr_Type zeroOp =
//                    Teuchos::rcp( new ZeroOp_Type( thyraRHS->getMultiVectorBlock(i)->range(),
//                                                   thyraX->getMultiVectorBlock(j)->range() ) ); //(range , domain)
//                thyraMatrixBlocksNonConst->getNonconstBlock(i,j) = zeroOp;
//            }
//        }
//    }
    
    
    
//    BlockMatrixPtr_Type system = problem->getSystem();
//    for (int i=0; i<system->size(); i++) {
//        for (int j=0; j<system->size(); j++) {
//            if (!system->blockExists(1,1)){
//                MatrixPtr_Type dummy = Teuchos::rcp( new Matrix_Type( system->getBlock(i,j)->getMap(), 0 ) );
//                dummy->fillComplete( problem->getSolution()->getBlock( j )->getMap(), problem->getRhs()->getBlock( i )->getMap() ); //(domain, range)
//                system->addBlock( dummy, i, j );
//            }
//        }
//    }
    
    ThyraLinOpConstPtr_Type thyraMatrix = problem->getSystem()->getThyraLinBlockOp();
    Thyra::initializePreconditionedOp<SC>(*lowsFactory, thyraMatrix, thyraPrec.getConst(), solver.ptr());
    {
        
        Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*solver, Thyra::NOTRANS, *thyraRHS, thyraX.ptr());
        if (verbose)
            std::cout << status << std::endl;
        
        its = status.extraParameters->get("Belos/Iteration Count",0);
        
        problem->getSolution()->fromThyraProdMultiVector( thyraX );
    }
    
    return its;
}
    
template<class SC,class LO,class GO,class NO>
int LinearSolver<SC,LO,GO,NO>::solveBlock(TimeProblem_Type* timeProblem, BlockMultiVectorPtr_Type rhs, std::string type ){
    typedef Thyra::DefaultZeroLinearOp<SC> ZeroOp_Type;
    typedef Teuchos::RCP<ZeroOp_Type> ZeroOpPtr_Type;

    int rank = timeProblem->getComm()->getRank();
    bool verbose(timeProblem->getVerbose());
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    int its=0;
    
    ProblemPtr_Type problem = timeProblem->getUnderlyingProblem();
    if (problem->getParameterList()->get("Zero Initial Guess",true)) {
        problem->getSolution()->putScalar(0.);
    }

    
    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraX = problem->getSolution()->getProdThyraMultiVector();
    
    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraRHS;
    if ( rhs.is_null() )
        thyraRHS = problem->getRhs()->getProdThyraMultiVector();
    else
        thyraRHS = rhs->getProdThyraMultiVector();
    
    ParameterListPtr_Type pListThyraSolver = sublist( problem->getParameterList(), "ThyraSolver" );

    
    problem->getLinearSolverBuilder()->setParameterList(pListThyraSolver);
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory = problem->getLinearSolverBuilder()->createLinearSolveStrategy("");
    
    problem->setupPreconditioner( type );

    lowsFactory->setOStream(out);
    lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);
    
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > solver = lowsFactory->createOp();
    
    ThyraPrecPtr_Type thyraPrec = problem->getPreconditioner()->getThyraPrec();

//    BlockMatrixPtr_Type system = timeProblem->getSystemCombined();
//    ThyraLinOpBlockConstPtr_Type thyraMatrixBlocks = system->getThyraLinBlockOp();
//    ThyraLinOpBlockPtr_Type thyraMatrixBlocksNonConst = Teuchos::rcp_const_cast<ThyraLinOpBlock_Type>( thyraMatrixBlocks );
//    
//    for (int i=0; i<system->size(); i++) {
//        for (int j=0; j<system->size(); j++) {
//            if ( !thyraMatrixBlocks->blockExists(i,j) ) {
//                ZeroOpPtr_Type zeroOp =
//                Teuchos::rcp( new ZeroOp_Type( thyraRHS->getMultiVectorBlock(i)->range(),
//                                              thyraX->getMultiVectorBlock(j)->range() ) ); //(range , domain)
//                thyraMatrixBlocksNonConst->getNonconstBlock(i,j) = zeroOp;
//            }
//        }
//    }

    BlockMatrixPtr_Type system = timeProblem->getSystemCombined();
//    for (int i=0; i<system->size(); i++) {
//        for (int j=0; j<system->size(); j++) {
//            if (!system->blockExists(1,1)){
//                MatrixPtr_Type dummy;
//                dummy.reset( new Matrix_Type( system->getBlock(i,j)->getMap(), 0 ) );
//                dummy->fillComplete( problem->getSolution()->getBlock( j )->getMap(), problem->getRhs()->getBlock( i )->getMap() ); //(domain, range)
//                system->addBlock( dummy, i, j );
//            }
//        }
//    }
    
    ThyraLinOpConstPtr_Type thyraMatrix = timeProblem->getSystemCombined()->getThyraLinBlockOp();
//    ThyraLinOpBlockConstPtr_Type thyraMatrixBlock = timeProblem->getSystemCombined()->getThyraLinBlockOp();
    Thyra::initializePreconditionedOp<SC>(*lowsFactory, thyraMatrix, thyraPrec.getConst(), solver.ptr());
    {
        
        Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*solver, Thyra::NOTRANS, *thyraRHS, thyraX.ptr());
        if (verbose)
            std::cout << status << std::endl;
        
        its = status.extraParameters->get("Belos/Iteration Count",0);
        
        problem->getSolution()->fromThyraProdMultiVector( thyraX );
    }
    
    return its;
}
}
#endif
