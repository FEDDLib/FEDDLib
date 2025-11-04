#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/General/HDF5Export.hpp"

#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NonLinLaplace.hpp"


void zeroDirichlet(double *x, double *res, double t, const double *parameters) { res[0] = 0.; }

void oneFunc2D(double *x, double *res, double *parameters) {
    res[0] = 1.;
    res[1] = 1.;
}

void oneFunc3D(double *x, double *res, double *parameters) {
    res[0] = 1.;
    res[1] = 1.;
    res[2] = 1.;
}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;
int main(int argc, char *argv[]) {

    typedef MeshUnstructured<SC, LO, GO, NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC, LO, GO, NO> Domain_Type;
    typedef RCP<Domain_Type> DomainPtr_Type;
    typedef ExporterParaView<SC, LO, GO, NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC, LO, GO, NO> MeshPartitioner_Type;

    typedef Map<LO, GO, NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef MultiVector<SC, LO, GO, NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC, LO, GO, NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    
    // ########################
    // Set default values for command line parameters
    // ########################
    Teuchos::CommandLineProcessor myCLP;
    string xmlProblemFile = "parametersProblemNonLinLaplace.xml";
    myCLP.setOption("problemfile", &xmlProblemFile, ".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrecNonLinLaplace.xml";
    myCLP.setOption("precfile", &xmlPrecFile, ".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile", &xmlSolverFile, ".xml file with Inputparameters.");
    double length = 4;
    myCLP.setOption("length", &length, "length of domain.");
    int dim = 2;
    myCLP.setOption("dim", &dim, "Dimension");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    bool verbose(comm->getRank() == 0); // Only first rank prints
    if (verbose) {
        cout << "###############################################################" << endl;
        cout << "############ Starting nonlinear Laplace ... ############" << endl;
        cout << "###############################################################" << endl;
    }

    ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
    ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
    ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

    // Set the dimension in the parameter list
    // Required because the AssembleFE class reads directly from the parameter file.
    parameterListProblem->sublist("Parameter").set("Dimension", dim);
    ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem));
    parameterListAll->setParameters(*parameterListPrec);
    parameterListAll->setParameters(*parameterListSolver);

    string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type", "structured");
    string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name", "");
    string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter", " ");
    int m = parameterListProblem->sublist("Parameter").get("H/h", 5);
    string FEType = parameterListProblem->sublist("Parameter").get("Discretization", "P1");

    int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse", 0);
    int size = comm->getSize() - numProcsCoarseSolve;

    Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
    Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
    Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

    int minNumberSubdomains, n;
    if (!meshType.compare("structured") || !meshType.compare("unstructured_struct")) {
        minNumberSubdomains = 1;
    } else if (!meshType.compare("structured_bfs") || !meshType.compare("unstructured_bfs")) {
        minNumberSubdomains = (int)2 * length + 1;
    }

    // ########################
    // Build mesh
    // ########################
    DomainPtr_Type domain;
    if (!meshType.compare("structured")) {
        TEUCHOS_TEST_FOR_EXCEPTION(size % minNumberSubdomains != 0, std::logic_error,
                                   "Wrong number of processors for structured mesh.");
        if (dim == 2) {
            n = (int)(std::pow(size, 1 / 2.) + 100. * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(2);
            x[0] = 0.0;
            x[1] = 0.0;
            domain = Teuchos::rcp(new Domain<SC, LO, GO, NO>(x, 1., 1., comm));
            domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        } else if (dim == 3) {
            n = (int)(std::pow(size, 1 / 3.) + 100. * Teuchos::ScalarTraits<SC>::eps()); // 1/H
            std::vector<double> x(3);
            x[0] = 0.0;
            x[1] = 0.0;
            x[2] = 0.0;
            domain = Teuchos::rcp(new Domain<SC, LO, GO, NO>(x, 1., 1., 1., comm));
            domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        }
    } else if (!meshType.compare("unstructured")) {
        domain.reset(new Domain<SC, LO, GO, NO>(comm, dim));
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
        domainP1Array[0] = domain;
        ParameterListPtr_Type pListPartitioner = sublist(parameterListAll, "Mesh Partitioner");
        MeshPartitioner<SC, LO, GO, NO> partitionerP1(domainP1Array, pListPartitioner, "P1", dim);
        partitionerP1.readAndPartition();
        domain = domain;
    }

    // ########################
    // Set flags for the boundary conditions
    // ########################

    Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(new BCBuilder<SC, LO, GO, NO>());
    bcFactory->addBC(zeroDirichlet, 1, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(zeroDirichlet, 2, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(zeroDirichlet, 3, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(zeroDirichlet, 4, 0, domain, "Dirichlet", 1);

    NonLinLaplace<SC, LO, GO, NO> nonLinLaplace(domain, FEType, parameterListAll);
    nonLinLaplace.addBoundaries(bcFactory);

    if (dim == 2)
        nonLinLaplace.addRhsFunction(oneFunc2D);
    else if (dim == 3)
        nonLinLaplace.addRhsFunction(oneFunc3D);

    // ######################
    // Assemble matrix, set boundary conditions and solve
    // ######################
    nonLinLaplace.initializeProblem();
    nonLinLaplace.assemble();
    nonLinLaplace.setBoundaries();
    nonLinLaplace.setBoundariesRHS();

    std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization", "NOX");
    NonLinearSolver<SC, LO, GO, NO> nlSolverAssFE(nlSolverType);
    nlSolverAssFE.solve(nonLinLaplace);
    comm->barrier();
    
    // For generating h5 solution files
    // HDF5Export<SC, LO, GO, NO> exporter(nonLinLaplace.getSolution()->getBlock(0)->getMap(),
    //                                     "ReferenceSolutions/solution_nonLinLaplace_" + std::to_string(dim) + "d_" + FEType +
    //                                         "_" + std::to_string(size) + "cores");     //  Map and file name
    // exporter.writeVariablesHDF5("solution", nonLinLaplace.getSolution()->getBlock(0)); // VariableName and Variable


    HDF5Import<SC, LO, GO, NO> importer(nonLinLaplace.getSolution()->getBlock(0)->getMap(),
                                        "ReferenceSolutions/solution_nonLinLaplace_" + std::to_string(dim) + "d_" + FEType + "_" +
                                            std::to_string(size) + "cores");
    Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImported = importer.readVariablesHDF5("solution");

    // We compare the imported solution to the current one
    Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionLaplace = nonLinLaplace.getSolution()->getBlock(0);
    // Calculating the error per node
    Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValues =
        Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(nonLinLaplace.getSolution()->getBlock(0)->getMap()));
    // this = alpha*A + beta*B + gamma*this
    errorValues->update(1., solutionLaplace, -1., solutionImported, 0.);
    // Computing norm
    Teuchos::Array<SC> norm(1);
    errorValues->norm2(norm);
    double normError = norm[0];

    // Output of error
    if (comm->getRank() == 0) {
        cout << " --------------------------------------------------" << endl;
        cout << "  Error Report " << endl;
        cout << "   || solution_current - solution_stored||_2 = " << normError << endl;
        cout << " --------------------------------------------------" << endl;
    }

    // Throwing exception, if error is too great.
    TEUCHOS_TEST_FOR_EXCEPTION(normError > 1.e-11, std::logic_error,
                               "Difference between current solution and stored solution greater than 1e-11.");

      // ########################
    // Export solution
    // ########################
    bool boolExportSolution = false;
    if (boolExportSolution) {
        Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exPara(new ExporterParaView<SC, LO, GO, NO>());
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolution = nonLinLaplace.getSolution()->getBlock(0);
        exPara->setup("solutionNonLinLaplace", domain->getMesh(), FEType);

        exPara->addVariable(exportSolution, "u", "Scalar", 1, domain->getMapUnique());

        exPara->save(0.0);
    }

    return EXIT_SUCCESS;
}
