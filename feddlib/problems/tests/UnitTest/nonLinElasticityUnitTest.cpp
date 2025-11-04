#include <Tpetra_Core.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_RCPBoostSharedPtrConversions.hpp>
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/General/HDF5Export.hpp"

#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"



/*!
main of nonlinear elasticity problem

 @brief NonLinElasticity main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;

void rhs2D(double *x, double *res, double *parameters) {
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];

    return;
}

void rhs(double *x, double *res, double *parameters) {
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    res[2] = 0.;
    return;
}

void rhsX(double *x, double *res, double *parameters) {
    // parameters[0] is the time, not needed here
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void zeroBC(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double *x, double *res, double t, const double *parameters) {

    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void oneBCVector(double *x, double *res, double t, const double *parameters) {

    res[0] = 1.;
    res[1] = 1.;
    res[2] = 1.;

    return;
}

void twoBCVector(double *x, double *res, double t, const double *parameters) {

    res[0] = 2.;
    res[1] = 2.;
    res[2] = 2.;

    return;
}

void threeBCVector(double *x, double *res, double t, const double *parameters) {

    res[0] = 3.;
    res[1] = 3.;
    res[2] = 3.;

    return;
}

void x1(double *x, double *res, double t, const double *parameters) {

    res[0] = -11.;

    return;
}
void x2(double *x, double *res, double t, const double *parameters) {

    res[0] = 2.;

    return;
}
void x3(double *x, double *res, double t, const double *parameters) {

    res[0] = 100.;

    return;
}
void x4(double *x, double *res, double t, const double *parameters) {

    res[0] = 4.;

    return;
}
void x5(double *x, double *res, double t, const double *parameters) {

    res[0] = 5.;

    return;
}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {

    typedef MeshPartitioner<SC, LO, GO, NO> MeshPartitioner_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    double length = 4.;
    myCLP.setOption("length", &length, "length of domain.");
    string xmlProblemFile = "parametersProblemNonLinElasticity.xml";
    myCLP.setOption("problemfile", &xmlProblemFile, ".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrecNonLinElasticity.xml";
    myCLP.setOption("precfile", &xmlPrecFile, ".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile", &xmlSolverFile, ".xml file with Inputparameters.");

    int dim = 2;
    myCLP.setOption("dim", &dim, "Dimension");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
    ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
    ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

    // Set the dimension in the parameter list
    // Required because the AssembleFE class reads directly from the parameter file.
    parameterListProblem->sublist("Parameter").set("Dimension", dim);
    parameterListPrec->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").set("DofsPerNode1", dim);

    std::string elasticityType = parameterListProblem->sublist("Parameter").get("Elasticity Type", "linear");
    string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type", "structured");
    string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name", "tetrahedron.mesh");
    string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter", " ");
    int m = parameterListProblem->sublist("Parameter").get("H/h", 5);

    int n;

    ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem));
    parameterListAll->setParameters(*parameterListPrec);
    parameterListAll->setParameters(*parameterListSolver);
    int minNumberSubdomains = (int)length;

    Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
    Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
    Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

    int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse", 0);

    int size = comm->getSize() - numProcsCoarseSolve;

    std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization", "P1");

    int partition = 0;
    bool verbose(comm->getRank() == 0);
    if (verbose) {
        cout << "######################################" << endl;
        cout << "########## Nonlinear Elasticity ######" << endl;
        cout << "######################################" << endl;
    }

    Teuchos::TimeMonitor totalTimeMonitor(*totalTime);
    Teuchos::RCP<Domain<SC, LO, GO, NO>> domain;
    if (!meshType.compare("structured")) {

        Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
        if (dim == 2) {
            n = (int)(std::pow(size, 1 / 2.) + 100. * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(2);
            x[0] = 0.0;
            x[1] = 0.0;
            domain = Teuchos::rcp(new Domain<SC, LO, GO, NO>(x, 1., 1., comm));
            domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        }
        if (dim == 3) {
            n = (int)(std::pow(size, 1 / 3.) + 100. * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(3);
            x[0] = 0.0;
            x[1] = 0.0;
            x[2] = 0.0;
            domain = Teuchos::rcp(new Domain<SC, LO, GO, NO>(x, 1., 1., 1., comm));
            domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        }
    } else if (!meshType.compare("unstructured")) {
        Teuchos::RCP<Domain<SC, LO, GO, NO>> domainP1;
        Teuchos::RCP<Domain<SC, LO, GO, NO>> domainP2;
        domainP1.reset(new Domain<SC, LO, GO, NO>(comm, dim));

        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
        domainP1Array[0] = domainP1;

        ParameterListPtr_Type pListPartitioner = sublist(parameterListAll, "Mesh Partitioner");
        MeshPartitioner<SC, LO, GO, NO> partitionerP1(domainP1Array, pListPartitioner, "P1", dim);

        partitionerP1.readAndPartition();

        if (FEType == "P2") {
            domainP2.reset(new Domain<SC, LO, GO, NO>(comm, dim));
            domainP2->buildP2ofP1Domain(domainP1);
            domain = domainP2;
        } else
            domain = domainP1;
    }

    Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(new BCBuilder<SC, LO, GO, NO>());
    if (meshType == "structured") {
        if (dim == 2)
            bcFactory->addBC(zeroDirichlet2D, 2, 0, domain, "Dirichlet", dim);
        else if (dim == 3)
            bcFactory->addBC(zeroDirichlet3D, 2, 0, domain, "Dirichlet", dim);
    } else if (meshType == "unstructured") {
        if (dim == 2)
            bcFactory->addBC(zeroDirichlet2D, 1, 0, domain, "Dirichlet", dim);
        else if (dim == 3)
            bcFactory->addBC(zeroDirichlet3D, 1, 0, domain, "Dirichlet", dim);
    }

    NonLinElasticity<SC, LO, GO, NO> elasticity(domain, FEType, parameterListAll);

    domain->info();
    elasticity.info();

    Teuchos::TimeMonitor solveTimeMonitor(*solveTime);

    elasticity.addBoundaries(bcFactory);

    if (dim == 2)
        elasticity.addRhsFunction(rhs2D);
    else if (dim == 3)
        elasticity.addRhsFunction(rhsX);

    double force = parameterListAll->sublist("Parameter").get("Volume force", 0.);
    double degree = 0;

    elasticity.addParemeterRhs(force);
    elasticity.addParemeterRhs(degree);
    elasticity.initializeProblem();
    elasticity.assemble();

    elasticity.setBoundaries();

    std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization", "Newton");
    NonLinearSolver<SC, LO, GO, NO> elasticitySolver(nlSolverType);
    elasticitySolver.solve(elasticity);
    comm->barrier();

    // For generating h5 solution files
    // HDF5Export<SC, LO, GO, NO> exporter(elasticity.getSolution()->getBlock(0)->getMap(),
    //                                     "ReferenceSolutions//solution_nonLinElasticity_" + std::to_string(dim) + "d_" + FEType +
    //                                         "_" + std::to_string(size) + "cores");     //  Map and file name
    // exporter.writeVariablesHDF5("solution", elasticity.getSolution()->getBlock(0)); // VariableName and Variable


    HDF5Import<SC, LO, GO, NO> importer(elasticity.getSolution()->getBlock(0)->getMap(),
                                        "ReferenceSolutions/solution_nonLinElasticity_" + std::to_string(dim) + "d_" + FEType + "_" +
                                            std::to_string(size) + "cores");
    Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImported = importer.readVariablesHDF5("solution");

    // We compare the imported solution to the current one
    Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionElasticity = elasticity.getSolution()->getBlock(0);
    // Calculating the error per node
    Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValues =
        Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(elasticity.getSolution()->getBlock(0)->getMap()));
    // this = alpha*A + beta*B + gamma*this
    errorValues->update(1., solutionElasticity, -1., solutionImported, 0.);
    // Computing norm
    Teuchos::Array<SC> norm(1);
    errorValues->normInf(norm);
    double normError = norm[0];

    // Output of error
    if (comm->getRank() == 0) {
        cout << " --------------------------------------------------" << endl;
        cout << "  Error Report " << endl;
        cout << "   || solution_current - solution_stored||_inf = " << normError << endl;
        cout << " --------------------------------------------------" << endl;
    }

    // Throwing exception, if error is too great.
    TEUCHOS_TEST_FOR_EXCEPTION(normError > 1.e-11, std::logic_error,
                               "Difference between current solution and stored solution greater than 1e-11.");


    // if (parameterListAll->sublist("General").get("ParaViewExport", false)) {
    if (true) {
        Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exParaDisp(new ExporterParaView<SC, LO, GO, NO>());
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolutionU = elasticity.getSolution()->getBlock(0);
        exParaDisp->setup("displacement", domain->getMesh(), domain->getFEType());
        UN dofsPerNode = dim;
        exParaDisp->addVariable(exportSolutionU, "u", "Vector", dofsPerNode, domain->getMapUnique());
        exParaDisp->save(0.0);
    }
    Teuchos::TimeMonitor::report(std::cout);

    return (EXIT_SUCCESS);
}
