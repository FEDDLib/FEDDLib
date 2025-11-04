#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/core/General/BCBuilder.hpp"


/*!
 main of linear elasticity unit test

 @brief Linear Elasticity main

 This laplace unit test compares the current solution of the linear elasticity problem in
 2D and 3D to older already stored solutions. The folling test constructs a
 structured mesh with H/h = 6. The stored solutions make sense for the following
 configurations:
 - 2D on 4 procs with P1 or P2 elements
 - 3D on 8 procs with P1 or P2 elements

 */

void zeroBC(double *x, double *res, double t, const double *parameters) { res[0] = 0.; }
void oneFunc(double *x, double *res, double *parameters) { res[0] = 1.; }

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

// Using a volume force in x direction
void rhsX2D(double* x, double* res, double* parameters){
    res[0] = parameters[1];
    res[1] = 0.;
    return;
}
void rhsX3D(double* x, double* res, double* parameters){
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
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

    int boolExportSolution = 0;
    myCLP.setOption("exportSolution", &boolExportSolution, "Export Solution");

    int dim = 2;
    myCLP.setOption("dim", &dim, "Dimension");

    std::string FEType = "P2";
    myCLP.setOption("FEType", &FEType, "Discretization");

    std::string xmlPrecFile;
    if(dim==2)
        xmlPrecFile = "parametersPrecNonLinElasticity.xml"; // We can use the stokes preconditioner list for now. Maybe we can use a more general name
    else if(dim==3)
        xmlPrecFile = "parametersPrecNonLinElasticity.xml";
    
    myCLP.setOption("precfile", &xmlPrecFile, ".xml file with Inputparameters.");
    std::string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile", &xmlSolverFile, ".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc, argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    {
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListPrec));
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int m = 4;
        std::string meshType = "structured";
        std::string meshDelimiter = " ";

        int n;
        int numProcsCoarseSolve =0;
        int size = comm->getSize();

        int minNumberSubdomains = 1;

        Teuchos::RCP<Domain<SC, LO, GO, NO>> domain;

        if (dim == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(std::floor(std::sqrt(size)) != std::sqrt(size), std::logic_error, "Wrong number of processors for structured squared mesh.");
            n = (int)(std::pow(size / minNumberSubdomains, 1 / 2.) + 100 * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(2);
            x[0] = 0.0;
            x[1] = 0.0;
            domain.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., comm));
        } else if (dim == 3) {
            TEUCHOS_TEST_FOR_EXCEPTION(std::floor(std::pow(size, 1 / 3)) != std::pow(size, 1 / 3), std::logic_error, "Wrong number of processors for structured squared mesh.");
            n = (int)(std::pow(size / minNumberSubdomains, 1 / 3.) + 100 * Teuchos::ScalarTraits<double>::eps()); // 1/H
            std::vector<double> x(3);
            x[0] = 0.0;
            x[1] = 0.0;
            x[2] = 0.0;
            domain.reset(new Domain<SC, LO, GO, NO>(x, 1., 1., 1., comm));
        }
        domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);

        // ####################         
        Teuchos::RCP<BCBuilder<SC, LO, GO, NO>> bcFactory(new BCBuilder<SC, LO, GO, NO>());
        if (dim == 2) { 
            bcFactory->addBC(zeroDirichlet2D, 1, 0, domain, "Dirichlet", dim);
        } else if (dim == 3) {
            bcFactory->addBC(zeroDirichlet3D, 1, 0, domain, "Dirichlet", dim);
        }

        // Setting parameters that are needed for linear elasticity. 
        // !! This avoids that errors occure when default values are changed along the way.
        parameterListAll->sublist("Parameter").set("Density",1.0);
        parameterListAll->sublist("Parameter").set("Poisson Ratio",0.4);
        parameterListAll->sublist("Parameter").set("Mu",3.0e0);
        parameterListAll->sublist("Parameter").set("Use AceGen Interface", false); // We don't want to use AceGen here

        LinElas<SC,LO,GO,NO> linElas( domain, FEType, parameterListAll );

        {       
            linElas.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

            if (dim==2)
                linElas.addRhsFunction( rhsX2D );
            else if(dim==3)
                linElas.addRhsFunction( rhsX3D );

            double force = 0.5;
            double degree = 0;

            linElas.addParemeterRhs( force );
            linElas.addParemeterRhs( degree );
            // ######################
            // Matrix assemblieren, RW setzen und System loesen
            // ######################
            linElas.initializeProblem();
            linElas.assemble();                
            linElas.setBoundaries(); // In der Klasse Problem
            linElas.solve();

        }

        // HDF5Export<SC, LO, GO, NO> exporter(linElas.getSolution()->getBlock(0)->getMap(),
        //                                     "ReferenceSolutions/solution_linElas_" + std::to_string(dim) + "d_" + FEType + "_" + std::to_string(size) + "cores"); //  Map and file name
        // exporter.writeVariablesHDF5("solution",
        //                             linElas.getSolution()->getBlock(0)); // VariableName and Variable

        // // We exclude any other tests, than the one prescribed
        if (dim == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION(!(size == 4 && m == 4), std::logic_error, "The 2D test solutions are only sensible for 4 processors.");
        } else if (dim == 3)
            TEUCHOS_TEST_FOR_EXCEPTION(!(size == 8 && m == 4), std::logic_error, "The 3D test solutions are only sensible for 8 processors.");

        HDF5Import<SC, LO, GO, NO> importer(linElas.getSolution()->getBlock(0)->getMap(),
                                            "ReferenceSolutions/solution_linElas_" + std::to_string(dim) + "d_" + FEType + "_" + std::to_string(size) + "cores");
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionImported = importer.readVariablesHDF5("solution");

        // We compare the imported solution to the current one
        Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> solutionlinElas = linElas.getSolution()->getBlock(0);
        // Calculating the error per node
        Teuchos::RCP<MultiVector<SC, LO, GO, NO>> errorValues = Teuchos::rcp(new MultiVector<SC, LO, GO, NO>(linElas.getSolution()->getBlock(0)->getMap()));
        // this = alpha*A + beta*B + gamma*this
        errorValues->update(1., solutionlinElas, -1., solutionImported, 0.);
        // Computing norm
        Teuchos::Array<SC> norm(1);
        errorValues->normInf(norm);
        double normError = norm[0];

        // Output of error
        if (comm->getRank() == 0) {
            std::cout << " --------------------------------------------------" << std::endl;
            std::cout << "  Error Report " << std::endl;
            std::cout << "   || solution_current - solution_stored||_inf = " << normError << std::endl;
            std::cout << " --------------------------------------------------" << std::endl;
        }
        // Throwing exception, if error is too great.
        TEUCHOS_TEST_FOR_EXCEPTION(normError > 1.e-11, std::logic_error,
                                    "Difference between current solution and "
                                    "stored solution greater than 1e-11.");
        if (boolExportSolution==1) {
            Teuchos::RCP<ExporterParaView<SC, LO, GO, NO>> exPara(new ExporterParaView<SC, LO, GO, NO>());

            Teuchos::RCP<const MultiVector<SC, LO, GO, NO>> exportSolution = linElas.getSolution()->getBlock(0);

            exPara->setup("solutionLinElas", domain->getMesh(), FEType);

            exPara->addVariable(exportSolution, "u_current", "Vector", dim, domain->getMapUnique());
            exPara->addVariable(solutionImported, "u_imported", "Vector", dim, domain->getMapUnique());

            exPara->save(0.0);
        }
    }
    return (EXIT_SUCCESS);
}
