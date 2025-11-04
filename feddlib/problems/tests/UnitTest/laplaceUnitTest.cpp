#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/BCBuilder.hpp"



/*!
 main of Laplace unit test

 @brief Laplace main

 This laplace unit test compares the current solution of the laplace problem in 2D and 3D to older already stored solutions. The folling test constructs a structured mesh with H/h = 6.
 The stored solutions make sense for the following configurations:
 - 2D on 4 procs with P1 or P2 elements
 - 3D on 8 procs with P1 or P2 elements
 
 */

void zeroBC(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
}
void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    // std::string xmlProblemFile = "parametersProblem.xml";
    // myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    std::string xmlPrecFile = "parametersPrec_Laplace.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    std::string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    int boolExportSolution = 0;
    myCLP.setOption("exportSolution", &boolExportSolution, "Export Solution");

    int dim = 2;
    myCLP.setOption("dim", &dim, "Dimension");

    std::string FEType = "P2";
    myCLP.setOption("FEType", &FEType, "Discretization");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    {
       // ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListPrec)) ;
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int m = 6;
        std::string meshType = "structured";
        std::string meshDelimiter = " ";

        int n;
        int size = comm->getSize();
        int numProcsCoarseSolve = 0;
        size -= numProcsCoarseSolve;

        int minNumberSubdomains = 1;
      
        Teuchos::RCP<Domain<SC,LO,GO,NO> > domain;
        if (dim == 2) {
            TEUCHOS_TEST_FOR_EXCEPTION( std::floor(std::sqrt(size)) != std::sqrt(size) , std::logic_error, "Wrong number of processors for structured squared mesh.");

            n = (int) (std::pow(size,1/2.) + 100.*Teuchos::ScalarTraits<double>::eps()); // 1/H

            std::vector<double> x(2);
            x[0]=0.0;    x[1]=0.0;
            domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., comm) ) ;
            domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        }
        else if (dim == 3) {
            TEUCHOS_TEST_FOR_EXCEPTION( std::floor(std::pow(size,1/3)) != std::pow(size,1/3) , std::logic_error, "Wrong number of processors for structured squared mesh.");

            n = (int) (std::pow(size,1/3.) + 100.*Teuchos::ScalarTraits< SC >::eps()); // 1/H
            std::vector<double> x(3);
            x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
            domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., 1., comm) ) ;
            domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
        }
        // ####################
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory(new BCBuilder<SC,LO,GO,NO>( ));

        bcFactory->addBC(zeroBC, 1, 0, domain, "Dirichlet", 1);
        bcFactory->addBC(zeroBC, 2, 0, domain, "Dirichlet", 1);
        bcFactory->addBC(zeroBC, 3, 0, domain, "Dirichlet", 1);
        

        Laplace<SC,LO,GO,NO> laplace(domain,FEType,parameterListAll,false);
        {
            laplace.addRhsFunction(oneFunc);
            laplace.addBoundaries(bcFactory);
            
            laplace.initializeProblem();
            laplace.assemble();
            laplace.setBoundaries();
            laplace.solve();
        }

        

        //  HDF5Export<SC,LO,GO,NO> exporter(laplace.getSolution()->getBlock(0)->getMap(), "ReferenceSolutions/solution_laplace_"+std::to_string(dim)+"d_"+FEType+"_"+std::to_string(size)+"cores"); //  Map and file name
        //  exporter.writeVariablesHDF5("solution",laplace.getSolution()->getBlock(0)); // VariableName and Variable

        // We exclude any other tests, than the one prescribed
        if(dim==2){
            TEUCHOS_TEST_FOR_EXCEPTION( !(size == 4 && m==6), std::logic_error, "The 2D test solutions are only sensible for 4 processors.");
        }
        else if(dim==3)
            TEUCHOS_TEST_FOR_EXCEPTION( !(size == 8 && m==6), std::logic_error, "The 3D test solutions are only sensible for 8 processors.");

        HDF5Import<SC,LO,GO,NO> importer(laplace.getSolution()->getBlock(0)->getMap(),"ReferenceSolutions/solution_laplace_"+std::to_string(dim)+"d_"+FEType+"_"+std::to_string(size)+"cores");
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > solutionImported = importer.readVariablesHDF5("solution");
    
        // We compare the imported solution to the current one
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > solutionLaplace = laplace.getSolution()->getBlock(0);
        // Calculating the error per node
        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > errorValues = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>( laplace.getSolution()->getBlock(0)->getMap() ) ); 
        //this = alpha*A + beta*B + gamma*this
        errorValues->update( 1., solutionLaplace, -1. ,solutionImported, 0.);
        // Computing norm
        Teuchos::Array<SC> norm(1); 
        errorValues->normInf(norm);
        double normError = norm[0];
                    
        // Output of error
        if(comm->getRank() ==0){
            std::cout << " --------------------------------------------------" << std::endl;
            std::cout << "  Error Report " << std::endl;
            std::cout << "   || solution_current - solution_stored||_inf = " << normError << std::endl;
            std::cout << " --------------------------------------------------" << std::endl;
        }
        // Throwing exception, if error is too great.

        TEUCHOS_TEST_FOR_EXCEPTION( normError > 1.e-11 , std::logic_error, "Difference between current solution and stored solution greater than 1e-11.");
        if (boolExportSolution==1) {
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolution = laplace.getSolution()->getBlock(0);

            exPara->setup("solutionLaplace", domain->getMesh(), FEType);
            
            exPara->addVariable(exportSolution, "u_current", "Scalar", 1, domain->getMapUnique());
            exPara->addVariable(solutionImported, "u_imported", "Scalar", 1, domain->getMapUnique());

            exPara->save(0.0);

        }
    }
    return(EXIT_SUCCESS);
}
