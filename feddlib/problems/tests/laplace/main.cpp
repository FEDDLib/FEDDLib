#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>


/*!
 main of Laplace problem

 @brief Laplace main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

void zeroBC(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
}
void oneBC(double* x, double* res, double t, const double* parameters){
    res[0] = 1.;
}
void twoBC(double* x, double* res, double t, const double* parameters){
    res[0] = 2.;
}
void threeBC(double* x, double* res, double t, const double* parameters){
    res[0] = 3.;
}
void zeroBC2D(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
    res[1] = 0.;
}
void zeroBC3D(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
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
using namespace std;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");


    std::string vectorLaplace = "false";
    myCLP.setOption("vectorLaplace",&vectorLaplace,"vectorLaplace");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
    double length = 4.;
    myCLP.setOption("length",&length,"length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    bool vL ( !vectorLaplace.compare("true") );
    
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",2);
        int m = parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization","P1");
        std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        std::string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","");
        std::string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");

        int n;
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;

        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("unstructured_struct")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_bfs") || !meshType.compare("unstructured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }


        Teuchos::RCP<Domain<SC,LO,GO,NO> > domain;
        if (!meshType.compare("structured")) {
            TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
            if (dim == 2) {
                n = (int) (std::pow(size,1/2.) + 100.*Teuchos::ScalarTraits<double>::eps()); // 1/H
                std::vector<double> x(2);
                x[0]=0.0;    x[1]=0.0;
                domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., comm) ) ;
                domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
            }
            else if (dim == 3){
                n = (int) (std::pow(size,1/3.) + 100.*Teuchos::ScalarTraits< SC >::eps()); // 1/H
                std::vector<double> x(3);
                x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
                domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., 1., comm) ) ;
                domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
            }
        }
        if (!meshType.compare("structured_bfs")) {
            TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured BFS mesh.");
            if (dim == 2) {
                n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                std::vector<double> x(2);
                x[0]=-1.0;    x[1]=-1.0;
                domain.reset(new Domain<SC,LO,GO,NO>( x, length+1., 2., comm ) );
            }
            else if (dim == 3){
                n = (int) (std::pow( size/minNumberSubdomains ,1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                std::vector<double> x(3);
                x[0]=-1.0;    x[1]=0.0;    x[2]=-1.0;
                domain.reset(new Domain<SC,LO,GO,NO>( x, length+1., 1., 2., comm));
            }
            domain->buildMesh( 2,"BFS", dim, FEType, n, m, numProcsCoarseSolve);
        }
        else if (!meshType.compare("unstructured")) {
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP1;
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP2;
            domainP1.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
            
            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
            domainP1Array[0] = domainP1;
            
            ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
            MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
            
            partitionerP1.readAndPartition();

            if (FEType=="P2") {
                domainP2.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
                domainP2->buildP2ofP1Domain( domainP1 );
                domain = domainP2;
            }
            else
                domain = domainP1;
        }



        // ####################
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory(new BCBuilder<SC,LO,GO,NO>( ));
        if (vL) {
            if (dim==2) {
                bcFactory->addBC(zeroBC2D, 1, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC2D, 2, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC2D, 3, 0, domain, "Dirichlet", dim);
            }
            else if(dim==3){
                bcFactory->addBC(zeroBC3D, 1, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC3D, 2, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroBC3D, 3, 0, domain, "Dirichlet", dim);
            }
        }
        else{
            bcFactory->addBC(zeroBC, 1, 0, domain, "Dirichlet", 1);
            bcFactory->addBC(zeroBC, 2, 0, domain, "Dirichlet", 1);
            bcFactory->addBC(zeroBC, 3, 0, domain, "Dirichlet", 1);
        }

        Laplace<SC,LO,GO,NO> laplace(domain,FEType,parameterListAll,vL);
        {
            laplace.addRhsFunction(oneFunc);
            laplace.addBoundaries(bcFactory);
            
            laplace.initializeProblem();
            laplace.assemble();
            laplace.setBoundaries();
            laplace.solve();
        }

        bool boolExportSolution = true;
        if (boolExportSolution) {
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolution = laplace.getSolution()->getBlock(0);

            exPara->setup("solutionLaplace", domain->getMesh(), FEType);
            
            if (vL)
                exPara->addVariable(exportSolution, "u", "Vector", dim, domain->getMapUnique(), domain->getMapUniqueP2());
            else
                exPara->addVariable(exportSolution, "u", "Scalar", 1, domain->getMapUnique(), domain->getMapUniqueP2());

            exPara->save(0.0);

        }
    }
    return(EXIT_SUCCESS);
}
