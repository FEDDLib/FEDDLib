#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/LaplaceBlocks.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>

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

void rhs2D(double* x, double* res, double* parameters){
    
    res[0] = 1.;
    
    return;
}


typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlPrecFile1 = "parametersPrec1.xml";
    myCLP.setOption("precfile1",&xmlPrecFile1,".xml file with Inputparameters.");
    string xmlPrecFile2 = "parametersPrec2.xml";
    myCLP.setOption("precfile2",&xmlPrecFile2,".xml file with Inputparameters.");

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
    
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListPrec1 = Teuchos::getParametersFromXmlFile(xmlPrecFile1);
        ParameterListPtr_Type parameterListPrec2 = Teuchos::getParametersFromXmlFile(xmlPrecFile2);
        
        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",2);
        TEUCHOS_TEST_FOR_EXCEPTION(dim!=2, std::runtime_error, "Only 2D.");
        int m = parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization","P1");
        std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        std::string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","");
        std::string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        TEUCHOS_TEST_FOR_EXCEPTION(meshType!="unstructured", std::runtime_error, "Only unstrutured.");
        std::string FEType1 = FEType;
        std::string FEType2 = FEType;
        
        int n;
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;

        Teuchos::RCP<Domain<SC,LO,GO,NO> > domain1;
        Teuchos::RCP<Domain<SC,LO,GO,NO> > domain2;
       
        if (!meshType.compare("unstructured")) {
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domain1P1;
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domain2P1;
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domain1P2;
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domain2P2;
            domain1P1.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
            domain2P1.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
            
            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
            domainP1Array[0] = domain1P1;
            domainP1Array[1] = domain2P1;
            ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
            MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
            
            partitionerP1.readAndPartition();

            if (FEType=="P2") {
                domain1P2.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
                domain1P2->buildP2ofP1Domain( domain1P1 );
                domain1 = domain1P2;
                domain2P2.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
                domain2P2->buildP2ofP1Domain( domain2P1 );
                domain2 = domain2P2;
            }
            else{
                domain1 = domain1P1;
                domain2 = domain2P1;
            }
        }



        // ####################
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory(new BCBuilder<SC,LO,GO,NO>( ));

        bcFactory->addBC(zeroBC, 1, 0, domain1, "Dirichlet", 1);
        bcFactory->addBC(zeroBC, 2, 0, domain1, "Dirichlet", 1);
        bcFactory->addBC(zeroBC, 3, 0, domain1, "Dirichlet", 1);

        bcFactory->addBC(zeroBC, 1, 1, domain2, "Dirichlet", 1);
        bcFactory->addBC(zeroBC, 2, 1, domain2, "Dirichlet", 1);
        bcFactory->addBC(zeroBC, 3, 1, domain2, "Dirichlet", 1);


        LaplaceBlocks<SC,LO,GO,NO> laplace(domain1,domain2,FEType1,FEType2,parameterListAll);
        
        laplace.setPrecLists(parameterListPrec1,parameterListPrec2);
        {
            laplace.addBoundaries(bcFactory);

            laplace.addRhsFunction( rhs2D );
            
            laplace.initializeProblem();
            
            laplace.assemble();
            
            laplace.addToRhs( laplace.getSourceTerm() );
            
            laplace.setBoundaries();
            
            laplace.solve();
        }

        bool boolExportSolution = true;
        if (boolExportSolution) {
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domain = domain1;
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolution = laplace.getSolution()->getBlock(0);

            exPara->setup( domain->getDimension(), domain->getNumElementsGlobal(), domain->getElements(), domain->getPointsUnique(), domain->getMapUnique(), domain->getMapRepeated(), FEType, "solutionLaplace1", 1, comm );


            exPara->addVariable(exportSolution, "u", "Scalar", 1, domain->getMapUnique());

            exPara->save(0.0);

        }
    }
    return(EXIT_SUCCESS);
}
