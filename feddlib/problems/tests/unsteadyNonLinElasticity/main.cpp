#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

void rhsY2D(double* x, double* res, double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    if (parameters[0]<=parameters[2])
        res[1] = parameters[1];
    
    std::cout << "res[1]:" << res[1] <<" parameters[0]:"<<parameters[0] << " parameters[2]:" << parameters[2]<< std::endl;
    return;
}

void rhsY(double* x, double* res, double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    if (parameters[0]<=parameters[2])
        res[1] = parameters[1];
    res[2] = 0.;
    return;
}

void rhsX2D(double* x, double* res, double* parameters){
    
    res[0] = 0.;
    if (parameters[0]<=parameters[2])
        res[0] = parameters[1];
    res[1] = 0.;

    return;
}

void rhsX(double* x, double* res, double* parameters){
    
    res[0] = 0.;
    if (parameters[0]<=parameters[2])
        res[0] = parameters[1];

    res[1] = 0.;
    res[2] = 0.;
    return;
}

void zeroDirichlet(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
}

void rhsYZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    double force = parameters[1];
    double TRamp = 10.0;
    if(parameters[0] <= TRamp+1e-06)
        force = parameters[0] * force * 1./(TRamp);
    else
        force = parameters[1];

    if(parameters[2] == 5)
        res[1] = force;
    else
        res[1] =0.;
        
    if (parameters[2] == 4)
        res[2] = force;
    else
        res[2] = 0.;
    
    return;
}


typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;
int main(int argc, char *argv[])
{

    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    // int dim = 2;
    // myCLP.setOption("dim",&dim,"dim");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    bool verbose (comm->getRank() == 0); // Print-Ausgaben nur auf rank = 0
    if (verbose)
    {
        cout << "###############################################################" <<endl;
        cout << "########## Starting unsteady nonlinear elasticity ... #########" <<endl;
        cout << "###############################################################" <<endl;
    }
    
    
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);
        
        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);
        
        // Mesh
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",3);
        int m = parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization","P1");
        std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        std::string bcType = parameterListProblem->sublist("Parameter").get("BC Type","volumeY");
        std::string bcPlace = parameterListProblem->sublist("Parameter").get("BC Placement","standard");
        
        int n;
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;
        
        int minNumberSubdomains = 1;
                
        Teuchos::RCP<Domain<SC,LO,GO,NO> > domain;
        if (!meshType.compare("structured")) {
            if (dim == 2) {
                n = (int) (std::pow(size,1/2.) + 100.*Teuchos::ScalarTraits<double>::eps()); // 1/H
                std::vector<double> x(2);
                x[0]=0.0;    x[1]=0.0;
                domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., comm) ) ;
                domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
            }
            else if (dim == 3) {
                n = (int) (std::pow(size,1/3.) + 100.*Teuchos::ScalarTraits<double>::eps()); // 1/H
                std::vector<double> x(3);
                x[0]=0.0;    x[1]=0.0; x[2]=0.0;
                domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., 1., comm) ) ;
                domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
            }
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
        if(dim==3){
			bcFactory->addBC(zeroDirichlet, 1, 0, domain, "Dirichlet_X", dim);
		    bcFactory->addBC(zeroDirichlet, 2, 0, domain, "Dirichlet_Y", dim);
		    bcFactory->addBC(zeroDirichlet, 3, 0, domain, "Dirichlet_Z", dim);
		    bcFactory->addBC(zeroDirichlet3D, 0, 0, domain, "Dirichlet", dim);
		    bcFactory->addBC(zeroDirichlet2D, 7, 0, domain, "Dirichlet_X_Y", dim);
		    bcFactory->addBC(zeroDirichlet2D, 8, 0, domain, "Dirichlet_Y_Z", dim);
		    bcFactory->addBC(zeroDirichlet2D, 9, 0, domain, "Dirichlet_X_Z", dim);
        }
        else{
        
        }
        
    
        {

            NonLinElasticity<SC,LO,GO,NO> nonLinElas(domain,FEType,parameterListAll);

            domain->info();
            nonLinElas.info();
           
            nonLinElas.addRhsFunction( rhsYZ );
            
            double force = parameterListAll->sublist("Parameter").get("Volume force",0.);
            double finalTimeRamp = parameterListAll->sublist("Timestepping Parameter").get("Final time force",0.1);
            double degree = 0;
            
            nonLinElas.addParemeterRhs( force );
            nonLinElas.addParemeterRhs( finalTimeRamp );
            nonLinElas.addParemeterRhs( degree );
            
            nonLinElas.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

            nonLinElas.initializeProblem();
            // initial assembly of system and vectors
            nonLinElas.assemble();

            // ######################
            // Zeitintegration
            // ######################
            DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

            // Only one block for structural problem
            SmallMatrix<int> defTS(1);
            defTS[0][0] = 1;

            // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
            // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
            daeTimeSolver.defineTimeStepping(defTS);

            // Uebergebe das (nicht) lineare Problem
            daeTimeSolver.setProblem(nonLinElas);

            // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
            // defTS definiert worden sind.
            daeTimeSolver.setupTimeStepping();

            // Fuehre die komplette Zeitintegration + Newton + Loesen + Exporter durch
            daeTimeSolver.advanceInTime();

        }
    }
    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
