// #define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
// #define VERBOSE
//
// #include "feddlib/core/Mesh/Mesh.hpp"
// #include "feddlib/core/Mesh/MyMeshConvert.hpp"
// #include "feddlib/core/General/ExporterParaView.hpp"
// #include "feddlib/problems/concrete/LinElas.hpp"
// // #include "feddlib/core/Solver/NonLinearSolver.hpp"
//
// #include "Teuchos_RCPDecl.hpp"
// #include "Teuchos_RCPBoostSharedPtrConversions.hpp"
// #include "Teuchos_ParameterList.hpp"
// #include "Teuchos_CommandLineProcessor.hpp"
// #include "Teuchos_XMLParameterListHelpers.hpp"
// #include "feddlib/core/Solver/DAESolverInTime.cpp"

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/LinElasFirstOrder.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
void rhs2D(double* x, double* res, double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    if (parameters[0]<=0.2)
        res[1] = parameters[1];
    return;
}

void rhs(double* x, double* res, double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    if (parameters[0]<=0.2)
        res[1] = parameters[1];
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
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;

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
        cout << "############ Starting Unsteady Linear Elasticity ... ############" <<endl;
        cout << "###############################################################" <<endl;
    }
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        string      precMethod      = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","dfg_fsi_solid.mesh");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int         n;
        int			zeroDirID       = parameterListProblem->sublist("Parameter").get("Homogeneous Dirichlet Flag",1); // Dirichlet-Flag in .mesh
        string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        std::string bcType = parameterListProblem->sublist("Parameter").get("BC Type","volumeY");
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

        DomainPtr_Type domainP1;
        DomainPtr_Type domainP2;

        // ########################
        // P1 und P2 Gitter bauen
        // ########################
        {
            Teuchos::TimeMonitor totalTimeMonitor(*totalTime);
            {
                Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
                if(verbose)
                {
                    cout << "-- Building Mesh ..." << flush;
                }
                
                // P1-Gitter bauen
                domainP1.reset( new Domain_Type( comm, dim ) );
                
                MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
                domainP1Array[0] = domainP1;
                
                ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
                MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                
                partitionerP1.readAndPartition();
                
                // P2-Giter bauen
                domainP2.reset( new Domain_Type( comm, dim ) );
                domainP2->buildP2ofP1Domain(domainP1);
                
                if(verbose)
                {
                    cout << " done! -- " << endl;
                }
            }

            // ########################
            // P1 oder P2 Gitter waehlen
            // ########################
            DomainPtr_Type domain;
            if(!discType.compare("P2"))
            {
                domain = domainP2;
                if(verbose)
                {
                    std::cout << "P2-mesh was chosen" << '\n';
                }
            }
            else if(!discType.compare("P1"))
            {
                domain = domainP1;
                if(verbose)
                {
                    std::cout << "P1-mesh was chosen" << '\n';
                }
            }

            // ######################
            // Setup fuer die RW
            // ######################
            // Die Null gibt an auf welchem Block die RW gesetzt werden sollen; hier gibt es nur einen
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
            if(dim == 2)
            {
                bcFactory->addBC(zeroDirichlet2D, zeroDirID, 0, domain, "Dirichlet", dim);
            }
            else if(dim == 3)
            {
                bcFactory->addBC(zeroDirichlet3D, zeroDirID, 0, domain, "Dirichlet", dim);
            }

            LinElasFirstOrder<SC,LO,GO,NO> lefo(domain,discType,parameterListAll);

            domain->info();
            lefo.info();
            if (bcType=="volumeY"){
                if (dim == 2)
                    lefo.addRhsFunction( rhs2D );
                else if (dim==3)
                    lefo.addRhsFunction( rhs );
            }
            else
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Unknown boundary function.");

            double force = this->parameterList_->sublist("Parameter").get("Volume force",0.);
            double finalTimeRamp = parameterListAll->sublist("Timestepping Parameter").get("Final time force",0.1);
            double degree = 0;
            
            LinElas.addParemeterRhs( force );
            LinElas.addParemeterRhs( finalTimeRamp );
            LinElas.addParemeterRhs( degree );

            lefo.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen
            lefo.initializeProblem();
            // Matrizen assemblieren
            lefo.assemble();

            // Wahrscheinlich nicht noetig
            // LinElas.SetBoundariesRHS();

            // ######################
            // Zeitintegration
            // ######################
            DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);


            // Fuer das Strukturproblem haben wir analog, da nur d_s als Varable vorhanden:
            SmallMatrix<int> defTS(2);
            defTS[0][0] = 0; defTS[0][1] = 2; defTS[1][0] = 2; defTS[1][1] = 0;
            // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
            // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
            daeTimeSolver.defineTimeStepping(defTS);

            // Uebergebe das (nicht) lineare Problem
            daeTimeSolver.setProblem(lefo);

            // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
            // defTS definiert worden sind.
            daeTimeSolver.setupTimeStepping();

            // Fuehre die komplette Zeitintegration + ggf. Newton/Fixpunkt + Loesen + Exporter durch
            daeTimeSolver.advanceInTime();


        }
    }
    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
