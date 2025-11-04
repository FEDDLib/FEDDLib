// #define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
// #define VERBOSE
//
// #include "feddlib/core/Mesh/Mesh.hpp"
// #include "feddlib/core/Mesh/MyMeshConvert.hpp"
// #include "feddlib/core/General/ExporterParaView.hpp"
// #include "feddlib/problems/concrete/LinElas.hpp"
// // #include "feddlib/core/Solver/NonLinearSolver.hpp"
//
// #include <Teuchos_RCPDecl.hpp>
// #include <Teuchos_RCPBoostSharedPtrConversions.hpp>
// #include <Teuchos_ParameterList.hpp>
// #include <Teuchos_CommandLineProcessor.hpp>
// #include <Teuchos_XMLParameterListHelpers.hpp>
// #include "feddlib/core/Solver/DAESolverInTime.cpp"

#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/core/General/BCBuilder.hpp"


void rhs2D(double* x, double* res, double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    if (parameters[0]<=1.)
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

void rhs3DX(double* x, double* res, double* parameters){
    
    res[2] = 0.;
    res[1] = 0.;
    res[0] = parameters[1];
    
    return;

}

void rhs2DX(double* x, double* res, double* parameters){
    
    res[1] = 0.;
    res[0] = parameters[1];
    
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

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

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
        return EXIT_SUCCESS;
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
                
                ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
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
                //bcFactory->addBC(threeBC, 1, 0, domain, "Dirichlet", 1);
                //bcFactory->addBC(zeroBC, 4, 0, domain, "Dirichlet", 1);
                bcFactory->addBC(zeroDirichlet2D, 4, 0, domain, "Dirichlet", dim);
                //bcFactory->addBC(zeroDirichlet, 1, 0, domain, "Dirichlet_X", dim);
                //bcFactory->addBC(zeroDirichlet, 1, 0, domain, "Dirichlet_X", dim);

            }
            else if(dim == 3)
            {
                bcFactory->addBC(zeroDirichlet, 1, 0, domain, "Dirichlet_X", dim);
                bcFactory->addBC(zeroDirichlet, 2, 0, domain, "Dirichlet_Y", dim);
                bcFactory->addBC(zeroDirichlet, 3, 0, domain, "Dirichlet_Z", dim);
                bcFactory->addBC(zeroDirichlet3D, 0, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(zeroDirichlet2D, 8, 0, domain, "Dirichlet_Y_Z", dim);
                bcFactory->addBC(zeroDirichlet2D, 9, 0, domain, "Dirichlet_X_Z", dim);
                bcFactory->addBC(zeroDirichlet2D, 7, 0, domain, "Dirichlet_X_Y", dim);
                //bcFactory->addBC(dummyFunc, 4, 0, domain, "Dirichlet", dim);
                //bcFactory->addBC(dummyFunc, 5, 0, domain, "Dirichlet", dim);
                bcFactory->addBC(dummyFunc, 6, 0, domain, "Neumann", dim);
            }

        
            LinElas<SC,LO,GO,NO> LinElas(domain,discType,parameterListAll);

        
            if (dim == 2)
                LinElas.addRhsFunction( rhs2DX );
            else if (dim==3)
                LinElas.addRhsFunction( rhs3DX );
            
            double force = parameterListAll->sublist("Parameter").get("Volume force",0.);
            double finalTimeRamp = parameterListAll->sublist("Timestepping Parameter").get("Final time force",0.1);
            double degree = 0;
            
            LinElas.addParemeterRhs( force );
            LinElas.addParemeterRhs( finalTimeRamp );
            LinElas.addParemeterRhs( degree );
            
            // LinElas Objekt erstellen

            // TODO: Keine Ahnung wieso auskommentieren?
            // {
            //     Teuchos::TimeMonitor solveTimeMonitor(*solveTime);

                LinElas.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen
                
                LinElas.initializeProblem();
                // Matrizen assemblieren
                LinElas.assemble();

                // Wahrscheinlich nicht noetig
                // LinElas.SetBoundariesRHS();

                // ######################
                // Zeitintegration
                // ######################
                DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

                // Das ist eigentlich fuer Stokes gedacht fuer ein System der Form
                // |u p| bzw. |K -B^T| gedacht.
                // |u p|      |-B   0|
                // Die 1 gibt dabei an, fuer welche Zeile und Spalte die Zeitintegration (Massematrix)
                // durchgefuehrt werden soll, hier also z.B. fuer die ganze erste Zeile.
                // Bei Stokes wird die Divergenz-Nebenbedingung, aufgrund der nicht-vorhandenen
                // Zeitabhaengigkeit naemlich nicht mit zeitintegriert.
                // SmallMatrix<int> defTS(2);
                // defTS[0][0] = 1;
                // defTS[0][1] = 1;
                // defTS[1][0] = 0;
                // defTS[1][1] = 0;

                // Fuer das Strukturproblem haben wir analog, da nur d_s als Varable vorhanden:
                SmallMatrix<int> defTS(1);
                defTS[0][0] = 1;

                // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
                // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
                daeTimeSolver.defineTimeStepping(defTS);

                // Uebergebe das (nicht) lineare Problem
                daeTimeSolver.setProblem(LinElas);

                // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
                // defTS definiert worden sind.
                daeTimeSolver.setupTimeStepping();

                // Fuehre die komplette Zeitintegration + ggf. Newton/Fixpunkt + Loesen + Exporter durch
                daeTimeSolver.advanceInTime();

            // TODO: Siehe oben
            // }

            // Exporter ist im Timestepping integriert

        }
    }
    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
