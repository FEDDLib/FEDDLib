#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/LinElas.hpp"

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

void rhs2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    
    return;
}

void rhs(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    res[2] = 0.;
    return;
}

void rhsX(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
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
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    bool verbose (comm->getRank() == 0); // Print-Ausgaben nur auf rank = 0
    if (verbose) {
        cout << "###############################################################" <<endl;
        cout << "############ Starting Steady Linear Elasticity ... ############" <<endl;
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
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","dfg_fsi_solid.mesh");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int         n;
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        string      FEType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

        DomainPtr_Type domain;

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


                if (!meshType.compare("structured")) {
                    
                    Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
                    if (dim == 2) {
                        n = (int) (std::pow(size,1/2.) + 100.*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(2);
                        x[0]=0.0;    x[1]=0.0;
                        domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., comm) ) ;
                        domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
                    }
                    if (dim == 3) {
                        n = (int) (std::pow(size,1/3.) + 100.*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(3);
                        x[0]=0.0;    x[1]=0.0; x[2]=0.0;
                        domain = Teuchos::rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., 1., comm) ) ;
                        domain->buildMesh(1, "Square", dim, FEType, n, m, numProcsCoarseSolve);
                    }
                }
                else if (!meshType.compare("unstructured")) {
                    domain.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
                    domainP1Array[0] = domain;
                    
                    ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
                    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                    
                    partitionerP1.readAndPartition();
                    if (FEType=="P2") {
                        Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP2;
                        domainP2.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                        domainP2->buildP2ofP1Domain( domain );
                        domain = domainP2;
                    }
                    
                }
                if(verbose)
                {
                    cout << " done! -- " << endl;
                }
            }

            // ######################
            // Setup fuer die RW
            // ######################
            // Die Null gibt an auf welchem Block die RW gesetzt werden sollen; hier gibt es nur einen
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
            
            if (meshType == "structured") {
                if (dim == 2)
                    bcFactory->addBC(zeroDirichlet2D, 2, 0, domain, "Dirichlet", dim);
                else if (dim == 3)
                    bcFactory->addBC(zeroDirichlet3D, 2, 0, domain, "Dirichlet", dim);
            }
            else if (meshType == "unstructured") {
                if (dim == 2)
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domain, "Dirichlet", dim);
                else if (dim == 3)
                    bcFactory->addBC(zeroDirichlet3D, 1, 0, domain, "Dirichlet", dim);
            }

            // LinElas Objekt erstellen
            LinElas<SC,LO,GO,NO> LinElas( domain, FEType, parameterListAll );

            {
                Teuchos::TimeMonitor solveTimeMonitor(*solveTime);

                LinElas.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

                if (dim==2)
                    LinElas.addRhsFunction( rhs2D );
                else if(dim==3)
                    LinElas.addRhsFunction( rhsX );

                double force = parameterListAll->sublist("Parameter").get("Volume force",0.);
                double degree = 0;
                
                LinElas.addParemeterRhs( force );
                LinElas.addParemeterRhs( degree );
                // ######################
                // Matrix assemblieren, RW setzen und System loesen
                // ######################
                LinElas.initializeProblem();
                LinElas.assemble();                
                LinElas.setBoundaries(); // In der Klasse Problem
                LinElas.solve();

            }

            // ######################
            // Exporter fuer die Loesung
            // ######################

            if ( parameterListAll->sublist("General").get("ParaViewExport",false) ) {

                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup( "displacements", domain->getMesh(), FEType );
                
                MultiVectorConstPtr_Type valuesSolidConst = LinElas.getSolution()->getBlock(0);
                exPara->addVariable( valuesSolidConst, "values", "Vector", dim, domain->getMapUnique());
                
                exPara->save(0.0);
                exPara->closeExporter();
                
            }

            // ######################
            // Mesh-Bewegung testen
            // ######################
            // Setze die aktuelle (nicht-deformierte) Konfiguration als Referenzkonfiguration
            domain->setReferenceConfiguration();

            typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
            typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
            typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

            MultiVectorConstPtr_Type displacementUniqueConst = LinElas.getSolution()->getBlock(0);
            MultiVectorPtr_Type displacementRepeated = rcp( new MultiVector_Type( LinElas.getDomain(0)->getMapVecFieldRepeated() ) );

            displacementRepeated->importFromVector( displacementUniqueConst );
            MultiVectorPtr_Type displacementUnique = rcp_const_cast<MultiVector_Type>(displacementUniqueConst);
            // Verschiebe das Gitter
            domain->moveMesh(displacementUnique, displacementRepeated);

            // Exportiere die berechnete Loesung, jedoch auf dem schon deformierten Gitter (die Loesung ist egal. Wir wollen nur das Gitter haben)
            if( parameterListAll->sublist("General").get("ParaViewExport",false) )
            {
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup( "moved_mesh", domain->getMesh(), FEType );
                
                MultiVectorConstPtr_Type valuesSolidConst = LinElas.getSolution()->getBlock(0);
                exPara->addVariable( valuesSolidConst, "values", "Vector", dim, domain->getMapUnique());
                
                exPara->save(0.0);
                exPara->closeExporter();
                
            }
        }
    }

    Teuchos::TimeMonitor::report(cout);
    return(EXIT_SUCCESS);
}
