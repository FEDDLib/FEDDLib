#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

#include "feddlib/problems/specific/Stokes.hpp"

#include <Teuchos_TestForException.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 main of Stokes problem
 
 @brief Stokes main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


void one(double* x, double* res, double t, const double* parameters){
    
    res[0] = 1.;
    res[1] = 1.;
    
    return;
}
void two(double* x, double* res, double t, const double* parameters){
    
    res[0] = 2.;
    res[1] = 2.;
    
    return;
}
void three(double* x, double* res, double t, const double* parameters){
    
    res[0] = 3.;
    res[1] = 3.;
    
    return;
}

// For Lid Driven Cavity Test
void ldcFunc2D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 1.;
    res[1] = 0.;
    
    return;
}

// For Lid Driven Cavity Test
void ldcFunc3D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 1.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void zeroDirichlet(double* x, double* res, double t, const double* parameters){
    
    res[0] = 0.;
    
    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    
    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters){
    
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}

void inflowParabolic2D(double* x, double* res, double t, const double* parameters){
    
    double H = parameters[1];
    res[0] = 4*parameters[0]*x[1]*(H-x[1])/(H*H);
    res[1] = 0.;
    
    return;
}

void inflowParabolic3D(double* x, double* res, double t, const double* parameters){
    
    double H = parameters[1];
    res[0] = 16*parameters[0]*x[1]*(H-x[1])*x[2]*(H-x[2])/(H*H*H*H);
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters){
    
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
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
    
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    bool verbose (comm->getRank() == 0);
    if (verbose) {
        cout << "#################################################" <<endl;
        cout << "################### Stokes ######################" <<endl;
        cout << "#################################################" <<endl;
    }

    
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlTekoPrecFile = "parametersTeko.xml";
    myCLP.setOption("tekoprecfile",&xmlTekoPrecFile,".xml file with Inputparameters.");
    string xmlBlockPrecFile = "parametersPrecBlock.xml";
    myCLP.setOption("blockprecfile",&xmlBlockPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    double length = 4.;
    myCLP.setOption("length",&length,"length of domain.");
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        MPI_Finalize();
        return 0;
    }
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);

        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);
        
        ParameterListPtr_Type parameterListPrecTeko = Teuchos::getParametersFromXmlFile(xmlTekoPrecFile);
        
        ParameterListPtr_Type parameterListPrecBlock = Teuchos::getParametersFromXmlFile(xmlBlockPrecFile);
        
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",3);
        
        std::string discVelocity = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        std::string discPressure = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");

        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","some_mesh_here");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int 		volumeID        = parameterListProblem->sublist("Parameter").get("Volume ID",0);
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);        
        int         n;
        int			inflowID        = parameterListProblem->sublist("Parameter").get("Inflow ID",2);
        int			inflowBCID      = parameterListProblem->sublist("Parameter").get("InflowBC ID",-1);
        string      bcType          = parameterListProblem->sublist("Parameter").get("BC Type","parabolic");
        string      precMethod      = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem));
        if (precMethod == "Monolithic")
            parameterListAll->setParameters(*parameterListPrec);
        else if(precMethod == "Teko")
            parameterListAll->setParameters(*parameterListPrecTeko);
        else if(precMethod == "Diagonal" || precMethod == "Triangular")
            parameterListAll->setParameters(*parameterListPrecBlock);
        
        parameterListAll->setParameters(*parameterListSolver);
        
        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("structured_ldc") || !meshType.compare("unstructured_struct")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_bfs") || !meshType.compare("unstructured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }
        
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;
        int numProcsProblem = size;

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));
        DomainPtr_Type domainPressure;
        DomainPtr_Type domainVelocity;
        {
            Teuchos::TimeMonitor totalTimeMonitor(*totalTime);
            {
                Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
                if (verbose) {
                    cout << "-- Building Mesh ..." << flush;
                }
                
                if (!meshType.compare("structured")) {
                    TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
                    if (dim == 2) {
                        n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(2);
                        x[0]=0.0;    x[1]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                    }
                    else if (dim == 3){
                        n = (int) (std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(3);
                        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                    }
                    domainPressure->buildMesh( 1,"Square", dim, discPressure, n, m, numProcsCoarseSolve);
                    domainVelocity->buildMesh( 1,"Square", dim, discVelocity, n, m, numProcsCoarseSolve);
                }
                if (!meshType.compare("structured_ldc")) {
                    // Structured Mesh for Lid-Driven Cavity Test
                    TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
                    if (dim == 2) {
                        n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(2);
                        x[0]=0.0;    x[1]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                    }
                    else if (dim == 3){
                        n = (int) (std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(3);
                        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                    }
                    domainPressure->buildMesh( 5,"Square", dim, discPressure, n, m, numProcsCoarseSolve);
                    domainVelocity->buildMesh( 5,"Square", dim, discVelocity, n, m, numProcsCoarseSolve);
                }
                if (!meshType.compare("structured_bfs")) {
                    TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured BFS mesh.");
                    if (dim == 2) {
                        n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(2);
                        x[0]=-1.0;    x[1]=-1.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, length+1., 2., comm ) );
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, length+1., 2., comm ) );
                    }
                    else if (dim == 3){
                        n = (int) (std::pow( size/minNumberSubdomains ,1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(3);
                        x[0]=-1.0;    x[1]=0.0;    x[2]=-1.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, length+1., 1., 2., comm));
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, length+1., 1., 2., comm));
                    }
                    domainPressure->buildMesh( 2,"BFS", dim, discPressure, n, m, numProcsCoarseSolve);
                    domainVelocity->buildMesh( 2,"BFS", dim, discVelocity, n, m, numProcsCoarseSolve);
                }
                else if (!meshType.compare("unstructured")) {
                    domainPressure.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                    domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                    
                    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
                    domainP1Array[0] = domainPressure;
                    
                    ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
                    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                    
                    partitionerP1.readAndPartition();
                    
                    domainVelocity->buildP2ofP1Domain( domainPressure );
                    
                    if (discVelocity=="P2")
                        domainVelocity->buildP2ofP1Domain( domainPressure );
                    else
                        domainVelocity = domainPressure;

                }
            }
            //domainVelocity->exportNodeFlags();
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("MaxVelocity",1.));
            
            // ####################
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
            
            if (!bcType.compare("parabolic"))
                parameter_vec.push_back(1.);//height of inflow region
            else if(!bcType.compare("parabolic_benchmark"))
                parameter_vec.push_back(.41);//height of inflow region
            else if(!bcType.compare("LDC")) // Lid Driven Cavity Test
                parameter_vec.push_back(0.);//Dummy
            else
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Select a valid boundary condition.");
            
            if (dim==2 && bcType.compare("LDC")) { // When it is not LDC ..
                bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
            }
            else if(dim==3 && bcType.compare("LDC")){
                bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
            }

            if (!bcType.compare("parabolic_benchmark")) {//flag of obstacle
                if (dim==2)
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                else if (dim==3)
                    bcFactory->addBC(zeroDirichlet3D, 4, 0, domainVelocity, "Dirichlet", dim);
            }
            if (!bcType.compare("LDC")) {
                if (dim==2){
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(ldcFunc2D, 2, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(ldcFunc3D, 2, 0, domainVelocity, "Dirichlet", dim);

                }
            }
            
            Stokes<SC,LO,GO,NO> stokes( domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll );

            domainVelocity->info();
            domainPressure->info();
            stokes.info();
            
            {
                Teuchos::TimeMonitor solveTimeMonitor(*solveTime);
                
                stokes.addBoundaries(bcFactory);
                
                stokes.initializeProblem();
                
                stokes.assemble();

                stokes.setBoundaries();
                                
                stokes.solve();

            }

            if ( parameterListAll->sublist("General").get("ParaViewExport",false) ) {
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaVelocity(new ExporterParaView<SC,LO,GO,NO>());
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaPressure(new ExporterParaView<SC,LO,GO,NO>());
                
                Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionV = stokes.getSolution()->getBlock(0);
                Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionP = stokes.getSolution()->getBlock(1);

                DomainPtr_Type dom = domainVelocity;
                
                exParaVelocity->setup("velocity", dom->getMesh(), dom->getFEType());
                
                UN dofsPerNode = dim;
                exParaVelocity->addVariable(exportSolutionV, "u", "Vector", dofsPerNode, dom->getMapUnique());
                
                dom = domainPressure;
                exParaPressure->setup("pressure", dom->getMesh(), dom->getFEType());
                
                if (dom->getFEType()=="P0")
                    exParaPressure->addVariable(exportSolutionP, "p", "Scalar", 1, dom->getElementMap());
                else
                    exParaPressure->addVariable(exportSolutionP, "p", "Scalar", 1, dom->getMapUnique());

                exParaVelocity->save(0.0);
                exParaPressure->save(0.0);
                
                exParaVelocity->closeExporter();
                exParaPressure->closeExporter();
                
            }
        }
    }
    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
