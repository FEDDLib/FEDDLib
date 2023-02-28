#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/specific/TPM.hpp"

#include <Teuchos_TestForException.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 main of unsteadyTPM problem
 
 @brief unsteadyTPM main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

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


void lineLoad2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;

    if (parameters[2] == 4. /*line at the top*/){
        if (parameters[0]<parameters[1]/*max time ramp*/)
            res[1] = - 1. * parameters[0] / parameters[1];
        else
            res[1] = - 1.;
    }
    
    return;
}

void lineLoad2DMini(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;

    if (parameters[2] == 3./*line at the right*/){
        if (parameters[0]<parameters[1]/*max time ramp*/){
            res[0] = - 10. * parameters[0] / parameters[1];
        }
        else
            res[0] = - 10.;
    }
    return;
}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;
using namespace std;

int main(int argc, char *argv[]) {
    
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
    
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    bool verbose (comm->getRank() == 0);
    if (verbose) {
        cout << "##############################################" <<endl;
        cout << "################ unsteady TPM ################" <<endl;
        cout << "##############################################" <<endl;
    }

    
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    string xmlTekoPrecFile = "parametersTeko.xml";
    myCLP.setOption("tekoprecfile",&xmlTekoPrecFile,".xml file with Inputparameters.");
    
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
        
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        
        std::string discVelocity = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        std::string discPressure = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");

        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","structured");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        int         n;
        string      bcType          = parameterListProblem->sublist("Parameter").get("BC Type","standard");
        string      precMethod      = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        
        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if (!precMethod.compare("Monolithic"))
            parameterListAll->setParameters(*parameterListPrec);
        
        parameterListAll->setParameters(*parameterListSolver);
        
        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("unstructured_struct") || !meshType.compare("structuredMiniTest")) {
            minNumberSubdomains = 1;
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
                    TEUCHOS_TEST_FOR_EXCEPTION( true , std::logic_error, "Flags and surface elements must be adjusted for tpm.");

                    domainPressure->buildMesh( 1, "Square", dim, discPressure, n, m, numProcsCoarseSolve);
                    domainVelocity->buildMesh( 1, "Square", dim, discVelocity, n, m, numProcsCoarseSolve);
                }
                else if (!meshType.compare("unstructured")) {
                    domainPressure.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                    domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                    
                    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
                    domainP1Array[0] = domainPressure;
                    
                    ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
                    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                    
                    partitionerP1.readAndPartition();
                    
                    domainVelocity->buildP2ofP1Domain( domainPressure );
                    if (dim==2) {
                        domainPressure->getMesh()->setElementFlags("TPM_square");
                        domainVelocity->getMesh()->setElementFlags("TPM_square");
                    }
                    else if(dim==3){
//                        domainPressure->getMesh()->setElementFlags("TPM_cube");
//                        domainVelocity->getMesh()->setElementFlags("TPM_cube");
                    }                
                }
                else if (!meshType.compare("structuredMiniTest")) {
                    TEUCHOS_TEST_FOR_EXCEPTION( size > 1 , std::logic_error, "We can only use this small test on a single core.");
                    if (dim == 2) {
                        n = (int) (std::pow( size/minNumberSubdomains ,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(2);
                        x[0]=0.0;    x[1]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, .2, .1, comm ) );
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, .2, .1, comm ) );
                    }
                    else if (dim == 3){
                        TEUCHOS_TEST_FOR_EXCEPTION( true , std::logic_error, "Only in 2d.");
                    }
                    domainPressure->buildMesh( 4, "structuredMiniTest", dim, discPressure, n, m, numProcsCoarseSolve);
                    domainVelocity->buildMesh( 4, "structuredMiniTest", dim, discVelocity, n, m, numProcsCoarseSolve);

                }
            }
                    
            
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("MaxVelocity",1.));
            
            // ####################
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
            
            if (dim==2 && meshType == "unstructured") {
                /*
                  --4--6--5--
                 |          |
                 |          |
                 2          3
                 |          |
                 |          |
                  ----1-----
                 */
                
                
                bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(zeroDirichlet2D, 2, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(zeroDirichlet2D, 3, 0, domainVelocity, "Dirichlet", dim);
//                bcFactory->addBC(dummyFunc, 4, 0, domainVelocity, "Neumann", dim);
                bcFactory->addBC(zeroDirichlet, 5, 1, domainPressure, "Dirichlet", 1);//part of top boundary
//                bcFactory->addBC(zeroDirichlet, 1, 1, domainPressure, "Dirichlet", 1);//part of top boundary
//                bcFactory->addBC(zeroDirichlet, 2, 1, domainPressure, "Dirichlet", 1);//part of top boundary
//                bcFactory->addBC(zeroDirichlet, 3, 1, domainPressure, "Dirichlet", 1);//part of top boundary
//                bcFactory->addBC(zeroDirichlet, 10, 1, domainPressure, "Dirichlet", 1);//part of top boundary
            }
            else if(dim==2 && meshType == "structuredMiniTest"){
                std::cout << "using structuredMiniTest BC" << std::endl;
                bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(zeroDirichlet2D, 2, 0, domainVelocity, "Dirichlet_Y", dim);
                bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet_Y", dim);
                bcFactory->addBC(zeroDirichlet2D, 3, 0, domainVelocity, "Dirichlet_X", dim);
//                bcFactory->addBC(zeroDirichlet, 1, 1, domainPressure, "Dirichlet", 1);
                bcFactory->addBC(zeroDirichlet, 2, 1, domainPressure, "Dirichlet", 1);
            }
            
            else if(dim==3){
                
                /*
                 --4--6--5--
                |          |
                |          |
                1          1
                |          |
                |          |
                 ----1-----
                */
                
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "No BC implemented in 3D yet.");
            }
            
            TPM<SC,LO,GO,NO> tpm( domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll );
            
            domainVelocity->info();
            domainPressure->info();

            {
                Teuchos::TimeMonitor solveTimeMonitor(*solveTime);
                
                if(meshType == "structuredMiniTest")
                    tpm.addRhsFunction( lineLoad2DMini );
                else
                    tpm.addRhsFunction( lineLoad2D );
                
                double finalTimeRamp = parameterListAll->sublist("Timestepping Parameter").get("Final time ramp",1.);
                double degree = 0;
                
                tpm.addParemeterRhs( finalTimeRamp );
                tpm.addParemeterRhs( degree );
                
                tpm.addBoundaries(bcFactory);
                
                tpm.initializeProblem();
                
                tpm.assemble("FirstAssemble");
                // ######################
                // Zeitintegration
                // ######################
                DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);
                
                // Only one block for structural problem
                SmallMatrix<int> defTS(2);//doesnt matter, we use the fully assemble TPM problem from AceGen Code
                
                // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
                // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
                daeTimeSolver.defineTimeStepping(defTS);
                
                // Uebergebe das (nicht) lineare Problem
                daeTimeSolver.setProblem(tpm);
                
                // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
                // defTS definiert worden sind.
                daeTimeSolver.setupTimeStepping();
                
                // Fuehre die komplette Zeitintegration + Newton + Loesen + Exporter durch
                daeTimeSolver.advanceInTime();

            }
        }
    }
    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
