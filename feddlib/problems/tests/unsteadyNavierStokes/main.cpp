#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"

#include <Xpetra_DefaultPlatform.hpp>

/*!
 main of time-dependent Navier-Stokes problem
 
 @brief time-dependent Navier-Stokes main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;

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

void inflowPartialCFD(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];

    if(t < 0.5)
    {
        res[0] = (4.0*1.5*parameters[0]*x[1]*(H-x[1])/(H*H))*((1 - cos(2.0*M_PI*t))/2.0);
        res[1] = 0.;
    }
    else
    {
        res[0] = 4.0*1.5*parameters[0]*x[1]*(H-x[1])/(H*H);
        res[1] = 0.;
    }

    return;
}

void inflowParabolic2D(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = 4.*parameters[0]*x[1]*(H-x[1])/(H*H);
    res[1] = 0.;
    
    return;
}

void inflowParabolic2DSin(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = sin(M_PI*t*0.125)*( 6*x[1]*(H-x[1]) ) / (H*H);
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

void inflow3DRichter(double* x, double* res, double t, const double* parameters)
{
    double H = parameters[1];
    
    if(t < 1.)
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) ) * ((1 - cos(2.0*M_PI*t))/2.0);
        res[1] = 0.;
        res[2] = 0.;
    }
    else
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) );
        res[1] = 0.;
        res[2] = 0.;
    }
    
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
        cout << "###############################################################" <<endl;
        cout << "################### Unsteady Navier-Stokes ####################" <<endl;
        cout << "###############################################################" <<endl;
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
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",3);
        std::string feTypeV = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        std::string feTypeP = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");
        string		meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","structured");
        string		meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int 		m = parameterListProblem->sublist("Parameter").get("H/h",5);
        string		linearization = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
        string		precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        bool computeInflow = parameterListProblem->sublist("Parameter").get("Compute Inflow",false);
        int         n;


        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if (!precMethod.compare("Monolithic"))
            parameterListAll->setParameters(*parameterListPrec);
        else
            parameterListAll->setParameters(*parameterListPrecTeko);

        parameterListAll->setParameters(*parameterListSolver);

        std::string bcType = parameterListProblem->sublist("Parameter").get("BC Type","parabolic");

        int minNumberSubdomains;
        if (!meshType.compare("structured")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_rec")){
            minNumberSubdomains = length;
        }
        else if(!meshType.compare("structured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

        double viscosity = parameterListProblem->sublist("Parameter").get("Viscosity",1.e-3);

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));
        
        {
            DomainPtr_Type domainPressure;
            DomainPtr_Type domainVelocity;

            Teuchos::TimeMonitor totalTimeMonitor(*totalTime);
            {
                Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
                if (verbose)
                    cout << "-- Building Mesh ..." << flush;
                
                if (!meshType.compare("structured")) {
                    if (dim == 2) {
                        n = (int) (std::pow(size,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(2);
                        x[0]=0.0;    x[1]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
                    }
                    else if (dim == 3){
                        n = (int) (std::pow(size,1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(3);
                        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
                        domainPressure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                        domainVelocity.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
                    }
                    domainPressure->buildMesh( 1,"Square", dim, feTypeP, n, m, numProcsCoarseSolve);
                    domainVelocity->buildMesh( 1,"Square", dim, feTypeV, n, m, numProcsCoarseSolve);
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
                    
                    if (feTypeV=="P2")
                        domainVelocity->buildP2ofP1Domain( domainPressure );
                    else
                        domainVelocity = domainPressure;
                    
                }
            }

            
//            ofstream myFile;
//            FILE * pFile;
//            std::cout << "FluidP2Elements..." << '\n';
//            myFile.open ("FluidP2ElementsH00.txt");
//            for(int i = 0; i < domainVelocity->getElements()->size(); i++)
//            {
//                for(int j= 0; j < domainVelocity->getElements()->at(i).size(); j++)
//                {
//                    myFile << domainVelocity->getElements()->at(i).at(j);
//                    myFile << " ";
//                }
//                myFile << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//
//            std::cout << "FluidP1Elements..." << '\n';
//            myFile.open ("FluidP1ElementsH00.txt");
//            for(int i = 0; i < domainPressure->getElements()->size(); i++)
//            {
//                for(int j= 0; j < domainPressure->getElements()->at(i).size(); j++)
//                {
//                    myFile << domainPressure->getElements()->at(i).at(j);
//                    myFile << " ";
//                }
//                myFile << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//
//            std::cout << "FluidP2Nodes..." << '\n';
////            myFile.open ("FluidP2NodesH00.txt");
//            pFile = fopen ("FluidP2NodesH00.txt","w");
//
//            for(int i = 0; i < domainVelocity->getPointsUnique()->size(); i++)
//            {
//                for(int j= 0; j < domainVelocity->getPointsUnique()->at(i).size(); j++)
//                {
//                    fprintf(pFile,"%4.10f ", domainVelocity->getPointsUnique()->at(i).at(j) );
//
////                    printf("%4.10f ", domainFluidVelocity->getPointsUnique()->at(i).at(j) );
////                    myFile << domainFluidVelocity->getPointsUnique()->at(i).at(j);
////                    myFile << " ";
//                }
//                fprintf(pFile,"\n");
////                myFile << endl;
//            }
//            fclose(pFile);
////            myFile.close();
//            std::cout << "done" << '\n';
//
//            std::cout << "FluidP1Nodes..." << '\n';
////            myFile.open ("FluidP1NodesH00.txt");
//            pFile = fopen ("FluidP1NodesH00.txt","w");
//
//            for(int i = 0; i < domainPressure->getPointsUnique()->size(); i++)
//            {
//                for(int j= 0; j < domainPressure->getPointsUnique()->at(i).size(); j++)
//                {
//                    fprintf(pFile,"%4.10f ", domainPressure->getPointsUnique()->at(i).at(j) );
////                    myFile << domainFluidPressure->getPointsUnique()->at(i).at(j);
////                    myFile << " ";
//                }
//                fprintf(pFile,"\n");
////                myFile << endl;
//            }
//            fclose(pFile);
////            myFile.close();
//            std::cout << "done" << '\n';
//            
//
//            std::cout << "FluidP2Flags..." << '\n';
//            myFile.open ("FluidP2FlagsH00.txt");
//            for(int i = 0; i < domainVelocity->getBCFlagUnique()->size(); i++)
//            {
//                myFile << domainVelocity->getBCFlagUnique()->at(i) << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//
//            std::cout << "FluidPFlags..." << '\n';
//            myFile.open ("FluidP1FlagsH00.txt");
//            for(int i = 0; i < domainPressure->getBCFlagUnique()->size(); i++)
//            {
//                myFile << domainPressure->getBCFlagUnique()->at(i) << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
            
            std::vector<double> parameter_vec(1);
            if ( !bcType.compare("parabolic") || !bcType.compare("parabolic_benchmark") || !bcType.compare("parabolic_benchmark_sin") )
                parameter_vec[0] = parameterListProblem->sublist("Parameter").get("MaxVelocity",1.5);
            else if ( !bcType.compare("partialCFD") ) //  Fuer CFD3
                parameter_vec[0] = parameterListProblem->sublist("Parameter").get("MeanVelocity",2.);
                        
            // ####################
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

            if (!bcType.compare("parabolic"))
                parameter_vec.push_back(1.);//height of inflow region
            else if(!bcType.compare("parabolic_benchmark_sin") || !bcType.compare("parabolic_benchmark") || !bcType.compare("partialCFD"))
                parameter_vec.push_back(.41);//height of inflow region
            else if(!bcType.compare("Richter3D"))
                parameter_vec.push_back(.4);
            else
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Select a valid boundary condition.");

            if (!bcType.compare("parabolic") || !bcType.compare("parabolic_benchmark")) {//flag of obstacle
                if (dim==2){
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet3D, 4, 0, domainVelocity, "Dirichlet", dim);
                    
                }
            }
            else if (!bcType.compare("parabolic_benchmark_sin")) {
                if (dim==2){
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic2DSin, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test for 2D only: parabolic_benchmark_sin");
                }
            }
            else if (!bcType.compare("partialCFD")) { // Fuer CFD3 Test
                if (dim==2)
                {
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowPartialCFD, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(zeroDirichlet2D, 5, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No partial CFD test in 3D.");
                }
            }
            else if (!bcType.compare("Richter3D")) {
                bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim); // wall
                bcFactory->addBC(inflow3DRichter, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // inflow
                bcFactory->addBC(zeroDirichlet3D, 3, 0, domainVelocity, "Dirichlet_Z", dim);
                bcFactory->addBC(zeroDirichlet3D, 5, 0, domainVelocity, "Dirichlet", dim);

//                bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);

            }
            
            int timeDisc = parameterListProblem->sublist("Timestepping Parameter").get("Butcher table",0);

            NavierStokes<SC,LO,GO,NO> navierStokes( domainVelocity, feTypeV, domainPressure, feTypeP, parameterListAll );

            navierStokes.addBoundaries(bcFactory);
            
            navierStokes.initializeProblem();
            
            navierStokes.assemble();

            navierStokes.setBoundariesRHS();

            DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);
            SmallMatrix<int> defTS(2);
            defTS[0][0] = 1;
            defTS[0][1] = 1;
            defTS[1][0] = 0;
            defTS[1][1] = 0;

            daeTimeSolver.defineTimeStepping(defTS);

            daeTimeSolver.setProblem(navierStokes);

            daeTimeSolver.setupTimeStepping();

            daeTimeSolver.advanceInTime();

        }
    }

    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
