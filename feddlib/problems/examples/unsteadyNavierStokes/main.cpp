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
#include <Teuchos_StackedTimer.hpp>

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

void couette2D(double* x, double* res, double t, const double* parameters){

    res[0] = 1.*parameters[0];
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

void inflowParabolic2DSinBenchmark(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = sin(M_PI*t*0.125)*( 6*x[1]*(H-x[1]) ) / (H*H);
    res[1] = 0.;

    return;
}


void inflowParabolic3DSinBenchmark(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];
    res[0] = 16*parameters[0]*x[1]*x[2]*sin(M_PI*t*0.125)*(H-x[1])*(H-x[2])/ (H*H*H*H);
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflowParabolic3D(double* x, double* res, double t, const double* parameters){

    double H = parameters[1];

    if(t < parameters[2]){
        res[0] = 16*parameters[0]*x[1]*(H-x[1])*x[2]*(H-x[2])/(H*H*H*H)* 0.5 *( 1 - cos( M_PI*t/parameters[2]) );
        res[1] = 0.;
        res[2] = 0.;
    }
    else{
        res[0] = 16*parameters[0]*x[1]*(H-x[1])*x[2]*(H-x[2])/(H*H*H*H);
        res[1] = 0.;
        res[2] = 0.;

    }
    

    return;
}
void inflowPoiseuille3D(double* x, double* res, double t, const double* parameters){

    double maxVelo = parameters[0];
    double iota = 16*maxVelo;

    res[0] = iota * x[1] * (1.-x[1]) * x[2]*(1-x[2]);
    res[1] = 0.;
    res[2] = 0.;

    return;
}
void inflowPoiseuille2D(double* x, double* res, double t, const double* parameters){

    double maxVelo = parameters[0];
    double iota = 4*maxVelo;

    res[0] = iota * x[1] * (1.-x[1]);
    res[1] = 0.;

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

void dummyFuncRhs(double* x, double* res, double* parameters){
    if(parameters[0]==2)
        res[0]=1;
    else
        res[0] = 0.;

    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters){

    return;
}
using namespace Teuchos;
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

    double maxVelocity = -999;
    myCLP.setOption("maxVelocity",&maxVelocity,"maximum inflow velocity");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        MPI_Finalize();
        return 0;
    }
    Teuchos::RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Unsteady Navier-Stokes",true));
    TimeMonitor::setStackedTimer(stackedTimer);
    {

        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);

        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);
        
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);

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

        if(maxVelocity < -998)
            maxVelocity= parameterListProblem->sublist("Parameter").get("MaxVelocity",2.);   

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
                            domainPressure->buildMesh( 2,"BFS", dim, feTypeP, n, m, numProcsCoarseSolve);
                            domainVelocity->buildMesh( 2,"BFS", dim, feTypeV, n, m, numProcsCoarseSolve);
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
            if(parameterListProblem->sublist("Parameter").get("Robin BC",false)==true)
            {
                // if(!meshType.compare("structured") || !meshType.compare("structured_bfs")){
                //     domainPressure->getMesh()->buildEdges(domainPressure->getElementsC());
                    
                //     domainPressure->setUnstructuredMesh(domainPressure->getMesh());
                //     domainVelocity->buildP2ofP1Domain( domainPressure );
                // }
                //domainPressure->exportMesh(true,false,"BFS_h_H_25_9_subdomains.mesh");
                //domainVelocity->exportNodeFlags();
                domainVelocity->preProcessMesh(true,false);
                domainVelocity->preProcessMesh(true,false);

                domainPressure->preProcessMesh(true,false);
            }
            if(verbose)
                cout << " Maximum Inflow Velocity " << maxVelocity << endl;
                
            std::vector<double> parameter_vec(1);
            if ( !bcType.compare("parabolic") || !bcType.compare("parabolic_benchmark") || !bcType.compare("parabolic_benchmark_sin") || !bcType.compare("poiseuille") )
                parameter_vec[0] = maxVelocity;
            else if ( !bcType.compare("partialCFD") ) //  Fuer CFD3
                parameter_vec[0] = parameterListProblem->sublist("Parameter").get("MeanVelocity",2.);
                        
            // ####################
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryPressureLaplace( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryPressureFp( new BCBuilder<SC,LO,GO,NO>( ) );

            if (!bcType.compare("parabolic") || !bcType.compare("Couette") || !bcType.compare("poiseuille"))
                parameter_vec.push_back(1.);//height of inflow region
            else if(!bcType.compare("parabolic_benchmark_sin") || !bcType.compare("parabolic_benchmark") || !bcType.compare("partialCFD"))
                parameter_vec.push_back(.41);//height of inflow region
            else if(!bcType.compare("Richter3D"))
                parameter_vec.push_back(.4);
            else
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Select a valid boundary condition.");

            parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Max Ramp Time",0.1) );

            string pcdBC = parameterListProblem->sublist("Parameter").get("PCD BC","Inlet");

            if (!bcType.compare("parabolic") || !bcType.compare("parabolic_benchmark")|| !bcType.compare("poiseuille")|| !bcType.compare("parabolic_benchmark_sin") ) {//flag of obstacle
                if (dim==2){
                    if(!bcType.compare("poiseuille"))
                        bcFactory->addBC(inflowPoiseuille2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                    else if(!bcType.compare("parabolic_benchmark_sin"))
                        bcFactory->addBC(inflowParabolic2DSinBenchmark, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                    else 
                        bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);

                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
//                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                  
                    // bcFactoryPressureLaplace->addBC(zeroDirichlet2D, 3, 0, domainPressure, "Dirichlet", 1);
                    if( !pcdBC.compare("Inlet")){

                        bcFactoryPressureLaplace->addBC(zeroDirichlet2D, 2, 0, domainPressure, "Dirichlet", 1);

                        bcFactoryPressureFp->addBC(zeroDirichlet2D, 2, 0, domainPressure, "Dirichlet", 1);
                    }
                    else if( !pcdBC.compare("Outlet")){
                        bcFactoryPressureLaplace->addBC(zeroDirichlet2D, 3, 0, domainPressure, "Dirichlet", 1);

                        bcFactoryPressureFp->addBC(zeroDirichlet2D, 3, 0, domainPressure, "Dirichlet", 1);
                    }
                    else if( !pcdBC.compare("Mixed")){
                        bcFactoryPressureLaplace->addBC(zeroDirichlet2D, 3, 0, domainPressure, "Dirichlet", 1);

                        bcFactoryPressureFp->addBC(zeroDirichlet2D, 2, 0, domainPressure, "Dirichlet", 1);
                    }
    

                }
                else if (dim==3){
                    bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
                    if(!bcType.compare("poiseuille"))
                        bcFactory->addBC(inflowPoiseuille3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                    else if(!bcType.compare("parabolic_benchmark_sin"))
                        bcFactory->addBC(inflowParabolic3DSinBenchmark, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                    else 
                        bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
                        
//                 //                    bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                    bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet3D, 4, 0, domainVelocity, "Dirichlet", dim);
                 
                    if( !pcdBC.compare("Inlet")){
                            if(verbose)
                                cout << " --------- PCD Info: Setting inlet of Laplace and Fp to Dirichlet ----------- " << endl;
                            bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 2, 0, domainPressure, "Dirichlet", 1);

                            bcFactoryPressureFp->addBC(zeroDirichlet3D, 2, 0, domainPressure, "Dirichlet", 1);
                        }
                        else if( !pcdBC.compare("BC0")){
                            if(verbose)
                                cout << " --------- PCD Info (BC-0): Setting outlet of Laplace and Fp to Dirichlet ----------- " << endl;   
                            bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                            bcFactoryPressureFp->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);
                        }
                        // else if( !pcdBC.compare("InletOutletWall")){
                        //      if(verbose)
                        //         cout << " --------- PCD Info: Setting outlet of Laplace to dirichlet. Setting inlet, outlet and wall of Fp to Dirichlet ----------- " << endl;
                        //     bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                        //     bcFactoryPressureFp->addBC(zeroDirichlet3D, 1, 0, domainPressure, "Dirichlet", 1);
                        //     bcFactoryPressureFp->addBC(zeroDirichlet3D, 2, 0, domainPressure, "Dirichlet", 1);
                        //     bcFactoryPressureFp->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                        // }
                        else if( !pcdBC.compare("BC3")){
                             if(verbose)
                                cout << " --------- PCD Info: Setting outlet of Laplace to dirichlet. Setting outlet and wall of Fp to Dirichlet ----------- " << endl;
                            bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                            bcFactoryPressureFp->addBC(zeroDirichlet3D, 1, 0, domainPressure, "Dirichlet", 1);
                            bcFactoryPressureFp->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                        }
                        // else if( !pcdBC.compare("OutletWall2")){
                        //     if(verbose)
                        //         cout << " --------- PCD Info: Setting outlet and wall of Laplace and Fp to Dirichlet ----------- " << endl;
                        //     bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 1, 0, domainPressure, "Dirichlet", 1);
                        //     bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                        //     bcFactoryPressureFp->addBC(zeroDirichlet3D, 1, 0, domainPressure, "Dirichlet", 1);
                        //     bcFactoryPressureFp->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                        // }
                        else if( !pcdBC.compare("BC1")){
                            if(verbose)
                                cout << " --------- PCD Info (BC-1): Setting outlet of Laplace and inlet and outlet of Fp to Dirichlet ----------- " << endl;
                            bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);

                            bcFactoryPressureFp->addBC(zeroDirichlet3D, 2, 0, domainPressure, "Dirichlet", 1);
                            bcFactoryPressureFp->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);
                        }
                        else if( !pcdBC.compare("BC4")){
                            if(verbose)
                                cout << " --------- PCD Info (BC-4): Setting outlet of Laplace ----------- " << endl;
                            bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);
                        }
                        else if( !pcdBC.compare("Mixed")){
                            if(verbose)
                                cout << " --------- PCD Info: Setting outlet of Laplace and inlet Fp to Dirichlet ----------- " << endl;
                            bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 3, 0, domainPressure, "Dirichlet", 1);
                            
                            bcFactoryPressureFp->addBC(zeroDirichlet3D, 2, 0, domainPressure, "Dirichlet", 1);
                        }
                       
    
                    
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
                
            }
            else if (!bcType.compare("Couette")){
              bcFactory->addBC(couette2D, 1, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // wall
              bcFactory->addBC(zeroDirichlet2D, 2, 0, domainVelocity, "Dirichlet", dim); // wall
              //bcFactory->addBC(inflow3DRichter, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // inflow
              //bcFactory->addBC(zeroDirichlet3D, 3, 0, domainVelocity, "Dirichlet_Z", dim);
              //bcFactory->addBC(zeroDirichlet3D, 5, 0, domainVelocity, "Dirichlet", dim);
              bcFactoryPressureLaplace->addBC(zeroDirichlet2D, 3, 0, domainPressure, "Dirichlet", 1);
              bcFactoryPressureFp->addBC(zeroDirichlet2D, 2, 0, domainPressure, "Dirichlet", 1);


            }

            
            int timeDisc = parameterListProblem->sublist("Timestepping Parameter").get("Butcher table",0);

            NavierStokes<SC,LO,GO,NO> navierStokes( domainVelocity, feTypeV, domainPressure, feTypeP, parameterListAll );

            navierStokes.addBoundaries(bcFactory);
            navierStokes.addBoundariesPressureLaplace(bcFactoryPressureLaplace);
            navierStokes.addBoundariesPressureFp(bcFactoryPressureFp);

            navierStokes.addRhsFunction( dummyFuncRhs );

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




            if (verbose) {
                cout << "###############################################################" <<endl;
                cout << "##################### Steady Navier-Stokes ####################" <<endl;
                cout << "Discretization: \t" << feTypeV << "-" << feTypeP  << endl;
                if (!precMethod.compare("Monolithic")){
                cout << "Coarse Opertor Type: \t" << parameterListPrec->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("CoarseOperator Type","NOTFOUND") << endl;
                cout << "IPOU Block 1: \t \t" << parameterListPrec->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").sublist("IPOUHarmonicCoarseOperator").sublist("Blocks").sublist("1").sublist("InterfacePartitionOfUnity").get("Type","NOTFOUND") << endl;
                cout << "IPOU Block 2: \t \t" << parameterListPrec->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").sublist("IPOUHarmonicCoarseOperator").sublist("Blocks").sublist("2").sublist("InterfacePartitionOfUnity").get("Type","NOTFOUND") << endl;
                }
                else if (!precMethod.compare("Teko")){
                    cout << "Block Preconditioner Type: \t" << parameterListAll->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").get("Inverse Type","SIMPLE") << endl;
                    cout << "Velocity Preconditioner: \t" << parameterListAll->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("FROSch-Velocity").get("CoarseOperator Type","GDSW#") << endl;
                    cout << "Pressure Preconditioner: \t" << parameterListAll->sublist("Teko Parameters").sublist("Preconditioner Types").sublist("Teko").sublist("Inverse Factory Library").sublist("FROSch-Pressure").get("CoarseOperator Type","GDSW#") << endl;

                }            
                cout << "###############################################################" <<endl;
            }


        }
    }

    Teuchos::TimeMonitor::report(cout);
    stackedTimer->stop("Unsteady Navier-Stokes");
	StackedTimer::OutputOptions options;
	options.output_fraction = options.output_histogram = options.output_minmax = true;
	stackedTimer->report((std::cout),comm,options);
    return(EXIT_SUCCESS);
}
