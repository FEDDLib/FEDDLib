#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_StackedTimer.hpp>

/*! Test case for specific artery geometrie or straight tube geometry. Inflow depends on inflow region
	-> straight tube: Inflow in (0,0,z)*laplaceInflow direction
	-> artery: Inflow scaled with normal vector on inflow (x,y,z) * laplaceInflow	

*/



void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void parabolicInflow3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];

    }

    return;
}

void parabolicInflow3DLin(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];
    }

    return;
}

void parabolicInflow3DArtery(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[1] = 0.;
        res[0] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];
    }

    return;
}

void parabolicInflow3DLinArtery(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];
    }

    return;
}
void parabolicInflow3DArteryHeartBeat(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // parameters[3] heart beat start
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    double heartBeatStart = parameters[3];
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else if(t > heartBeatStart)
    {
    
        double a0    = 11.693284502463376;
        double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
            0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
            0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
        double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
            -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
            0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                    
        double Q = 0.5*a0;
        

        double t_min = t - fmod(t,1.0)+heartBeatStart-std::floor(t)+0.25; ; //FlowConditions::t_start_unsteady;
        double t_max = t_min + 1.0; // One heartbeat lasts 1.0 second    
        double y = M_PI * ( 2.0*( t-t_min ) / ( t_max - t_min ) -1.0)  ;
        
        for(int i=0; i< 20; i++)
            Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
        
        
        // Remove initial offset due to FFT
        Q -= 0.026039341343493;
        Q = (Q - 2.85489)/(7.96908-2.85489);
        double lambda = 1.;

        if( t+1.0e-12 < heartBeatStart + 0.25)
		    lambda = 0.90 + 0.1*cos(4*M_PI*(t-heartBeatStart));
        else 
    	    lambda= 0.8 + 1.2*Q;

        res[0] = 0.;
        res[1] = 0.;
        res[2] = (parameters[0] / parameters[2]) * (x[0] * lambda) ;
        
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0];

    }

    return;
}

void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}


void dummyFunc(double* x, double* res, double* parameters){
    if(parameters[0]==4)
        res[0]=1;
    else
        res[0] = 0.;

    return;
}


typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;

int main(int argc, char *argv[])
{


    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
     string xmlSolverFile = "parametersSolver.xml"; // GI
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
  
    string xmlPrecFileFluidMono = "parametersPrec.xml";
    string xmlPrecFileFluidTeko = "parametersPrecTeko.xml";
    myCLP.setOption("precfile",&xmlPrecFileFluidMono,".xml file with Inputparameters.");
    myCLP.setOption("precfileTeko",&xmlPrecFileFluidTeko,".xml file with Inputparameters.");
      
    string xmlProbL = "plistProblemLaplace.xml";
    myCLP.setOption("probLaplace",&xmlProbL,".xml file with Inputparameters.");
    string xmlPrecL = "plistPrecLaplace.xml";
    myCLP.setOption("precLaplace",&xmlPrecL,".xml file with Inputparameters.");
    string xmlSolverL = "plistSolverLaplace.xml";
    myCLP.setOption("solverLaplace",&xmlSolverL,".xml file with Inputparameters.");
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    Teuchos::RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Unsteady Navier-Stokes",true));
    TimeMonitor::setStackedTimer(stackedTimer);

    bool verbose (comm->getRank() == 0);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);
        ParameterListPtr_Type parameterListPrecFluidMono = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidMono);
        ParameterListPtr_Type parameterListPrecFluidTeko = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidTeko);
        
        string		precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if (!precMethod.compare("Monolithic"))
            parameterListAll->setParameters(*parameterListPrecFluidMono);
        else
            parameterListAll->setParameters(*parameterListPrecFluidTeko);

        parameterListAll->setParameters(*parameterListSolver);

   
        // Fuer das Geometrieproblem, falls GE       
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");

        string feTypeV = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        string feTypeP = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");
        string preconditionerMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);        

        TimePtr_Type totalTime(TimeMonitor_Type::getNewCounter("FEDD - main - Total Time"));
        TimePtr_Type buildMesh(TimeMonitor_Type::getNewCounter("FEDD - main - Build Mesh"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        // #####################
        // Mesh bauen und wahlen
        // #####################
        {
            if (verbose)
            {
                cout << "###############################################" <<endl;
                cout << "############ Starting Unsteady NS... ##########" <<endl;
                cout << "###############################################" <<endl;
            }

            DomainPtr_Type domainP1fluid;
            DomainPtr_Type domainP2fluid;
                    
            DomainPtr_Type domainFluidVelocity;
            DomainPtr_Type domainFluidPressure;

            
            std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","Compute Inflow");
            std::string geometryType = parameterListAll->sublist("Parameter").get("Geometry Type","Artery");
            
            {
                TimeMonitor_Type totalTimeMonitor(*totalTime);
                {
                    TimeMonitor_Type buildMeshMonitor(*buildMesh);
                    if (verbose)
                    {
                        cout << " -- Building Mesh ... " << flush;
                    }

                    domainP1fluid.reset( new Domain_Type( comm, dim ) );
                    domainP2fluid.reset( new Domain_Type( comm, dim ) );
                    //                    
                    if (!meshType.compare("structured")) {
                        int minNumberSubdomains = 1;
                        TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
                        
                        n = (int) (std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                        std::vector<double> x(3);
                        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
                        // setting length, width and depth in cm. Approximating the realistic geometry with a rectangular channel with similar size. 
                        domainFluidPressure.reset(new Domain<SC,LO,GO,NO>( x, 0.27, 0.27, 0.54, comm)); // 5 Subcubes 1.35
                        domainFluidVelocity.reset(new Domain<SC,LO,GO,NO>( x, 0.27, 0.27, 0.54, comm));
                        
                        domainFluidPressure->buildMesh( 4,"Tube", dim, feTypeP, n, m, numProcsCoarseSolve);
                        domainFluidVelocity->buildMesh( 4,"Tube", dim, feTypeV, n, m, numProcsCoarseSolve);
                    }
                    else if (!meshType.compare("unstructured")) {
                                                
                        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
                        domainP1Array[0] = domainP1fluid;
                        
                        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
                        if (!feTypeV.compare("P2")){
                            pListPartitioner->set("Build Edge List",true);
                            pListPartitioner->set("Build Surface List",true);
                        }
                        else{
                            pListPartitioner->set("Build Edge List",false);
                            pListPartitioner->set("Build Surface List",false);
                        }
                        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                        
                        partitionerP1.readAndPartition(15, "mm",true); // converting mesh from mm unit to cm unit
                        
                        if(parameterListProblem->sublist("General").get("ParaViewCoarse",false)){
                            domainP1fluid->exportElementFlags("Fluid");
                            domainP1fluid->exportNodeFlags("Fluid");
                            domainP1fluid->exportProcessor("Distribution")
                        }



                        if (!feTypeV.compare("P2")){
                            domainP2fluid->buildP2ofP1Domain( domainP1fluid );
                        }
                        
                        
						if (!feTypeV.compare("P2"))
						{
							domainFluidVelocity = domainP2fluid;
							domainFluidPressure = domainP1fluid;
						}
						else
						{
							domainFluidVelocity = domainP1fluid;
							domainFluidPressure = domainP1fluid;
						}
                    }
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
                domainFluidVelocity->preProcessMesh(true,false);

                domainFluidPressure->preProcessMesh(true,false);
            }
            //domainFluidPressure->setUnstructuredMesh(domainFluidPressure->getMesh());
            //domainFluidPressure->exportMesh(" ");
            //domainFluidVelocity->exportProcessor("Fluid");

                     
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("Max Velocity",1.));
            parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Max Ramp Time",0.1) );

            TEUCHOS_TEST_FOR_EXCEPTION(bcType != "Compute Inflow", std::logic_error, "Select a valid boundary condition. Only Compute Inflow available.");

            //#############################################
            //#############################################
            //#### Compute parabolic inflow with laplacian
            //#############################################
            //#############################################
            MultiVectorConstPtr_Type solutionLaplace;
            {
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryLaplace(new BCBuilder<SC,LO,GO,NO>( ));
                
                bcFactoryLaplace->addBC(zeroBC, 9, 0, domainFluidVelocity, "Dirichlet", 1); //inflow ring
                bcFactoryLaplace->addBC(zeroBC, 10, 0, domainFluidVelocity, "Dirichlet", 1); //outflow ring
                bcFactoryLaplace->addBC(zeroBC, 6, 0, domainFluidVelocity, "Dirichlet", 1); //surface
                
                ParameterListPtr_Type parameterListProblemL = Teuchos::getParametersFromXmlFile(xmlProbL);
                ParameterListPtr_Type parameterListPrecL = Teuchos::getParametersFromXmlFile(xmlPrecL);
                ParameterListPtr_Type parameterListSolverL = Teuchos::getParametersFromXmlFile(xmlSolverL);

                ParameterListPtr_Type parameterListLaplace(new Teuchos::ParameterList(*parameterListProblemL)) ;
                parameterListLaplace->setParameters(*parameterListPrecL);
                parameterListLaplace->setParameters(*parameterListSolverL);
                
                Laplace<SC,LO,GO,NO> laplace( domainFluidVelocity, feTypeV, parameterListLaplace, false );
                {
                    laplace.addRhsFunction(oneFunc);
                    laplace.addBoundaries(bcFactoryLaplace);
                    
                    laplace.initializeProblem();
                    laplace.assemble();
                    laplace.setBoundaries();
                    laplace.solve();
                }
                
                //We need the values in the inflow area. Therefore, we use the above bcFactory and the volume flag 10 and the outlet flag 5 and set zero Dirichlet boundary values
                bcFactoryLaplace->addBC(zeroBC, 5, 0, domainFluidVelocity, "Dirichlet", 1);
                bcFactoryLaplace->addBC(zeroBC, 15, 0, domainFluidVelocity, "Dirichlet", 1);
                bcFactoryLaplace->setRHS( laplace.getSolution(), 0./*time; does not matter here*/ );
                solutionLaplace = laplace.getSolution()->getBlock(0);
            
                SC maxValue = solutionLaplace->getMax();
                
                parameter_vec.push_back(maxValue);



                /*Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup("parabolicInflow", domainFluidVelocity->getMesh(), feTypeV);
                
                MultiVectorConstPtr_Type valuesConst = laplace.getSolution()->getBlock(0);
                exPara->addVariable( valuesConst, "values", "Scalar", 1, domainFluidVelocity->getMapUnique() );

                exPara->save(0.0);
                exPara->closeExporter();*/

            }
            parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Heart Beat Start",0.2) ); // Adding the heart beat start last

            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryPressureLaplace( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryPressureFp( new BCBuilder<SC,LO,GO,NO>( ) );

            // TODO: Vermutlich braucht man keine bcFactoryFluid und bcFactoryStructure,
            // da die RW sowieso auf dem FSI-Problem gesetzt werden.

            // Fluid-RW
            string pcdBC = parameterListProblem->sublist("Parameter").get("PCD BC","Inlet");

                                           
            //bcFactory->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
            string rampType = parameterListProblem->sublist("Parameter Fluid").get("Ramp type","cos");
            
            bcFactory->addBC(zeroDirichlet3D, 9, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow ring
            bcFactory->addBC(parabolicInflow3DArteryHeartBeat, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow
            bcFactory->addBC(zeroDirichlet3D, 6, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // Wall
            bcFactory->addBC(zeroDirichlet3D, 10, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // outflow ring
            
            if( !pcdBC.compare("Inlet")){
                if(verbose)
                    cout << " --------- PCD Info: Setting inlet of Laplace and Fp to Dirichlet ----------- " << endl;
                bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 4, 0, domainFluidPressure, "Dirichlet", 1);

                bcFactoryPressureFp->addBC(zeroDirichlet3D, 4, 0, domainFluidPressure, "Dirichlet", 1);
            }
            else if( !pcdBC.compare("BC0")){
                if(verbose)
                    cout << " --------- PCD Info (BC-0): Setting outlet of Laplace and Fp to Dirichlet ----------- " << endl;   
                bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 5, 0, domainFluidPressure, "Dirichlet", 1);

                bcFactoryPressureFp->addBC(zeroDirichlet3D, 5, 0, domainFluidPressure, "Dirichlet", 1);
            }
            else if( !pcdBC.compare("BC4")){
                if(verbose)
                    cout << " --------- PCD Info (BC-4): Setting outlet of Laplace ----------- " << endl;
                bcFactoryPressureLaplace->addBC(zeroDirichlet3D, 5, 0, domainFluidPressure, "Dirichlet", 1);
            }
                        
                
            

            int timeDisc = parameterListProblem->sublist("Timestepping Parameter").get("Butcher table",0);

            NavierStokes<SC,LO,GO,NO> navierStokes( domainFluidVelocity, feTypeV, domainFluidPressure, feTypeP, parameterListAll );

            navierStokes.addBoundaries(bcFactory);
            navierStokes.addBoundariesPressureLaplace(bcFactoryPressureLaplace);
            navierStokes.addBoundariesPressureFp(bcFactoryPressureFp);

            navierStokes.addRhsFunction( dummyFunc );

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
                cout << "Coarse Opertor Type: \t" << parameterListPrecFluidMono->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").get("CoarseOperator Type","NOTFOUND") << endl;
                cout << "IPOU Block 1: \t \t" << parameterListPrecFluidMono->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").sublist("IPOUHarmonicCoarseOperator").sublist("Blocks").sublist("1").sublist("InterfacePartitionOfUnity").get("Type","NOTFOUND") << endl;
                cout << "IPOU Block 2: \t \t" << parameterListPrecFluidMono->sublist("ThyraPreconditioner").sublist("Preconditioner Types").sublist("FROSch").sublist("IPOUHarmonicCoarseOperator").sublist("Blocks").sublist("2").sublist("InterfacePartitionOfUnity").get("Type","NOTFOUND") << endl;
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
