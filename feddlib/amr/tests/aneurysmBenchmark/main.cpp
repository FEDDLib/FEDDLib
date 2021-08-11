#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include "feddlib/amr/AdaptiveMeshRefinement.hpp"

#include <Teuchos_TestForException.hpp>

#include "feddlib/problems/abstract/Problem.hpp"

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"

#include "feddlib/problems/specific/Laplace.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 main of time-dependent Navier-Stokes problem
 
 @brief time-dependent Navier-Stokes main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

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
void rhs0( double* p, double* res, const double* parameters){

	res[0] =0;
    res[1] =0;
	res[2] =0;

	//cout << " res[0] " << res[0] << " res[1] " << res[1] << endl;
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


void dummyFunc(double* x, double* res, double t, const double* parameters){

    return;
}
// Function for the rhs
void rhs3D(double* x, double* res, double* parameters){
    
    res[0] = 1.;


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
        res[2] = parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = parameters[0] / parameters[2] * x[0];
        res[1] = 0.;
        res[2] = 0.;
    }

    return;
}
void parabolicInflow3DStokes(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    

    res[0] = parameters[0] / parameters[2] * x[0];
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
int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;
typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
	typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;

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

    string xmlProblemFileLaplace = "parametersProblem_Laplace.xml";
    myCLP.setOption("problemfile",&xmlProblemFileLaplace,".xml file with Inputparameters.");
    string xmlPrecFileLaplace = "parametersPrec_Laplace.xml";
    myCLP.setOption("precfile",&xmlPrecFileLaplace,".xml file with Inputparameters.");
    string xmlSolverFileLaplace = "parametersSolver_Laplace.xml";
    myCLP.setOption("solverfile",&xmlSolverFileLaplace,".xml file with Inputparameters.");

 	myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }


    ParameterListPtr_Type parameterListProblemL = Teuchos::getParametersFromXmlFile(xmlProblemFileLaplace);
    ParameterListPtr_Type parameterListPrecL = Teuchos::getParametersFromXmlFile(xmlPrecFileLaplace);
    ParameterListPtr_Type parameterListSolverL = Teuchos::getParametersFromXmlFile(xmlSolverFileLaplace);

    ParameterListPtr_Type parameterListAllL(new Teuchos::ParameterList(*parameterListProblemL)) ;
    parameterListAllL->setParameters(*parameterListPrecL);
    parameterListAllL->setParameters(*parameterListSolverL);



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
		int 		maxIter 		= parameterListProblem->sublist("Mesh Refinement").get("MaxIter",5);
		double maxVel 				= parameterListProblem->sublist("Parameter").get("MaxVelocity",2.);


        std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("Max Velocity",30));
        parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Max Ramp Time",3.) );
        //parameter_vec.push_back(1.);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if (!precMethod.compare("Monolithic"))
            parameterListAll->setParameters(*parameterListPrec);
        else
            parameterListAll->setParameters(*parameterListPrecTeko);

        parameterListAll->setParameters(*parameterListSolver);

        std::string bcType = parameterListProblem->sublist("Parameter").get("BC Type","parabolic");

        int minNumberSubdomains;

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

        double viscosity = parameterListProblem->sublist("Parameter").get("Viscosity",1.e-3);

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));
        
        {
            DomainPtr_Type domainPressure;
            DomainPtr_Type domainVelocity;

            DomainPtr_Type domainRefined;

            Teuchos::TimeMonitor totalTimeMonitor(*totalTime);
            
            Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
            if (verbose)
                cout << "-- Building Mesh ..." << flush;
            

	        domainPressure.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
	        domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
	        
	        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
	        domainP1Array[0] = domainPressure;
	        
	        ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
	        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
	        
	        partitionerP1.readAndPartition();
		


			AdaptiveMeshRefinement<SC,LO,GO,NO> meshRefiner("NavierStokes",parameterListAll);

		    int j=0;
			MAIN_TIMER_START(Total," Step 4:	 Total RefinementAlgorithm");
			while(j<maxIter+1 ){
				MAIN_TIMER_START(buildP2," Step 0:	 buildP2Mesh");
				if (feTypeV=="P2" ) {
					domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
				    domainVelocity->buildP2ofP1Domain( domainPressure );
					}
				else
				    domainVelocity = domainPressure;
				MAIN_TIMER_STOP(buildP2);		
				// ############################################################
				// Genereating Vector for inlet 
				// ############################################################
				MAIN_TIMER_START(LaplaceInlet," Step 1:	 laplaceVector for Inlet");
				Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryL(new BCBuilder<SC,LO,GO,NO>( ));

				// Apply Boundary conditions - laplace at inlet
				bcFactoryL->addBC(zeroDirichlet, 1, 0, domainVelocity, "Dirichlet", 1);
				//bcFactoryL->addBC(zeroDirichlet, 2, 0, domainVelocity, "Neumann", 1);

      			Laplace<SC,LO,GO,NO> laplace(domainVelocity,feTypeV,parameterListAllL,false);

				{
					laplace.addBoundaries(bcFactoryL);
					laplace.addRhsFunction(rhs3D);
					laplace.initializeProblem();
			   		laplace.assemble();
			   		laplace.setBoundaries();
			   		laplace.solve();
				}

				bcFactoryL->addBC(zeroBC, 3, 0, domainVelocity, "Dirichlet", 1);
		        bcFactoryL->addBC(zeroBC, 10, 0, domainVelocity, "Dirichlet", 1);

		        bcFactoryL->setRHS( laplace.getSolution(), 0.);


				Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

				Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolution = laplace.getSolution()->getBlock(0);

				/*exPara->setup("solutionLaplace", domainVelocity->getMesh(), feTypeV);
			    
			    exPara->addVariable(exportSolution, "u", "Scalar", 1, domainVelocity->getMapUnique(), domainVelocity->getMapUniqueP2());

			    exPara->save(0.0);*/
		        

				// #################################################
				// Setting up and solving steady Navier-Stokes Problem 
				// #################################################		    
				    
				MAIN_TIMER_STOP(LaplaceInlet);	
							
				MAIN_TIMER_START(Bounds," Step 2:	 bcFactory");
				
				SC maxValue = exportSolution->getMax();
		            
		        parameter_vec.push_back(maxValue);

		        //parameter_vec[0] = parameterListProblem->sublist("Parameter").get("MaxVelocity",1.5);
				MultiVectorConstPtr_Type inletSol = laplace.getSolution()->getBlock(0);               

				Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory(new BCBuilder<SC,LO,GO,NO>( ));

		        bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim, parameter_vec);

		        bcFactory->addBC(parabolicInflow3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec, inletSol);

				//bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
		        //bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
		        
				MAIN_TIMER_STOP(Bounds);	
				MAIN_TIMER_START(Solver," Step 3:	 solving PDE");

				Teuchos::RCP<NavierStokes<SC,LO,GO,NO>> navierStokes( new NavierStokes<SC,LO,GO,NO> (domainVelocity, feTypeV, domainPressure, feTypeP, parameterListAll ));

			    navierStokes->info();

			    {
			        Teuchos::TimeMonitor solveTimeMonitor(*solveTime);

			        navierStokes->addBoundaries(bcFactory);
					navierStokes->addRhsFunction(rhs0);						    
			        navierStokes->initializeProblem();
			        navierStokes->assemble();
			        navierStokes->setBoundariesRHS();

			        std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
			        NonLinearSolver<SC,LO,GO,NO> nlSolver( nlSolverType );
			        nlSolver.solve( *navierStokes );
			        comm->barrier();
			    }

				MAIN_TIMER_STOP(Solver);	


				MAIN_TIMER_START(Refinement," Step 3:	 meshRefinement");

				// Refinement
				domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );

				{

					ProblemPtr_Type problem = Teuchos::rcp_dynamic_cast<Problem_Type>( navierStokes , true);
					domainRefined = meshRefiner.globalAlgorithm( domainPressure,  domainVelocity, navierStokes->getSolution(), problem, rhs0 );
				}

				domainPressure = domainRefined;
				domainVelocity = domainPressure;
				
				j++;
				MAIN_TIMER_STOP(Refinement);	
		    
		    // ####################
		   
		        
		    }

        }
    }

    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
