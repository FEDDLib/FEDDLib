#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

#include "feddlib/problems/specific/Stokes.hpp"
#include "feddlib/amr/AdaptiveMeshRefinement.hpp"

#include <Teuchos_TestForException.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

#include "feddlib/problems/abstract/Problem.hpp"

/*!
 main of Stokes problem
 
 @brief Stokes main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */



// ######################
// Paper
// ######################
void rhsPaper1( double* p, double* res, const double* parameters){

	double x = p[0];
	double y = p[1];
	res[0] =-(-4*y*pow((1-x),2)+16*x*y*(1-x)-4*pow(x,2)*y)*(1-3*y+2*pow(y,2))-(-4*pow(x,2)*(-3+4*y)-8*pow(x,2)*y)*pow((1-x),2)+1;
    res[1] =-(4*pow(y,2)*(-3+4*x)+8*pow(y,2)*x)*pow((1-y),2)-(4*x*pow((1-y),2)-16*x*y*(1-y)+4*pow(y,2)*x)*(1-3*x+2*pow(x,2))+1;

	//cout << " res[0] " << res[0] << " res[1] " << res[1] << endl;
}


void rhsPaper2( double* p, double* res, const double* parameters){

	double x = p[0];
	double y = p[1];

    res[0] =-(-4*y*pow((1-x),2)+16*x*y*(1-x)-4*pow(x,2)*y)*(1-3*y+2*pow(y,2))-(-4*pow(x,2)*(-3+4*y)-8*pow(x,2)*y)*pow((1-x),2)-sin(M_PI*x)*cos(M_PI*y)*M_PI;
    res[1] = -(4*pow(y,2)*(-3+4*x)+8*pow(y,2)*x)*pow((1-y),2)-(4*x*pow((1-y),2)-16*x*y*(1-y)+4*pow(y,2)*x)*(1-3*x+2*pow(x,2))-sin(M_PI*y)*cos(M_PI*x)*M_PI;
}

void exactSolutionPaperU1( double* p, double* res){

	double x = p[0];
	double y = p[1];

	res[0] =  -2*pow(x,2)*y*pow((1-x),2)*(1-3*y+2*pow(y,2));
	res[1] =  2*x*pow(y,2)*pow((1-y),2)*(1-3*x+2*pow(x,2));
}

void exactSolutionPaperP1( double* p, double* res){

	double x = p[0];
	double y = p[1];

	res[0] =  x+y-1;
}

void exactSolutionPaperP2( double* p, double* res){

	double x = p[0];
	double y = p[1];

	res[0] =  cos(x*M_PI) *cos (M_PI *y);	
}

// ####################################
// ####################################

void rhs0( double* p, double* res, const double* parameters){

	res[0] =0;
    res[1] =0;
	res[2] =0;

	//cout << " res[0] " << res[0] << " res[1] " << res[1] << endl;
}

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

	typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;


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
		int 		maxIter 		= parameterListProblem->sublist("Mesh Refinement").get("MaxIter",5);
		double maxVel 				= parameterListProblem->sublist("Parameter").get("MaxVelocity",2.);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem));
        if (precMethod == "Monolithic")
            parameterListAll->setParameters(*parameterListPrec);
        else if(precMethod == "Teko")
            parameterListAll->setParameters(*parameterListPrecTeko);
        else if(precMethod == "Diagonal" || precMethod == "Triangular")
            parameterListAll->setParameters(*parameterListPrecBlock);
        
        parameterListAll->setParameters(*parameterListSolver);
        
        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("unstructured_struct")) {
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
  
		domainPressure.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
		domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
		
		MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
		domainP1Array[0] = domainPressure;
		
		ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
		MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
		
		partitionerP1.readAndPartition();
		
		Teuchos::RCP<Domain<SC,LO,GO,NO> > domainRefined;

		AdaptiveMeshRefinement<SC,LO,GO,NO> meshRefiner("Stokes",parameterListProblem,exactSolutionPaperU1); 
		
		std::vector<double> parameter_vec(0);
		parameter_vec.push_back(maxVel);//height of inflow region
		parameter_vec.push_back(0.41);//height of inflow region

		int j=0;
		MAIN_TIMER_START(Total," Step 4:	 Total RefinementAlgorithm");
		while(j<maxIter+1 ){

			MAIN_TIMER_START(buildP2," Step 0:	 buildP2Mesh");
			if (discVelocity=="P2" ) {
				domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
		        domainVelocity->buildP2ofP1Domain( domainPressure );
				}
		    else
		        domainVelocity = domainPressure;
			
			MAIN_TIMER_STOP(buildP2);		

			MAIN_TIMER_START(Bounds," Step 1:	 bcFactory");
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

			if (dim==2) {
				bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
				bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
				bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
	//                bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
	//                bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
			}
			else if(dim==3){
				bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
				bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
				bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
	//                bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
	//                bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
			}

			MAIN_TIMER_STOP(Bounds);	
			MAIN_TIMER_START(Solver," Step 2:	 solving PDE");

			
            Teuchos::RCP<Stokes<SC,LO,GO,NO> > stokes( new Stokes<SC,LO,GO,NO>(domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll ));

			//domainVelocity->info();
			//domainPressure->info();
			//stokes->info();
			
			{
				Teuchos::TimeMonitor solveTimeMonitor(*solveTime);
				
				stokes->addBoundaries(bcFactory);
				stokes->addRhsFunction(rhs0);						    
				stokes->initializeProblem();						    
				stokes->assemble();
				stokes->setBoundaries();             
				stokes->solve();

			}
			MAIN_TIMER_STOP(Solver);	

						
			
	
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaVelocity(new ExporterParaView<SC,LO,GO,NO>());
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaPressure(new ExporterParaView<SC,LO,GO,NO>());
            
            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionV = stokes->getSolution()->getBlock(0);
            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionP = stokes->getSolution()->getBlock(1);

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

			MAIN_TIMER_START(Refinement," Step 3:	 meshRefinement");

			// Refinement
			domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );

			{

				ProblemPtr_Type problem = Teuchos::rcp_dynamic_cast<Problem_Type>( stokes , true);
				domainRefined = meshRefiner.globalAlgorithm( domainPressure,  domainVelocity, stokes->getSolution(), problem, rhs0 );
			}

			domainPressure = domainRefined;
			domainVelocity = domainPressure;
			
			j++;
			MAIN_TIMER_STOP(Refinement);	
        
        // ####################
       
            
        }

		MAIN_TIMER_STOP(Total);	
		Teuchos::TimeMonitor::report(cout,"Main");
   
    }

    return(EXIT_SUCCESS);
}
