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

void bcPaper( double* p, double* res, double t, const double* parameters){

	double x = p[0];
	double y = p[1];

	res[0] =  -2*pow(x,2)*y*pow((1-x),2)*(1-3*y+2*pow(y,2));
	res[1] =  2*x*pow(y,2)*pow((1-y),2)*(1-3*x+2*pow(x,2));
}

// ####################################
// ######################
// Paper
// ######################
void rhsPaperSin( double* p, double* res, const double* parameters){

	double x = M_PI*p[0];
	double y = M_PI*p[1];

	res[0] = -1*(4*pow(M_PI,3)*cos(y)*sin(y)*(pow(cos(x),2)-pow(sin(x),2)) - 8*pow(M_PI,3)*pow(sin(x),2)*cos(y)*sin(y) )- M_PI*sin(x)*cos(y);
    res[1] = -1*(-4*pow(M_PI,3)*cos(x)*sin(x)*(pow(cos(y),2)-pow(sin(y),2)) + 8*pow(M_PI,3)*pow(sin(y),2)*cos(x)*sin(x)) - M_PI*cos(x)*sin(y);

}


void exactSolutionPaperSinU( double* p, double* res){

	double x = M_PI*p[0];
	double y = M_PI*p[1];

	res[0] =  2*M_PI*sin(x)*sin(x)*cos(y)*sin(y);
	res[1] =  -2*M_PI*sin(x)*sin(y)*cos(x)*sin(y);
}

void exactSolutionPaperSinP( double* p, double* res){

	double x = M_PI*p[0];
	double y = M_PI*p[1];

	res[0] =  cos(x)*cos(y);
}

void bcPaperSin( double* p, double* res, double t, const double* parameters){

	double x = p[0];
	double y = p[1];

	res[0] =  0.; 
	res[1] =  0.; 
}

// ####################################
// ######################
// Paper
// ######################
void rhsVer( double* x, double* res, const double* parameters){

	double r = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi=0.;
    if(x[1] < 0.0)
		phi = 2.0*M_PI+atan2(x[1],x[0]);
    else
		phi = atan2(x[1],x[0]);


	double psi = 3*sin(0.5*phi)-sin(1.5*phi);
	double dPsi = 1.5*cos(0.5*phi)-1.5*cos(1.5*phi);
	double ddPsi = -0.75*sin(0.5*phi)+2.25*sin(1.5*phi);
	double dddPsi = -0.375*cos(0.5*phi)+3.375*cos(1.5*phi);
	double ddddPsi = 0.1875*sin(0.5*phi)-5.0625*sin(1.5*phi);

	double alpha= 0.5;
	double omega= 2*M_PI;

	res[0] =0;
	res[1] =0;
	
}


void exactSolutionUVer( double* x, double* res){

	double r = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi;
    if(x[1] < 0.0)
		phi = 2.0*M_PI+atan2(x[1],x[0]);
    else
		phi = atan2(x[1],x[0]);

	double psi = 3*sin(0.5*phi)-sin(1.5*phi);
	double deriPsi = 1.5*cos(0.5*phi)-1.5*cos(1.5*phi);

	double alpha= 0.5;
	double omega= 2*M_PI;

	res[0] = 	pow(r,alpha)*((1+alpha)*sin(phi)*psi + cos(phi)*deriPsi);
	res[1] = 	pow(r,alpha)*(sin(phi)*deriPsi -(1+alpha)* cos(phi)*psi);
}

void exactSolutionPVer( double* x, double* res){

	double r = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi;
    if(x[1] < 0.0)
		phi = 2.0*M_PI+atan2(x[1],x[0]);
    else
		phi = atan2(x[1],x[0]);

	double psi = 3*sin(0.5*phi)-sin(1.5*phi);
	double dPsi = 1.5*cos(0.5*phi)-1.5*cos(1.5*phi);
	double ddPsi = -0.75*sin(0.5*phi)+2.25*sin(1.5*phi);
	double dddPsi = -0.375*cos(0.5*phi)+3.375*cos(1.5*phi);
	double ddddPsi = 0.1875*sin(0.5*phi)-5.0625*sin(1.5*phi);

	if(r==0)
		res[0] = 0;
	else 	
		res[0] = 	-2*pow(r,-0.5)*(2.25*dPsi+dddPsi);
}
void bcSolutionPVer( double* x, double* res, double t, const double* parameters){

	double r = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi;
    if(x[1] < 0.0)
		phi = 2.0*M_PI+atan2(x[1],x[0]);
    else
		phi = atan2(x[1],x[0]);

	double psi = 3*sin(0.5*phi)-sin(1.5*phi);
	double dPsi = 1.5*cos(0.5*phi)-1.5*cos(1.5*phi);
	double ddPsi = -0.75*sin(0.5*phi)+2.25*sin(1.5*phi);
	double dddPsi = -0.375*cos(0.5*phi)+3.375*cos(1.5*phi);
	double ddddPsi = 0.1875*sin(0.5*phi)-5.0625*sin(1.5*phi);

	res[0] = 	-2*pow(r,-0.5)*(2.25*dPsi+dddPsi);
}

void bcVer( double* x, double* res, double t, const double* parameters){

	double r = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi;
    if(x[1] < 0.0)
		phi = 2.0*M_PI+atan2(x[1],x[0]);
    else
		phi = atan2(x[1],x[0]);

	double psi = 3*sin(0.5*phi)-sin(1.5*phi);
	double deriPsi = 1.5*cos(0.5*phi)-1.5*cos(1.5*phi);

	double alpha= 0.5;
	double omega= 2*M_PI;

	res[0] = 	pow(r,alpha)*((1+alpha)*sin(phi)*psi + cos(phi)*deriPsi);
	res[1] = 	pow(r,alpha)*(sin(phi)*deriPsi -(1+alpha)* cos(phi)*psi);

}

// ####################################
// ######################
// Paper Han
// ######################
void rhsHan( double* p, double* res, const double* parameters){

	double x= p[0];
	double y=p[1];
	res[0] = -(12*pow(x,2)-12*x+2)*(4*pow(y,3)-6*pow(y,2)+2*y)+1;
	res[1] = -(12*pow(y,2)-12*y+2)*(4*pow(x,3)-6*pow(x,2)+2*x)+1;
}


void exactSolutionUHan( double* p, double* res){

	double x= p[0];
	double y=p[1];

	res[0] = pow(x-1,2)*pow(x,2)*2*(y-1)*pow(y,2)+pow(x-1,2)*x*x*pow((y-1),2)*2*y;
	res[1] = pow(y-1,2)*pow(y,2)*2*(x-1)*pow(x,2)+pow(y-1,2)*y*y*pow((x-1),2)*2*x;
}

void exactSolutionPHan( double* p, double* res){

	double x= p[0];
	double y=p[1];
	res[0] = x+y-1;

}
void bcHan( double* x, double* res, double t, const double* parameters){

	res[0] = 0;
	res[1] = 0;
}
// ####################################
// ######################
// Paper
// ######################
void rhsPaperExp( double* p, double* res, const double* parameters){

	double x = 2*M_PI*p[0];
	double y = 2*M_PI*p[1];

	double R = 1000;

	double lambda = R/2. - sqrt(pow(R,2)/4.+4.*M_PI*M_PI);

	res[0] = -1*(-lambda*lambda*exp(lambda*p[0]) *cos(y) + exp(lambda*p[0])*4*M_PI*M_PI*cos(y)) +lambda*exp(2*lambda*p[0]);
    res[1] = -1*(pow(lambda,3)/(2*M_PI)*exp(lambda*p[0]) *sin(y)-lambda*exp(lambda*p[0])*2*M_PI*sin(y));

}


void exactSolutionPaperExp( double* p, double* res){

	double x = 2*M_PI*p[0];
	double y = 2*M_PI*p[1];


	double R = 1000;

	double lambda = R/2. - sqrt(pow(R,2)/4.+4.*M_PI*M_PI);

	res[0] =  1-exp(lambda*p[0])*cos(y);
	res[1] =  lambda/(2*M_PI) * exp(lambda*p[0])*sin(y);
}

void exactSolutionPaperExpP( double* p, double* res){

	double R = 10;

	double lambda = R/2. - sqrt(pow(R,2)/4.+4.*M_PI*M_PI);

	res[0] =  0.5*exp(2*lambda*p[0]);

}


void bcPaperExp( double* p, double* res, double t, const double* parameters){

	double x = 2*M_PI*p[0];
	double y = 2*M_PI*p[1];

	double R = 1000;

	double lambda = R/2. - sqrt(pow(R,2)/4.+4.*M_PI*M_PI);

	res[0] =  1-exp(lambda*p[0])*cos(y);
	res[1] =  lambda/(2*M_PI) * exp(lambda*p[0])*sin(y);
}
// ####################################

void rhs0( double* p, double* res, const double* parameters){

	res[0] =0;
    res[1] =0;
	res[2] =0;
}

void one(double* x, double* res, double t, const double* parameters){
    
    res[0] = 1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void two(double* x, double* res, double t, const double* parameters){
    
    res[0] = 2.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void three(double* x, double* res, double t, const double* parameters){
    
    res[0] = 3.;
    res[1] = 0.;
	res[2] = 0.;
    
    return;
}
void four(double* x, double* res, double t, const double* parameters){
    
    res[0] =parameters[0]*1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void five(double* x, double* res, double t, const double* parameters){
    
    res[0] = parameters[0]*1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}
void six(double* x, double* res, double t, const double* parameters){
    
    res[0] = 6.;
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
    
    return ;
}

void dummyFuncSol(double* x, double* res){
    
	res[0] = 0.;

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

	typedef boost::function<void(double* x, double* res, double t, const double* parameters)>   BCFunc_Type;  

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

		string modellProblem = parameterListProblem->sublist("Mesh Refinement").get("Modell Problem","Seminar1");

		std::vector<double> parameter_vec(0);
		parameter_vec.push_back(maxVel);//height of inflow region


		// Parameters determined by Modell Problem
		RhsFunc_Type rhs;
		Func_Type exactSolU;
		Func_Type exactSolP;
		BCFunc_Type flag1Func;
		BCFunc_Type flag2Func;
		BCFunc_Type flag3Func;
		BCFunc_Type flag4Func;
		BCFunc_Type flag5Func;
		BCFunc_Type flag6Func;
		BCFunc_Type flag7Func;

		if(modellProblem == "Seminar1" && dim ==2){
			rhs = rhsPaper1;
			exactSolU= exactSolutionPaperU1;
			exactSolP = exactSolutionPaperP1;
			flag1Func = bcPaper;
		}
		else if(modellProblem == "Seminar2" && dim == 2){
			rhs = rhsPaper2;
			exactSolU= exactSolutionPaperU1;
			exactSolP = exactSolutionPaperP2;
			flag1Func = bcPaper;
		}
		else if(modellProblem == "Paper" && dim == 2){
			rhs = rhsPaperExp;
			exactSolU= exactSolutionPaperExp;
			exactSolP = exactSolutionPaperExpP;
			flag1Func = bcPaperExp;
		}
		else if(modellProblem == "Verfuerth" && dim == 2){
			rhs = rhsVer;
			exactSolU= exactSolutionUVer;
			exactSolP= exactSolutionPVer;
			flag1Func = bcVer;
			flag5Func = bcSolutionPVer;
		}
		else if(modellProblem == "Han" && dim == 2){
			rhs = rhsHan;
			exactSolU= exactSolutionUHan;
			exactSolP= exactSolutionPHan;
			flag1Func = bcHan;
		}
		else if(modellProblem == "PaperSin" && dim == 2){
			rhs = rhsPaperSin;
			exactSolU= exactSolutionPaperSinU;
			exactSolP = exactSolutionPaperSinP;
			flag1Func = bcPaperSin;
		}
		else if(modellProblem == "LDC" && dim ==2){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet2D;
			flag2Func = zeroDirichlet2D; //two;
			flag3Func = zeroDirichlet2D; //three;
			flag4Func = four;
			flag5Func = five;
			//flag6Func = six;

		}
		else if(modellProblem == "LDC" && dim ==3){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet3D;
			flag2Func = zeroDirichlet2D; //two;
			flag3Func = zeroDirichlet2D; //three;
			flag4Func = four;
			flag5Func = five;
			//flag6Func = six;
		}
		else if(modellProblem == "BFS" && dim ==2){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet2D;
			flag2Func = inflowParabolic2D;
			parameter_vec.push_back(1);//height of inflow region	
		}
		else if(modellProblem == "BFS" && dim == 3){
			rhs = rhs0;
			exactSolU = dummyFuncSol;
			exactSolP = dummyFuncSol;
			flag1Func = zeroDirichlet3D;
			flag2Func = inflowParabolic3D;
			parameter_vec.push_back(1);//height of inflow region				
		}


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

		AdaptiveMeshRefinement<SC,LO,GO,NO> meshRefiner("Stokes",parameterListProblem,exactSolU,exactSolP); 
		
		

		int j=0;
		vec_int_Type iterations(0);
		MAIN_TIMER_START(Total," Step 4:	 Total RefinementAlgorithm");
		while(j<maxIter+1 ){

			MAIN_TIMER_START(buildP2," Step 0:	 buildP2Mesh");
			if (discVelocity=="P2" ) {
				domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
		        domainVelocity->buildP2ofP1Domain( domainPressure );
				}
		    else
		        domainVelocity = domainPressure;
			
			domainVelocity->exportNodeFlags();
			domainPressure->preProcessMesh(true,true); // Preprocessing pressure mesh
			domainVelocity->preProcessMesh(true,true); // Preprocessing velocity mesh

   			domainVelocity->exportSurfaceNormals("domain"); // exporting to check if correct
    		domainVelocity->exportElementOrientation("domain"); // exporting to check if correct

			MAIN_TIMER_STOP(buildP2);		

			MAIN_TIMER_START(Bounds," Step 1:	 bcFactory");
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

			bcFactory->addBC(flag1Func, 1, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			bcFactory->addBC(flag2Func, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			//bcFactory->addBC(flag3Func, 3, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			bcFactory->addBC(flag4Func, 4, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			bcFactory->addBC(flag5Func, 5, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
			bcFactory->addBC(zeroDirichlet2D, 7, 1, domainPressure, "Dirichlet", dim, parameter_vec);

      
			MAIN_TIMER_STOP(Bounds);	
			MAIN_TIMER_START(Solver," Step 2:	 solving PDE");

			
            Teuchos::RCP<Stokes<SC,LO,GO,NO> > stokes( new Stokes<SC,LO,GO,NO>(domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll ));

			{
				Teuchos::TimeMonitor solveTimeMonitor(*solveTime);
				
				stokes->addBoundaries(bcFactory);
				stokes->addRhsFunction(rhs);						    
				stokes->initializeProblem();						    
				stokes->assemble();
				stokes->setBoundaries();          
				stokes->setBoundariesRHS();                
				int iter = stokes->solve();
				iterations.push_back(iter);

			}
			MAIN_TIMER_STOP(Solver);	

						


			MAIN_TIMER_START(Refinement," Step 3:	 meshRefinement");

			// Refinement
			domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );

			{

				ProblemPtr_Type problem = Teuchos::rcp_dynamic_cast<Problem_Type>( stokes , true);
				domainRefined = meshRefiner.globalAlgorithm( domainPressure,  domainVelocity, stokes->getSolution(), problem, rhs );
			}

			domainPressure = domainRefined;
			domainVelocity = domainPressure;
			
			
			//domainRefined->exportNodeFlags("Pressure");

			j++;
			MAIN_TIMER_STOP(Refinement);	
        
        // ####################
       
            
        }
		if(verbose)
			for(int i=0; i< iterations.size(); i++)
				cout << " Iteration in step " << i << ": " << iterations[i] << endl;

		MAIN_TIMER_STOP(Total);	
		Teuchos::TimeMonitor::report(cout,"Main");
   
    }

    return(EXIT_SUCCESS);
}
