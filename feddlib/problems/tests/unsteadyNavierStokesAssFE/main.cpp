#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include "feddlib/problems/specific/NavierStokesAssFE.hpp"
#include "feddlib/core/General/BCBuilder.hpp"

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

using namespace Teuchos;
using namespace FEDD;
int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

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
        return EXIT_SUCCESS;
    }

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);

        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);

        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListPrecTeko = Teuchos::getParametersFromXmlFile(xmlTekoPrecFile);
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",3);
        std::string discVelocity = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        std::string discPressure = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");
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
                              
            domainPressure.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
            domainVelocity.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
            
            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
            domainP1Array[0] = domainPressure;
            
            ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
            MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
            
            partitionerP1.readAndPartition();

            if (discVelocity=="P2")
                domainVelocity->buildP2ofP1Domain( domainPressure );
            else
                domainVelocity = domainPressure;
      
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("MaxVelocity",1.));

            // ####################
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );


          	 parameter_vec.push_back(0.41);//height of inflow region

            if (dim==2){
                bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                       bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                        bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                //bcFactory->addBC(zeroDirichlet2D, 3, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(zeroDirichlet2D, 5, 0, domainVelocity, "Dirichlet", dim);
            }
            else if (dim==3){
                bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                        bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                        bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                bcFactory->addBC(zeroDirichlet3D, 2, 0, domainVelocity, "Dirichlet", dim);
                
            }

            

            
//          
                     
            int timeDisc = parameterListProblem->sublist("Timestepping Parameter").get("Butcher table",0);

			DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);
            SmallMatrix<int> defTS(2);
            defTS[0][0] = 1;
            defTS[0][1] = 1;
            defTS[1][0] = 0;
            defTS[1][1] = 0;

            daeTimeSolver.defineTimeStepping(defTS);

			// ###########################################################################################################
			// OLD ASSEMBLY
			MAIN_TIMER_START(FE," Old: Solve equation");
            NavierStokes<SC,LO,GO,NO> navierStokes( domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll );

            navierStokes.addBoundaries(bcFactory);
            
            navierStokes.initializeProblem();
            
            navierStokes.assemble();

            navierStokes.setBoundariesRHS();

            daeTimeSolver.setProblem(navierStokes);

            daeTimeSolver.setupTimeStepping();

            daeTimeSolver.advanceInTime();
			MAIN_TIMER_STOP(FE);	

			// ###########################################################################################################

			DAESolverInTime<SC,LO,GO,NO> daeTimeSolverAssFE(parameterListAll, comm);


            daeTimeSolverAssFE.defineTimeStepping(defTS);

			// ###########################################################################################################
			// New Assembly
			MAIN_TIMER_START(FE_test," New: Solve equation");

  			NavierStokesAssFE<SC,LO,GO,NO> navierStokesAssFE( domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll );
			navierStokesAssFE.addBoundaries(bcFactory);
            
            navierStokesAssFE.initializeProblem();
            
            navierStokesAssFE.assemble();

            navierStokesAssFE.setBoundariesRHS();

            daeTimeSolverAssFE.setProblem(navierStokesAssFE);

            daeTimeSolverAssFE.setupTimeStepping();

            daeTimeSolverAssFE.advanceInTime();
			MAIN_TIMER_STOP(FE_test);	
			Teuchos::TimeMonitor::report(cout,"Main");
			// ###########################################################################################################


			Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaVelocity(new ExporterParaView<SC,LO,GO,NO>());
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaPressure(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionV = navierStokes.getSolution()->getBlock(0);
            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionP = navierStokes.getSolution()->getBlock(1);

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionVAssFE = navierStokesAssFE.getSolution()->getBlock(0);
            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionPAssFE = navierStokesAssFE.getSolution()->getBlock(1);




			// Calculating the error per node
			Teuchos::RCP<MultiVector<SC,LO,GO,NO> > errorValues = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>( navierStokes.getSolution()->getBlock(0)->getMap() ) ); 
			//this = alpha*A + beta*B + gamma*this
			errorValues->update( 1., exportSolutionV, -1. ,exportSolutionVAssFE, 0.);

			// Taking abs norm
			Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > errorValuesAbs = errorValues;

			errorValues->abs(errorValuesAbs);

 			Teuchos::Array<SC> norm(1); 
    		errorValues->norm2(norm);//const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms);
			double res = norm[0];
			if(comm->getRank() ==0)
				cout << " 2 Norm of Error of Solution Velocity " << res << endl;
			double NormError = res;
		
			navierStokes.getSolution()->norm2(norm);
			res = norm[0];
			if(comm->getRank() ==0)
				cout << " 2 rel. Norm to solution navier stokes " << NormError/res << endl;

			navierStokesAssFE.getSolution()->norm2(norm);
			res = norm[0];
			if(comm->getRank() ==0)
				cout << " 2 rel. Norm to solutions navier stokes assemFE " << NormError/res << endl;

            DomainPtr_Type dom = domainVelocity;

            exParaVelocity->setup("velocity", dom->getMesh(), dom->getFEType());
                                
            UN dofsPerNode = dim;
            exParaVelocity->addVariable(exportSolutionV, "u", "Vector", dofsPerNode, dom->getMapUnique());
            exParaVelocity->addVariable(exportSolutionVAssFE, "uAssFE", "Vector", dofsPerNode, dom->getMapUnique());
            exParaVelocity->addVariable(errorValuesAbs, "u-uAssFE", "Vector", dofsPerNode, dom->getMapUnique());

            dom = domainPressure;
            exParaPressure->setup("pressure", dom->getMesh(), dom->getFEType());

            exParaPressure->addVariable(exportSolutionP, "p", "Scalar", 1, dom->getMapUnique());
            exParaPressure->addVariable(exportSolutionPAssFE, "pAssFE", "Scalar", 1, dom->getMapUnique());


            exParaVelocity->save(0.0);
            exParaPressure->save(0.0); 	


           //TEUCHOS_TEST_FOR_EXCEPTION( infNormError > 1e-11 , std::logic_error, "Inf Norm of Error between calculated solutions is too great. Exceeded 1e-11. ");


        }
    }

    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
