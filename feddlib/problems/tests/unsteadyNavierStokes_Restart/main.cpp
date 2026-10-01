#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif

#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include "feddlib/problems/specific/NavierStokesAssFE.hpp"
#include <Xpetra_DefaultPlatform.hpp>

/*!
 main of time-dependent Navier-Stokes problem with test for restart function  
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

using namespace Teuchos;
using namespace FEDD;
int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    bool verbose (comm->getRank() == 0);
    if (verbose) {
        std::cout << "###############################################################" <<std::endl;
        std::cout << "################### Unsteady Navier-Stokes ####################" <<std::endl;
        std::cout << "###############################################################" <<std::endl;
    }

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    std::string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    std::string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    std::string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    std::string xmlTekoPrecFile = "parametersTeko.xml";
    myCLP.setOption("tekoprecfile",&xmlTekoPrecFile,".xml file with Inputparameters.");

    double length = 4.;
    myCLP.setOption("length",&length,"length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    bool testPassed = true;
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
        int 		dim				= 3;
        std::string discVelocity = "P2";
        std::string discPressure = "P1";
        std::string		meshType ="unstructured";
        std::string		meshName = "BFS3dCC.mesh";
        std::string		meshDelimiter = " ";
        std::string		linearization = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
        std::string		precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;


        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if (!precMethod.compare("Monolithic"))
            parameterListAll->setParameters(*parameterListPrec);
        else
            parameterListAll->setParameters(*parameterListPrecTeko);

        parameterListAll->setParameters(*parameterListSolver);

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


          	parameter_vec.push_back(1.0);//height of inflow region

        
            bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
            bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
            bcFactory->addBC(zeroDirichlet3D, 2, 0, domainVelocity, "Dirichlet", dim);
  
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
            // Comparing computed solution from restart to the hdf5 saved solution.
            double finalTime = parameterListProblem->sublist("Timestepping Parameter").get("Final time compare", 0.0);
            double dt = parameterListProblem->sublist("Timestepping Parameter").get("dt", 0.01);
            HDF5Import<SC,LO,GO,NO> importer(navierStokes.getSolution()->getBlock(0)->getMap(),restartFile(parameterListProblem, "Solutionu"));
            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > solutionImported = importer.readVariablesHDF5( std::to_string(finalTime));
		
			Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaVelocity(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionV = navierStokes.getSolution()->getBlock(0);

			// Calculating the error per node
			Teuchos::RCP<MultiVector<SC,LO,GO,NO> > errorValues = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>( navierStokes.getSolution()->getBlock(0)->getMap() ) ); 
			//this = alpha*A + beta*B + gamma*this
			errorValues->update( 1., exportSolutionV, -1. ,solutionImported, 0.);

			// Taking abs norm
			Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > errorValuesAbs = errorValues;

			errorValues->abs(errorValuesAbs);

 			Teuchos::Array<SC> norm(1); 
    		errorValues->norm2(norm);//const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms);
			double res = norm[0];
			if(comm->getRank() ==0)
				std::cout << " 2 Norm of Error of Solution Velocity " << res << std::endl;
			double NormError = res;
		
			navierStokes.getSolution()->norm2(norm);
			res = norm[0];
			if(comm->getRank() ==0)
				std::cout << " 2 rel. Norm to solution navier stokes " << NormError/res << std::endl;

			double tolerance = parameterListProblem->sublist("Timestepping Parameter").get("Restart tolerance", 1.e-6);
			if(!(NormError/res <= tolerance))
				testPassed = false;

			navierStokes.getSolution()->norm2(norm);
			res = norm[0];
			if(comm->getRank() ==0)
				std::cout << " 2 rel. Norm to solutions navier stokes assemFE " << NormError/res << std::endl;

            DomainPtr_Type dom = domainVelocity;

            exParaVelocity->setup("velocity", dom->getMesh(), dom->getFEType());
                                
            UN dofsPerNode = dim;
            exParaVelocity->addVariable(exportSolutionV, "u", "Vector", dofsPerNode, dom->getMapUnique());
            exParaVelocity->addVariable(errorValuesAbs, "u-uImport", "Vector", dofsPerNode, dom->getMapUnique());

            exParaVelocity->save(0.0);

           //TEUCHOS_TEST_FOR_EXCEPTION( infNormError > 1e-11 , std::logic_error, "Inf Norm of Error between calculated solutions is too great. Exceeded 1e-11. ");


        }
    }

    Teuchos::TimeMonitor::report(std::cout);

    return testPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}
