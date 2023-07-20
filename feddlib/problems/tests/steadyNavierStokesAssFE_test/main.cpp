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

#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/specific/NavierStokesAssFE.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>


/*!
 main of steady-state Navier-Stokes problem

 @brief steady-state Navier-Stokes main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
using namespace Teuchos;
using namespace FEDD;

void zeroDirichlet(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;
    res[1] = 0.;

    return;
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
void four(double* x, double* res, double t, const double* parameters){

    res[0] = 4.;
    res[1] = 4.;

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
void inflow3DRichter(double* x, double* res, double t, const double* parameters)
{

    double H = parameters[1];
    
    res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) );
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

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    bool verbose (comm->getRank() == 0);

//    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    if (verbose) {
        cout << "###############################################################" <<endl;
        cout << "##################### Steady Navier-Stokes ####################" <<endl;
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

        std::string discVelocity = parameterListProblem->sublist("Parameter").get("Discretization Velocity","P2");
        std::string discPressure = parameterListProblem->sublist("Parameter").get("Discretization Pressure","P1");


        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","circle2D_1800.mesh");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        string		linearization	= parameterListProblem->sublist("General").get("Linearization","FixedPoint");
        string		precMethod      = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         mixedFPIts		= parameterListProblem->sublist("General").get("MixedFPIts",1);        
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
        else if(!meshType.compare("structured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

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
			// ---------------------
			// Old Assembly Rouine
        
          	NavierStokes<SC,LO,GO,NO> navierStokes( domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll );

            {
				MAIN_TIMER_START(NavierStokes," Assemble System and solve");
                navierStokes.addBoundaries(bcFactory);
                navierStokes.initializeProblem();
                navierStokes.assemble();

                navierStokes.setBoundariesRHS();

                std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
                NonLinearSolver<SC,LO,GO,NO> nlSolver( nlSolverType );
                nlSolver.solve( navierStokes );
				MAIN_TIMER_STOP(NavierStokes);	
                comm->barrier();
            }

           	// ---------------------
			// New Assembly Rouine
        
            NavierStokesAssFE<SC,LO,GO,NO> navierStokesAssFE( domainVelocity, discVelocity, domainPressure, discPressure, parameterListAll );

            {
				MAIN_TIMER_START(NavierStokesAssFE," AssFE:   Assemble System and solve");
                navierStokesAssFE.addBoundaries(bcFactory);
                navierStokesAssFE.initializeProblem();
                navierStokesAssFE.assemble();

                navierStokesAssFE.setBoundariesRHS();

                std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
                NonLinearSolver<SC,LO,GO,NO> nlSolverAssFE( nlSolverType );
                nlSolverAssFE.solve( navierStokesAssFE );
				MAIN_TIMER_STOP(NavierStokesAssFE);	
                comm->barrier();
            }
 
			Teuchos::TimeMonitor::report(cout,"Main");	
       
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
				cout << " Inf Norm of Error of Solutions " << res << endl;
			double twoNormError = res;

			navierStokes.getSolution()->norm2(norm);
			res = norm[0];
			if(comm->getRank() ==0)
				cout << " 2 rel. Norm of solution navier stokes " << twoNormError/res << endl;

			navierStokesAssFE.getSolution()->norm2(norm);
			res = norm[0];
			if(comm->getRank() ==0)
				cout << " 2 rel. Norm of solutions navier stokes assemFE " << twoNormError/res << endl;

			MatrixPtr_Type Sum2= Teuchos::rcp(new Matrix_Type( domainVelocity->getMapVecFieldUnique(), domainVelocity->getDimension() * domainVelocity->getApproxEntriesPerRow() )  );
			navierStokes.getSystem()->getBlock(0,0)->addMatrix(1, Sum2, 1);
			navierStokesAssFE.getSystem()->getBlock(0,0)->addMatrix(-1, Sum2, 1);

			Teuchos::ArrayView<const GO> indices;
			Teuchos::ArrayView<const SC> values;
			res=0.;
			for (UN i=0; i < domainVelocity->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
				for(int d=0; d< dim ; d++){
					GO row = dim*domainVelocity->getMapUnique()->getGlobalElement( i )+d;
					Sum2->getGlobalRowView(row, indices,values);
					
					for(int j=0; j< values.size() ; j++){
						if(fabs(values[j])>res)
							res = fabs(values[j]);			
					}	
				}	
			}
			res = fabs(res);
			reduceAll<int, double> (*comm, REDUCE_MAX, res, outArg (res));
           
			if(comm->getRank() == 0)
				cout << "Inf Norm of Difference between Block A: " << res << endl;

			MatrixPtr_Type Sum1= Teuchos::rcp(new Matrix_Type( domainPressure->getMapUnique(), domainVelocity->getDimension() * domainVelocity->getApproxEntriesPerRow() )  );
			navierStokes.getSystem()->getBlock(1,0)->addMatrix(1, Sum1, 1);
			navierStokesAssFE.getSystem()->getBlock(1,0)->addMatrix(-1, Sum1, 1);

			res=0.;
			for (UN i=0; i < domainPressure->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
				GO row = domainPressure->getMapUnique()->getGlobalElement( i );
				Sum1->getGlobalRowView(row, indices,values);
				
				for(int j=0; j< values.size() ; j++){
					res += fabs(values[j]);			
				}	
			}	
			
			res = fabs(res);
			reduceAll<int, double> (*comm, REDUCE_SUM, res, outArg (res));
		
			if(comm->getRank() == 0)
				cout << " Norm of Difference between Block B: " << res << endl;

		  			

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

    return(EXIT_SUCCESS);
}
