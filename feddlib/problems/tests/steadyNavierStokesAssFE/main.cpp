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

void couette2D(double* x, double* res, double t, const double* parameters){

    res[0] = parameters[0]; // da in parameters 0 maxvelocity
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
    res[0] = 4.*parameters[0]*x[1]*(H-x[1])/(H*H);
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

      //  bool        newtonian       = parameterListProblem->sublist("Parameter").get("Newtonian","true");
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
                domainPressure->buildMesh( 2,"BFS", dim, discPressure, n, m, numProcsCoarseSolve);
                domainVelocity->buildMesh( 2,"BFS", dim, discVelocity, n, m, numProcsCoarseSolve);
            }
	        else if (!meshType.compare("unstructured")) {
	        
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
	        }
        

      
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("MaxVelocity",1.));

            // ####################
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
			if (!bcType.compare("parabolic"))
                parameter_vec.push_back(1.);//height of inflow region
            else if(!bcType.compare("parabolic_benchmark") || !bcType.compare("partialCFD"))
                parameter_vec.push_back(.41);//height of inflow region
            else if(!bcType.compare("Richter3D"))
                parameter_vec.push_back(.4);
            else
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Select a valid boundary condition.");


            if ( !bcType.compare("parabolic") || !bcType.compare("parabolic_benchmark") ) {//flag of obstacle
                if (dim==2){
                    bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                        bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                        bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                }
                else if (dim==3){
                    bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim);
                    bcFactory->addBC(inflowParabolic3D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec);
//                        bcFactory->addBC(dummyFunc, 3, 0, domainVelocity, "Neumann", dim);
//                        bcFactory->addBC(dummyFunc, 666, 1, domainPressure, "Neumann", 1);
                    bcFactory->addBC(zeroDirichlet3D, 4, 0, domainVelocity, "Dirichlet", dim);
                    
                }
            }
            else if (!bcType.compare("partialCFD")) {
                bcFactory->addBC(zeroDirichlet2D, 1, 0, domainVelocity, "Dirichlet", dim); // wall
                bcFactory->addBC(inflowParabolic2D, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // inflow
                bcFactory->addBC(zeroDirichlet2D, 4, 0, domainVelocity, "Dirichlet", dim);
                bcFactory->addBC(zeroDirichlet2D, 5, 0, domainVelocity, "Dirichlet", dim);
            }
            
            else if (!bcType.compare("Richter3D")) {
                bcFactory->addBC(zeroDirichlet3D, 1, 0, domainVelocity, "Dirichlet", dim); // wall
                bcFactory->addBC(inflow3DRichter, 2, 0, domainVelocity, "Dirichlet", dim, parameter_vec); // inflow
                bcFactory->addBC(zeroDirichlet3D, 3, 0, domainVelocity, "Dirichlet_Z", dim);
                bcFactory->addBC(zeroDirichlet3D, 5, 0, domainVelocity, "Dirichlet", dim);
            }
         
             // Flag Check
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

			Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domainVelocity->getMapUnique()));
			vec_int_ptr_Type BCFlags = domainVelocity->getBCFlagUnique();

			Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
			for(int i=0; i< entries.size(); i++){
				entries[i] = BCFlags->at(i);
			}

			Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

			exPara->setup("FlagsFluid",domainVelocity->getMesh(), discVelocity);

			exPara->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domainVelocity->getMapUnique(), domainVelocity->getMapUniqueP2());

			exPara->save(0.0);
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
			// New Assembly Routine
       // /*

        if (verbose) {
        cout << "###############################################################" <<endl;
        cout << "##################### Start NEW ASSEMBLY ROUTINE ####################" <<endl;
        cout << "###############################################################" <<endl;
         }

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
 //*/
			Teuchos::TimeMonitor::report(cout,"Main");	
       
			Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaVelocity(new ExporterParaView<SC,LO,GO,NO>());
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaPressure(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionV = navierStokes.getSolution()->getBlock(0);
            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionP = navierStokes.getSolution()->getBlock(1);
///*
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
    		errorValues->normInf(norm);//const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms);
			double res = norm[0];
			if(comm->getRank() ==0)
				cout << " Inf Norm of Error of Solutions " << res << endl;
			double twoNormError = res;

			navierStokes.getSolution()->norm2(norm);
			res = norm[0];
			if(comm->getRank() ==0)
				cout << " 2 rel. Norm of solution navier stokes " << twoNormError/res << endl;
            double relativeError=twoNormError/res;

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
				cout << " Inf Norm of Difference between Block A: " << res << endl;

            if(discVelocity=="P1"){
            MatrixPtr_Type Sum3= Teuchos::rcp(new Matrix_Type( domainPressure->getMapUnique(), domainPressure->getDimension() * domainPressure->getApproxEntriesPerRow() )  );
			navierStokes.getSystem()->getBlock(1,1)->addMatrix(1, Sum3, 1);
			navierStokesAssFE.getSystem()->getBlock(1,1)->addMatrix(-1, Sum3, 1);
			
			res=0.;
			for (UN i=0; i < domainPressure->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
				
                GO row = domainPressure->getMapUnique()->getGlobalElement( i );
                Sum3->getGlobalRowView(row, indices,values);
                
                for(int j=0; j< values.size() ; j++){
                    if(fabs(values[j])>res)
                        res = fabs(values[j]);			
                }	
            
			}
			res = fabs(res);
			reduceAll<int, double> (*comm, REDUCE_MAX, res, outArg (res));
           
			if(comm->getRank() == 0)
				cout << " Inf Norm of Difference between Block C: " << res << endl;
            }
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

            MatrixPtr_Type Sum4= Teuchos::rcp(new Matrix_Type( domainVelocity->getMapVecFieldUnique(), domainVelocity->getDimension() * domainVelocity->getApproxEntriesPerRow() )  );
			navierStokes.getSystem()->getBlock(0,1)->addMatrix(1, Sum4, 1);
			navierStokesAssFE.getSystem()->getBlock(0,1)->addMatrix(-1, Sum4, 1);

			res=0.;
			for (UN i=0; i < domainVelocity->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
				for(int d=0; d< dim ; d++){
					GO row = dim*domainVelocity->getMapUnique()->getGlobalElement( i )+d;
					Sum4->getGlobalRowView(row, indices,values);
				
                    for(int j=0; j< values.size() ; j++){
                        res += fabs(values[j]);			
                    }	
                }
			}	
			
			res = fabs(res);
			reduceAll<int, double> (*comm, REDUCE_SUM, res, outArg (res));
		
			if(comm->getRank() == 0)
				cout << " Norm of Difference between Block BT: " << res << endl;
		  			
//*/
            DomainPtr_Type dom = domainVelocity;

            exParaVelocity->setup("velocity", dom->getMesh(), dom->getFEType());
                                
            UN dofsPerNode = dim;
            exParaVelocity->addVariable(exportSolutionV, "u", "Vector", dofsPerNode, dom->getMapUnique());
///*       
            exParaVelocity->addVariable(exportSolutionVAssFE, "uAssFE", "Vector", dofsPerNode, dom->getMapUnique());
            exParaVelocity->addVariable(errorValuesAbs, "u-uAssFE", "Vector", dofsPerNode, dom->getMapUnique());
//*/
            dom = domainPressure;
            exParaPressure->setup("pressure", dom->getMesh(), dom->getFEType());

            exParaPressure->addVariable(exportSolutionP, "p", "Scalar", 1, dom->getMapUnique());
 //*/
            exParaPressure->addVariable(exportSolutionPAssFE, "pAssFE", "Scalar", 1, dom->getMapUnique());


            exParaVelocity->save(0.0);
            exParaPressure->save(0.0); 


           TEUCHOS_TEST_FOR_EXCEPTION( relativeError > 1e-12 , std::logic_error, "Relative Error between calculated solutions is too great. Exceeded 1e-12. ");

			
            
        }
    }

    return(EXIT_SUCCESS);
}
