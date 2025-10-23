#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/amr/AdaptiveMeshRefinement.hpp"

void zeroDirichlet(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void zeroDirichletX(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = x[1];
    res[2] = x[2];

    return;
}

void zeroDirichletY(double* x, double* res, double t, const double* parameters)
{
    res[0] = x[0];
    res[1] = 0.;
    res[2] = x[2];


    return;
}

void zeroDirichletZ(double* x, double* res, double t, const double* parameters)
{
    res[0] = x[0];
    res[1] = x[1];
    res[2] = 0.;

    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
    return;
}

void rhs2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    
    return;
}

void rhsY(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    res[2] = 0.;
    return;
}

void rhsX(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void rhsYZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    double force = parameters[1];

    if(parameters[2] == 5)
        res[1] = force;
    else
        res[1] =0.;
        
    if (parameters[2] == 4)
        res[2] = force;
    else
        res[2] = 0.;
    
    return;
}


typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;
int main(int argc, char *argv[])
{

    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    // int dim = 2;
    // myCLP.setOption("dim",&dim,"dim");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    bool verbose (comm->getRank() == 0); // Print-Ausgaben nur auf rank = 0
    if (verbose) {
        cout << "###############################################################" <<endl;
        cout << "############ Starting Steady Nonlinear Elasticity ... ############" <<endl;
        cout << "###############################################################" <<endl;
    }

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",3);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        string		meshName    	= parameterListProblem->sublist("Parameter").get("Mesh Name","cube_0_1.mesh");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        int         n;
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        string      FEType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        int size = comm->getSize() - numProcsCoarseSolve;

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

        DomainPtr_Type domain;

        // ########################
        // P1 und P2 Gitter bauen
        // ########################

        domain.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
        domainP1Array[0] = domain;
        
        ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
        
        partitionerP1.readAndPartition();
 
            
		// Refining mesh uniformly or area wise
		// ----------
		int level = parameterListProblem->sublist("Parameter").get("Refinement Level",1);
		string type = parameterListProblem->sublist("Parameter").get("Refinement Type", "uniform");
		AdaptiveMeshRefinement<SC,LO,GO,NO> meshRefiner(parameterListProblem); // exactSolLShape
		Teuchos::RCP<Domain<SC,LO,GO,NO> > domainRefined;
		domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
		vec2D_dbl_Type area(dim,vec_dbl_Type(2));
	    area[0][0] = 0.0;
		area[0][1] = 0.15;
		area[1][0] = 0.0;
		area[1][1] = 0.15;
		if(dim==3){
			area[2][0] = 0.0;
			area[2][1] = 0.15;
		}
		
		if(level>0){
			if(type=="area"){
				domainRefined = meshRefiner.refineArea(domain,area,level );
				domain=domainRefined;
			}
			else if(type=="uniform"){
				domainRefined = meshRefiner.refineUniform(domain,level );
				domain=domainRefined;
			}
				
		}
		
		
       if (FEType=="P2") {
            Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP2;
            domainP2.reset( new Domain_Type( comm, dim ) );
            domainP2->buildP2ofP1Domain( domain );
            domain = domainP2;
        }
        

       // ########################
        // Flags setzen
        // ########################

		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaF(new ExporterParaView<SC,LO,GO,NO>());

		Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domain->getMapUnique()));
		vec_int_ptr_Type BCFlags = domain->getBCFlagUnique();

		Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
		for(int i=0; i< entries.size(); i++){
			entries[i] = BCFlags->at(i);
		}

		Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

		exParaF->setup("Flags", domain->getMesh(), FEType);

		exParaF->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domain->getMapUnique());

		exParaF->save(0.0);

		Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
        if (dim == 2)
            bcFactory->addBC(zeroDirichlet2D, 1, 0, domain, "Dirichlet", dim);
        else if (dim == 3){
                      
            bcFactory->addBC(zeroDirichlet, 1, 0, domain, "Dirichlet_X", dim);
            bcFactory->addBC(zeroDirichlet, 2, 0, domain, "Dirichlet_Y", dim);
            bcFactory->addBC(zeroDirichlet, 3, 0, domain, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet3D, 0, 0, domain, "Dirichlet", dim);
            bcFactory->addBC(zeroDirichlet2D, 7, 0, domain, "Dirichlet_X_Y", dim);
            bcFactory->addBC(zeroDirichlet2D, 8, 0, domain, "Dirichlet_Y_Z", dim);
            bcFactory->addBC(zeroDirichlet2D, 9, 0, domain, "Dirichlet_X_Z", dim);
        
		}

        // LinElas Objekt erstellen
        NonLinElasticity<SC,LO,GO,NO> NonLinElasAssFE( domain, FEType, parameterListAll );
        NonLinElasAssFE.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        if (dim==2)
            NonLinElasAssFE.addRhsFunction( rhs2D );
        else if(dim==3)
            NonLinElasAssFE.addRhsFunction( rhsYZ );
        
        double force = parameterListAll->sublist("Parameter").get("Volume force",0.);
        double degree = 0;      
  
        NonLinElasAssFE.addParemeterRhs( force );
        NonLinElasAssFE.addParemeterRhs( degree );
        // ######################
        // Matrix assemblieren, RW setzen und System loesen
        // ######################
        NonLinElasAssFE.initializeProblem();
        NonLinElasAssFE.assemble();                
        NonLinElasAssFE.setBoundaries(); // In der Klasse Problem
        NonLinElasAssFE.setBoundariesRHS();

		std::string nlSolverType = parameterListProblem->sublist("General").get("Linearization","FixedPoint");
        NonLinearSolver<SC,LO,GO,NO> nlSolverAssFE( nlSolverType );
        nlSolverAssFE.solve( NonLinElasAssFE );
        comm->barrier(); 

        // LinElas Objekt erstellen
        parameterListAll->sublist("Parameter").set("Use AceGen Interface",false); // Not use interface
        NonLinElasticity<SC,LO,GO,NO> NonLinElas( domain, FEType, parameterListAll );

        NonLinElas.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        if (dim==2)
            NonLinElas.addRhsFunction( rhs2D );
        else if(dim==3)
            NonLinElas.addRhsFunction( rhsYZ );

        NonLinElas.addParemeterRhs( force );
        NonLinElas.addParemeterRhs( degree );
        // ######################
        // Matrix assemblieren, RW setzen und System loesen
        // ######################
        NonLinElas.initializeProblem();
        NonLinElas.assemble();                
        NonLinElas.setBoundaries(); // In der Klasse Problem
		NonLinElas.setBoundariesRHS();

        NonLinearSolver<SC,LO,GO,NO> nlSolver( nlSolverType );
        nlSolver.solve( NonLinElas );
        comm->barrier();	


        
        
		if(comm->getRank() ==0){
			cout << " ############################################### " << endl;
			cout << " Nonlinear Iterations FEDDLib Assembly : " << nlSolver.getNonLinIts() << endl; 
			cout << " Nonlinear Iterations AceGEN Assembly  : " << nlSolverAssFE.getNonLinIts() << endl;
			cout << " ############################################### " << endl;

		}

        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        exPara->setup( "displacements", domain->getMesh(), FEType );

        MultiVectorConstPtr_Type valuesSolidConst1 = NonLinElas.getSolution()->getBlock(0);
        exPara->addVariable( valuesSolidConst1, "valuesNonLinElas", "Vector", dim, domain->getMapUnique());

        MultiVectorConstPtr_Type valuesSolidConst2 = NonLinElasAssFE.getSolution()->getBlock(0);
        exPara->addVariable( valuesSolidConst2, "valuesNonLinElasAssFE", "Vector", dim, domain->getMapUnique());

        // Calculating the error per node
        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > errorValues = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>( valuesSolidConst1->getMap() ) ); 
        //this = alpha*A + beta*B + gamma*this
        errorValues->update( 1., valuesSolidConst2, -1. ,valuesSolidConst1, 0.);

        // Taking abs norm
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > errorValuesAbs = errorValues;

        errorValues->abs(errorValuesAbs);

        exPara->addVariable( errorValuesAbs, "erroeValues", "Vector", dim, domain->getMapUnique());
        exPara->save(0.0);

        Teuchos::Array<SC> norm(1); 
        errorValues->normInf(norm);//const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms);
        double res = norm[0];
        if(comm->getRank() ==0)
            cout << " Inf Norm of Error of Solutions " << res << endl;
        double infNormError = res;
    
        NonLinElas.getSolution()->getBlock(0)->normInf(norm);
        res = norm[0];
        if(comm->getRank() ==0)
            cout << " Relative error Inf-Norm of solution nonlinear elasticity " << infNormError/res << endl;

        NonLinElasAssFE.getSolution()->getBlock(0)->normInf(norm);
        res = norm[0];
        if(comm->getRank() ==0)
            cout << " Relative error Inf-Norm of solutions nonlinear elasticity assemFE " << infNormError/res << endl;
    

          // TEUCHOS_TEST_FOR_EXCEPTION( infNormError > 1e-11 , std::logic_error, "Inf Norm of Error between calculated solutions is too great. Exceeded 1e-11. ");

    }

    return(EXIT_SUCCESS);
}
