#ifndef MAIN_TIMER_START
#define MAIN_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Main") + std::string(S))));
#endif

#ifndef MAIN_TIMER_STOP
#define MAIN_TIMER_STOP(A) A.reset();
#endif


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"


void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}

using namespace std;
using namespace Teuchos;
using namespace FEDD;

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

int main(int argc, char *argv[]) {

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

 	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

  
 	typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;


    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type ;
    typedef RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    //myCLP.setOption("dim",&dim,"dim");
    //int m = 2;
    //myCLP.setOption("m",&m,"H/h");
    
    // Mesh
	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemNavierStokes.xml");
    string filename = params->sublist("Parameter").get("Mesh Name","square.mesh");
    std::string FETypeV=params->sublist("Parameter").get("DiscretizationV","P2");
    std::string FETypeP=params->sublist("Parameter").get("DiscretizationP","P1");
    int dofsV =params->sublist("Parameter").get("DofsV",2);
    int dofsP =params->sublist("Parameter").get("DofsP",1);
    int dim = params->sublist("Parameter").get("Dimension",3);
    int numProcsCoarseSolve = 0;
    int n;
    int size = comm->getSize();


    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    

	DomainPtr_Type domainP1;
    DomainPtr_Type domainP2;
    DomainPtr_Type domain;
    
    ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
    pListPartitioner->set( "Mesh 1 Name", filename );
    
    domainP1.reset( new Domain_Type( comm, dim ) );
    domainP2.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
    domainP1Array[0] = domainP1;

    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
    
    partitionerP1.readAndPartition();

    if (FETypeV == "P2") {
        domainP2->buildP2ofP1Domain( domainP1 );
        domain = domainP2;
    }
    else
        domain = domainP1;

	// Here we check the different blocks of the NavierStokes system.

    FE<SC,LO,GO,NO> fe;
 	FE<SC,LO,GO,NO> fe_test;
    fe.addFE(domain);
    fe.addFE(domainP1);
	fe.doSetZeros(std::pow(10,-13));
	// Solution
	MultiVectorPtr_Type u_rep = Teuchos::rcp(new MultiVector_Type(domain->getMapVecFieldRepeated(),1));
	MultiVectorPtr_Type p_rep = Teuchos::rcp(new MultiVector_Type(domain->getMapRepeated(),1));
	u_rep->putScalar(100.);
	p_rep->putScalar(0.);

	// Checking first Block:

	MAIN_TIMER_START(FE," FE Routine:   Assemble System");
	BlockMatrixPtr_Type systemFE= Teuchos::rcp(new BlockMatrix_Type(2 ) );  
	MatrixPtr_Type ANW = Teuchos::rcp(new Matrix_Type(domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() ) );  

	MatrixPtr_Type A = Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) );
    fe.assemblyLaplaceVecField(dim, FETypeV, 2, A, true/*call fillComplete*/);

	MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type(domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() ) );      
    fe.assemblyAdvectionVecField( dim,"P2", N, u_rep, true );

	N->resumeFill();
	//N->scale(density);
	N->fillComplete(domain->getMapVecFieldUnique(), domain->getMapVecFieldUnique());

	A->addMatrix(1.,ANW,0.); // ANW = A
	N->addMatrix(1.,ANW,1.); // ANW = ANW + N

	MatrixPtr_Type W = Teuchos::rcp(new Matrix_Type(domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() ) );
    fe.assemblyAdvectionInUVecField( dim,"P2", W, u_rep, true );

	W->resumeFill();
    //W->scale(density);
    W->fillComplete( domain->getMapVecFieldUnique(), domain->getMapVecFieldUnique());
    W->addMatrix(1.,ANW,1.);

	ANW->fillComplete(domain->getMapVecFieldUnique(), domain->getMapVecFieldUnique());

	systemFE->addBlock(ANW,0,0);
	MatrixPtr_Type BT= Teuchos::rcp(new Matrix_Type(domain->getMapVecFieldUnique(), domainP1->getDimension() * domainP1->getApproxEntriesPerRow() )  );
    MatrixPtr_Type B= Teuchos::rcp(new Matrix_Type(domainP1->getMapUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() )  );

    fe.assemblyDivAndDivTFast(dim, "P2", "P1", 2, B, BT, domain->getMapVecFieldUnique(), domainP1->getMapUnique(), true );

	B->resumeFill();
	B->scale(-1.);	
	B->fillComplete();
	BT->resumeFill();
	BT->scale(-1.);
	BT->fillComplete();

	MAIN_TIMER_STOP(FE);	
	cout << " Done for FE routine " << endl;
	// ANW is first block 
	// --------------------------------------------------------------
	MAIN_TIMER_START(FE_test," Assemble_FE_Routine: Assemble System");

    fe_test.addFE(domain);
    fe_test.addFE(domainP1);
    BlockMatrixPtr_Type systemFETest= Teuchos::rcp(new BlockMatrix_Type(2 ) );  

  	MatrixPtr_Type A_test= Teuchos::rcp(new Matrix_Type(domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() )  );
    MatrixPtr_Type BT_test= Teuchos::rcp(new Matrix_Type(domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() )  );
    MatrixPtr_Type B_test= Teuchos::rcp(new Matrix_Type(domainP1->getMapUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() )  );
	MatrixPtr_Type dummy = Teuchos::rcp( new Matrix_Type( domainP1->getMapUnique(), 1 ) );
	
	systemFETest->addBlock(A_test,0,0);
	systemFETest->addBlock(BT_test,0,1);
	systemFETest->addBlock(B_test,1,0);
	systemFETest->addBlock(dummy,1,1);

	BlockMultiVectorPtr_Type resVec = Teuchos::rcp( new BlockMultiVector_Type( 2 ) ); // RHS vector
    MultiVectorPtr_Type u_res = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldUnique(), 1 ) ); // d_vec
    MultiVectorPtr_Type p_res = Teuchos::rcp( new MultiVector_Type( domainP1->getMapUnique(), 1 ) ); // d_vec
    resVec->addBlock(u_res,0);
    resVec->addBlock(p_res,1);
    
    SmallMatrix<SC> defTS(2);
    defTS[0][0] = 1;
    defTS[0][1] = 1;
    defTS[1][0] = 0;
    defTS[1][1] = 0;


    fe_test.assemblyNavierStokes(dim, FETypeV, FETypeP, 2,dofsV,dofsP,u_rep, p_rep, systemFETest, resVec, defTS, params, false, "Jacobian", true/*call fillComplete*/);
    
	//B->print();
	//B_test->print();

	cout << " Done for FE Test" << endl;
	MAIN_TIMER_STOP(FE_test);	
	Teuchos::TimeMonitor::report(cout,"Main");
	MatrixPtr_Type Sum= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() )  );
	ANW->addMatrix(1, Sum, 0);
	A_test->addMatrix(-1, Sum, 1);

	int maxRank = std::get<1>(domain->getMesh()->rankRange_);

	double res=0.;
	Teuchos::ArrayView<const GO> indices;
	Teuchos::ArrayView<const SC> values;

	for (UN i=0; i < domain->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
		for(int d=0; d< dofsV ; d++){
			GO row = dofsV*domain->getMapUnique()->getGlobalElement( i )+d;
			Sum->getGlobalRowView(row, indices,values);
			
			for(int j=0; j< values.size() ; j++){
				if(std::fabs(values[j])>res)
					res = std::fabs(values[j]);			
			}	
		}	
	}
	Sum->fillComplete();
	//Sum->print();
	res = std::fabs(res);
	reduceAll<int, double> (*comm, REDUCE_SUM, res, outArg (res));

	if(comm->getRank() == 0)
		cout << " Max Norm of Difference between Block A: " << res << endl;
	

	MatrixPtr_Type Sum2= Teuchos::rcp(new Matrix_Type( domainP1->getMapUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() )  );
	B->addMatrix(1, Sum2, 0);
	B_test->addMatrix(-1, Sum2, 1);

	res=0.;
	for (UN i=0; i < domainP1->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
		for(int d=0; d< dofsP ; d++){
			GO row = dofsP*domainP1->getMapUnique()->getGlobalElement( i )+d;
			Sum2->getGlobalRowView(row, indices,values);
			
			for(int j=0; j< values.size() ; j++){
				res += std::fabs(values[j]);			
			}	
		}	
	}
	res = std::fabs(res);
	reduceAll<int, double> (*comm, REDUCE_SUM, res, outArg (res));
	//Sum2->fillComplete();
	//Sum2->print();
	if(comm->getRank() == 0)
		cout << " Norm of Difference between Block B: " << res << endl;



    return(EXIT_SUCCESS);
}
