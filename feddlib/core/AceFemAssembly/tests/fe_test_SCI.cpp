#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/AceFemAssembly/TestFE/FE_Test.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/LinearAlgebra/BlockMultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"

// #######################################################
// Execute programm: make, then ./core_fe_test_linElas.exe 
// #######################################################


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

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemDeformDiffu.xml");
    Teuchos::CommandLineProcessor myCLP;
    int dim = params->sublist("Parameter").get("Dimension",3);
    myCLP.setOption("dim",&dim,"dim");
    int m = 2;
    myCLP.setOption("m",&m,"H/h");
    
    string filename = params->sublist("Mesh Partitioner").get("Mesh 1 Name","cube.mesh");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    

    // Mesh
    std::string FEType=params->sublist("Parameter").get("Discretization","P2");
    double poissonRatio=params->sublist("Parameter").get("Poisson Ratio",0.4e-0);

    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstanten \lambda
    double youngModulus=params->sublist("Parameter").get("E",1000.0);
    double solConst=params->sublist("Parameter").get("Solution",1.0);
    double lambda = (poissonRatio*youngModulus)/((1. + poissonRatio)*(1. - 2.*poissonRatio));
    double mu = youngModulus/(2.*(1.+poissonRatio));

    int dofs =params->sublist("Parameter").get("Dofs",3);
    int numProcsCoarseSolve = 0;
    int n;
    int size = comm->getSize();

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

	// Building a P2 Domain if requested
    if (FEType == "P2") {
        domainP2->buildP2ofP1Domain( domainP1 );
        domain = domainP2;
    }
    else
        domain = domainP1;

	// Class for assembling linear Elasticity in FEDDLib
    FE<SC,LO,GO,NO> fe;
    fe.addFE(domain);
    fe.addFE(domain);
     
    MultiVectorPtr_Type d_rep = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldRepeated(), 1 ) ); // d_vec
	d_rep->putScalar(solConst);
	MultiVectorPtr_Type c_rep = Teuchos::rcp( new MultiVector_Type( domain->getMapRepeated(), 1 ) ); // c_vec
	c_rep->putScalar(solConst);
	
	MatrixPtr_Type A(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() ) );                     
    MatrixPtr_Type B(new Matrix_Type( domain->getMapUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() ) );
    MatrixPtr_Type BT(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow() ) );
    MatrixPtr_Type C(new Matrix_Type( domain->getMapUnique(),domain->getDimension() * domain->getApproxEntriesPerRow() ));

            // For implicit the system is ordered differently with solid block in 0,0 and diffusion in 1,1
    BlockMatrixPtr_Type system = Teuchos::rcp( new BlockMatrix_Type( 2 ) ); // 
    
	system->addBlock(A,0,0);
	system->addBlock(BT,0,1);
	system->addBlock(B,1,0);
	system->addBlock(C,1,1);

    BlockMultiVectorPtr_Type resVec = Teuchos::rcp( new BlockMultiVector_Type( 2 ) ); // RHS vector
    MultiVectorPtr_Type d_res = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldUnique(), 1 ) ); // d_vec
    MultiVectorPtr_Type c_res = Teuchos::rcp( new MultiVector_Type( domain->getMapUnique(), 1 ) ); // d_vec
    resVec->addBlock(d_res,0);
    resVec->addBlock(c_res,1);
        
    cout << " ACEFem Implementation .... " << endl;
    {
        fe.assemblyAceDeformDiffu(dim,  FEType, FEType, 2 ,  1 , 3, c_rep, d_rep, system, resVec, params, "Jacobian" , true);
        fe.assemblyAceDeformDiffu(dim,  FEType, FEType, 2 ,  1 , 3, c_rep, d_rep, system, resVec, params, "Rhs" , true);
                                                  
    }
    cout << " ... done " << endl;
   
  
		    


    return(EXIT_SUCCESS);
}
