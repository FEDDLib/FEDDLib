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

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;  
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    int dim = 2;
    myCLP.setOption("dim",&dim,"dim");
    int m = 2;
    myCLP.setOption("m",&m,"H/h");
    
    string filename = "square.mesh";

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    

    // Mesh
	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemLaplace.xml");
    std::string FEType=params->sublist("Parameter").get("Discretization1","P1");
    int dofs =params->sublist("Parameter").get("Dofs",1);
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

    if (FEType == "P2") {
        domainP2->buildP2ofP1Domain( domainP1 );
        domain = domainP2;
    }
    else
        domain = domainP1;


    std::vector<double> funcParameter(1);

    FE<SC,LO,GO,NO> fe;
    fe.addFE(domain);

    MatrixPtr_Type A= Teuchos::rcp(new Matrix_Type( domain->getMapUnique(), 30  ) ); // with approx entries per row

	if(dofs == dim)
		A.reset(new Matrix_Type( domain->getMapVecFieldUnique(), 30  ) );

    {
		if(dofs == dim)
       		 fe.assemblyLaplaceVecField(dim, FEType, 2, A, true/*call fillComplete*/);
		else
			 fe.assemblyLaplace(dim, FEType, 2, A, true/*call fillComplete*/);

    }
    
 	FE<SC,LO,GO,NO> fe_test;
    fe_test.addFE(domain);
    
    MatrixPtr_Type A_test= Teuchos::rcp(new Matrix_Type( domain->getMapUnique(),domain->getDimension() * domain->getApproxEntriesPerRow()   ) );
    BlockMatrixPtr_Type system = Teuchos::rcp( new BlockMatrix_Type( 2 ) ); // 
    
	if(dofs == dim)
		A_test.reset(new Matrix_Type( domain->getMapVecFieldUnique(), 30  ) );
		
	system->addBlock(A_test,0,0);
    {
        fe_test.assemblyLaplaceAssFE(dim, FEType, 2,dofs, system, true/*call fillComplete*/);
    }
   
	MatrixPtr_Type Sum= Teuchos::rcp(new Matrix_Type(domain->getMapUnique(), 30  ) ); // with approx entries per row
	if(dofs == dim)
		Sum.reset(new Matrix_Type( domain->getMapVecFieldUnique(), 30  ) );
	A->addMatrix(-1, Sum, 1);
	A_test->addMatrix(1, Sum, -1);


	int maxRank = std::get<1>(domain->getMesh()->rankRange_);

	double res=0.;
	Teuchos::ArrayView<const GO> indices;
	Teuchos::ArrayView<const SC> values;

	for (UN i=0; i < domain->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
		for(int d=0; d< dofs ; d++){
			GO row = dofs*domain->getMapUnique()->getGlobalElement( i )+d;
			Sum->getGlobalRowView(row, indices,values);
			
			for(int j=0; j< values.size() ; j++){
				res += fabs(values[j]);			
			}	
		}	
	}
	res = fabs(res);
	reduceAll<int, double> (*comm, REDUCE_SUM, res, outArg (res));

	if(comm->getRank() == 0)
		cout << " Norm of Difference between StiffnessMatrices: " << res << endl;
	

    return(EXIT_SUCCESS);
}
