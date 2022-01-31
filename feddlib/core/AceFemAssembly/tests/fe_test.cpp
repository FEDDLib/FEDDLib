#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/AceFemAssembly/TestFE/FE_Test.hpp"
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
typedef KokkosClassic::DefaultNode::DefaultNodeType NO;

int main(int argc, char *argv[]) {

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

 	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
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
    
//    Xpetra::UnderlyingLib ulib;
//    if (!ulib_str.compare("UseTpetra"))
//        ulib = Xpetra::UseTpetra;
//    else if (!ulib_str.compare("UseEpetra"))
//        ulib = Xpetra::UseEpetra;

    // Mesh
	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblem.xml");
    std::string FEType=params->sublist("Parameter").get("Discretization1","P1");
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

    MultiVectorPtr_Type a= Teuchos::rcp(new MultiVector_Type(domain->getMapRepeated()  ) ); 

    {
        fe.assemblyLaplace(dim, FEType, 2, A, true/*call fillComplete*/);
		fe.assemblyRHS(dim, FEType, a, "Scalar", oneFunc,funcParameter);
    }
    
 	FE_Test<SC,LO,GO,NO> fe_test;
    fe_test.addFE(domain);
    
    MatrixPtr_Type A_test= Teuchos::rcp(new Matrix_Type(domain->getMapUnique(), 30  ) ); // with approx entries per row

    MultiVectorPtr_Type a_test= Teuchos::rcp(new MultiVector_Type(domain->getMapRepeated()  ) ); 
    {
        fe_test.assemblyLaplace(dim, FEType, 2, A_test, true/*call fillComplete*/);
		fe_test.assemblyRHS(dim, FEType, a_test, "Scalar" , oneFunc,funcParameter);
    }

   
	MatrixPtr_Type Sum= Teuchos::rcp(new Matrix_Type(domain->getMapUnique(), 30  ) ); // with approx entries per row
	A->addMatrix(-1, Sum, 1);
	A_test->addMatrix(1, Sum, -1);

	//this = alpha*A + beta*this
	MultiVectorConstPtr_Type tmp = a_test;
    a->update( -1., tmp,1. );

	vec_dbl_Type resNorm(1);
	const Teuchos::ArrayView<SC> normRhs = Teuchos::arrayViewFromVector( resNorm);
	a->norm2(normRhs);
	if(comm->getRank() == 0)
		cout << " Norm of Difference between right hand sides: " << normRhs[0] << endl;

	
	int maxRank = std::get<1>(domain->getMesh()->rankRange_);

	double res=0.;
	Teuchos::ArrayView<const GO> indices;
	Teuchos::ArrayView<const SC> values;

	for (UN i=0; i < domain->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
		GO row = domain->getMapUnique()->getGlobalElement( i );
		Sum->getGlobalRowView(row, indices,values);
		
		for(int j=0; j< values.size() ; j++){
			res += values[j];			
		}		
	}
	res = fabs(res);
	reduceAll<int, double> (*comm, REDUCE_SUM, res, outArg (res));

	if(comm->getRank() == 0)
		cout << " Norm of Difference between StiffnessMatrices: " << res << endl;
	

    return(EXIT_SUCCESS);
}
