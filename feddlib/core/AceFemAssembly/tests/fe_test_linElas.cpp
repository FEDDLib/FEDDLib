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
	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemLinElas.xml");
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
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
    std::string FEType=params->sublist("Parameter").get("Discretization","P1");
    double mu=params->sublist("Parameter").get("Mu",0.3571);
    double poissonRatio=params->sublist("Parameter").get("Poisson Ratio",0.4e-0);

    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstanten \lambda
    double youngModulus = mu*2.*(1 + poissonRatio);
    double lambda = (poissonRatio*youngModulus)/((1 + poissonRatio)*(1 - 2*poissonRatio));

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

    MatrixPtr_Type A= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) );

    {
		fe.assemblyLinElasXDim( dim, domain->getFEType(), A, lambda, mu );
    }
    
	// Class for assembling linear Elasticity via Acefem implementation
 	FE_Test<SC,LO,GO,NO> fe_test;
    fe_test.addFE(domain);
    
    MatrixPtr_Type A_test= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(),domain->getDimension() * domain->getApproxEntriesPerRow()   ) );

    {
        fe_test.assemblyLinElas(dim, FEType, 2,dofs, A_test, true/*call fillComplete*/);
    }

   	// Comparing matrices
	MatrixPtr_Type Sum= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) );
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
