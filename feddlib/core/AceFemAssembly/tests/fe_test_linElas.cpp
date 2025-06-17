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
	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemLinElas.xml");
    Teuchos::CommandLineProcessor myCLP;
    int dim = params->sublist("Parameter").get("Dimension",3);
    myCLP.setOption("dim",&dim,"dim");
    int m = 2;
    myCLP.setOption("m",&m,"H/h");
    
    string filename = params->sublist("Mesh Partitioner").get("Mesh 1 Name","tetrahedron.mesh");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    

    // Mesh
    std::string FEType=params->sublist("Parameter").get("Discretization","P12");
    double poissonRatio=params->sublist("Parameter").get("Poisson Ratio",0.4e-0);

    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstanten \lambda
    double youngModulus=params->sublist("Parameter").get("E",1000.0);
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

    MatrixPtr_Type A= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) );

    {
		fe.assemblyLinElasXDim( dim, domain->getFEType(), A, lambda, mu );
    }
    A->writeMM("A");
	// Class for assembling linear Elasticity via Acefem implementation
 	FE<SC,LO,GO,NO> fe_test;
    fe_test.addFE(domain);
    
    MatrixPtr_Type A_test= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(),domain->getDimension() * domain->getApproxEntriesPerRow()   ) );
     BlockMatrixPtr_Type system = Teuchos::rcp( new BlockMatrix_Type( 1 ) ); //
	system->addBlock(A_test,0,0);
    MultiVectorPtr_Type d_rep = Teuchos::rcp(new MultiVector_Type(domain->getMapVecFieldRepeated(),1));
	d_rep->putScalar(0.);
    BlockMultiVectorPtr_Type resVec = Teuchos::rcp( new BlockMultiVector_Type( 2 ) ); // RHS vector
    MultiVectorPtr_Type d_res = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldUnique(), 1 ) ); // d_vec
    resVec->addBlock(d_res,0);
    {
        fe_test.assemblyLinearElasticity(dim, FEType, 2,dofs, d_rep, system, resVec, params,true, " ",true/*call fillComplete*/);
    }
	//A_test->writeMM("A_ACE");


    // Comparing matrices
	MatrixPtr_Type Sum= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) );
	A->addMatrix(1, Sum, 0);
	A_test->addMatrix(-1, Sum, 1);

	// Build A_test as an 'unfilled' Matrix in order to access its entries
	MatrixPtr_Type A_AssFE= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) );
	A_test->addMatrix(1, A_AssFE, 0);

	// Build A_test as an 'unfilled' Matrix in order to access its entries
	MatrixPtr_Type A_Ass= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) );
	A->addMatrix(1, A_Ass, 0);

	int maxRank = std::get<1>(domain->getMesh()->rankRange_);
	double res=0.;
	double relRes =0.;
	Teuchos::ArrayView<const GO> indices;
	Teuchos::ArrayView<const SC> values;

	Teuchos::ArrayView<const GO> indices2;
	Teuchos::ArrayView<const SC> values2;

	Teuchos::ArrayView<const GO> indices3;
	Teuchos::ArrayView<const SC> values3;
	double epsilon = 1e-13;
	int approxEqual = 1;
	int essentEqual = 1;
	for (UN i=0; i < domain->getMapUnique()->getMaxLocalIndex()+1 ; i++) {
		for(int d=0; d< dofs ; d++){
			GO row = dofs*domain->getMapUnique()->getGlobalElement( i )+d;

			A_Ass->getGlobalRowView(row, indices,values);
			
			A_AssFE->getGlobalRowView(row,indices2,values2);

			Sum->getGlobalRowView(row,indices3,values3);
			// We assume it is correct enough and the same indices are filled
			for(int j=0; j< values.size() ; j++){
				if(std::abs(values3[j])>res){
					res = std::abs(values3[j]);	
					cout << "Values1 " << values[j] << " Value2 " << values2[j] << endl;
				}
				if(!(std::fabs(values[j] - values2[j]) <= ( (std::fabs(values[j]) > std::fabs(values2[j]) ? std::fabs(values2[j]) : std::fabs(values[j])) * epsilon))){
					if(std::abs(values[j] - values2[j]) >= 1e-12) // Ignore differences smaller than 1e-12
						essentEqual=0;
						
				}
				if(!(std::fabs(values[j] - values2[j]) <= ( (std::fabs(values[j]) < std::fabs(values2[j]) ? std::fabs(values2[j]) : std::fabs(values[j])) * epsilon))){
					if(std::abs(values[j] - values2[j]) >= 1e-12)
						approxEqual=0;
				}
				
			}	
		}	
	}

	reduceAll<int, double> (*comm, REDUCE_MAX, res, outArg (res));

	reduceAll<int, int> (*comm, REDUCE_MIN, approxEqual, outArg (approxEqual));

	reduceAll<int, int> (*comm, REDUCE_MIN, essentEqual, outArg (essentEqual));

	if(comm->getRank() == 0){
		cout << "Max Difference between StiffnessMatrices: " << res << endl;
		if( essentEqual >0)
			cout << "Matrices are essentially equal" << endl;
		if( approxEqual >0)
			cout << "Matrices are approximately equal" <<  endl;
	}

    return(EXIT_SUCCESS);
}
