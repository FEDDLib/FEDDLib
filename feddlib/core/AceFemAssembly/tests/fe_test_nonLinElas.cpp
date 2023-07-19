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

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemNonLinElas.xml");
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
    
    MultiVectorPtr_Type d_rep = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldRepeated(), 1 ) ); // RHS vector
	d_rep->putScalar(solConst);
    MatrixPtr_Type A= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(), domain->getDimension() * domain->getApproxEntriesPerRow()  ) ); // Jacobi Matrix
    MultiVectorPtr_Type f = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldRepeated(), 1 ) ); // RHS vector
    cout << " FEDDLib Implementation .... " << endl;
    {
        //fe.assemblyElasticityJacobianAndStressAceFEM(dim, domain->getFEType(), A, f, d_rep, params, 1);
        fe.assemblyElasticityJacobianAceFEM(dim, domain->getFEType(), A,d_rep, "Neo-Hooke",youngModulus,poissonRatio,1.,true);
        fe.assemblyElasticityStressesAceFEM(dim, domain->getFEType(), f,d_rep, "Neo-Hooke",youngModulus,poissonRatio,1.,true);
                                                  

    }
    cout << " ... done " << endl;
    //A->print();
	// Class for assembling linear Elasticity via Acefem implementation
 	FE_Test<SC,LO,GO,NO> fe_test;
    fe_test.addFE(domain);
    
    MatrixPtr_Type A_test= Teuchos::rcp(new Matrix_Type( domain->getMapVecFieldUnique(),domain->getDimension() * domain->getApproxEntriesPerRow()   ) );
    MultiVectorPtr_Type f_test = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldRepeated(), 1 ) ); // RHS vector
    
    cout << " ACEGen Implementation .... " << endl;

    {
        fe_test.assemblyNonLinElas(dim, FEType, 2,dofs,d_rep, A_test,f_test,params,false,"Jacobian",true);
        fe_test.assemblyNonLinElas(dim, FEType, 2,dofs,d_rep, A_test,f_test,params,false,"Rhs",true);
        
    }
    cout << " ... done " << endl;

	//A_test->print();

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
				}
				if(!(fabs(values[j] - values2[j]) <= ( (fabs(values[j]) > fabs(values2[j]) ? fabs(values2[j]) : fabs(values[j])) * epsilon))){
					if(std::abs(values[j] - values2[j]) >= 1e-12) // Ignore differences smaller than 1e-12
						essentEqual=0;
						
				}
				if(!(fabs(values[j] - values2[j]) <= ( (fabs(values[j]) < fabs(values2[j]) ? fabs(values2[j]) : fabs(values[j])) * epsilon))){
					if(std::abs(values[j] - values2[j]) >= 1e-12)
						approxEqual=0;
				}
				
			}	
		}	
	}
		    
	//MultiVectorPtr_Type A_AssMV;
	//MultiVectorPtr_Type A_AssFEMV;
	
	//A->toMV( A_AssMV );
	//A_test->toMV( A_AssFEMV ); // A_test is Acegen

	//Teuchos::Array<SC> normMatrix(1); 
   // A_AssMV->normInf(normMatrix);

	//Teuchos::Array<SC> normMatrix2(1); 
    //A_AssFEMV->normInf(normMatrix2);

    
	reduceAll<int, double> (*comm, REDUCE_MAX, res, outArg (res));

	reduceAll<int, int> (*comm, REDUCE_MIN, approxEqual, outArg (approxEqual));

	reduceAll<int, int> (*comm, REDUCE_MIN, essentEqual, outArg (essentEqual));

	if(comm->getRank() == 0){
		cout << "Max Difference between StiffnessMatrices:\t \t \t \t" << res << endl;
		//cout << "Max Difference between StiffnessMatrices relative to A: \t \t " << res/normMatrix[0] << endl;
		//cout << "Max Difference between StiffnessMatrices relative to A from AceGen: \t " << res/normMatrix[0] << endl;

		if( essentEqual >0)
			cout << "Matrices are essentially equal" << endl;
		if( approxEqual >0)
			cout << "Matrices are approximately equal" <<  endl;
	}

    return(EXIT_SUCCESS);
}
