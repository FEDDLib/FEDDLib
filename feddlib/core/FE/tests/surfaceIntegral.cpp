#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/abstract/Problem.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>


void fScalar(double* x, double* res, double* parameters){

    if (x[0]<0.000001) {
        res[0] = -1.;
    }
    if (x[0]>0.999999) {
        res[0] = 2.;
    }
    if (x[1]<0.000001) {
        res[0] = -3.;
    }
    if (x[1]>0.999999) {
        res[0] = 4.;
    }
    return;
}

void fVector(double* x, double* res, double* parameters){
    
    if (x[0]<0.000001) {
        res[0] = -1.;
        res[1] = -1.*2.;
    }
    if (x[0]>0.999999) {
        res[0] = 2.;
        res[1] = 2.*2.;
    }
    if (x[1]<0.000001) {
        res[0] = -3.;
        res[1] = -3.*2.;
    }
    if (x[1]>0.999999) {
        res[0] = 4.;
        res[1] = 4.*2.;
    }
    return;

}

void fScalar3D(double* x, double* res, double* parameters){
    // parameters[0]: desired flag, parameters[1]: current flag of element
    if (parameters[1]==parameters[0]) {
        res[0] = -1.;
    }
    else{
        res[0] = 0.;
    }

    return;
}

void fVector3D(double* x, double* res, double* parameters){
    
    // parameters[0]: desired flag, parameters[1]: current flag of element
    if (parameters[1]==parameters[0]) {
        res[0] = 0.;
        res[1] = -1.;
        res[2] = 0.;
    }
    else {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = 0.;
    }
    return;

}

using namespace std;
using namespace Teuchos;
using namespace FEDD;

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

int main(int argc, char *argv[]) {

    
    typedef Domain<> Domain_Type;
    typedef RCP<Domain_Type> DomainPtr_Type;
    typedef FE<SC,LO,GO,NO> FEFac_Type;
    typedef Teuchos::RCP<FEFac_Type> FEFacPtr_Type;
    typedef MultiVector<> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string type = "unstructured";
    myCLP.setOption("type",&type,"structured or unstructured.");
    string FEType = "P1";
    myCLP.setOption("FEType",&FEType,"P1 or P2.");
    int dim = 2;
    myCLP.setOption("dim",&dim,"2 or 3.");
    string filename="square.mesh";
    myCLP.setOption("mesh",&filename,"Name of mesh-file.");
    
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    int n;
    int size = comm->getSize();
    int m = 2;
    int numProcsCoarseSolve = 0;
    bool verbose(comm->getRank()==0);


    {
        // Mesh
        bool exportParaView = true;
        ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
        pListPartitioner->set( "Mesh 1 Name", filename );
        DomainPtr_Type domain;
        if ( type == "structured" ) {
            if (dim == 2) {
                n = (int) (std::pow(size,1/2.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
                std::vector<double> x(2);
                x[0]=0.0;    x[1]=0.0;
                domain.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
            }
            
            domain->buildMesh( 3, "SquareTPM", dim, FEType, n, m, numProcsCoarseSolve );
        }
        else if ( type == "unstructured" ) {
            DomainPtr_Type domainP1 = Teuchos::rcp( new Domain_Type( comm, dim ) );
            
            MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
            domainP1Array[0] = domainP1;

            MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
            
            partitionerP1.readAndPartition();
            
            if (FEType == "P2") {
                domain.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
                domain->buildP2ofP1Domain( domainP1 );
            }
            else
                domain = domainP1;
        }

        FEFacPtr_Type feFactory = Teuchos::rcp( new FEFac_Type() );
        
        feFactory->addFE(domain);
        
        MultiVectorPtr_Type  a = Teuchos::rcp( new MultiVector_Type( domain->getMapRepeated(), 1 ) );
        MultiVectorPtr_Type  aUnique = Teuchos::rcp( new MultiVector_Type( domain->getMapUnique(), 1 ) );
        
        std::vector<SC> funcParameter(3 , 0.);//0: surface flag (no specific flag is used in function above), 1:placeholder for surface flag of element, 2:order
        if (dim==3) {
            funcParameter[0] = 2.; // specific flag 2 is used in the 3D case
            feFactory->assemblySurfaceIntegral( dim, FEType, a, "Scalar", fScalar3D, funcParameter);
        }
        else{
            feFactory->assemblySurfaceIntegral( dim, FEType, a, "Scalar", fScalar, funcParameter);
        }
        aUnique->exportFromVector( a, false, "Add" );
        
        MultiVectorPtr_Type  b = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldRepeated(), 1 ) );
        MultiVectorPtr_Type  bUnique = Teuchos::rcp( new MultiVector_Type( domain->getMapVecFieldUnique(), 1 ) );
        
        if (dim==3) {
            funcParameter[0] = 2.; // specific flag 2 is used in the 3D case
            feFactory->assemblySurfaceIntegral( dim, FEType, b, "Vector", fVector3D, funcParameter);
        }
        else{
            feFactory->assemblySurfaceIntegral( dim, FEType, b, "Vector", fVector, funcParameter);
        }
        bUnique->exportFromVector( b, false, "Add" );
        
        if (exportParaView) {
            
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
            std::string filename = "surfaceIntegrals";
            
            exPara->setup(filename, domain->getMesh(), FEType);
            
            MultiVectorConstPtr_Type aUniqueConst = aUnique;
            MultiVectorConstPtr_Type bUniqueConst = bUnique;
            exPara->addVariable( aUniqueConst, "scalar surface", "Scalar", 1, domain->getMapUnique() );
            UN dofsPerNode = dim;
            exPara->addVariable( bUniqueConst, "vector surface", "Vector", dofsPerNode, domain->getMapUnique() );
            exPara->save(0.0);
            exPara->closeExporter();
            
        }
        
    }
    
    
    return(EXIT_SUCCESS);
}
