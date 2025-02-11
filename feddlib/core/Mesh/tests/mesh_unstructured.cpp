#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/BCBuilder.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 MeshUnstructured test

 @brief  MeshUnstructured test
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


void x1(double* x, double* res, double t, const double* parameters){
    res[0] = 1.;
}
void x2(double* x, double* res, double t, const double* parameters){
    res[0] = 2.;
}
void x3(double* x, double* res, double t, const double* parameters){
    res[0] = 3.;
}
void x4(double* x, double* res, double t, const double* parameters){
    res[0] = 4.;
}
void x5(double* x, double* res, double t, const double* parameters){
    res[0] = 5.;
}
void x6(double* x, double* res, double t, const double* parameters){
    res[0] = 6.;
}
void x7(double* x, double* res, double t, const double* parameters){
    res[0] = 7.;
}
void x8(double* x, double* res, double t, const double* parameters){
    res[0] = 8.;
}
void x9(double* x, double* res, double t, const double* parameters){
    res[0] = 9.;
}
void x10(double* x, double* res, double t, const double* parameters){
    res[0] = 10.;
}
void x11(double* x, double* res, double t, const double* parameters){
    res[0] = 11.;
}


using namespace std;
using namespace Teuchos;

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
using namespace FEDD;
int main(int argc, char *argv[]) {

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra"; //this does nothing atm
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    string filename = "dfg_fsi_fluid_h004.mesh";
    myCLP.setOption("file",&filename,"Mesh filename");
    int dim = 2;
    myCLP.setOption("dim",&dim,"Dimension");
    string delimiter = " ";
    myCLP.setOption("delimiter",&delimiter,"Delimiter in mesh-file");
    string FEType="P1";
    myCLP.setOption("FEType",&FEType,"FEType");
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    // Mesh
    
    int numProcsCoarseSolve = 0;
    bool boolExportMesh = true;
    bool boolExportSubdomains = false;
    int volumeID = 10;
    if (filename=="some_tetrahedron.mesh")
        volumeID = 12;

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

    
    Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

    bcFactory->addBC(x1, 1, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x2, 2, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x3, 3, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x4, 4, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x5, 5, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x6, 6, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x7, 7, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x8, 8, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x9, 9, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x10, 10, 0, domain, "Dirichlet", 1);
    bcFactory->addBC(x11, 11, 0, domain, "Dirichlet", 1);

    MultiVectorPtr_Type values = rcp(new MultiVector_Type( domain->getMapUnique() ) );
    BlockMultiVectorPtr_Type valuesBlock = rcp(new BlockMultiVector_Type( 1 ) );

    valuesBlock->addBlock( values, 0 );

    bcFactory->setRHS( valuesBlock );

    if (boolExportMesh) {
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        std::string filename = "unstructuredMesh";
        
        exPara->setup(filename, domain->getMesh(), FEType);

        MultiVectorConstPtr_Type valuesConst = valuesBlock->getBlock( 0 );
        exPara->addVariable( valuesConst, "values", "Scalar", 1, domain->getMapUnique() );

        exPara->save(0.0);
        exPara->closeExporter();

    }

    if (boolExportSubdomains) {

        MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domainP1->getElementMap() ) );
        MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
        vecDecomposition->putScalar(comm->getRank()+1.);

        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        exPara->setup( "subdomains", domainP1->getMesh(), "P0" );

        exPara->addVariable( vecDecompositionConst, "subdomain", "Scalar", 1, domainP1->getElementMap());

        exPara->save(0.0);
        exPara->closeExporter();
    }


    return(EXIT_SUCCESS);
}
