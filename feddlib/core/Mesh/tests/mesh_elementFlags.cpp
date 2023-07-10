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
 Mesh Element Flags test

 @brief  Mesh Element Flags test
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


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
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
    
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra"; //this does nothing atm
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    string filename = "testFoam2.mesh";
    myCLP.setOption("file",&filename,"Mesh filename");
    int dim = 3;
    myCLP.setOption("dim",&dim,"Dimension");
    string delimiter = " ";
    myCLP.setOption("delimiter",&delimiter,"Delimiter in mesh-file");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    // Mesh
    std::string FEType="P1";
    int numProcsCoarseSolve = 0;
    bool boolExportMesh = true;
    int volumeID = 10;

    DomainPtr_Type domainP1;
    DomainPtr_Type domain;

    
    ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
    pListPartitioner->set( "Mesh 1 Name", filename );
    pListPartitioner->set( "Build Edge List", false );
    pListPartitioner->set( "Build Surface List", false );
    
    //Paralleles P1 mesh (Partitionierung mit metis)
    domainP1.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
    domainP1Array[0] = domainP1;

    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
    
    partitionerP1.readAndPartition();

    if (boolExportMesh) {

        MultiVectorPtr_Type mvFlag = rcp(new MultiVector_Type( domainP1->getElementMap() ) );
        ElementsPtr_Type elements = domainP1->getElementsC();
        Teuchos::ArrayRCP< SC > values = mvFlag->getDataNonConst(0);
        TEUCHOS_TEST_FOR_EXCEPTION( values.size()!= elements->numberElements(), std::runtime_error, "MultiVector and element list have different sizes.");
        for (int i=0; i<values.size(); i++)
            values[i] = elements->getElement(i).getFlag();
        
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        
        exPara->setup( "elementFlags", domainP1->getMesh(), "P0" );
        MultiVectorConstPtr_Type mvFlagConst = mvFlag;
        exPara->addVariable( mvFlagConst, "flag", "Scalar", 1, domainP1->getElementMap());
        
        exPara->save(0.0);
        exPara->closeExporter();

    }

    return(EXIT_SUCCESS);
}
