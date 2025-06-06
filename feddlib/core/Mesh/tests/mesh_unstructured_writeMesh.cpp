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



using namespace std;
using namespace Teuchos;

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
    string filename = "square.mesh";
    myCLP.setOption("file",&filename,"Mesh filename");
    string exportfilename = "export.mesh";
    myCLP.setOption("fileExport",&exportfilename,"Export Mesh filename");
    int dim = 2;
    myCLP.setOption("dim",&dim,"Dimension");
    string delimiter = " ";
    myCLP.setOption("delimiter",&delimiter,"Delimiter in mesh-file");
    string FEType="P1";
    myCLP.setOption("FEType",&FEType,"FEType");
   /* bool exportEdges= true;
    myCLP.setOption("exportEdges",&exportEdges,"Exporting Edges");
    bool exportSurfaces= true;
    myCLP.setOption("exportSurfaces",&exportSurfaces,"Exporting Surfaces");
    */
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
    int volumeID = 16;
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
    
    partitionerP1.readAndPartition(volumeID);

    if (FEType == "P2") {
        domainP2->buildP2ofP1Domain( domainP1 );
        domain = domainP2;
    }
    else
        domain = domainP1;

    domain->exportMesh(true,true,exportfilename);
    comm->barrier();
    if(FEType == "P1"){
        ParameterListPtr_Type pListPartitionerTest = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
        pListPartitionerTest->set( "Mesh 1 Name", exportfilename );
        
        DomainPtr_Type domainP1Test;
        DomainPtr_Type domainTest;

        domainP1Test.reset( new Domain_Type( comm, dim ) );
        MeshPartitioner_Type::DomainPtrArray_Type domainP1ArrayTest(1);
        domainP1ArrayTest[0] = domainP1Test;

        MeshPartitioner<SC,LO,GO,NO> partitionerP1Test ( domainP1ArrayTest, pListPartitionerTest, "P1", dim );
        comm->barrier();

        partitionerP1Test.readAndPartition(volumeID);

        domainTest = domainP1Test;

        if(comm->getRank() == 0){
            cout << " ------------------------------------------------------------------- " << endl;
            cout << " Comparing nodes between originally read mesh and exported mesh ..." << endl;
        }
        vec2D_dbl_ptr_Type points1 = domain->getPointsUnique();
        vec_int_ptr_Type bc1 = domain->getBCFlagUnique();

        vec2D_dbl_ptr_Type points2 = domainTest->getPointsUnique();
        vec_int_ptr_Type bc2 = domainTest->getBCFlagUnique();

        if(points1->size() != points2->size())
            cout << " ERROR !! the meshes seemed to have been partitioned differently " << endl;

        double error=0.;
        for(int i=0; i< points1->size() ; i++){
            if(bc1->at(i) != bc2->at(i))
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error ,"Different Flag for two points..");

            for (int j=0; j< dim; j++){
                double errorTmp = std::abs(points1->at(i).at(j) - points2->at(i).at(j));
                if(errorTmp > error){
                    error = errorTmp/std::abs(points1->at(i).at(j));
                }
            }
        }
        reduceAll<int, double> (*comm, REDUCE_MAX, error, outArg (error));
           
        if(comm->getRank() == 0){
            cout << " Maximal relative difference between nodes: " << error << endl;
            cout << " ------------------------------------------------------------------- " << endl;
            TEUCHOS_TEST_FOR_EXCEPTION( error > 1.e-5, std::runtime_error ,"Error between nodes to great. Either read or write went wrong.");

        }
        MultiVectorPtr_Type vecDecomposition =  rcp(new MultiVector_Type( domain->getMapUnique() ) );
        MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
        vecDecomposition->putScalar(comm->getRank()+1.);

        MultiVectorPtr_Type vecDecompositionTest = rcp(new MultiVector_Type( domainTest->getMapUnique() ) );
        MultiVectorConstPtr_Type vecDecompositionConstTest = vecDecompositionTest;
        vecDecompositionTest->putScalar(comm->getRank()+1.);


        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        exPara->setup( "Mesh_Export_Test_Compare", domainTest->getMesh(), "P1" );

        exPara->addVariable( vecDecompositionConst, "orignal", "Scalar", 1, domain->getMapUnique());
        exPara->addVariable( vecDecompositionConstTest, "exportMesh", "Scalar", 1, domainTest->getMapUnique());

        exPara->save(0.0);
        exPara->closeExporter();
    }
    
   

    

    


    return(EXIT_SUCCESS);
}
