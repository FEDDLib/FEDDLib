#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/problems/abstract/Problem.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/General/HDF5Export.hpp"


void scalarFunc(double* x, double* res, double t, const double* parameters){
   
    res[0] = 1.;
    return;
}
void zeroDirichlet(double* x, double* res, double t, const double* parameters){
   
    res[0] = 0.;
    res[1] = 0.;

    return;
}


void vector2DFunc(double* x, double* res, double t, const double* parameters){
    
    res[0] = -1.;
    res[1] = -2.;
    
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
    
    typedef BlockMultiVector<> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
    
    typedef Map<> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef RCP<const Map_Type> MapConstPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    typedef HDF5Export<SC,LO,GO,NO> Export_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    int intFEType = 1;
    myCLP.setOption("FEType",&intFEType,"P1 or P2. Choose 1 or 2.");

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
    std::string FEType;
    if (intFEType==1)
        FEType = "P1";
    else
        FEType = "P2";
    

    {
        // Reading Mesh
        std::string filename="square.mesh";
        int dim = 2;
        ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
        pListPartitioner->set( "Mesh 1 Name", filename );
        DomainPtr_Type domain;
        
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
        
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );
    
        bcFactory->addBC(zeroDirichlet, 3, 0, domain, "Dirichlet", 1); // Setting some sort of boundary conditions
        
        std::vector<MapConstPtr_Type> maps(1);
        maps[0] = domain->getMapUnique();
        BlockMultiVectorPtr_Type  vector = Teuchos::rcp( new BlockMultiVector_Type( maps ) ) ;
        vector->getBlockNonConst(0)->putScalar(1.);
        bcFactory->setRHS(vector);
        
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        std::string filenameExport = "exportValues";

        exPara->setup(filenameExport, domain->getMesh(), FEType);

        MultiVectorConstPtr_Type aUniqueConst = vector->getBlock(0);
        exPara->addVariable( aUniqueConst, "Export", "Scalar", 1, domain->getMapUnique() );
        exPara->save(0.0);
        exPara->closeExporter();

        // ----------------------
        // Exporting as HDF5 Type
        HDF5Export<SC,LO,GO,NO> exporter(domain->getMapUnique(), "exportVector"); //  Map and file name
        exporter.writeVariablesHDF5("Test",aUniqueConst); // VariableName and Variable

        HDF5Import<SC,LO,GO,NO> importer(domain->getMapUnique(),"exportVector"); // Import Map and  file name to read
        MultiVectorPtr_Type aImported = importer.readVariablesHDF5("Test"); // VariableName 
        
        TEUCHOS_TEST_FOR_EXCEPTION(std::fabs(aUniqueConst->getLocalLength() - aImported->getLocalLength())>1e-13  , std::runtime_error, "Local length of written and read vector is differnt"); // Checking of vectors have the same length

        // Testing difference
        double error= 0.;
        for(int i=0; i<aUniqueConst->getLocalLength(); i++)
            error += std::abs((aUniqueConst->getData(0))[i] - (aImported->getData(0))[i]);
        
        TEUCHOS_TEST_FOR_EXCEPTION( error>1e-13, std::runtime_error, "Error between written and read vector is to great: " << error); // Checking if error between solutions is to great

    }
    
    
    return(EXIT_SUCCESS);
}
