#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>

/*!
 MeshStructured test
 
 @brief  MeshStructured test
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

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
   
    string FEType = "P1";
    myCLP.setOption("FEType",&FEType,"FEType");
    int dim = 2;
    myCLP.setOption("dim",&dim,"dim.");
    int m=1;
    myCLP.setOption("m",&m,"H/h.");
    std::string meshType = "structured";
    myCLP.setOption("meshType",&meshType,"Mesh type");
    int length = 2;
    myCLP.setOption("length",&length,"length of bfs.");
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    
    // Mesh
    int minNumberSubdomains;
    std::string meshGeometry;
    if (!meshType.compare("structured")) {
        minNumberSubdomains = 1;
        meshGeometry="cube";
    }
    else if(!meshType.compare("structured_bfs")){
        minNumberSubdomains = (int) 2*length+1;
        meshGeometry="bfs";
    }
    

    int numProcsCoarseSolve = 0;
    int n;
    int size = comm->getSize();
    bool boolExportMesh = true;
    bool boolExportMeshSubdomains = true;
    RCP<MeshStructured<SC,LO,GO,NO> > meshStr;
    if (dim == 2) {
        n = (int) (std::pow(size/minNumberSubdomains,1/2.) + 100.*ScalarTraits< SC >::eps()); // 1/H
        std::vector<double> x(2);
        x[0]=0.0;    x[1]=0.0;
        meshStr = rcp(new MeshStructured<SC,LO,GO,NO>(comm));
        meshStr->setGeometry2DRectangle(x, 1., 1.);
        meshStr->buildMesh2D( FEType,n,m,numProcsCoarseSolve);
    }
    else if (dim == 3){
        n = (int) (std::pow(size/minNumberSubdomains,1/3.) + 100.*ScalarTraits< SC >::eps()); // 1/H
        std::vector<double> x(3);
        if (meshGeometry=="cube") {
            x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
            meshStr = rcp(new MeshStructured<SC,LO,GO,NO>(comm));
            meshStr->setGeometry3DBox(x, 1., 1., 1.);
            meshStr->buildMesh3D( FEType,n,m,numProcsCoarseSolve);
        }
        else if(meshGeometry=="bfs") {
            n = (int) (std::pow(size,1/3.) + 100.*ScalarTraits< SC >::eps()); // 1/H
            x[0]=-1.0;    x[1]=0.0;    x[2]=-1.0;
            meshStr = rcp(new MeshStructured<SC,LO,GO,NO>(comm));
            meshStr->setGeometry3DBox(x, length+1., 1., 2.);
            meshStr->buildMesh3DBFS( FEType,n,m,numProcsCoarseSolve);
        }

    }
    meshStr->buildElementMap();
    
    if (boolExportMesh) {
        RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        
        typedef Tpetra::MultiVector<SC,LO,GO,NO> TpetraMultiVector_Type;
        typedef RCP<TpetraMultiVector_Type> TpetraMultiVectorPtr_Type;
        
        TpetraMultiVectorPtr_Type tMV = RCP( new TpetraMultiVector_Type( meshStr->getMapUnique()->getTpetraMap(), 1));

        RCP<const MultiVector<SC,LO,GO,NO> > exportDummy = rcp(new MultiVector<SC,LO,GO,NO>(tMV));
        
        exPara->setup("u", meshStr, FEType);

        exPara->addVariable(exportDummy, "u", "Scalar", 1, meshStr->getMapUnique());
        
        exPara->save(0.0);
        
    }
        
    
    return(EXIT_SUCCESS);
}
