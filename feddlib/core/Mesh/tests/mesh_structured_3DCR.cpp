#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>

/*!
 MeshStructured 3D CR
 
 @brief  MeshStructured 3D Crouzeix-Raviart
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

    int dim = 3;
    myCLP.setOption("dim",&dim,"dim.");
    int M = 2;
    myCLP.setOption("M",&M,"H/h.");
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    
    // Mesh
    std::string FEType="P2-CR";
    int numProcsCoarseSolve = 0;
    int n;
    int size = comm->getSize();
    bool boolExportMesh = true;
    RCP<MeshStructured<SC,LO,GO,NO> > meshStr;
    if (dim == 3) {
        n = (int) (std::pow(size,1/3.) + 100.*ScalarTraits< SC >::eps()); // 1/H
        std::vector<double> x(3);
        x[0]=0.0; x[1]=0.0; x[2]=0.0;
        meshStr = rcp(new MeshStructured<SC,LO,GO,NO>(comm));
        meshStr->setGeometry3DBox(x, 1., 1., 1.);
        meshStr->buildMesh3D( FEType,n,M,numProcsCoarseSolve);
    }
    
    if (boolExportMesh) {
        RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        
        typedef Tpetra::MultiVector<SC,LO,GO,NO> TpetraMultiVector_Type;
        typedef RCP<TpetraMultiVector_Type> TpetraMultiVectorPtr_Type;
        
        TpetraMultiVectorPtr_Type tMV = RCP( new TpetraMultiVector_Type( meshStr->getMapUnique()->getTpetraMap(), 1));

        RCP<const MultiVector<SC,LO,GO,NO> > exportDummy = rcp(new MultiVector<SC,LO,GO,NO>(tMV));
        
        exPara->setup("u_p2_CR", meshStr, FEType);

        exPara->addVariable(exportDummy, "u", "Scalar", 1, meshStr->getMapUnique());
        
        exPara->save(0.0);
        
    }
    
    return(EXIT_SUCCESS);
}
