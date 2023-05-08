#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

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

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    int dim = 3;
    myCLP.setOption("dim",&dim,"dim.");
    int M = 2;
    myCLP.setOption("M",&M,"H/h.");
    int length = 1;
    myCLP.setOption("length",&length,"length");
    
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
    int minNumberSubdomains = (int) 2*length+1;
    bool boolExportMesh = false;
    RCP<MeshStructured<SC,LO,GO,NO> > meshStr;
    if (dim == 3) {
        n = (int) (std::pow(size/minNumberSubdomains,1/3.) + 100.*ScalarTraits< SC >::eps()); // 1/H
        std::vector<double> x(3);
        x[0]=-1.0;    x[1]=0.0;    x[2]=-1.0;
        meshStr = rcp(new MeshStructured<SC,LO,GO,NO>(comm));
        meshStr->setGeometry3DBox(x, length+1., 1., 2.);
        meshStr->buildMesh3DBFS( FEType,n,M,numProcsCoarseSolve);
    }
    
    if (boolExportMesh) {
        RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> XpetraMV_Type;
        typedef RCP<XpetraMV_Type> XpetraMVPtr_Type;
        
        XpetraMVPtr_Type xMV = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(meshStr->getMapUnique()->getXpetraMap(), 1);

        RCP<const MultiVector<SC,LO,GO,NO> > exportDummy = rcp(new MultiVector<SC,LO,GO,NO>(xMV));
        
        exPara->setup(meshStr->getDimension(), meshStr->getNumElementsGlobal(), meshStr->getElements(), meshStr->getPointsUnique(), meshStr->getMapUnique(), meshStr->getMapRepeated(), FEType, "u",1,comm);
        
        exPara->addVariable(exportDummy, "u", "Scalar", 1, meshStr->getMapUnique());
        
        exPara->save(0.0);
        
    }
    
    return(EXIT_SUCCESS);
}
