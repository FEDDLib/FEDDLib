#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 MeshStructured 3D 
 
 @brief  MeshStructured test for cube with 5 element subcube
 @version 1.0
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
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    int dim = 3;
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
    std::string FEType="P2";
    int numProcsCoarseSolve = 0;
    int n;
    int size = comm->getSize();
    bool boolExportMesh = true;
   // RCP<MeshStructured<SC,LO,GO,NO> > meshStr;

    DomainPtr_Type domain;
    if (dim == 3) {
        n = (int) (std::pow(size,1/3.) + 100.*ScalarTraits< SC >::eps()); // 1/H
        std::vector<double> x(3);
        x[0]=0.0; x[1]=0.0; x[2]=0.0;
        domain.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
		
        domain->buildMesh( 3,"Square5Element", dim, FEType, n, M, numProcsCoarseSolve);

    }
    
    if (boolExportMesh) {
        MultiVectorPtr_Type mvFlag = rcp(new MultiVector_Type(domain->getMapUnique() ) );
        vec_int_ptr_Type BCFlags = domain->getBCFlagUnique();
        Teuchos::ArrayRCP< SC > values = mvFlag->getDataNonConst(0);
        TEUCHOS_TEST_FOR_EXCEPTION( values.size()!= BCFlags->size(), std::runtime_error, "MultiVector and node list have different sizes.");
        for (int i=0; i<values.size(); i++)
            values[i] = BCFlags->at(i);
        
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        
        exPara->setup( "Element5SubcubeMeshWithFlagsOption3", domain->getMesh(), FEType );
        MultiVectorConstPtr_Type mvFlagConst = mvFlag;
        exPara->addVariable( mvFlagConst, "Flags", "Scalar", 1,domain->getMapUnique());
        
        exPara->save(0.0);
        exPara->closeExporter();
    }
    
    return(EXIT_SUCCESS);
}
