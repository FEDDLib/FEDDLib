#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>

using namespace std;
using namespace Teuchos;
using namespace FEDD;

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

int main(int argc, char *argv[]) {

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();    
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    int dim = 2;
    myCLP.setOption("dim",&dim,"dim");
    int m = 2;
    myCLP.setOption("m",&m,"H/h");
    
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
    int n;
    int size = comm->getSize();
    
    bool exportMesh = true;
    bool exportSubdomains = true;
    RCP<Domain<SC,LO,GO,NO> > domain;
    if (dim == 2) {
        n = (int) (std::pow(size,1/2.) + 100.*ScalarTraits< SC >::eps()); // 1/H
        std::vector<double> x(2);
        x[0]=0.0;    x[1]=0.0;
        domain = rcp( new Domain<SC,LO,GO,NO>(x, 1., 1., comm) ) ;
        domain->buildMesh(0, "Square", dim, FEType, n, m, 0);
    }
    else if (dim == 3){
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Implement 3D test.");
    }
    
    
    
    if (exportMesh) {
        RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        RCP<const MultiVector<SC,LO,GO,NO> > exportDummy = rcp( new MultiVector<SC,LO,GO,NO>( domain->getMapUnique() ) );
        
        exPara->setup("u", domain->getMesh(), FEType);
        
        exPara->addVariable(exportDummy, "u", "Scalar", 1, domain->getMapUnique());
        
        exPara->save(0.0);
    }
    
    if ( exportSubdomains ){
    
        typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
        typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
        typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
        typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
        typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

        {
            MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domain->getElementMap() ) );
            MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
            vecDecomposition->putScalar(comm->getRank()+1.);
            
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
            
            exPara->setup( "structured_subdomains", domain->getMesh(), "P0" );
            
            exPara->addVariable( vecDecompositionConst, "subdomains", "Scalar", 1, domain->getElementMap());
            exPara->save(0.0);
            exPara->closeExporter();
        }
    }
    return(EXIT_SUCCESS);
}
