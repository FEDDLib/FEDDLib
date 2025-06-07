#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

using namespace std;
using namespace Teuchos;
using namespace FEDD;

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

int main(int argc, char *argv[]) {

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

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

    RCP<Domain<SC,LO,GO,NO> > domain;
    if (dim == 2) {
        n = (int) (std::pow(size,1/2.) + 100.*ScalarTraits< SC >::eps()); // 1/H
        std::vector<double> x(2);
        x[0]=0.0;    x[1]=0.0;
        domain.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., comm ) );
    }
    else if (dim == 3){
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Implement 3D test.");
    }
    domain->buildMesh( 1,"Square", dim, FEType, n, m, numProcsCoarseSolve);
    
    FE<SC,LO,GO,NO> fe;
    fe.addFE(domain);
    
    Teuchos::RCP<Teuchos::Time> timer1(Teuchos::TimeMonitor::getNewCounter("ASSEMBLY: scalar procedure"));
    Teuchos::RCP<Teuchos::Time> timer2(Teuchos::TimeMonitor::getNewCounter("ASSEMBLY: vector blas"));
    
    MatrixPtr_Type A1(new Matrix_Type( domain->getMapVecFieldUnique(), 10 ) );
    {
        Teuchos::TimeMonitor startTM(*timer1);
        fe.assemblyLaplaceVecField(dim, FEType, 2, A1, true/*call fillComplete*/);
    }
    MatrixPtr_Type A2(new Matrix_Type( domain->getMapVecFieldUnique(), 20 ) );
    {
        Teuchos::TimeMonitor startTM(*timer2);
        fe.assemblyLaplaceVecFieldV2(dim, FEType, 2, A2, true/*call fillComplete*/);
    }
    
    
    bool exportMesh = false;
    if (exportMesh) {
        RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        RCP<const MultiVector<SC,LO,GO,NO> > exportDummy = rcp( new MultiVector<SC,LO,GO,NO>( domain->getMapUnique() ) );
        
        exPara->setup("u", domain->getMesh(), FEType);

        exPara->addVariable(exportDummy, "u", "Scalar", 1, domain->getMapUnique());
        
        exPara->save(0.0);
    }
    Teuchos::TimeMonitor::report(std::cout);
    return(EXIT_SUCCESS);
}
