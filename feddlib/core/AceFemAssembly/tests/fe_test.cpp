#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/AceFemAssembly/TestFE/FE_Test.hpp"
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
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
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
    
//    Xpetra::UnderlyingLib ulib;
//    if (!ulib_str.compare("UseTpetra"))
//        ulib = Xpetra::UseTpetra;
//    else if (!ulib_str.compare("UseEpetra"))
//        ulib = Xpetra::UseEpetra;

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
   

    FE_Test<SC,LO,GO,NO> fe_test;
    fe_test.addFE(domain);
    
    MatrixPtr_Type A1(new Matrix_Type( domain->getMapUnique(), 10 ) );
    {
        fe_test.assemblyLaplace(dim, FEType, 2, A1, true/*call fillComplete*/);
    }

    FE<SC,LO,GO,NO> fe;
    fe.addFE(domain);
    MatrixPtr_Type A2(new Matrix_Type( domain->getMapUnique(), 10 ) );
    {
        fe.assemblyLaplace(dim, FEType, 2, A2, true/*call fillComplete*/);
    }
    
    if(comm->getRank() == 0){
		cout << " Compare Assemblys " << endl;
	
	}
    bool exportMesh = false;
    if (exportMesh) {
        RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        RCP<const MultiVector<SC,LO,GO,NO> > exportDummy = rcp( new MultiVector<SC,LO,GO,NO>( domain->getMapUnique() ) );
        
        exPara->setup("u", domain->getMesh(), FEType);

        exPara->addVariable(exportDummy, "u", "Scalar", 1, domain->getMapUnique());
        
        exPara->save(0.0);
    }

    return(EXIT_SUCCESS);
}
