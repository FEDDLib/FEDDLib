#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

void constFunc(double* x, double* res, double* parameters){
    
    res[0] = 1.;

}

void linFunc(double* x, double* res, double* parameters){
    
    res[0] = x[0] + 2. * x[1];

}

void quadFunc(double* x, double* res, double* parameters){
    
    res[0] = x[0]+ 2. * x[0]*x[1];
}

void cubicFunc(double* x, double* res, double* parameters){
    
    res[0] = x[0]+ 2. * x[0]*x[1]*x[1];
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

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

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
    MultiVectorPtr_Type FERhs;
    MultiVectorPtr_Type feRhsUnique;
    for (int i=1; i<6; i++) {
        FERhs = Teuchos::rcp(new MultiVector_Type( domain->getMapRepeated() ) );
        vec_dbl_Type funcParameter(2,0.); // second value does not matter for our functions used
        fe.assemblyRHSDegTest(domain->getDimension(),
                                             domain->getFEType(),
                                             FERhs,
                                             "Scalar",
                                             constFunc,
                                             funcParameter,
                                             i);
        feRhsUnique = Teuchos::rcp(new MultiVector_Type( domain->getMapUnique() ) );
        feRhsUnique->exportFromVector( FERhs, false, "Add" );
        feRhsUnique->print();
    }
    for (int i=1; i<6; i++) {
        FERhs = Teuchos::rcp(new MultiVector_Type( domain->getMapRepeated() ) );
        vec_dbl_Type funcParameter(2,1.); // second value does not matter for our functions used
        fe.assemblyRHSDegTest(domain->getDimension(),
                                             domain->getFEType(),
                                             FERhs,
                                             "Scalar",
                                             linFunc,
                                             funcParameter,
                                             i);
        feRhsUnique = Teuchos::rcp(new MultiVector_Type( domain->getMapUnique() ) );
        feRhsUnique->exportFromVector( FERhs, false, "Add" );
        feRhsUnique->print();
    }
    for (int i=1; i<7; i++) {
        FERhs = Teuchos::rcp(new MultiVector_Type( domain->getMapRepeated() ) );
        vec_dbl_Type funcParameter(2,2.); // second value does not matter for our functions used
        fe.assemblyRHSDegTest(domain->getDimension(),
                                             domain->getFEType(),
                                             FERhs,
                                             "Scalar",
                                             quadFunc,
                                             funcParameter,
                                             i);
        feRhsUnique = Teuchos::rcp(new MultiVector_Type( domain->getMapUnique() ) );
        feRhsUnique->exportFromVector( FERhs, false, "Add" );
        feRhsUnique->print();
    }
    for (int i=1; i<7; i++) {
        FERhs = Teuchos::rcp(new MultiVector_Type( domain->getMapRepeated() ) );
        vec_dbl_Type funcParameter(2,3.); // second value does not matter for our functions used
        fe.assemblyRHSDegTest(domain->getDimension(),
                                             domain->getFEType(),
                                             FERhs,
                                             "Scalar",
                                             cubicFunc,
                                             funcParameter,
                                             i);
        Teuchos::rcp(new MultiVector_Type( domain->getMapUnique() ) );
        feRhsUnique->exportFromVector( FERhs, false, "Add" );
        feRhsUnique->print();
    }
    
    
    return(EXIT_SUCCESS);
}
