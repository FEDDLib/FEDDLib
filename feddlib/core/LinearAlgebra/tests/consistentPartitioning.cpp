#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

#include <unistd.h>

/*!
 Conistency test for different partitioning types. The linear algebra behaviour (Trilinos Tpetra) is generally not consistent for different partitions.
 Using a distributed object, which is only distributed to a real subset available ranks, might lead to inconsistent behaviour.

 @brief  Consistent partitioning
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

    RCP<const Comm<int> > comm  = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();
    int rank = comm->getRank();
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    GO numGlobalElements1 = 10;
    myCLP.setOption("nge1",&numGlobalElements1,"numGlobalElements1.");
    GO numGlobalElements2 = 20;
    myCLP.setOption("nge2",&numGlobalElements2,"numGlobalElements2.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef RCP<Matrix_Type> MatrixPtr_Type;

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef RCP<const Map_Type> MapConstPtr_Type;
    
    typedef BlockMap<LO,GO,NO> BlockMap_Type;
    typedef RCP<BlockMap_Type> BlockMapPtr_Type;

    TEUCHOS_TEST_FOR_EXCEPTION( ulib_str != "Tpetra", std::runtime_error,"Only Tpetra available.");

    MapConstPtr_Type map1;
    MapConstPtr_Type map2;
    MapConstPtr_Type map3;
    MapConstPtr_Type map4dotRes;
    BlockMultiVectorPtr_Type a = rcp(new BlockMultiVector_Type( 2 ));
    BlockMultiVectorPtr_Type b = rcp(new BlockMultiVector_Type( 2 ));
    {
        Teuchos::Array<GO> indices(1);
        indices[0] = comm->getRank();
        map1 = rcp( new Map_Type( ulib_str, (GO) -1, indices(), 0, comm ) );
        MultiVectorPtr_Type part1 = rcp( new MultiVector_Type( map1 ) );
        part1->putScalar(2.);
        MultiVectorPtr_Type part2 = rcp( new MultiVector_Type( map1 ) );
        part2->putScalar(3.);
        a->addBlock( part1, 0);
        a->addBlock( part2, 1);
    }
    {
        Teuchos::Array<GO> indices;
        if (rank==0){
            indices.push_back(0); indices.push_back(1);
        }
        
        map2 = rcp( new Map_Type( ulib_str, (GO) -1, indices(), 0, comm ) );
        MultiVectorPtr_Type part1 = rcp( new MultiVector_Type( map2 ) );
        part1->putScalar(2.);
        b->addBlock( part1, 0);
    }
    {
        Teuchos::Array<GO> indices;
        if (rank==1){
            indices.push_back(0); indices.push_back(1);
        }
        map3 = rcp( new Map_Type( ulib_str, (GO) -1, indices(), 0, comm ) );
        MultiVectorPtr_Type part2 = rcp( new MultiVector_Type( map3 ) );
        part2->putScalar(3.);
        b->addBlock( part2, 1);
    }

    {
        Teuchos::Array<GO> indices;

        indices.push_back(0);
        Teuchos::RCP<const Teuchos::Comm<LO> > commSerial = rcp(new Teuchos::MpiComm<LO>(MPI_COMM_SELF));
        map4dotRes = rcp( new Map_Type( ulib_str, (GO) -1, indices(), 0, commSerial ) );
    }

    Array< ScalarTraits<SC>::magnitudeType> dots(1);
    
    map1->print();
    map2->print();
    map3->print();
    a->print();
    b->print();
    
    a->dot(a, dots());
    std::cout << "rank: " << rank << " a dot a:" << dots[0] << std::endl;

    b->dot(b, dots());
    std::cout << "rank: " << rank << " b dot b:" << dots[0] << std::endl;
    
    MultiVectorPtr_Type c = rcp( new MultiVector_Type( map4dotRes ) );
    map4dotRes->print();
    
    BlockMultiVectorConstPtr_Type aConst = a.getConst();
    BlockMultiVectorConstPtr_Type bConst = b.getConst();
    c->multiply( TRANS, NO_TRANS, 1., aConst, aConst, 0. );
    std::cout << "version 1"<< std::endl;
    c->print();
    usleep(1e6);
    
    c->multiply( TRANS, NO_TRANS, 1., bConst, bConst, 0. );
    std::cout << "version 2"<< std::endl;
    c->print();
    usleep(1e6);
    
//    c->putScalar(2.);
//    c->multiply( TRANS, NO_TRANS, 1., aConst, aConst, 1. );
//    c->print();
//    c->putScalar(2.);
//    c->multiply( TRANS, NO_TRANS, 1., bConst, bConst, 1. );
//    c->print();
    
    return(EXIT_SUCCESS);
}
