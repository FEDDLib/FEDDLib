#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/BlockMultiVector.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/problems/Solver/PrecOpFaCSI.hpp"

using namespace std;
using namespace Teuchos;
using namespace FEDD;

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

int main(int argc, char *argv[]) {

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
   
    int option = 1;
    myCLP.setOption("option",&option,"option");
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }

    TEUCHOS_TEST_FOR_EXCEPTION( option>3 || option<1, std::logic_error, "Chooes option beteween 1 and 3 for FaCSI test.");
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef RCP<const Map_Type> MapConstPtr_Type;
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
    typedef Teuchos::RCP<const Matrix_Type> MatrixConstPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
    typedef Teuchos::RCP<const Thyra::LinearOpBase<SC> > ThyraConstLinOpPtr_Type;
    typedef Teuchos::RCP<Thyra::LinearOpBase<SC> > ThyraLinOpPtr_Type;
    
    
    int rank = comm->getRank();
    
    Teuchos::Array< GO> elementList1( 0 );
    if ( comm->getRank() == 0)
        elementList1.push_back(0);
    Teuchos::Array< GO> elementList2( 0 );
    if ( comm->getRank() == 1)
        elementList2.push_back(0);

    
    MapConstPtr_Type map1 = Teuchos::rcp( new Map_Type( -1, elementList1(), 0, comm) );
    MapConstPtr_Type map2 = Teuchos::rcp( new Map_Type( -1, elementList2(), 0, comm) );
    
    MatrixPtr_Type mm1 = Teuchos::rcp( new Matrix_Type( map1, 1 )  );
    Array<GO> col;
    Array<SC> val;
    if (rank==0) {
        col.push_back(0);
        val.push_back(1.);
    }
    mm1->insertGlobalValues( 0, col(), val() );
    mm1->fillComplete();
    ThyraConstLinOpPtr_Type mc = mm1->getThyraLinOp();
    ThyraLinOpPtr_Type m = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> > (mc);

    MatrixPtr_Type mm2 = Teuchos::rcp( new Matrix_Type( map2, 1 )  );
    Array<GO> col2;
    Array<SC> val2;
    if (rank==1) {
        col2.push_back(0);
        val2.push_back(2.);
    }
    mm2->insertGlobalValues( 0, col2(), val2() );
    mm2->fillComplete();
    ThyraConstLinOpPtr_Type mc2 = mm2->getThyraLinOp();
    ThyraLinOpPtr_Type m2 = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> > (mc2);

    MatrixPtr_Type mm3 = Teuchos::rcp( new Matrix_Type( map1, 1 )  );
    Array<GO> col3;
    Array<SC> val3;
    if (rank==1) {
        col3.push_back(0);
        val3.push_back(3.);
    }
    mm3->insertGlobalValues( 0, col3(), val3() );
    mm3->fillComplete( map2, map1 );
    ThyraConstLinOpPtr_Type mc3 = mm3->getThyraLinOp();
    ThyraLinOpPtr_Type m3 = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> > (mc3);
    
    mm1->print();
    mm2->print();
    mm3->print();
    
    
    BlockMatrixPtr_Type bmm1 = Teuchos::rcp( new BlockMatrix_Type( 2 ) );
    bmm1->addBlock(mm1,0,0);
    bmm1->addBlock(mm1,1,1);
    ThyraConstLinOpPtr_Type bmc = bmm1->getThyraLinOp();
    ThyraLinOpPtr_Type bm = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> > (bmc);

    PrecOpFaCSI<SC,LO,GO,NO> facsi( comm, true );

    // This is just a simple test so we use the same matrix in every block
    if (option==1)
        facsi.setGE( m, m, m3, m2, bm, m, m );
//    else if(option==2)
//        facsi.setGI( m, m, m, m, m, m, m, m );
//    else if(option==3)
//        facsi.setGIShape( m, m, m, m, m, m, m, m );

    BlockMultiVectorPtr_Type X = Teuchos::rcp( new BlockMultiVector_Type( 4 ) );
    BlockMultiVectorPtr_Type Y = Teuchos::rcp( new BlockMultiVector_Type( 4 ) );
    if (option==1) {
        MultiVectorPtr_Type X1 = Teuchos::rcp( new MultiVector_Type( map1 ) ); X1->putScalar(1.);
        MultiVectorPtr_Type X2 = Teuchos::rcp( new MultiVector_Type( map1 ) ); X2->putScalar(2.);
        MultiVectorPtr_Type X3 = Teuchos::rcp( new MultiVector_Type( map2 ) ); X3->putScalar(3.);
        MultiVectorPtr_Type X4 = Teuchos::rcp( new MultiVector_Type( map1 ) ); X4->putScalar(4.);
        MultiVectorPtr_Type Y1 = Teuchos::rcp( new MultiVector_Type( map1 ) );
        MultiVectorPtr_Type Y2 = Teuchos::rcp( new MultiVector_Type( map1 ) );
        MultiVectorPtr_Type Y3 = Teuchos::rcp( new MultiVector_Type( map2 ) );
        MultiVectorPtr_Type Y4 = Teuchos::rcp( new MultiVector_Type( map1 ) );
        
        X->addBlock(X1,0); X->addBlock(X2,1); X->addBlock(X3,2); X->addBlock(X4,3);
        Y->addBlock(Y1,0); Y->addBlock(Y2,1); Y->addBlock(Y3,2); Y->addBlock(Y4,3);
    }
    
    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraX = X->getProdThyraMultiVector();
    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraY = Y->getProdThyraMultiVector();

    facsi.applyIt(Thyra::NOTRANS, *thyraX, thyraY.ptr(), 1., 0.);

    Y->getBlockNonConst(0)->fromThyraMultiVector( thyraY->getNonconstMultiVectorBlock( 0 ) );
    Y->getBlockNonConst(1)->fromThyraMultiVector( thyraY->getNonconstMultiVectorBlock( 1 ) );
    Y->getBlockNonConst(2)->fromThyraMultiVector( thyraY->getNonconstMultiVectorBlock( 2 ) );
    Y->getBlockNonConst(3)->fromThyraMultiVector( thyraY->getNonconstMultiVectorBlock( 3 ) );
    
    Y->print();
    return(EXIT_SUCCESS);
}
