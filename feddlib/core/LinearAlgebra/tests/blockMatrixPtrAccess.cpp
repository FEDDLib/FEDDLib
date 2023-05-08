#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 Matrix test

 @brief  Matrix test
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

    RCP<const Comm<int> > commWorld = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    GO numGlobalElements = 4;
    myCLP.setOption("nge",&numGlobalElements,"numGlobalElements.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;

    typedef Xpetra::Matrix<SC,LO,GO,NO> XpetraMatrix_Type;
    typedef RCP<XpetraMatrix_Type> XpetraMatrixPtr_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef RCP<Matrix_Type> MatrixPtr_Type;
    
    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    Xpetra::UnderlyingLib ulib;
    if (!ulib_str.compare("Tpetra"))
        ulib = Xpetra::UseTpetra;
    else if (!ulib_str.compare("Epetra"))
        ulib = Xpetra::UseEpetra;


    XpetraMapConstPtr_Type xmap = Xpetra::MapFactory<LO,GO,NO>::Build(ulib, numGlobalElements, 0, commWorld);
    XpetraMatrixPtr_Type xmatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(xmap, 1);
    
    MatrixPtr_Type matrix = rcp( new Matrix_Type( xmatrix ) );

    for (int i=0; i<matrix->getMap()->getNodeNumElements(); i++) {
        Array<SC> values( 1 , 1.);
        Array<GO> indicesCol( 1 , matrix->getMap()->getGlobalElement( i ) );
        matrix->insertGlobalValues( indicesCol[0], indicesCol(), values() );
    }
    matrix->fillComplete();
    matrix->print();

    BlockMatrixPtr_Type bmatrix = Teuchos::rcp(new BlockMatrix_Type( 1 ) );
    bmatrix->addBlock(matrix,0,0);
    
    bmatrix->print();
    
    matrix->resumeFill();
    matrix->scale(2.);
    matrix->fillComplete();
    
    matrix->print();
    bmatrix->print();
    
    return(EXIT_SUCCESS);
}
