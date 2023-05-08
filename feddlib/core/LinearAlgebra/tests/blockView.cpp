#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 blockView test

 @brief  blockView
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
    GO numGlobalElements = 2;
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

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef RCP<const Map_Type> MapConstPtr_Type;

    Xpetra::UnderlyingLib ulib;
    if (!ulib_str.compare("Tpetra"))
        ulib = Xpetra::UseTpetra;
    else if (!ulib_str.compare("Epetra"))
        ulib = Xpetra::UseEpetra;


    XpetraMapConstPtr_Type xmap = Xpetra::MapFactory<LO,GO,NO>::Build(ulib, numGlobalElements, 0, comm);
    XpetraMatrixPtr_Type xmatrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(xmap, 1);

    MatrixPtr_Type matrix1 = rcp( new Matrix_Type( xmatrix ) );
    MatrixPtr_Type matrix2 = rcp( new Matrix_Type( xmatrix ) );
    MapConstPtr_Type map = matrix1->getMap();
    for (UN i=0 ; i<matrix1->getNodeNumRows(); i++) {
        Array<GO> indices(1);
        Array<SC> values1(1,1.);
        Array<SC> values2(1,-1.);
        indices[0] = map->getGlobalElement( i );
        matrix1->insertGlobalValues( indices[0], indices(), values1() );
        matrix2->insertGlobalValues( indices[0], indices(), values2() );
    }

    matrix1->fillComplete();
    matrix2->fillComplete();

    BlockMatrixPtr_Type system = rcp(new BlockMatrix_Type(2));

    system->addBlock( matrix1, 0, 0 );
    system->addBlock( matrix2, 1, 1 );

    system->print(VERB_EXTREME);

    system->merge();

    system->printMerge(VERB_EXTREME);

    matrix1->resumeFill();
    MapConstPtr_Type colMap = matrix1->getMap("col");
    for (UN i=0 ; i<matrix1->getNodeNumRows(); i++) {
        Array<LO> indices(1);
        Array<SC> values1(1,2.);
        indices[0] = colMap->getLocalElement( map->getGlobalElement( i ) );
//        matrix1->insertLocalValues( i, indices(), values1() );
    }

    matrix1->fillComplete();
    system->print(VERB_EXTREME);
    system->printMerge(VERB_EXTREME);


    return(EXIT_SUCCESS);
}
