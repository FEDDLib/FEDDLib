#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>

/*!
 Matrix assembly test

 @brief  Matrix assembly
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

    Teuchos::RCP<const Teuchos::Comm<int> > commWorld = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    
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

    typedef Tpetra::Map<LO,GO,NO> TpetraMap_Type;
    typedef RCP<TpetraMap_Type> TpetraMapPtr_Type;
    typedef RCP<const TpetraMap_Type> TpetraMapConstPtr_Type;

    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
    typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef RCP<Matrix_Type> MatrixPtr_Type;

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;

    Teuchos::Array<GO> indices(2);
    for (int i=0; i<indices.size(); i++) {
        indices[i] = i + commWorld->getRank();
    }

    MapPtr_Type map1rep = rcp( new Map_Type( (GO) -1, indices(), 0, commWorld ) );
    MapPtr_Type map1unique = map1rep->buildUniqueMap();

    Teuchos::Array<GO> indices2(3);
    for (int i=0; i<indices2.size(); i++) {
        indices2[i] = i + commWorld->getRank();
    }
    MapPtr_Type map2rep = rcp( new Map_Type( (GO) -1, indices2(), 0, commWorld ) );
    MapPtr_Type map2unique = map2rep->buildUniqueMap();

    MatrixPtr_Type matrix = rcp( new Matrix_Type( map2unique, 3 ) );
    for (int i=0; i<map2rep->getNodeNumElements(); i++) {
        Array<SC> values( map1unique->getGlobalNumElements() , 1.);
        Array<GO> indicesCol( map1unique->getGlobalNumElements() , 0);
        for (int j=0; j<indicesCol.size(); j++) {
            indicesCol[j]  = j;
        }
        matrix->insertGlobalValues( map2rep->getGlobalElement(i), indicesCol(), values() );
    }

    matrix->fillComplete( map1unique, map2unique );
    matrix->print(VERB_EXTREME);

    return(EXIT_SUCCESS);
}
