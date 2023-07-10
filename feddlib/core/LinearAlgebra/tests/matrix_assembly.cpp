#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

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

    RCP<const Comm<int> > commWorld = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

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

    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;

    typedef Xpetra::Matrix<SC,LO,GO,NO> XpetraMatrix_Type;
    typedef RCP<XpetraMatrix_Type> XpetraMatrixPtr_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef RCP<Matrix_Type> MatrixPtr_Type;

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;

    TEUCHOS_TEST_FOR_EXCEPTION(!(!ulib_str.compare("Tpetra") || !ulib_str.compare("Epetra") ) , std::runtime_error,"Unknown algebra type");


    Teuchos::Array<GO> indices(2);
    for (int i=0; i<indices.size(); i++) {
        indices[i] = i + commWorld->getRank();
    }

    MapPtr_Type map1rep = rcp( new Map_Type( ulib_str, (GO) -1, indices(), 0, commWorld ) );
    MapPtr_Type map1unique = map1rep->buildUniqueMap();

    Teuchos::Array<GO> indices2(3);
    for (int i=0; i<indices2.size(); i++) {
        indices2[i] = i + commWorld->getRank();
    }
    MapPtr_Type map2rep = rcp( new Map_Type( ulib_str, (GO) -1, indices2(), 0, commWorld ) );
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
