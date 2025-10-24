#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>

/*!
 Map test

 @brief  Map test
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


using namespace std;
using namespace Teuchos;

using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::outArg;

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
using namespace FEDD;
int main(int argc, char *argv[]) {

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > commWorld = Tpetra::getDefaultComm();
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    LO numLocalElements = 10;
   // LO numGlobalElements = 10*(commWorld->getSize()-2);

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;


    UN lowerOffset = 10 * commWorld->getRank() - commWorld->getRank()*2 ;
    UN upperOffset = 10*(commWorld->getRank()+1);

    LO numLocalIDs = upperOffset-lowerOffset;

    Array<GO> indices(numLocalIDs);

    for (UN i=0; i<numLocalIDs; i++) {
        indices[i] = lowerOffset+i;
    }
    GO numGlobalIDs = 0;

	reduceAll<int, GO> (*commWorld, REDUCE_SUM, numLocalIDs, outArg (numGlobalIDs));

    MapPtr_Type map = rcp( new Map_Type(numGlobalIDs, indices(), 0, commWorld) );

    map->print();

    MapPtr_Type mapUnique =  map->buildUniqueMap();

    mapUnique->print();


    return(EXIT_SUCCESS);
}
