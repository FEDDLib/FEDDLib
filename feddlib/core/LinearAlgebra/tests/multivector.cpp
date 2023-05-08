#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 MultiVector test

 @brief  MultiVector test
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
    int rank = commWorld->getRank();
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    GO numGlobalElements = 3;
    myCLP.setOption("nge",&numGlobalElements,"numGlobalElements.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }


    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef RCP<const Map_Type> MapConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MV_Type;
    typedef RCP<MV_Type> MVPtr_Type;

    TEUCHOS_TEST_FOR_EXCEPTION(!(!ulib_str.compare("Tpetra") || !ulib_str.compare("Epetra") ) , std::runtime_error, "Unknown algebra type");

    Array<GO> indices(numGlobalElements);
    for (UN i=0; i<indices.size(); i++) {
        indices[i] = i;
    }

    MapConstPtr_Type mapRepeated = rcp( new Map_Type(ulib_str, commWorld->getSize()*numGlobalElements, indices(), 0, commWorld) );

    MapConstPtr_Type mapUnique = mapRepeated->buildUniqueMap();

    
    MVPtr_Type mvRep = rcp( new MV_Type( mapRepeated ) );
    MVPtr_Type mvUni = rcp( new MV_Type( mapUnique ) );

    mvUni->putScalar( rank + 1. );
    mvUni->print();

    mvRep->importFromVector(mvUni);
    mvRep->print();
    mvRep->putScalar( 0. );
    mvRep->exportFromVector(mvUni);
    mvRep->print();
    
    return(EXIT_SUCCESS);
}
