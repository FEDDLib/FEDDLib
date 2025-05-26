#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 Map test

 @brief  Map test
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
    GO numGlobalElements = 10;
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

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;

    TEUCHOS_TEST_FOR_EXCEPTION(!(!ulib_str.compare("Tpetra") || !ulib_str.compare("Epetra") ) , std::runtime_error,"Unknown algebra type");


    Array<GO> indices(numGlobalElements);
    for (UN i=0; i<indices.size(); i++) {
        indices[i] = i;
    }

    MapPtr_Type map = rcp( new Map_Type(ulib_str, commWorld->getSize()*numGlobalElements, indices(), 0, commWorld) );

    map->print();
    
    
    // Determine the globalInterfaceIDs of tagged edges
    vec_GO_Type globalInterfaceIDs;
    if(commWorld->getRank()== 0)
        globalInterfaceIDs = {2,7,10,11,12};

    if(commWorld->getRank()== 1)
        globalInterfaceIDs = {2,7,12,13,14};

    if(commWorld->getRank()== 2)
        globalInterfaceIDs = {10,11,12,17,22};

    if(commWorld->getRank()== 3)
        globalInterfaceIDs = {12,13,14,17,22};
   
    tuple_intint_Type rankRange_;
    get<0>(rankRange_) = 0;
    get<1>(rankRange_) = commWorld->getSize() - 1;

    MapPtr_Type mapInterfaceNodes = rcp( new Map_Type( ulib_str, Teuchos::OrdinalTraits<GO>::invalid(), globalInterfaceIDs, 0, commWorld) );
    MapPtr_Type mapInterfaceNodesUnique = mapInterfaceNodes->buildUniqueMap( rankRange_ );

    // Multivector containing only Uniquely distributed entries
    MultiVectorPtr_Type nodesUnique = Teuchos::rcp( new MultiVector_Type(mapInterfaceNodesUnique, 1 ) );
    nodesUnique->putScalar(0);

    // Multivector containing only Repeatedly distributed entries
    MultiVectorPtr_Type nodesRepeated = Teuchos::rcp( new MultiVector_Type( mapInterfaceNodes, 1 ) );
    nodesRepeated->putScalar(1); // INSERT REAL NODE VALUES HERE

    nodesUnique->exportFromVector(nodesRepeated, true, "Add"); // All repeated values are added in unique vector
    nodesUnique->print();

    // Now the uniquely distributed result is send back to the repeated distribution
    MultiVectorPtr_Type resultDist = Teuchos::rcp( new MultiVector_Type( mapInterfaceNodes, 1 ) );
    resultDist->putScalar(0);

    resultDist->importFromVector(nodesUnique, true, "Insert");
    resultDist->print();





    return(EXIT_SUCCESS);
}
