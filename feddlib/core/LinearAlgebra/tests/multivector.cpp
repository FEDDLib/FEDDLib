#define MAIN_ASSERT(A,S) if(!(A)) { cerr<<"Assertion failed. "<<S<<endl; cout.flush(); throw out_of_range("Assertion.");};
#define VERBOSE

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_Core.hpp>
#include <Thyra_DefaultZeroLinearOp_decl.hpp>

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

    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > commWorld = Tpetra::getDefaultComm();

    int rank = commWorld->getRank();
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    GO numGlobalElements = 4;
    myCLP.setOption("nge",&numGlobalElements,"numGlobalElements.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    std::cout << "SC = " << typeid(SC).name() << std::endl;
    std::cout << "LO = " << typeid(LO).name() << std::endl;
    std::cout << "GO = " << typeid(GO).name() << std::endl;
    std::cout << "NO = " << typeid(NO).name() << std::endl;

    auto demangle = [](const char* name){
      int st;
      std::unique_ptr<char, void(*)(void*)> res{
        abi::__cxa_demangle(name, nullptr, nullptr, &st),
        std::free
      };
      return std::string(st == 0 ? res.get() : name);
    };

    std::cout << "Default SC = " << demangle(typeid(SC).name()) << "\n";
    std::cout << "Default LO = " << demangle(typeid(LO).name()) << "\n";
    std::cout << "Default GO = " << demangle(typeid(GO).name()) << "\n";
    std::cout << "Default NO = " << demangle(typeid(NO).name()) << "\n";

    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef RCP<const Map_Type> MapConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MV_Type;
    typedef RCP<MV_Type> MVPtr_Type;

    Array<GO> indices(numGlobalElements);
    for (UN i=0; i<indices.size(); i++) {
        indices[i] = i;
    }

    MapConstPtr_Type mapRepeated = rcp( new Map_Type( commWorld->getSize()*numGlobalElements, indices(), 0, commWorld) );

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
    
    // Scalar Type tests
    typedef MultiVector<LO,LO,GO,NO> MVLO_Type;
    typedef RCP<MVLO_Type> MVLOPtr_Type;
    
    MVLOPtr_Type mvLO = rcp( new MVLO_Type( mapUnique,1 ) );
    mvLO->putScalar( rank + 1 );
    mvLO->print();
	Teuchos::ArrayRCP< SC > flagExportEntries  = mvLO->getDataNonConst(0);
    // 



    return(EXIT_SUCCESS);
}
