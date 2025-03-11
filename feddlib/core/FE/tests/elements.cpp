#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/EntitiesOfElements.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>


#include "feddlib/core/FE/Elements.hpp"

using namespace std;
using namespace Teuchos;
using namespace FEDD;

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

int main(int argc, char *argv[]) {

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }
    
    bool verbose(comm->getRank()==0);
    
    {
        if (verbose)
            std::cout << "-- Testing elements --" << std::endl;
        Elements element;
        vec_int_Type localNodeList1(3,1);
        localNodeList1[2] = 2;
        vec_int_Type localNodeList2(3,2);
        vec_int_Type localNodeList3(3,1);
        localNodeList3[2] = 2;
        vec_int_Type localNodeList4(3,1);
        localNodeList4[1] = 2;
        vec_int_Type localNodeList5(3,1);
        localNodeList5[2] = 2;
        FiniteElement fe1( localNodeList1 );
        FiniteElement fe2( localNodeList2 );
        FiniteElement fe3( localNodeList3 );
        FiniteElement fe4( localNodeList4 );
        FiniteElement fe5( localNodeList5 );
        element.addElement( fe1 );
        element.addElement( fe2 );
        element.addElement( fe3 );
        element.addElement( fe4 );
        element.addElement( fe5 );
        if (verbose)
            std::cout << "-- Pre unique --" << std::endl;
            
        element.print();
        
        element.sortUnique();

        if (verbose)
            std::cout << "-- After unique --" << std::endl;
        
        element.print();
    }

    if (verbose)
        std::cout << "\n\n";

    {
        if (verbose)
            std::cout << "-- Testing sub elements --" << std::endl;
    
        Elements elements("P1",3);
        
        vec_LO_Type localNodeIDs(4,0);
        localNodeIDs[1] = 1; localNodeIDs[2] = 2; localNodeIDs[3] = 3;
        FiniteElement fe1( localNodeIDs );
        
        elements.addElement(fe1);
        elements.getElement(0).initializeSubElements("P1",2);
        
        elements.print();
        vec_LO_Type subLocalNodeIDs(3,0);
        subLocalNodeIDs[1] = 1; subLocalNodeIDs[2] = 2;
        FiniteElement feSub1( subLocalNodeIDs );
        vec_LO_Type subLocalNodeIDs2(3,0);
        subLocalNodeIDs2[1] = 1; subLocalNodeIDs2[2] = 3;
        FiniteElement feSub2( subLocalNodeIDs2 );
        
        elements.getElement(0).addSubElement( feSub1 );
        elements.getElement(0).addSubElement( feSub2 );
        
    }
    
    
    
    return(EXIT_SUCCESS);
}
