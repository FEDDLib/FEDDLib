#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include <Tpetra_Core.hpp>

#include <Teuchos_GlobalMPISession.hpp>

/*!
 Mesh Surface Orientation

 @brief  Mesh Element Flags test
 @version 1.0
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

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
    
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string filename = "cube_tetr6_struct_h=0_2.mesh";
    myCLP.setOption("file",&filename,"Mesh filename");
    int dim = 3;
    myCLP.setOption("dim",&dim,"Dimension");
    string delimiter = " ";
    myCLP.setOption("delimiter",&delimiter,"Delimiter in mesh-file");
	int volumeID = 10;
    myCLP.setOption("volumeID",&volumeID,"Volume ID");
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    // Mesh
    std::string FEType="P2";
    int numProcsCoarseSolve = 0;
    bool boolExportMesh = true;


    DomainPtr_Type domainP1;
    DomainPtr_Type domain;

    
    ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
    pListPartitioner->set( "Mesh 1 Name", filename );
    
    //Paralleles P1 mesh (Partitionierung mit metis)
    domainP1.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
    domainP1Array[0] = domainP1;

    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
    
    partitionerP1.readAndPartition(volumeID);

    domain.reset( new Domain_Type( comm, dim ) );
    domain->buildP2ofP1Domain( domainP1 );

    if (boolExportMesh) {

        MultiVectorPtr_Type mvFlag = rcp(new MultiVector_Type( domainP1->getElementMap() ) );
        ElementsPtr_Type elements = domainP1->getElementsC();
        vec2D_dbl_ptr_Type pointsRep = domainP1->getPointsRepeated();

        Teuchos::ArrayRCP< SC > value = mvFlag->getDataNonConst(0);
        TEUCHOS_TEST_FOR_EXCEPTION( value.size()!= elements->numberElements(), std::runtime_error, "MultiVector and element list have different sizes.");
       
       
        for (UN T=0; T<elements->numberElements(); T++) {
            FiniteElement fe = elements->getElement( T );
            ElementsPtr_Type subEl = fe.getSubElements(); // might be null
            for (int surface=0; surface<fe.numSubElements(); surface++) {
                FiniteElement feSub = subEl->getElement( surface  );
                if(subEl->getDimension() == dim-1){
                    // Setting flag to the placeholder (second last entry). The last entry at (funcParameter.size() - 1) should always be the degree of the surface function
                
                    vec_int_Type nodeList = feSub.getVectorNodeListNonConst ();

            
                    vec_dbl_Type p1(dim),p2(dim),v_E(dim,1.);
                    double sum_v_E = 1.;
                    double norm_v_E = 1.;
                    
                    p1[0] = pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[1]).at(0);
					p1[1] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[1]).at(1);
					p1[2] = pointsRep->at(nodeList[0]).at(2) - pointsRep->at(nodeList[1]).at(2);

					p2[0] = pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[2]).at(0);
					p2[1] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[2]).at(1);
					p2[2] = pointsRep->at(nodeList[0]).at(2) - pointsRep->at(nodeList[2]).at(2);

					v_E[0] = p1[1]*p2[2] - p1[2]*p2[1];
					v_E[1] = p1[2]*p2[0] - p1[0]*p2[2];
					v_E[2] = p1[0]*p2[1] - p1[1]*p2[0];
		            
				    norm_v_E = std::sqrt(std::pow(v_E[0],2)+std::pow(v_E[1],2)+std::pow(v_E[2],2));

                    norm_v_E = std::sqrt(std::pow(v_E[0],2) + std::pow(v_E[1],2) + std::pow(v_E[2],2));
                    sum_v_E = v_E[0] + v_E[1] + v_E[2]; // as we are looking at a cube, all surface normals have on entry, either positive or negative
                   
                    value[T] = sum_v_E/norm_v_E; 
                }
            }
        }
        
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
        
        exPara->setup( "surfaceNormalDirectionOnSurface", domainP1->getMesh(), "P0" );
        MultiVectorConstPtr_Type mvFlagConst = mvFlag;
        exPara->addVariable( mvFlagConst, "flag", "Scalar", 1, domainP1->getElementMap());
        
        exPara->save(0.0);
        exPara->closeExporter();

    }

    return(EXIT_SUCCESS);
}
