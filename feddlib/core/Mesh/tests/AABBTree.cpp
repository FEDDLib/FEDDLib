#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/Elements.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/AABBTree.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
/*!
AABBTree test

@brief  This tests demonstrates how to create an instance of AABBTree and
query it. It also verifies the basis functions of the class AABBTree.
@author Viktor Grimm
@version 1.0
@copyright VG
*/
using namespace std;
using namespace Teuchos;
using namespace FEDD;

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
typedef Elements Elements_Type;
typedef RCP<Elements_Type> ElementsPtr_Type;
typedef Domain<SC,LO,GO,NO> Domain_Type;
typedef RCP<Domain_Type > DomainPtr_Type;

typedef AABBTree<SC,LO,GO,NO> AABBTree_Type;
typedef RCP<AABBTree_Type > AABBTreePtr_Type;

int main(int argc, char *argv[]) {

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();


    //*************************************************************************//
    //
    // Load simple mesh, e.g. rectangle.mesh
    //
    //*************************************************************************//
    int dim = 2;
    DomainPtr_Type domain;
    domain.reset( new Domain_Type( comm, dim ) );
    string filename = "rectangle.mesh";

    //Paralleles P1 mesh (Partitionierung mit metis)
    domain.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner<SC,LO,GO,NO>::DomainPtrArray_Type domainArray(1);
    domainArray[0] = domain;
    ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
    pListPartitioner->set( "Mesh 1 Name", filename );
    MeshPartitioner<SC,LO,GO,NO> partitioner ( domainArray, pListPartitioner, "P1", dim );
    
    partitioner.readAndPartition();

    //*************************************************************************//
    //
    // Create Tree
    //
    //*************************************************************************//
    ElementsPtr_Type elementsMesh = domain->getMesh()->getElementsC();
    vec2D_dbl_ptr_Type nodesRepeated = domain->getMesh()->getPointsRepeated();


    int numberObjects = 32;
    double longTol = 0.75;
    double volumeTol = 0.55;
    bool verbose = true;
    
    AABBTreePtr_Type testTree(new AABBTree_Type());
    testTree->createTreeFromElements(
        elementsMesh,
        nodesRepeated,
        numberObjects,
        longTol,
        volumeTol,
        verbose
    );

    assert(testTree->getNumNodes() == 17);

    //*************************************************************************//
    //
    // Query tree
    //
    //*************************************************************************//
    vec2D_dbl_ptr_Type queryPoints(
        new vec2D_dbl_Type(
            5,
            vec_dbl_Type(2, 0.0)
        )
    );
    for (int point = 0; point < 5; point++){
        queryPoints->at(point).at(0) = point * 0.2;
        queryPoints->at(point).at(1) = point * 0.2;
    }
    map<int, list<int> > treeToItem;
    map<int, list<int> > itemToTree;
    tie(treeToItem, itemToTree) = testTree->scanTree(queryPoints, false);

    for (auto keyValue: treeToItem){
        std::cout << "  node " << keyValue.first << " contains points = " ;
        for (auto value: keyValue.second){
            std::cout << value <<" " ;
        }
        std::cout << " " << '\n';
    }

    for (auto keyValue: itemToTree){
        std::cout << "  point " << keyValue.first << " is in nodes = " ;
        for (auto value: keyValue.second){
            std::cout << value <<" " ;
        }
        std::cout << " " << '\n';
    }

}
