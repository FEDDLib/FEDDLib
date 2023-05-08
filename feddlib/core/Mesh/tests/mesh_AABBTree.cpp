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
mesh_AABB_Ttest test

@brief  This test demonstrates how to create and query an instance of AABBTree
for a given mesh.
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

int main(int argc, char *argv[]) {
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;

    typedef AABBTree<SC,LO,GO,NO> AABBTree_Type;
    typedef RCP<AABBTree_Type > AABBTreePtr_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();


    //*************************************************************************//
    //
    // Load mesh_A and mesh_B, e.g. rectangle.mesh and big_rectangle.mesh
    //
    //*************************************************************************//
    int dim = 2;
    DomainPtr_Type domainA, domainB;
    domainA.reset( new Domain_Type( comm, dim ) );
    string filename = "rectangle.mesh";
    
    {
        //Paralleles P1 mesh (Partitionierung mit metis)
        MeshPartitioner<SC,LO,GO,NO>::DomainPtrArray_Type domainArray(1);
        domainArray[0] = domainA;
        ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
        pListPartitioner->set( "Mesh 1 Name", filename );
        MeshPartitioner<SC,LO,GO,NO> partitioner ( domainArray, pListPartitioner, "P1", dim );
        
        partitioner.readAndPartition();
    }
    domainA->getMesh()->create_AABBTree();

    domainB.reset( new Domain_Type( comm, dim ) );
    filename = "big_rectangle.mesh";
    {
        //Paralleles P1 mesh (Partitionierung mit metis)
        MeshPartitioner<SC,LO,GO,NO>::DomainPtrArray_Type domainArray(1);
        domainArray[0] = domainB;
        ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
        pListPartitioner->set( "Mesh 1 Name", filename );
        MeshPartitioner<SC,LO,GO,NO> partitioner ( domainArray, pListPartitioner, "P1", dim );
        
        partitioner.readAndPartition();
    }
    
    // Find for each point of mesh_B in which element of mesh_A it is
    vec2D_dbl_ptr_Type queryPoints = domainB->getMesh()->getPointsRepeated();
    vec_int_ptr_Type elemsForPoints;
    elemsForPoints = domainA->getMesh()->findElemsForPoints(queryPoints);


    // Print out for each point in which element it is
    ElementsPtr_Type elementsMesh = domainA->getMesh()->getElementsC();
    vec2D_dbl_ptr_Type nodesRepeated = domainA->getMesh()->getPointsRepeated();


    int counterNotInElement = 0;
    for (int point = 0; point < queryPoints->size(); point++){
        if(queryPoints->at(point).at(0) >= 0 && queryPoints->at(point).at(0) <= 1
        && queryPoints->at(point).at(1) >= 0 && queryPoints->at(point).at(1) <= 1){
            // Code will fail here if findElemsForPoints stopped working
            assert(elemsForPoints->at(point) != -1);
        }
        else {
            assert(elemsForPoints->at(point) == -1);
        }

        if (elemsForPoints->at(point) != -1){
            std::cout << "   Query Point " << point << " = (" << queryPoints->at(point).at(0) << ", " << queryPoints->at(point).at(1) << ") is in element " << elemsForPoints->at(point) << '\n';
            vec_int_Type localNodes = elementsMesh->getElement(elemsForPoints->at(point)).getVectorNodeList();
            vec2D_dbl_Type coords(
                localNodes.size(),
                vec_dbl_Type(2, 0.0) // FIXME: this depends on the dimension
            );
            int node = 0;
            for (int localNode=0; localNode<localNodes.size(); localNode++){
                node = localNodes.at(localNode); // get global node_id
                coords.at(localNode) = nodesRepeated->at(node); //get global coordinates
            }

            cout << "      with node " << localNodes.at(0) << " (" <<coords.at(0).at(0) << "," << coords.at(0).at(1) << "), ";
            cout << localNodes.at(1) <<" ("<<coords.at(1).at(0) << "," << coords.at(1).at(1) << "), ";
            cout << localNodes.at(2) << " ("<<coords.at(2).at(0) << "," << coords.at(2).at(1) << ")"  << '\n';
        }
        else {
            counterNotInElement++;
        }
    }
    std::cout << "   Count of points not in any element = " << counterNotInElement << '\n';

    assert(counterNotInElement == 2085);
}
