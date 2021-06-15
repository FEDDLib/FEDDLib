#ifndef MESHUNSTRUCTURED_decl_hpp
#define MESHUNSTRUCTURED_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "Mesh.hpp"
#include "MeshInterface.hpp"
#include "MeshFileReader.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/FE/TriangleElements.hpp"

/*!
 Declaration of MeshUnstructured
 
 @brief  MeshUnstructured
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshUnstructured : public Mesh<SC,LO,GO,NO> {
    
public:
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshUnstrPtr_Type;

    typedef std::vector<MeshUnstrPtr_Type> MeshUnstrPtrArray_Type;

    typedef typename Mesh_Type::CommPtr_Type CommPtr_Type;
    typedef typename Mesh_Type::CommConstPtr_Type CommConstPtr_Type;
    typedef typename Mesh_Type::Elements_Type Elements_Type;
    typedef typename Mesh_Type::ElementsPtr_Type ElementsPtr_Type;    
    
    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;

    typedef SurfaceElements SurfaceElements_Type;
    typedef Teuchos::RCP<SurfaceElements_Type> SurfaceElementsPtr_Type;
    
    typedef MeshInterface<SC,LO,GO,NO> MeshInterface_Type;
    typedef Teuchos::RCP<MeshInterface_Type> MeshInterfacePtr_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef Teuchos::OrdinalTraits<LO> OTLO;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef MultiVector<LO,LO,GO,NO> MultiVectorLO_Type;
	typedef Teuchos::RCP<MultiVectorLO_Type> MultiVectorLOPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;

    MeshUnstructured();
    
    MeshUnstructured( CommConstPtr_Type comm, int volumeID=10 );
    
    ~MeshUnstructured();
    
//    virtual vec2D_int_ptr_Type getElements();
    
    virtual void dummy() {};
     
    void buildP2ofP1MeshEdge( MeshUnstrPtr_Type meshP1 );

    void setP2SurfaceElements( MeshUnstrPtr_Type meshP1 );
    
    void setSurfaceP2( FiniteElement &feP2, const FiniteElement &surfFeP1, const vec2D_int_Type &surfacePermutation, int dim );
    
    vec_int_Type reorderP2SurfaceIndices( vec_int_Type& additionalP2IDs, vec_int_Type& index, bool track=false);
    
    void getLocalSurfaceIndices( vec2D_int_Type& surfacePermutation, int surfaceElementOrder );
    
    void getEdgeCombinations( vec2D_int_Type& edgeCombinations );
        
    void determinePositionInElementP2( vec_int_Type& positions, vec_GO_Type& elementsGlobalOfEdge, LO p1ID, LO p2ID, MeshUnstrPtr_Type meshP1 );
    
    int determineFlagP2( FiniteElement& fe, LO p1ID, LO p2ID, vec2D_int_Type& permutation );
    
    int determineFlagP2( MeshUnstrPtr_Type meshP1, LO p1ID, LO p2ID,  LO localEdgeID, vec2D_LO_Type& markedPoint );
    
    void getTriangles(int vertex1ID, int vertex2ID, vec_int_Type &vertices3ID);

    SurfaceElementsPtr_Type getSurfaceTriangleElements(){return surfaceTriangleElements_;};
    
    void findSurfaces( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localSurfaceNodeList_vec, vec_int_Type& locSurfaces, bool critical = false );

    void findEdges( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localEdgeNodeList_vec, vec_int_Type& locEdges);
    
    MeshInterfacePtr_Type getMeshInterface();        
    
    void buildMeshInterfaceParallelAndDistance( MeshUnstrPtr_Type mesh, vec_int_Type flag_vec, vec_dbl_ptr_Type &distancesToInterface );
    
    void partitionInterface();
    
    void setEdgeElements( EdgeElementsPtr_Type edgeElements ){ edgeElements_ = edgeElements; };
    
    EdgeElementsPtr_Type getEdgeElements( ){ return edgeElements_; };
    
    ElementsPtr_Type getSurfaceEdgeElements(){return surfaceEdgeElements_;};
    
    void readMeshSize();
    
    void readMeshEntity(string entityType);
    
    void setMeshFileName(string meshFileName, string delimiter);
    
    int getSurfaceElementOrder(){return surfaceElementOrder_;};
    
    int getEdgeElementOrder(){return edgesElementOrder_;};
    
    int getNumGlobalNodes(){return numNodes_;};
    
    /* ###################################################################### */
    
    MeshInterfacePtr_Type meshInterface_;
        
    int volumeID_;
    /* ###################################################################### */


 	EdgeElementsPtr_Type edgeElements_;    
    ElementsPtr_Type surfaceEdgeElements_;
	SurfaceElementsPtr_Type surfaceTriangleElements_;

 	string meshFileName_;
    string delimiter_;

    int elementOrder_;
    int surfaceElementOrder_;
    int edgesElementOrder_;
    int numElements_;
    int numSurfaces_;
    int numEdges_;
    int numNodes_;

private:
    
    void readSurfaces();
    
    void readLines();
    
    void readElements();
    
    void readNodes();
     
};
}
#endif
