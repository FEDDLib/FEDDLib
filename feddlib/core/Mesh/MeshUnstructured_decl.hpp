#ifndef MESHUNSTRUCTURED_decl_hpp
#define MESHUNSTRUCTURED_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "Mesh.hpp"
#include "MeshInterface.hpp"
#include "MeshFileReader.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/core/FE/TriangleElements.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
/*!
 Declaration of MeshUnstructured
 
 @brief  MeshUnstructured
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

/*! 
	Extension of mesh class. MeshUnstructured represents meshes that are read through .mesh input file. In contrast to structured meshes which are 
	meshes we build internally. Features functions to read a mesh and extend it to a P2 mesh.
*/



namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshUnstructured : public Mesh<SC,LO,GO,NO> {
    
public:
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshUnstrPtr_Type;
    typedef Teuchos::RCP<Mesh_Type> MeshPtr_Type;

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

	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    MeshUnstructured();
    
    MeshUnstructured( CommConstPtr_Type comm, int volumeID=10 );
    
    ~MeshUnstructured();
       
	/*! 
		\brief Function to build a P2 mesh of a P1 mesh
		@param[in] meshP1 The p1 mesh we use for building P2 mesh
	*/
    void buildP2ofP1MeshEdge( MeshUnstrPtr_Type meshP1 );

	/*! 
		\brief Adding the correct surface subelement to the new P2 Elements based on the P1 subelements
		@param[in] meshP1 The p1 mesh we use for building P2 mesh
	*/
    void setP2SurfaceElements( MeshUnstrPtr_Type meshP1 );
    
	/*! 
		\brief Helper function for setP2SurfaceElements. Adds the correct nodes to the meshP1 subelements. Based on sorted Elements.
		@param[in] feP2 P2 element 
		@param[in] surfFeP1 P1 surface element that need new P2 nodes
		@param[in] surfacePermutation Surface permutations of element
		@param[in] dim Dimension
	*/
    void setSurfaceP2( FiniteElement &feP2, const FiniteElement &surfFeP1, const vec2D_int_Type &surfacePermutation, int dim );
    
	/*! 
		\brief Helper function for setP2SurfaceElements. Adds the correct nodes to the meshP1 subelements. Based on unsorted elements. Different approach than above. Based on edge Midpoints.

		@param[in] feP2 P2 element 
		@param[in] surfFeP1 P1 surface element that need new P2 nodes
		@param[in] surfacePermutation Surface permutations of element
		@param[in] dim Dimension
	*/
	void addSurfaceP2Nodes( FiniteElement &feP2, const FiniteElement &surfFeP1, const vec2D_int_Type &surfacePermutation, int dim );
	/*! 
		\brief Depending on the sorting of P1 surface nodes we have to adjust the new ordering of P2 edge midpoints for surfaces in 3D
	*/	
    vec_int_Type reorderP2SurfaceIndices( vec_int_Type& additionalP2IDs, vec_int_Type& index, bool track=false);
    
	/*! 
		\brief Get local Surface Indices
	*/
    void getLocalSurfaceIndices( vec2D_int_Type& surfacePermutation, int surfaceElementOrder );
    
	/*! 
		\brief Get edge combinations
	*/
    void getEdgeCombinations( vec2D_int_Type& edgeCombinations );
        
	/*! 
		\brief Determine position of new P2 node in element, as all elements should follow the same structure
		@param[in] positions
		@param[in] elementsGlobalOfEdge Global IDs of elements connected to the edge (numbering needs to be consistent)
		@param[in] p1ID Local ID of first node of edge
		@param[in] p2ID Local ID of second node of edge
		@param[in] meshP1 The p1 mesh we use for building P2 mesh
	*/
    void determinePositionInElementP2( vec_int_Type& positions, vec_GO_Type& elementsGlobalOfEdge, LO p1ID, LO p2ID, MeshUnstrPtr_Type meshP1 );
    
	/*! 
		\brief Essentially determine flag of an edge. Thus determine P2 Flag for building P2 mesh.
		@param[in] fe Element the edge belongs to. Use to look through subelements
		@param[in] p1ID Local ID of first node of edge
		@param[in] p2ID Local ID of second node of edge
		@param[in] permutation Edge combinations of element
	*/
    int determineFlagP2( FiniteElement& fe, LO p1ID, LO p2ID, vec2D_int_Type& permutation );
    
	/*! 
		\brief Essentially determine flag of an edge. Thus determine P2 Flag for building P2 mesh. Longer version of other determineFlagP2 function. 
				Also informs whether flag could be found or not (not always the case in parallel)
		@param[in] p1ID Local ID of first node of edge
		@param[in] p2ID Local ID of second node of edge
		@param[in] localEdgeID LocalEdge ID of second node of edge
		@param[in] markedPoint Edge information, if no flag was found
	*/
    int determineFlagP2( LO p1ID, LO p2ID,  LO localEdgeID, vec2D_LO_Type& markedPoint );
    
    void getTriangles(int vertex1ID, int vertex2ID, vec_int_Type &vertices3ID);

    SurfaceElementsPtr_Type getSurfaceTriangleElements(){return surfaceTriangleElements_;}
    
    void findSurfaces( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localSurfaceNodeList_vec, vec_int_Type& locSurfaces, bool critical = false );

	/*! 
		\brief Determine which edges belong to an element
		@param[in] elementNodeList local node IDs of one element 
		@param[in] numbering
		@param[in] localEdgeNodeList_vec local node IDs of edges
		@param[in] locEdges vector that stores the local IDs of edges belonging to element of elementNodeList
	*/
    void findEdges( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localEdgeNodeList_vec, vec_int_Type& locEdges);
    
	/*! 
		\brief Get mesh interface
		\return meshInterface_
	*/
    MeshInterfacePtr_Type getMeshInterface();        
    
    void buildMeshInterfaceParallelAndDistance( MeshUnstrPtr_Type mesh, vec_int_Type flag_vec, vec_dbl_ptr_Type &distancesToInterface );
    
    void partitionInterface();

    /*!
            \brief setEdgeElements with external edges
            @param[in] edgeElements
    */
    void setEdgeElements(EdgeElementsPtr_Type edgeElements) { edgeElements_ = edgeElements; }

    /*!
            \brief Get EdgeElements
            \return edgeElements_
    */
    EdgeElementsPtr_Type getEdgeElements() { return edgeElements_; }

    /*!
            \brief Get SurfaceEdgeElements. Edges as only surface elements (i.e. when reading .mesh file). Used in mesh
       partitioner \return surfaceEdgeElements_
    */
    ElementsPtr_Type getSurfaceEdgeElements() { return surfaceEdgeElements_; }

    /*!
            \brief Reading mesh size
    */
    void readMeshSize();

    /*!
            \brief Reading the .mesh files entities
            @param[in] entityType i.e. nodes, edges, elements...
    */
    void readMeshEntity(std::string entityType);

    /*!
            \brief Set the .mesh file name
    */
    void setMeshFileName(std::string meshFileName, std::string delimiter);

    /*!
            \brief Get global number of nodes
            \return numNodes_
    */
    int getNumGlobalNodes() { return numNodes_; }

    /*!
            \brief Assigning flags to all edges
    */
    void assignEdgeFlags();

	/*! 
		\brief Building an edgemap from scratch when edges are already distributed parallel
	*/
	void buildEdgeMap();

	/*! 
		\brief Exporting Mesh as .mesh file. For most detailed export we also write surfaces and edges. This can be useful/necessary
		@param[in] meshName for export
		@param[in] exportEdges wether to export edges or not
		@param[in] exportSurfaces whether to export surfaces or not

	*/
	void exportMesh(MapConstPtr_Type mapUnique, MapConstPtr_Type mapRep, bool exportEdges=false, bool exportSurface=false, std::string meshName="export.mesh");

	/*!
	
	
	
	
	*/
	void exportNodeFlags();
	
	/*!
	
	
	
	
	*/
	void exportElementFlags();

    /* ###################################################################### */
    
    MeshInterfacePtr_Type meshInterface_;
        
    int volumeID_;
    /* ###################################################################### */


 	EdgeElementsPtr_Type edgeElements_;    
    ElementsPtr_Type surfaceEdgeElements_;
	SurfaceElementsPtr_Type surfaceTriangleElements_;

 	std::string meshFileName_;
    std::string delimiter_;

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
