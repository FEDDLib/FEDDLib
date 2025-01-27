#ifndef Mesh_decl_hpp
#define Mesh_decl_hpp

#define FULL_Mesh_TIMER

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Elements.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/FE/TriangleElements.hpp"
#include "feddlib/core/FE/FiniteElement.hpp"
#include "feddlib/core/FE/Helper.hpp"
#include "feddlib/core/Mesh/AABBTree.hpp"

/*!
Defintion of Mesh

@brief  Mesh
@author Christian Hochmuth
@version 1.0
@copyright CH
*/

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Mesh {
    
public:
    typedef Elements Elements_Type;
    typedef FiniteElement FiniteElement_Type;
    typedef Teuchos::RCP<FiniteElement_Type> FiniteElementPtr_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;

    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;

    typedef SurfaceElements SurfaceElements_Type;
    typedef Teuchos::RCP<SurfaceElements_Type> SurfaceElementsPtr_Type;
    
    typedef Teuchos::RCP<Mesh> Mesh_ptr_Type;
    
    typedef Teuchos::RCP<Teuchos::Comm<int> > CommPtr_Type;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > CommConstPtr_Type;
    typedef const CommConstPtr_Type CommConstPtrConst_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtrConst_Type;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef MultiVector<SC,LO,GO,NO> MultiVectorLO_Type;
	typedef Teuchos::RCP<MultiVectorLO_Type> MultiVectorLOPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVectorGO_Type;
    typedef Teuchos::RCP<MultiVectorGO_Type> MultiVectorGOPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;
    typedef Teuchos::OrdinalTraits<LO> OTLO;

	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef AABBTree<SC,LO,GO,NO> AABBTree_Type;
    typedef Teuchos::RCP<AABBTree_Type > AABBTreePtr_Type;
    
    /* ###################################################################### */
    //
    Mesh();

    Mesh(CommConstPtrConst_Type& comm);
    
    ~Mesh();
    
    /*!
     Delete all member variables
     */
    void deleteData(){};
    
    /// @brief Setting input parameterlist to be parameterlist here
    /// @param pL 
    void setParameterList( ParameterListPtr_Type& pL );

    /// @brief Getter for paramaeterlist that is set here
    /// @return pL
    ParameterListConstPtr_Type getParameterList( ) const;
    
    /// @brief Getter for element flags
    /// @return elements flag vector -- not sure if this is used anywhere
    vec_int_ptr_Type getElementsFlag() const;
    
    /// @brief Getter for unique node map
    /// @return mapUnique_
    MapConstPtr_Type getMapUnique() const;
    
    /// @brief Getter for repeated node mal
    /// @return mapRepeated_
    MapConstPtr_Type getMapRepeated() const;

    /// @brief Getter for unique node P2 map. Dont know what this is for exactly. Think this object is empty
    /// @return 
    MapConstPtr_Type getMapUniqueP2() const;
    
    /// @brief Getter for repeated node P2 map. Dont know what this is for exactly. Think this object is empty
    /// @return 
    MapConstPtr_Type getMapRepeatedP2() const;
    
    /// @brief Getter for element map 
    /// @return 
    MapConstPtr_Type getElementMap();
	
    /// @brief Getter for edge map
    /// @return 
    MapConstPtr_Type getEdgeMap(); // Edge Map
    
    /// @brief getter for list of repeated points with x,y,z coordinates in each row
    /// @return pointsRep_
    vec2D_dbl_ptr_Type getPointsRepeated() const;

    /// @brief getter for list of unique points with x,y,z coordinates in each row
    /// @return pointsUni_
    vec2D_dbl_ptr_Type getPointsUnique() const;
    
    /// @brief Getter for flags corresponting to repeated points
    /// @return bcFlagRep_
    vec_int_ptr_Type getBCFlagRepeated() const;
    
    /// @brief Getter for flags corresponting to unique points
    /// @return bcFlagUni_
    vec_int_ptr_Type getBCFlagUnique() const;
        
    virtual void dummy() = 0;
    
    /// @brief Returns element list as c-object
    /// @return elementsC_
    ElementsPtr_Type getElementsC();

    /// @brief Getter for surface elements. Probably set in mesh partitioner. They are generally the dim-1 surface elements
    /// @return surfaceElements_
    ElementsPtr_Type getSurfaceElements();
    
    /// @brief 
    /// @return dim_ 
    int getDimension();
    
    /// @brief Global number of elements
    /// @return numElementsGlob_
    GO getNumElementsGlobal();
    
    /// @brief Local number of elements
    /// @return 
    LO getNumElements();
    
    /// @brief Get local (LO) number of points either in unique or repeated version
    /// @param type LO (local ordinal)
    /// @return numer of points
    LO getNumPoints(std::string type="Unique");
    
    /// @brief 
    /// @return 
    int getOrderElement();

    /// @brief Communicator object
    /// @return 
    CommConstPtrConst_Type getComm(){return comm_;};
    
    /// @brief This is done in meshStructured. Maybe we should move it here or delete this.
    /// @param flags 
    /// @return 
    int setStructuredMeshFlags(int flags){return 0;};
    
    /// @brief Something for TPM. Can be deprecated soon.
    /// @param type 
    void setElementFlags(std::string type="");
    
    /// @brief Setting current coordinates as reference configuration. Should only be called once :D
    void setReferenceConfiguration();
    
    /// @brief Moving mesh according to displacement based on reference configuration.
    /// @param displacementUnique displacement in unqiue dist.
    /// @param displacementRepeated displacement in repeated dist.
    void moveMesh( MultiVectorPtr_Type displacementUnique, MultiVectorPtr_Type displacementRepeated );
    
    // Creates an AABBTree from own vertice- and elementlist.
    void create_AABBTree();
    
    vec_int_ptr_Type  findElemsForPoints(vec2D_dbl_ptr_Type query_points);
    
    vec_dbl_Type getBaryCoords(vec_dbl_Type point, int element);
    
    bool isPointInElem(vec_dbl_Type point, int element);

    /// @brief 
    /// @return 
    tuple_intint_Type getRankRange() const {return rankRange_;};
    
    /// @brief Deleting surface elements, called after reading input mesh. 
    void deleteSurfaceElements(){ surfaceElements_.reset(); };

    /// @brief Correcting the normal direction of all surface normals set as subelements of volume elements to be outward normals
    void correctNormalDirections();

    /// @brief Correct the element orientation of all elements to have positive volume / det when doint transformation
    void correctElementOrientation();    

    /// @brief Generating elements list corresponding to element list 
    void buildEdges(ElementsPtr_Type elements);

    /// Rebuilding the information of elementsOfEdgeLocal and elementsOfEdgeGlobal of edges
    void updateElementsOfEdgesLocalAndGlobal(int maxRank);

    /// @brief This just gives us the element locations of nodes that determine an edge
    void setLocalEdgeIndices(vec2D_int_Type &localEdgeIndices);

    /// @brief Compute the volume of a tetrahedra
    vec_dbl_Type determineVolTet(ElementsPtr_Type elements, vec2D_dbl_ptr_Type points);

    /// @brief Compute the circum diameter of tetrahedra
    void calcDiamTetraeder();

    /*!
            \brief Returns elements as a vector type contrary to the C-object list.
    */
    vec2D_int_ptr_Type getElements();

    /*!
            \brief Building Edge Map
    */
    void buildEdgeMap();
    /* ###################################################################### */
    
    int                     dim_; // Dimension
    long long               numElementsGlob_; // global number of elements

    std::string 			FEType_; // Finite Element Type
    MapPtr_Type             mapUnique_; // Unique node map
    MapPtr_Type 			mapRepeated_; // Repeated node map
    vec2D_dbl_ptr_Type		pointsRep_; // repeated points
    vec2D_dbl_ptr_Type 		pointsUni_; // unique points
    vec_int_ptr_Type 		bcFlagRep_; // repeated flags
    vec_int_ptr_Type		bcFlagUni_; // unique flags

    ElementsPtr_Type        surfaceElements_; // Surface elements as elements type

    ElementsPtr_Type        elementsC_;
    MapPtr_Type				elementMap_;
    MapPtr_Type				edgeMap_;

    // April 2024: EdgeElements and SurfaceTriangleElements were move from unstructuredMeshes to Meshes, as they can be a feature for structured meshes now.
    EdgeElementsPtr_Type edgeElements_;    
	SurfaceElementsPtr_Type surfaceTriangleElements_; 

    CommConstPtrConst_Type  comm_;
    
	vec2D_int_ptr_Type  elementsVec_; // Elements as a vector not a c-object
    
    vec2D_dbl_ptr_Type		pointsRepRef_; // Repeated Referenzkonfiguration
    vec2D_dbl_ptr_Type		pointsUniRef_; // Unique Referenzkonfiguration
    
    //vec_int_ptr_Type 		elementFlag_;

    MapPtr_Type mapUniqueP2Map_;
    MapPtr_Type mapRepeatedP2Map_;
    
    ParameterListPtr_Type pList_;

    int volumeID_;

    int elementOrder_;
    int surfaceElementOrder_;
    int edgesElementOrder_;

    AABBTreePtr_Type AABBTree_;
    
    tuple_intint_Type rankRange_;
    
    /* ###################################################################### */
private:

    void flipSurface(ElementsPtr_Type subEl, int surfaceNumber);
    void flipElement(ElementsPtr_Type elements, int elementNumber);
};
}

#endif
