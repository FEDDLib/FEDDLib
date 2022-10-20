#ifndef MESHUNSTRUCTUREDREFINEMENT_decl_hpp
#define MESHUNSTRUCTUREDREFINEMENT_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "Mesh.hpp"
#include "MeshUnstructured.hpp"
#include "MeshInterface.hpp"
#include "MeshFileReader.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"

/*!
 Declaration of MeshUnstructuredRefinement
 
 @brief  MeshUnstructuredRefinement
 @author Lea Saßmannshausen
 @version 1.0
 @copyright CH
 */

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshUnstructuredRefinement : public MeshUnstructured<SC,LO,GO,NO> {
    
public:
	typedef Mesh<SC,LO,GO,NO> Mesh_Type;
	typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshUnstrPtr_Type;

	typedef std::vector<MeshUnstrPtr_Type> MeshUnstrPtrArray_Type;

    typedef MeshUnstructuredRefinement<SC,LO,GO,NO> MeshUnstrRef_Type;
    typedef Teuchos::RCP<MeshUnstrRef_Type> MeshUnstrRefPtr_Type;
    typedef std::vector<MeshUnstrRefPtr_Type> MeshUnstrRefPtrArray_Type; // Array of meshUnstr for meshRefinement
    
    typedef typename Mesh_Type::CommPtr_Type CommPtr_Type;
    typedef typename Mesh_Type::CommConstPtr_Type CommConstPtr_Type;
    typedef typename Mesh_Type::Elements_Type Elements_Type;
    typedef typename Mesh_Type::ElementsPtr_Type ElementsPtr_Type;
    
    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;
    
    typedef MeshInterface<SC,LO,GO,NO> MeshInterface_Type;
    typedef Teuchos::RCP<MeshInterface_Type> MeshInterfacePtr_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef MultiVector<LO,LO,GO,NO> MultiVectorLO_Type;
	typedef Teuchos::RCP<MultiVectorLO_Type> MultiVectorLOPtr_Type;
	typedef MultiVector<GO,LO,GO,NO> MultiVectorGO_Type;
	typedef Teuchos::RCP<MultiVectorGO_Type> MultiVectorGOPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;
	typedef Teuchos::OrdinalTraits<LO> OTLO;

    MeshUnstructuredRefinement();
    
    MeshUnstructuredRefinement( CommConstPtr_Type comm, int volumeID=10 );

    MeshUnstructuredRefinement( CommConstPtr_Type comm, int volumeID, MeshUnstrPtr_Type meshP1);
    
    ~MeshUnstructuredRefinement();
    
    void refineMesh(  MeshUnstrRefPtrArray_Type meshP1, int iteration, bool checkRestrictions, string restriction); // MeshRefinement
		
	vec_dbl_Type errorEstimation(MultiVectorPtrConst_Type valuesSolution, double theta, string strategy);

protected: 



private:   	
	void refineRegular( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i); // aka red refinement
	
	void refineGreen( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i); // green refinement

	void refineBlue( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i); // blue refinement

	void refineRed( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i); // red refinement

	void addMidpoint(MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, int i); // Adding midpoint on edge

	int determineLongestEdge( EdgeElementsPtr_Type edgeElements, vec_int_Type edgeVec, vec2D_dbl_ptr_Type points); // Determines longest edge in triangle

	void addRepeatedIDsEdges(vec_GO_Type& vecGlobalIDsEdges, vec_GO_Type edgesInterfaceUntagged, vec_GO_Type edgesInterfaceTagged ,int procOffsetEdgesUniqueSum,vec2D_int_Type interfaceEdgesLocalId);

	void updateElementsOfEdgesLocalAndGlobal(vec2D_int_Type interfaceEdgesLocalId, int maxRank, vec_GO_Type edgesInterface, MapConstPtr_Type edgeMap);

	void checkGreenTags(MeshUnstrRefPtr_Type meshP1, ElementsPtr_Type elements ,EdgeElementsPtr_Type edgeElements, int iteration, int& newPoints, int& newPointsCommon, vec_GO_Type& globalTaggedInterfaceIDs, MapConstPtr_Type mapInterfaceEdges, string restriction ); // check if Element that is tagged to be refined green has previously been refined green
	vec2D_int_Type determineInterfaceEdgesNewLocalIds(EdgeElementsPtr_Type edgeElements,vec_GO_Type edgesInterface, MapConstPtr_Type edgeMap);

	vec2D_int_Type determineNewInterfaceEdgesLocalIds(vec_GO_Type& globalNodeIdsNewInterfaceEdges) ;

	void distributeTaggedAndUntaggedEdges(MapConstPtr_Type mapProc, MapConstPtr_Type mapGlobalProc, MapConstPtr_Type mapInterfaceEdgesTaggedUnique, MapConstPtr_Type mapInterfaceEdgesUntaggedUnique, vec_GO_Type& edgesInterfaceUntagged, vec_GO_Type& edgesInterfaceTagged) ;

    
};
}
#endif
