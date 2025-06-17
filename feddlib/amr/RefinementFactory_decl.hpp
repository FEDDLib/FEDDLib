#ifndef RefinementFactory_decl_hpp
#define RefinementFactory_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "feddlib/core/Mesh/Mesh.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/Mesh/MeshInterface.hpp"
#include "feddlib/core/Mesh/MeshFileReader.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/FE/TriangleElements.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include <Tpetra_CrsMatrix.hpp>
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/amr/RefinementFactory.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include  <boost/function.hpp>

/*!
 Declaration of RefinementFactory
 
 @brief  RefinementFactory
 @author Lea Saßmannshausen
 @version 1.0
 @copyright CH
 */

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class RefinementFactory : public MeshUnstructured<SC,LO,GO,NO> {
    
public:
   typedef Mesh<SC,LO,GO,NO> Mesh_Type;
   typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshUnstrPtr_Type;

    typedef std::vector<MeshUnstrPtr_Type> MeshUnstrPtrArray_Type;

    /*typedef MeshUnstructuredRefinement<SC,LO,GO,NO> MeshUnstrRef_Type;
    typedef Teuchos::RCP<MeshUnstrRef_Type> MeshUnstrRefPtr_Type;
    typedef std::vector<MeshUnstrRefPtr_Type> MeshUnstrRefPtrArray_Type; // Array of meshUnstr for meshRefinement*/

	typedef typename Mesh_Type::CommPtr_Type CommPtr_Type;
    typedef typename Mesh_Type::CommConstPtr_Type CommConstPtr_Type;

    
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type>  ElementsPtr_Type;
    typedef SurfaceElements SurfaceElements_Type;
    typedef Teuchos::RCP<SurfaceElements_Type> SurfaceElementsPtr_Type;
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

	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    RefinementFactory();
    
    RefinementFactory( CommConstPtr_Type comm, int volumeID=10 );

    RefinementFactory( CommConstPtr_Type comm, int volumeID,  ParameterListPtr_Type parameterListAll);
    
    ~RefinementFactory();
    

    void refineMesh( MeshUnstrPtr_Type meshP1, int iteration, MeshUnstrPtr_Type outputMesh, std::string refinementMode); // MeshRefinement
    	
	void refineRegular(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i, SurfaceElementsPtr_Type surfaceTriangleElements); // aka red refinement
	
	void refineGreen( EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i); // green refinement

	void refineBlue(  EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i); // blue refinement

	void refineRed(  EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements,  int i); // red refinement

	void refineType1( EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements);
	void refineType2( EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements);
	void refineType3( EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements);
	void refineType4(  EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements);


	void addMidpoint(EdgeElementsPtr_Type edgeElements, int i); // Adding midpoint on edge

	int determineLongestEdge( EdgeElementsPtr_Type edgeElements, vec_int_Type edgeVec, vec2D_dbl_ptr_Type points); // Determines longest edge in triangle

	void buildEdgeMap(MapConstPtr_Type mapGlobalProc,MapConstPtr_Type mapProc);
	void buildNodeMap(EdgeElementsPtr_Type edgeElements, MapConstPtr_Type mapGlobalProc, MapConstPtr_Type mapProc, int newPoints, int newPointsRepeated);

	
	void updateElementsOfEdgesLocalAndGlobal(int maxRank, MapConstPtr_Type edgeMap);

	void updateElementsOfSurfaceLocalAndGlobal(EdgeElementsPtr_Type edgeElements);

	vec_bool_Type checkInterfaceSurface( EdgeElementsPtr_Type edgeElements,vec_int_Type originFlag, vec_int_Type edgeNumbers, int indexElement);

	void refinementRestrictions(MeshUnstrPtr_Type meshP1, ElementsPtr_Type elements ,EdgeElementsPtr_Type edgeElements,SurfaceElementsPtr_Type surfaceTriangleElements, int& newPoints, int& newPointsCommon, vec_GO_Type& globalInterfaceIDsTagged, MapConstPtr_Type mapInterfaceEdges, int& newElements); // check if Element that is tagged to be refined green has previously been refined green

	void refineMeshRegIreg(ElementsPtr_Type elements, EdgeElementsPtr_Type edgeElements, int& newElements, MapConstPtr_Type edgeMap, SurfaceElementsPtr_Type surfaceTriangleElements);

	void buildSurfaceTriangleElements(ElementsPtr_Type elements, EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceTriangleElements, MapConstPtr_Type edgeMap, MapConstPtr_Type elementMap );

	void setErrorEstimate(vec_dbl_Type errorElements) { errorEstimation_ = errorElements;};	

	vec_dbl_Type getErrorEstimate() { return errorEstimation_ ; };	

	void bisectEdges(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements, std::string mode = "default");

	void bisectElement3(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElementp);

	std::string refinementRestriction_ = "none";

	bool writeRefinementTime_ = "true";

	int refinement3DDiagonal_ = 0; // 0 beeing the shortest interior Diagonal, 1 the second shortest and 2 the longest interior Diagonal 

	int currentIter_ = 0;

	std::string refinementMode_ = "Regular";

protected: 
	vec_GO_Type globalInterfaceIDs_;
	vec_dbl_Type errorEstimation_; // error estimated according to A-posteriori Error Estimator. Sorted according to loca Element IDs

	vec_dbl_Type areaTriangles_;
	vec_dbl_Type volTetraeders_;
	vec_dbl_Type h_T_diam_E_; // Diameter of 2D Triangular Elements
	vec_dbl_Type h_T_min_;	// Diameter in the 3D Tetraeder Sense
	
	MapConstPtr_Type surfaceTriangleMap_;



private:
		

    
 
    
};
}
#endif
