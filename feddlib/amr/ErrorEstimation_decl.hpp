#ifndef ErrorEstimation_decl_hpp
#define ErrorEstimation_decl_hpp

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
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/amr/ErrorEstimation.hpp"
#include "feddlib/amr/RefinementFactory.hpp"
#include "feddlib/amr/AdaptiveMeshRefinement.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include  <boost/function.hpp>

/*!
 Declaration of ErrorEstimation
 
 @brief  ErrorEstimation
 @author Lea Saßmannshausen
 @version 1.0
 @copyright CH
 */

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class ErrorEstimation {
    
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

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef MultiVector<LO,LO,GO,NO> MultiVectorLO_Type;
	typedef Teuchos::RCP<MultiVectorLO_Type> MultiVectorLOPtr_Type;
    typedef MultiVector<GO,LO,GO,NO> MultiVectorGO_Type;
    typedef Teuchos::RCP<MultiVectorGO_Type> MultiVectorGOPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef Teuchos::OrdinalTraits<LO> OTLO;

	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

 typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
  


    ErrorEstimation();
    
    ErrorEstimation(int dim, string problemType);
    
    ~ErrorEstimation();
    
	// Error Estimation Functions that will be reallocted soon
	MultiVectorPtr_Type estimateError(MeshUnstrPtr_Type inputMeshP12, MeshUnstrPtr_Type inputMeshP1, BlockMultiVectorConstPtr_Type valuesSolution, RhsFunc_Type rhsFunc, string FEType);

	void identifyProblem(BlockMultiVectorConstPtr_Type valuesSolution);

	void makeRepeatedSolution(BlockMultiVectorConstPtr_Type valuesSolution);

	vec3D_dbl_Type calcNPhi(string phiDerivative, int dofsSol, string FEType);

	vec_dbl_Type calculateJump();

	vec2D_dbl_Type gradPhi(int dim,int intFE,vec_dbl_Type &p);
	vec_dbl_Type phi(int dim,int intFE,vec_dbl_Type &p);
	vec_dbl_Type divPhi(int dim,int intFE,vec_dbl_Type &p);

	MultiVectorPtr_Type determineCoarseningError(MeshUnstrPtr_Type mesh_k, MeshUnstrPtr_Type mesh_k_m, MultiVectorPtr_Type errorElementMv_k,  string distribution, string markingStrategy, double theta); // Mesh Coarsening

	double determineResElement(FiniteElement element, RhsFunc_Type rhsFunc);

	double determineDivU(FiniteElement element);

	vec2D_dbl_Type getQuadValues(int dim, string FEType, string Type, vec_dbl_Type &QuadW, FiniteElement surface);

	void markElements(MultiVectorPtr_Type errorElementMv, double theta, string strategy,  MeshUnstrPtr_Type meshUnstr);
	
	vec_dbl_Type determineVolTet(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points);

	vec_dbl_Type calcDiamTriangles(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points, vec_dbl_Type& areaTriangles, vec_dbl_Type& rho_T,vec_dbl_Type& C_T);
	vec_dbl_Type calcDiamTriangles3D(SurfaceElementsPtr_Type surfaceTriangleElements,vec2D_dbl_ptr_Type points,vec_dbl_Type& areaTriangles, vec_dbl_Type& rho_T,vec_dbl_Type& C_T);

	vec_dbl_Type calcDiamTetraeder(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points, vec_dbl_Type volTet);

	vec_dbl_Type calcRhoTetraeder(ElementsPtr_Type elements,SurfaceElementsPtr_Type surfaceTriangleElements, vec_dbl_Type volTet, vec_dbl_Type areaTriangles);

	vec_dbl_Type determineAreaTriangles(ElementsPtr_Type elements,EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceElements, vec2D_dbl_ptr_Type points);

	void buildTriangleMap();

	void updateElementsOfSurfaceLocalAndGlobal(EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceTriangleElements);

	void setErrorEstimate(MultiVectorPtr_Type errorElements) { errorEstimation_ = errorElements;};	

	MultiVectorPtr_Type getErrorEstimate() { return errorEstimation_ ; };	

	void tagArea(MeshUnstrPtr_Type meshUnstr,vec2D_dbl_Type area);

	string refinementRestriction_ = "none";
	string markingStrategy_ = "Maximum";
	double theta_ = 0.5;

	bool meshQualityPrint_ = "false";
	bool timeTablePrint_ = "false";
	int refinement3DDiagonal_ = 0; // 0 beeing the shortest interior Diagonal, 1 the second shortest and 2 the longest interior Diagonal 

	int dim_;
	string problemType_;

protected: 
	vec_GO_Type globalInterfaceIDs_;
	MultiVectorPtr_Type errorEstimation_; // error estimated according to A-posteriori Error Estimator. Sorted according to loca Element IDs

	vec_dbl_Type areaTriangles_;
	vec_dbl_Type volTetraeders_;
	vec_dbl_Type h_T_diam_E_; // Diameter of 2D Triangular Elements
	vec_dbl_Type h_T_min_;	// Diameter in the 3D Tetraeder Sense
	
	MapConstPtr_Type surfaceTriangleMap_;
	SurfaceElementsPtr_Type surfaceElements_;

	int dofs_;
	int dofsP_;

	bool calculatePressure_ = false;

	BlockMultiVectorConstPtr_Type valuesSolutionRepVel_;
	BlockMultiVectorConstPtr_Type valuesSolutionRepPre_;


private:
	MeshUnstrPtr_Type inputMesh_;
	MeshUnstrPtr_Type inputMeshP1_;
	string FEType1_;
	string FEType2_;

    
 
    
};
}
#endif
