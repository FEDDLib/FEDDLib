#ifndef AdaptiveMeshRefinement_decl_hpp
#define AdaptiveMeshRefinement_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "feddlib/core/Mesh/Mesh.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/Mesh/MeshInterface.hpp"
#include "feddlib/core/Mesh/MeshFileReader.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/FE/TriangleElements.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include <Tpetra_CrsMatrix.hpp>
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include  <boost/function.hpp>
#include "feddlib/problems/abstract/Problem.hpp"
#include "feddlib/amr/ExporterParaViewAMR.hpp"
#include "feddlib/amr/ErrorEstimation.hpp"
#include "feddlib/amr/RefinementFactory.hpp"

/*!
 Declaration of Adaptive Mesh Refinement
 
 @brief  Adaptive Mesh Refinement
 @author Lea Saßmannshausen

 */

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AdaptiveMeshRefinement {
    
public:
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef MeshUnstructured<SC,LO,GO,NO>  MeshUnstr_Type;
    typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshUnstrPtr_Type;

    typedef std::vector<MeshUnstrPtr_Type> MeshUnstrPtrArray_Type;

 	typedef typename Mesh_Type::CommPtr_Type CommPtr_Type;
    typedef typename Mesh_Type::CommConstPtr_Type CommConstPtr_Type;
    
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type>  ElementsPtr_Type;
    typedef SurfaceElements SurfaceElements_Type;
    typedef Teuchos::RCP<SurfaceElements_Type> SurfaceElementsPtr_Type;
    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;
     
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

 	typedef ExporterParaViewAMR<SC,LO,GO,NO> Exporter_Type;
    typedef Teuchos::RCP<Exporter_Type> ExporterPtr_Type;
    typedef Teuchos::RCP<ExporterTxt> ExporterTxtPtr_Type;

	typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;

	typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef Teuchos::RCP<Domain_Type> DomainPtr_Type;
    typedef std::vector<DomainPtr_Type> DomainPtrArray_Type;

    typedef std::vector<MultiVectorPtr_Type> MultiVectorPtrArray_Type;

	typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
  
    AdaptiveMeshRefinement();

    AdaptiveMeshRefinement(ParameterListPtr_Type parameterListAll);
    
    AdaptiveMeshRefinement(std::string problemType, ParameterListPtr_Type parameterListAll);

	AdaptiveMeshRefinement(std::string problemType, ParameterListPtr_Type parameterListAll , Func_Type exactSolFunc );
	AdaptiveMeshRefinement(std::string problemType, ParameterListPtr_Type parameterListAll , Func_Type exactSolFuncU,Func_Type exactSolFuncP );
    
    ~AdaptiveMeshRefinement();

	DomainPtr_Type globalAlgorithm(DomainPtr_Type domainP1, DomainPtr_Type domainP12, BlockMultiVectorConstPtr_Type solution,ProblemPtr_Type problem, RhsFunc_Type rhsFunc );

	DomainPtr_Type refineArea(DomainPtr_Type domainP1, vec2D_dbl_Type area, int level);
    
    DomainPtr_Type refineUniform(DomainPtr_Type domainP1, int level);
    
	DomainPtr_Type refineFlag(DomainPtr_Type domainP1, int level ,int flag);
	
	MultiVectorConstPtr_Type calcExactSolution();
	MultiVectorConstPtr_Type calcExactSolutionP();
	//void determineCoarsening();

	void identifyProblem(BlockMultiVectorConstPtr_Type valuesSolution);

	void calcErrorNorms(MultiVectorConstPtr_Type exactSolution, MultiVectorConstPtr_Type solutionP12,MultiVectorConstPtr_Type exactSolutionP);

	void initExporter( ParameterListPtr_Type parameterListAll);

	void exportSolution(MeshUnstrPtr_Type mesh, MultiVectorConstPtr_Type exportSolutionMv, MultiVectorConstPtr_Type errorValues, MultiVectorConstPtr_Type exactSolutionMv,MultiVectorConstPtr_Type exportSolutionPMv, MultiVectorConstPtr_Type exactSolutionPMv);

	void exportError(MeshUnstrPtr_Type mesh, MultiVectorConstPtr_Type errorElConst, MultiVectorConstPtr_Type errorElConstH1 , MultiVectorConstPtr_Type difH1Eta ,MultiVectorConstPtr_Type vecDecompositionConst );

	void writeRefinementInfo();
	
	void buildSurfaceTriangleElements(ElementsPtr_Type elements, EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceTriangleElements );

	vec_bool_Type checkInterfaceSurface( EdgeElementsPtr_Type edgeElements,vec_int_Type originFlag, vec_int_Type edgeNumbers, int indexElement);

protected: 
	


private:

	RhsFunc_Type rhsFunc_;
	Func_Type exactSolFunc_;

	Func_Type exactSolPFunc_;
	
	MeshUnstrPtr_Type inputMeshP1_;
	MeshUnstrPtr_Type inputMeshP12_;
	MeshUnstrPtr_Type outputMesh_;

	MultiVectorPtrArray_Type errorEstimationMv_;

	MultiVectorPtr_Type errorElementsMv_;
	MultiVectorPtr_Type errorH1ElementsMv_;
	MultiVectorPtr_Type difH1EtaElementsMv_;

	MultiVectorConstPtr_Type errorNodesMv_;
	MultiVectorConstPtr_Type errorNodesPMv_;

	BlockMultiVectorConstPtr_Type solution_;

   	CommConstPtr_Type comm_;

	bool exportWithParaview_ = true;
	bool initExporter_=false;

	ExporterPtr_Type exporterSol_;
	ExporterPtr_Type exporterSolP_;
	ExporterPtr_Type exporterError_;

	DomainPtrArray_Type domainsP1_;
	DomainPtrArray_Type domainsP12_;

	DomainPtr_Type domainP1_;
	DomainPtr_Type domainP12_;

	ProblemPtr_Type problem_;
	
	std::string refinementRestriction_ = "keepRegularity";
	std::string markingStrategy_ = "Maximum";

	double theta_ = 0.5;
	double tol_= 0.001 ;

	bool meshQualityPrint_ = "false";
	bool timeTablePrint_ = "false";
	int refinement3DDiagonal_ = 0; // 0 beeing the shortest interior Diagonal, 1 the second shortest and 2 the longest interior Diagonal 

	std::string problemType_;
	int dim_;

	int currentIter_;
	int maxIter_ = 5;
	int maxRank_;

	std::string FEType1_;
	std::string FEType2_;

	vec_dbl_Type maxErrorEl;
	vec_dbl_Type maxErrorKn;
	vec_int_Type numElements;
	vec_int_Type numElementsProc;
	vec_dbl_Type relError;
	vec_dbl_Type eRelError;
	vec_dbl_Type errorH1;
	vec_dbl_Type errorL2;
	vec_dbl_Type errorL2P;
	vec_int_Type numNodes;

	bool writeRefinementTime_ = true ;
	bool writeRefinementInfo_ = true ;
	bool writeMeshQuality_ = true ;

	bool hasProblemType_=true;

	ParameterListPtr_Type parameterListAll_ ;

	int dofs_;
	int dofsP_;

	bool exactSolInput_ = false ;
	bool exactSolPInput_ = false ;

	bool calculatePressure_=false;

	int restrictionLayer_=2;
	
	int coarseningCycle_=0 ;
	int coarseningM_ =  1;
	int coarseningN_  = 1;
		
	std::string refinementMode_ = "Regular";

  
    
};
}
#endif
