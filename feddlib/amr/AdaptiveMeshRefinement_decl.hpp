#ifndef ADAPTIVEMESHREFINEMENT_decl_hpp
#define ADAPTIVEMESHREFINEMENT_decl_hpp

#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "feddlic/core/Mesh/Mesh.hpp"
#include "feddlic/core/Mesh/MeshUnstructured.hpp"
#include "feddlic/core/Mesh/MeshInterface.hpp"
#include "feddlic/core/Mesh/MeshFileReader.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include "feddlib/core/FE/TriangleElements.hpp"
#include "feddlib/core/FE/EdgeElements.hpp"
#include <Tpetra_CrsMatrix.hpp>
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructuredRefinement.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include  <boost/function.hpp>

/*!
 Declaration of Adaptive Mesh Refinement
 
 @brief  MeshUnstructuredRefinement
 @author Lea Saßmannshausen
 @version 1.0
 @copyright CH
 */

namespace FEDD {
    
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AdaptiveMeshRefinement {
    
public:
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<MeshUnstructured<SC,LO,GO,NO> > MeshUnstrPtr_Type;

    typedef std::vector<MeshUnstrPtr_Type> MeshUnstrPtrArray_Type;

    typedef MeshUnstructuredRefinement<SC,LO,GO,NO> MeshUnstrRef_Type;
    typedef Teuchos::RCP<MeshUnstrRef_Type> MeshUnstrRefPtr_Type;
    typedef std::vector<MeshUnstrRefPtr_Type> MeshUnstrRefPtrArray_Type; // Array of meshUnstr for meshRefinement

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

 	typedef ExporterParaView<SC,LO,GO,NO> Exporter_Type;
    typedef Teuchos::RCP<Exporter_Type> ExporterPtr_Type;
    typedef Teuchos::RCP<ExporterTxt> ExporterTxtPtr_Type;


    AdaptiveMeshRefinement();
    
    AdaptiveMeshRefinement( CommConstPtr_Type comm, int volumeID=10 );

	AdaptiveMeshRefinement(string problemType, int dim);

	AdaptiveMeshRefinement(string problemType, int dim, vec_dbl_Type parasDbl, vec_string_type parasString, , ExporterPtr_Type exportSol, ExporterPtr_Type exportError )
    
    ~AdaptiveMeshRefinement();

	void globalAlgorithm(DomainPtr_Type domainP1, DomainPtr_Type domainP2, vec_dbl_Type parasDbl, vec_string_type parasString, ExporterPtr_Type exportSol, ExporterPtr_Type exportError, MultiVectorPtr_Type solutionP1, MultiVectorPtr_Type solutionP2 );
    
	void exportSolution();
	void exportError();


protected: 
	


private:

	bool exportWithParaview_ = true;

	ExporterPtr_Type exportSol_;
	ExporterPtr_Type exportError_;

	DomainPtr_Type domainsP1_;
	DomainPtr_Type domainsP12_;

	string refinementRestrictions_ = "none";
	string markingStrategy_ = "Maximum";

	double theta_ = 0.5;

	bool meshQualityPrint_ = "false";
	bool timeTablePrint_ = "false";
	int refinement3DDiagonal_ = 0; // 0 beeing the shortest interior Diagonal, 1 the second shortest and 2 the longest interior Diagonal 



	
	
		

    
 
    
};
}
#endif
