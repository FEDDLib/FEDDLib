#ifndef AdaptiveMeshRefinement_def_hpp
#define AdaptiveMeshRefinement_def_hpp

#ifndef MESH_TIMER_START
#define MESH_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Mesh Refinement") + std::string(S))));
#endif

#ifndef MESH_TIMER_STOP
#define MESH_TIMER_STOP(A) A.reset();
#endif

#include "AdaptiveMeshRefinement_decl.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector_def.hpp"
#include <chrono> 
/*!
 Definition of AdaptiveMeshRefinement
 
 @brief  AdaptiveMeshRefinement
 @version 1.0

 */

using namespace std;
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::outArg;

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement():
MeshUnstructured<SC,LO,GO,NO>()
{

  
}

template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement(CommConstPtr_Type comm, int volumeID):
MeshUnstructured<SC,LO,GO,NO>(comm,volumeID)
{

}
/// 
/// Initializing problem with the kind of problem (e.g. Laplace, Stokes) for determining the correct error estimation and the dimension
///
template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement(string problemType, int dim)
{
	this->dim_ = dim;
	this->problemType_ = problemType;

	this->exportWithParaview_ = false;
	
}
/// 
/// Initializing problem with the kind of problem, dimension and refinement spectific parameters
///
template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement(string problemType, int dim, vec_dbl_Type parasDbl, vec_string_type parasString, , ExporterPtr_Type exportSol, ExporterPtr_Type exportError )
{
	this->dim_ = dim;
	this->problemType_ = problemType;

	// set refinement specific parameters with parasDbl and parasString that are usually set with default values


	
}

template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::~AdaptiveMeshRefinement(){

}

///
/// Global Algorithm of Mesh Refinement
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::globalAlgorithmglobalAlgorithm(DomainPtr_Type domainP1, DomainPtr_Type domainP2 MultiVectorPtr_Type solutionP1, MultiVectorPtr_Type solutionP2 ){

	// The P1 Mesh is always used for refinement while the P1 or P2 Mesh is used for error Estimation depending on Discretisation

	// Reading Mesh from domainP1 as we always refine the P1 Mesh
	MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp( new MeshUnstrRef_Type( comm_,(Teuchos::rcp_dynamic_cast<MeshUnstrRef_Type>( domainsP1[iter]->mesh_)->volumeID_) ) );


	/*
	// Set Domain Values
	n_ = domainsP1[iter]->n_;
	m_ = domainsP1[iter]->m_;
    dim_ = domainsP1[iter]->dim_;
	FEType_ = domainsP1[iter]->FEType_;

    numProcsCoarseSolve_ = 0;
    flagsOption_ = -1;

    meshType_ = "unstructured";


	meshUnstructuredRefined->refinementRestriction_ = restriction;
	meshUnstructuredRefined->meshQualityPrint_ = writeMeshQuality;
	meshUnstructuredRefined->timeTablePrint_ =writeTime;
	meshUnstructuredRefined->refinement3DDiagonal_ = diagonal;

	MeshUnstrRefPtr_Type meshUnstructuredRefinedP2 = Teuchos::rcp_dynamic_cast<MeshUnstrRef_Type>( domainP2->mesh_ , true);


	// If coarsen Mesh is false, so consequently we refine the Mesh we go about as folows:
	if(coarsen == false){
		// We pass on all previously refined meshes 
		// Useful for later coarsening of elements or determining restrictions of refinement strategies 
		// i.e. no previously green refined element is refined green again
		MeshUnstrRefPtrArray_Type meshUnstructuredP1(iter+1);

		for(int i=0; i<iter+1; i++)
			meshUnstructuredP1[i] = Teuchos::rcp_dynamic_cast<MeshUnstrRef_Type>( domainsP1[i]->mesh_ , true);
			

		// Estimating the error with the Discretizations Mesh.
		vec_dbl_Type errorElements = meshUnstructuredRefinedP2->estimateError(valuesSolution, theta, strategy, rhsFunc);

		meshUnstructuredRefined->setErrorEstimate(errorElements); 
		meshUnstructuredP1[iter]->markElements(errorElements,strategy);

   		meshUnstructuredRefined->refineMesh(meshUnstructuredP1,iter);
	}
	else {
		int m=1;
		int n=m+1;
		int k = iter;
		int iterC;
		MeshUnstrRefPtrArray_Type meshUnstructuredP1(iter+1-m);

		for(int i=0; i<iter+1-m; i++)
			meshUnstructuredP1[i] = Teuchos::rcp_dynamic_cast<MeshUnstrRef_Type>( domainsP1[i]->mesh_ , true);
			
		// We extract the error estimation of the mesh iter

		vec_dbl_Type errorElements; // = meshUnstructuredRefined_k->mesh_->getErrorEstimate();
		MeshUnstrRefPtr_Type meshUnstructuredRefined_k ;
		MeshUnstrRefPtr_Type meshUnstructuredRefined_k_1;
		MeshUnstrRefPtr_Type meshUnstructuredRefined_k_m_1;
		for(int i=0; i<m-1 ; i++){
			meshUnstructuredRefined_k = Teuchos::rcp_dynamic_cast<MeshUnstrRef_Type>( domainsP1[iter-i]->mesh_ , true); 
			meshUnstructuredRefined_k_1 = Teuchos::rcp_dynamic_cast<MeshUnstrRef_Type>( domainsP1[iter-1-i]->mesh_ , true); 

			errorElements = meshUnstructuredRefined_k_1->determineCoarseningError(meshUnstructuredRefined_k,"backwards");

			meshUnstructuredRefined_k_1->setErrorEstimate(errorElements);

		}

		// Setting error of: Mesh_(k-m+1) with the previous error ->downscaling errors
		if(m>1)
			meshUnstructuredRefined_k_m_1 = meshUnstructuredRefined_k_1; 
		else
			meshUnstructuredRefined_k_m_1 =Teuchos::rcp_dynamic_cast<MeshUnstrRef_Type>( domainsP1[iter]->mesh_ , true); 


		for(int i=0; i< n; i++){
			iterC = k-m+i;
			if(i==0){
				errorElements = meshUnstructuredP1[k-m]->determineCoarseningError(meshUnstructuredRefined_k_m_1,"backwards");				
			}
			else
				errorElements = meshUnstructuredP1[k-m+i]->determineCoarseningError(meshUnstructuredP1[k-m+i-1],"forwards");

			meshUnstructuredRefined->setErrorEstimate(errorElements); 
			meshUnstructuredRefined->refineMesh(meshUnstructuredP1, iterC);

			meshUnstructuredP1.push_back(meshUnstructuredRefined);

			meshUnstructuredRefined.reset( new MeshUnstrRef_Type( comm_,  meshUnstructuredP1[0]->volumeID_) );

			meshUnstructuredRefined->refinementRestriction_ = restriction;
			meshUnstructuredRefined->meshQualityPrint_ = writeMeshQuality;
			meshUnstructuredRefined->timeTablePrint_ =writeTime;
			meshUnstructuredRefined->refinement3DDiagonal_ = diagonal;

		}
		meshUnstructuredRefined = meshUnstructuredP1[iter+1];
	}

    mesh_ = meshUnstructuredRefined;
	*/

}




}

template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportSolution(){



}

template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportError(){



}
}















