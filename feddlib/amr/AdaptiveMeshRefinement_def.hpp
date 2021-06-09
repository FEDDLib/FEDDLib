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

/*!
Initializing problem with the kind of problem (e.g. Laplace, Stokes) for determining the correct error estimation and the dimension
@param[in] problemType, dim
*/
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
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement(string problemType, int dim, vec_dbl_Type parasDbl, vec_string_type parasString, ExporterPtr_Type exportSol, ExporterPtr_Type exportError )
{
	this->dim_ = dim;
	this->problemType_ = problemType;

	/* Still need some way to set the parameters efficiently
	refinementRestriction_ = restriction;
	meshQualityPrint_ = writeMeshQuality;
	timeTablePrint_ =writeTime;
	refinement3DDiagonal_ = diagonal;
	*/

	// set refinement specific parameters with parasDbl and parasString that are usually set with default values


	
}

template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::~AdaptiveMeshRefinement(){

}

/*!
\brief Global Algorithm of Mesh Refinement 

\brief Given domains and solutions depending on problem global mesh refinement algorithm and error estimation is performed

\brief i.e. if to solve simple laplace problem, we have only one solution to put in, if to estimate error for Navier-Stokes equation we need pressure and velocity solution

@param[in] domainP1  domain with P1 discretization, always neccesary as refinement is performed on P1 Mesh
@param[in] domainP12 domain with P1 or P2 discretization if available, otherwise input domainP1
@param[in] solution1 solution of problem on P1 or P2 discretization
@param[in] solution2 solution of problem on P1 or P2 discretization if available, otherwise input solutionP1
*/
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::globalAlgorithm(DomainPtr_Type domainP1, DomainPtr_Type domainP12 MultiVectorPtr_Type solutionP1, MultiVectorPtr_Type solutionP2 ){

	int iter = domainsP1_->size();
	// We save the domains of each step
	// The P1 Mesh is always used for refinement while the P1 or P2 Mesh is used for error Estimation depending on Discretisation
	domainsP1_.push_back(domainP1);
	domainsP12_.push_back(domainP12);

	// Reading Mesh from domainP1 as we always refine the P1 Mesh, here defined as inputMesh_
	inputMeshP1_ = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP1->mesh_ , true);
	
	// With the global Algorithm we create a new P1 domain with a new mesh
	DomainPtr_Type domainRefined;

	// !!!!! Not yes existing function 
	domainRefined->initWithDomain(domainP1);
	/*
	n_ = domainsP1[iter]->n_;
	m_ = domainsP1[iter]->m_;
    dim_ = domainsP1[iter]->dim_;
	FEType_ = domainsP1[iter]->FEType_;

    numProcsCoarseSolve_ = 0;
    flagsOption_ = -1;

    meshType_ = "unstructured";*/

	inputMeshP12_ = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP12->mesh_ , true);


	// Error Estimation object
    ErrorEstimation<SC,LO,GO,NO> errorEstimator (dim_, problemType_ );

	// Refinement Factory object
	RefinementFactory<SC,LO,GO,NO> refinementFactory( domainP1->getComm(), inputMeshP1_->volumeID, inputMeshP1_);

	if(dim_ ==3){
		SurfaceElementsPtr_Type surfaceTriangleElements = inputMeshP12->getSurfaceTriangleElements(); // Surfaces
		if(surfaceTriangleElements.is_null()){
			surfaceTriangleElements.reset(new SurfaceElements()); // Surface
			cout << " Building surfaceTriangleElemenets ... " << endl;
			surfaceTriangleElements = this->buildSurfaceTriangleElements(elements,edgeElements);
			cout << " ... done " << endl;
		}
		else if(surfaceTriangleElements->numberElements() ==0){
			cout << " Building surfaceTriangleElemenets " << endl;
			surfaceTriangleElements = this->buildSurfaceTriangleElements(elements,edgeElements);
			cout << " ... done " << endl;
		}
	}

	// MESH COARSENING
	coarsen = false;

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
	

}

template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::determineCoarsening(){



}



/*!
\brief Calculating exact solution if possible
*/
template <class SC, class LO, class GO, class NO>
MultiVectorPtrConst_Type exactSolution AdaptiveMeshRefinement<SC,LO,GO,NO>::calcExactSolution(){
	
	MultiVectorPtrConst_Type exactSolution = Teuchos::rcp(new MultiVector_Type( domain->getMapUnique() ) ); 
	Teuchos::ArrayRCP<SC> exactSolA = exactSolution->getDataNonConst(0);

	vec2D_dbl_ptr_Type points = inputMesh_->getPointsUnique();

	for(int i=0; i< points.size(); i++)
		exactSol_(points->at(i),exactSolA[i]);


	return exactSolution;
}


/*!
Calculating error norms. If the exact solution is unknown we use approxmated errorNorm and error indicators
@param[in] exact solution if known
@param[in] FE solution
@param[in] error estimation
*/

template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::calcErrorNorms(){

	errorH1[j] = sqrt(laplace.calculateH1Norm(errorValues) + laplace.calculateL2Norm(errorValues));
	errorL2[j] = sqrt(laplace.calculateL2Norm(errorValues));

	double solElementH1=sqrt(laplace.calculateH1Norm(exactSolutionMv) + laplace.calculateL2Norm(exactSolutionMv));
	double solhElementH1=sqrt(laplace.calculateH1Norm(valuesSolution) + laplace.calculateL2Norm(valuesSolution));





}



///
/// ParaView exporter setup
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::initExporter(){

	exParaSol->setup( "Refinement" , domain->getMesh(),  FEType, parameterListAll );

	exParaSol->addVariable( exportSolution, "u", "Scalar", 1, domain->getMapUnique() );
	exParaSol->addVariable( errorValues, "Error |u-u_h|", "Scalar", 1, domain->getMapUnique() );
	exParaSol->addVariable( flagValues, "Flags", "Scalar", 1, domain->getMapUnique() );


	exParaEr->setup("Error_and_Dist", domainP1->getMesh(), "P0",parameterListAll );
	exParaEr->addVariable( errorElConst, "ErrorEstimate", "Scalar", 1, domainP1->getElementMap());
	exParaEr->addVariable( vecDecompositionConst, "Proc", "Scalar", 1, domainP1->getElementMap());

}

///
/// ParaView exporter export of solution on current mesh
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportSolution(){

	exParaSol->reSetup(domain->getMesh());
	exParaSol->updateVariables(exportSolution, "u");
	exParaSol->updateVariables(errorValues, "Error |u-u_h|");
	exParaSol->updateVariables(flagValues, "Flags");
			
	exPara->save( (double) j);

}


///
/// ParaView exporter export of error values and element distribution on current mesh
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportError(){

	exParaEr->reSetup(domain->getMesh());

	exParaEr->updateVariables(errorElConst,"ErrorEstimate");
	exParaEr->updateVariables(vecDecompositionConst,"Proc");

    exParaEr->save((double) j);


}

/*!
Writing refinement information
*/
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::writeRefinementInfo(){


	if(comm->getRank() == 0){	
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Summary Mesh Refinement" << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Marking Strategy:	" << strategy << endl;
			cout << " Theta:			" << theta << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Tolerance:			" << tol << endl;
			cout << " Max number of Iterations:	" <<  maxIter << endl;
			cout << " Number of Processors:		" << maxRank+1 << endl;
			cout << " Number of Refinements:		" << j << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Number of elements after Refinement.... " << endl;
			for(int i=1; i<=j ; i++)
				cout <<" "<< i << ":	" << numElements[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Number of Nodes after Refinement.... " << endl;
			for(int i=1; i<=j ; i++)
				cout <<" "<< i << ":	" << numNodes[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Errorestimation: max error in Elements according to error Estimator after Refinement.... " << endl;
			for (int i=1; i<=j ; i++)
				cout <<" "<< i << ":	" << maxErrorEl[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Maximal error in nodes after Refinement. " << endl;
			for (int i=1; i<=j ; i++)
				cout <<" "<< i << ":	" << maxErrorKn[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;

			cout << " ||u-u_h||_H1 / ||u ||_H1 	||  eta / ||u_h ||_H1	||	|| u-u_h ||_H1	||	|| u-u_h ||_L2	...." << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			for (int i=1; i<=j ; i++)
				cout <<" "<< i << ":	" << relError[i] << " 		||	" << eRelError[i] << "  	||	" << errorH1[i]<< "	||	" << errorL2[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << "Distribution of elements on .. " << endl;
			for(int l=0; l<maxRank+1 ; l++)
				cout <<" Processor "<< l << " carries " << elementProcList[l] << " Elements "<< endl; 
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
		}
		
}

// Mesh Refinement 
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::buildSurfaceTriangleElements(ElementsPtr_Type elements, EdgeElementsPtr_Type edgeElements){
    TEUCHOS_TEST_FOR_EXCEPTION( elements.is_null(), std::runtime_error, "Elements not initialized!");

	vec_LO_Type nodeInd(4);
	vec2D_int_Type newTriangles(4,vec_int_Type(0)); // new Triangles
	
	vec_GO_Type globalInterfaceIDs;
	if(globalInterfaceIDs_.size() == 0 )
		globalInterfaceIDs = edgeElements->determineInterfaceEdges(this->edgeMap_);
	//if(edgeElements->getEdgesOfElement(0) ) here we need some sort of test if the function was already used
	edgeElements->matchEdgesToElements(this->elementMap_);


	for(int i=0; i<elements->numberElements(); i++){
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(i); // indeces of edges belonging to element

		// Extract the four points of tetraeder
		vec_int_Type nodeInd(0);
		for(int i=0; i<6; i++)	{
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
		}
		sort( nodeInd.begin(), nodeInd.end() );
		nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );

		// With our sorted Nodes we construct edges as follows

		// Edge_0 = [x_0,x_1]
		// Edge_1 = [x_0,x_2]	 
		// Edge_2 = [x_0,x_3]	 
		// Edge_3 = [x_1,x_2]
		// Edge_4 = [x_1,x_3]	 
		// Edge_5 = [x_2,x_3]

		vec_int_Type edgeNumbersTmp = edgeNumbers;
		for(int i=0; i<6; i++){
			if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[0] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[0]){
				if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[1] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[1])
					edgeNumbers[0] = edgeNumbersTmp[i];
				else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[2] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[2])
					edgeNumbers[1] = edgeNumbersTmp[i];
				else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[3] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[3])
					edgeNumbers[2] = edgeNumbersTmp[i];	
			}
			else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[1] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[1]){
				if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[2] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[2])
					edgeNumbers[3] = edgeNumbersTmp[i];
				else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[3] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[3])
					edgeNumbers[4] = edgeNumbersTmp[i];
			}
			 else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[2] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[2]){
				if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == nodeInd[3] || edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == nodeInd[3])
					edgeNumbers[5] = edgeNumbersTmp[i];
			}
		}

		// We hace 4 Triangles in our Tetraedron
		// If one or more of those Triangles are Part of the domains' boundaries, they are added to the element in question as subelement
		// We extract them in follwing pattern:

		// Tri_0 = [x_0,x_1,x_2] -> Edge 0,1,3
		// Tri_1 = [x_0,x_1,x_3] -> Edge 0,2,4
		// Tri_2 = [x_0,x_2,x_3] -> Edge 1,2,5
		// Tri_3 = [x_1,x_2,x_3] -> Edge 3,4,5

		/*vec2D_int_Type subElements(4,vec_int_Type(3));
		subElements[0] = {edgeNumbers[0],edgeNumbers[1],edgeNumbers[3]};
		subElements[1] = {edgeNumbers[0],edgeNumbers[2],edgeNumbers[4]};
		subElements[2] = {edgeNumbers[1],edgeNumbers[2],edgeNumbers[5]};
		subElements[3] = {edgeNumbers[3],edgeNumbers[4],edgeNumbers[5]};*/

		// We check if one or more of these triangles are part of the boundary surface and determine there flag

		vec2D_int_Type originTriangles(4,vec_int_Type(3));
		originTriangles[0] = {nodeInd[0],nodeInd[1],nodeInd[2]};
		originTriangles[1] = {nodeInd[0],nodeInd[1],nodeInd[3]};
		originTriangles[2] = {nodeInd[0],nodeInd[2],nodeInd[3]};
		originTriangles[3] = {nodeInd[1],nodeInd[2],nodeInd[3]};
		
		
		vec_int_Type originFlag(4,10); // Triangle Flag

		int numberSubElSurf=0;
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);
		int entry; 
		if (elements->getElement(i).subElementsInitialized() ){
			numberSubElSurf = elements->getElement(i).getSubElements()->numberElements();
			for(int k=0; k< numberSubElSurf ; k++){
				triTmp =elements->getElement(i).getSubElements()->getElement(k).getVectorNodeList();
				for(int j=0; j<4 ; j++){
					originTriangleTmp = originTriangles[j];
					sort(originTriangleTmp.begin(),originTriangleTmp.end());
					sort(triTmp.begin(),triTmp.end());
					if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] &&  triTmp[2] == originTriangleTmp[2] ) 
						originFlag[j] = elements->getElement(i).getSubElements()->getElement(k).getFlag();
					//auto it1 = find( originTriangles.begin(), originTriangles.end() ,triTmp );
               		//entry = distance( originTriangles.begin() , it1 );
				}
			}
		}

		// Furthermore we have to determine whether the triangles are part of the interface between processors, as we need this information to determine if edges
		// that emerge on the triangles are part of the interface
		// A triangle is part of the interface if all of its edges are part of the interface (the information if edges are part of the interface was determined
		// in the beginning of the Mesh Refinement by 'determineInterfaceEdges')

		vec_bool_Type interfaceSurface = checkInterfaceSurface(edgeElements,originFlag, edgeNumbers,i);
		
		for(int j=0; j<4; j++){	
			sort( newTriangles.at(j).begin(), newTriangles.at(j).end() );
			FiniteElementNew feNew(originTriangles[j],originFlag[j]);
			feNew.setInterfaceElement(interfaceSurface[j]);
			this->surfaceTriangleElements_->addSurface(feNew, i);
		}
	}
	vec2D_GO_Type combinedSurfaceElements;
	this->surfaceTriangleElements_->sortUniqueAndSetGlobalIDsParallel(this->elementMap_,combinedSurfaceElements);
	
	this->surfaceTriangleElements_->setElementsSurface( combinedSurfaceElements );

	this->surfaceTriangleElements_->setUpElementsOfSurface( this->elementMap_, this->edgeMap_, edgeElements);

	this->updateElementsOfSurfaceLocalAndGlobal(edgeElements);

	vec2D_GO_Type elementsOfSurfaceGlobal = this->surfaceTriangleElements_->getElementsOfSurfaceGlobal();

	/*for(int i=0; i< this->surfaceTriangleElements_->numberElements() ;i++){
		cout << " Surface i="<< i << " liegt an Element " << elementsOfSurfaceGlobal[i][0] ;
		if(elementsOfSurfaceGlobal[i].size()>1)
			cout << " " << elementsOfSurfaceGlobal[i][1] ;
		cout << endl;
	}*/

}



}














