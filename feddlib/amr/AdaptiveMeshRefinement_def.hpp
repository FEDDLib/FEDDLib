#ifndef AdaptiveMeshRefinement_def_hpp
#define AdaptiveMeshRefinement_def_hpp

#ifndef MESH_TIMER_START
#define MESH_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Mesh Refinement") + std::string(S))));
#endif

#ifndef MESH_TIMER_STOP
#define MESH_TIMER_STOP(A) A.reset();
#endif

#include "AdaptiveMeshRefinement_decl.hpp"
#include <chrono> 
/*!
 Definition of AdaptiveMeshRefinement
 
 @brief  AdaptiveMeshRefinement
 @author Lea Saßmannshausen

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
inputMeshP1_(),
inputMeshP12_(),
outputMesh_(),
errorElementsMv_(),
errorEstimationMv_(0),
domainsP1_(0)
{

}


/*!
Initializing problem with the kind of problem (e.g. Laplace, Stokes) for determining the correct error estimation and the dimension
@param[in] problemType, dim
*/
template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement(string problemType, ParameterListPtr_Type parameterListAll ):
inputMeshP1_(),
inputMeshP12_(),
outputMesh_(),
errorElementsMv_(),
errorEstimationMv_(0),
domainsP1_(0)
{
	parameterListAll_ = parameterListAll;
	this->dim_ = parameterListAll->sublist("Parameter").get("Dimension",2);;
	this->problemType_ = problemType;

	this->FEType1_ = "P1";
	this->FEType2_ = parameterListAll->sublist("Parameter").get("Discretization","P1");

	this->exportWithParaview_ = false;

	tol_= parameterListAll->sublist("Mesh Refinement").get("Toleranz",0.001);
	theta_ = parameterListAll->sublist("Mesh Refinement").get("Theta",0.35);
	markingStrategy_ = parameterListAll->sublist("Mesh Refinement").get("RefinementType","Uniform");
	maxIter_ = parameterListAll->sublist("Mesh Refinement").get("MaxIter",3);
	refinementRestriction_ = parameterListAll->sublist("Mesh Refinement").get("Refinement Restriction","Bisection");
	refinement3DDiagonal_ = parameterListAll->sublist("Mesh Refinement").get("3D regular Refinement Diagonal Pick",0);

	writeRefinementTime_ = parameterListAll->sublist("Mesh Refinement").get("Write Refinement Time",true);
	writeMeshQuality_ = parameterListAll->sublist("Mesh Refinement").get("Write Mesh Quality",true);

	restrictionLayer_ =  parameterListAll->sublist("Mesh Refinement").get("Restriction Layer",2);

	coarseningCycle_ =  parameterListAll->sublist("Mesh Refinement").get("Coarsening Cycle",0);
	coarseningM_ =  parameterListAll->sublist("Mesh Refinement").get("Coarsening m",1);
	coarseningN_  = parameterListAll->sublist("Mesh Refinement").get("Coarsening n" ,1);

	refinementMode_  = parameterListAll->sublist("Mesh Refinement").get("Refinement Mode" ,"Regular");
	// If no exact solution is given, we use a dummy function == 0 !!!


}
/// 
/// Initializing problem with the kind of problem, dimension and refinement spectific parameters
///
template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement(string problemType, ParameterListPtr_Type parameterListAll , Func_Type exactSolFunc ):
inputMeshP1_(),
inputMeshP12_(),
outputMesh_(),
errorElementsMv_(),
errorEstimationMv_(0),
domainsP1_(0)
{
	parameterListAll_ = parameterListAll;
	this->dim_ = parameterListAll->sublist("Parameter").get("Dimension",2);;
	this->problemType_ = problemType;

	this->FEType1_ = "P1";
	this->FEType2_ = parameterListAll->sublist("Parameter").get("Discretization","P1");

	exactSolFunc_ = exactSolFunc;
	exactSolInput_ = true;

	exactSolPInput_ = false;

	tol_= parameterListAll->sublist("Mesh Refinement").get("Toleranz",0.001);
	theta_ = parameterListAll->sublist("Mesh Refinement").get("Theta",0.35);
	markingStrategy_ = parameterListAll->sublist("Mesh Refinement").get("RefinementType","Uniform");
	maxIter_ = parameterListAll->sublist("Mesh Refinement").get("MaxIter",3);
	refinementRestriction_ = parameterListAll->sublist("Mesh Refinement").get("Refinement Restriction","Bisection");
	refinement3DDiagonal_ = parameterListAll->sublist("Mesh Refinement").get("3D regular Refinement Diagonal Pick",0);

	writeRefinementTime_ = parameterListAll->sublist("Mesh Refinement").get("Write Refinement Time",true);
	writeMeshQuality_ = parameterListAll->sublist("Mesh Refinement").get("Write Mesh Quality",true);

	restrictionLayer_ =  parameterListAll->sublist("Mesh Refinement").get("Restriction Layer",2);

	coarseningCycle_ =  parameterListAll->sublist("Mesh Refinement").get("Coarsening Cycle",0);
	coarseningM_ =  parameterListAll->sublist("Mesh Refinement").get("Coarsening m",1);
	coarseningN_  = parameterListAll->sublist("Mesh Refinement").get("Coarsening n" ,1);

	refinementMode_  = parameterListAll->sublist("Mesh Refinement").get("Refinement Mode" ,"Regular");

}

/// 
/// Initializing problem with the kind of problem, dimension and refinement spectific parameters
///
template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::AdaptiveMeshRefinement(string problemType, ParameterListPtr_Type parameterListAll , Func_Type exactSolFuncU ,Func_Type exactSolFuncP ):
inputMeshP1_(),
inputMeshP12_(),
outputMesh_(),
errorElementsMv_(),
errorEstimationMv_(0),
domainsP1_(0)
{
	parameterListAll_ = parameterListAll;
	this->dim_ = parameterListAll->sublist("Parameter").get("Dimension",2);;
	this->problemType_ = problemType;

	this->FEType1_ = "P1";
	this->FEType2_ = parameterListAll->sublist("Parameter").get("Discretization","P1");

	exactSolFunc_ = exactSolFuncU;
	exactSolPFunc_ = exactSolFuncP;

	exactSolInput_ = true;
	exactSolPInput_ = true;
	calculatePressure_ = true;

	tol_= parameterListAll->sublist("Mesh Refinement").get("Toleranz",0.001);
	theta_ = parameterListAll->sublist("Mesh Refinement").get("Theta",0.35);
	markingStrategy_ = parameterListAll->sublist("Mesh Refinement").get("RefinementType","Uniform");
	maxIter_ = parameterListAll->sublist("Mesh Refinement").get("MaxIter",3);
	refinementRestriction_ = parameterListAll->sublist("Mesh Refinement").get("Refinement Restriction","Bisection");
	refinement3DDiagonal_ = parameterListAll->sublist("Mesh Refinement").get("3D regular Refinement Diagonal Pick",0);

	writeRefinementTime_ = parameterListAll->sublist("Mesh Refinement").get("Write Refinement Time",true);
	writeMeshQuality_ = parameterListAll->sublist("Mesh Refinement").get("Write Mesh Quality",true);

	restrictionLayer_ =  parameterListAll->sublist("Mesh Refinement").get("Restriction Layer",2);

	coarseningCycle_ =  parameterListAll->sublist("Mesh Refinement").get("Coarsening Cycle",0);
	coarseningM_ =  parameterListAll->sublist("Mesh Refinement").get("Coarsening m",1);
	coarseningN_  = parameterListAll->sublist("Mesh Refinement").get("Coarsening n" ,1);

	refinementMode_  = parameterListAll->sublist("Mesh Refinement").get("Refinement Mode" ,"Regular");
		
	
}
template <class SC, class LO, class GO, class NO>
AdaptiveMeshRefinement<SC,LO,GO,NO>::~AdaptiveMeshRefinement(){

}

template <class SC, class LO, class GO, class NO>
typename AdaptiveMeshRefinement<SC,LO,GO,NO>::DomainPtr_Type AdaptiveMeshRefinement<SC,LO,GO,NO>:: refineArea(DomainPtr_Type domainP1, vec2D_dbl_Type area, int level ){

	DomainPtr_Type domainRefined(new Domain<SC,LO,GO,NO>( domainP1->getComm() , dim_ ));

	inputMeshP1_ = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP1->getMesh() , true);
	inputMeshP1_->FEType_ = domainP1->getFEType();

	MeshUnstrPtr_Type outputMesh(new MeshUnstr_Type(domainP1->getComm(),  inputMeshP1_->volumeID_));

	// !!!!! Not yes existing function 
	domainRefined->initWithDomain(domainP1);

	inputMeshP1_ = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP1->getMesh() , true);
	inputMeshP1_->FEType_ = domainP1->getFEType();

	// Error Estimation object
    ErrorEstimation<SC,LO,GO,NO> errorEstimator (dim_, problemType_ );

	// Refinement Factory object

	RefinementFactory<SC,LO,GO,NO> refinementFactory( domainP1->getComm(), inputMeshP1_->volumeID_, refinementRestriction_, refinement3DDiagonal_, restrictionLayer_); 

	// Estimating the error with the Discretizations Mesh.
	int currentLevel =0;
	while(currentLevel < level){		
		errorEstimator.tagArea(inputMeshP1_,area);
		refinementFactory.refineMesh(inputMeshP1_,currentLevel, outputMesh, refinementMode_);

		inputMeshP1_ = outputMesh;
		currentLevel++;
	}

    domainRefined->setMesh(outputMesh);
	
	return domainRefined;

}

template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::identifyProblem(BlockMultiVectorConstPtr_Type valuesSolution){

	// dofs can be determined by checking the solutions' blocks an comparing the length of the vectors to the number of unique nodes. 
	// If the length is the same: dofs =1 , if it's twice as long: dofs =2 .. and so on

	// P_12 generally represents the velocity domain:
	if(inputMeshP12_->getMapUnique()->getNodeNumElements() == valuesSolution->getBlock(0)->getDataNonConst(0).size())
		dofs_ = 1;
	else if(2*(inputMeshP12_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(0)->getDataNonConst(0).size())
		dofs_ = 2;
	else if(3*(inputMeshP12_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(0)->getDataNonConst(0).size())
		dofs_ = 3;
	
	if(valuesSolution->size() > 1){
		if(inputMeshP1_->getMapUnique()->getNodeNumElements() == valuesSolution->getBlock(1)->getDataNonConst(0).size())
			dofsP_ = 1;
		else if(2*(inputMeshP1_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(1)->getDataNonConst(0).size())
			dofsP_ = 2;
		else if(3*(inputMeshP1_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(1)->getDataNonConst(0).size())
			dofsP_ = 3;
		calculatePressure_=true;
	}

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
typename AdaptiveMeshRefinement<SC,LO,GO,NO>::DomainPtr_Type AdaptiveMeshRefinement<SC,LO,GO,NO>::globalAlgorithm(DomainPtr_Type domainP1, DomainPtr_Type domainP12, BlockMultiVectorConstPtr_Type solution,ProblemPtr_Type problem, RhsFunc_Type rhsFunc ){

	solution_ = solution;

	currentIter_ = domainsP1_.size() ;

	rhsFunc_ = rhsFunc;

	problem_ = problem;

	comm_ = domainP1 ->getComm();

	maxRank_ = std::get<1>(domainP1->getMesh()->rankRange_);
	// We save the domains of each step
	// The P1 Mesh is always used for refinement while the P1 or P2 Mesh is used for error Estimation depending on Discretisation
	domainsP1_.push_back(domainP1);
	domainsP12_.push_back(domainP12);

	domainP1_ = domainP1;
	domainP12_ = domainP12;


	// Reading Mesh from domainP1 as we always refine the P1 Mesh, here defined as inputMesh_
	inputMeshP1_ = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP1->getMesh() , true);
	inputMeshP1_->FEType_ = domainP1->getFEType();
	
	// With the global Algorithm we create a new P1 domain with a new mesh
	DomainPtr_Type domainRefined(new Domain<SC,LO,GO,NO>( domainP1->getComm() , dim_ ));

	// Output Mesh
	MeshUnstrPtr_Type outputMesh(new MeshUnstr_Type(domainP1->getComm(),  inputMeshP1_->volumeID_));

	// Init to be refined domain with inputDomain
	domainRefined->initWithDomain(domainP1);

	inputMeshP12_ = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP12->getMesh() , true);
	inputMeshP12_->FEType_ = domainP12->getFEType();

	this->identifyProblem(solution);

	// Error Estimation object
    ErrorEstimation<SC,LO,GO,NO> errorEstimator (dim_, problemType_ );

	// Refinement Factory object
	RefinementFactory<SC,LO,GO,NO> refinementFactory( domainP1->getComm(), inputMeshP1_->volumeID_, refinementRestriction_, refinement3DDiagonal_,restrictionLayer_); 

	if(dim_ ==3){
		SurfaceElementsPtr_Type surfaceTriangleElements = inputMeshP12_->getSurfaceTriangleElements(); // Surfaces
		if(surfaceTriangleElements.is_null()){
			surfaceTriangleElements.reset(new SurfaceElements()); // Surface
			refinementFactory.buildSurfaceTriangleElements( inputMeshP12_->getElementsC(),inputMeshP12_->getEdgeElements(),surfaceTriangleElements, inputMeshP12_->getEdgeMap(),inputMeshP12_->getElementMap() );
			inputMeshP12_->surfaceTriangleElements_ = surfaceTriangleElements;
			inputMeshP1_->surfaceTriangleElements_ = surfaceTriangleElements;
		}
		else if(surfaceTriangleElements->numberElements() ==0){
			refinementFactory.buildSurfaceTriangleElements( inputMeshP12_->getElementsC(),inputMeshP12_->getEdgeElements(),inputMeshP12_->getSurfaceTriangleElements() , inputMeshP12_->getEdgeMap(),inputMeshP12_->getElementMap() );
			inputMeshP12_->surfaceTriangleElements_ = surfaceTriangleElements;
			inputMeshP1_->surfaceTriangleElements_ = surfaceTriangleElements;
		}
	}
	if(currentIter_ == 0 && dim_ == 2){
		refinementFactory.assignEdgeFlags( inputMeshP1_,inputMeshP12_->getEdgeElements());
		refinementFactory.assignEdgeFlags( inputMeshP1_,inputMeshP1_->getEdgeElements());
	}
	// If coarsen Mesh is false, so consequently we refine the Mesh we go about as folows:	
	bool coarsening= false;
	if(coarseningCycle_ > 0 && currentIter_>0){
		if(currentIter_ % coarseningCycle_ == 0)
			coarsening = true;
	}
	errorElementsMv_ =Teuchos::rcp( new MultiVector_Type( domainP12_ ->getElementMap(), 1 ) );
	errorElementsMv_->putScalar(0.);
	if( coarsening== true &&  currentIter_ < maxIter_ ){

		// We start by calculating the error of the current mesh. As this is out starting point for mesh coarsening. In the previous iteration we calculated the error estimation beforehand.
		errorElementsMv_ = errorEstimator.estimateError(inputMeshP12_, inputMeshP1_, solution, rhsFunc_, domainP12->getFEType());

		errorEstimationMv_.push_back(errorElementsMv_);

		int m= coarseningM_;
		int n= coarseningN_;
		int k = currentIter_;
		int iterC;
		MeshUnstrPtrArray_Type meshUnstructuredP1(currentIter_+n);

		for(int i=0; i<currentIter_+1-m; i++)
			meshUnstructuredP1[i] = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainsP1_[i]->getMesh() , true);
			
		// We extract the error estimation of the mesh iter

		MultiVectorPtr_Type errorElements; // = meshUnstructuredRefined_k->mesh_->getErrorEstimate();
		MeshUnstrPtr_Type meshUnstructuredRefined_k ;
		MeshUnstrPtr_Type meshUnstructuredRefined_k_1;
		MeshUnstrPtr_Type meshUnstructuredRefined_k_m_1;
		for(int i=0; i<m-1 ; i++){
			meshUnstructuredRefined_k = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainsP1_[currentIter_-i]->getMesh() , true); 
			meshUnstructuredRefined_k_1 = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainsP1_[currentIter_-1-i]->getMesh() , true); 

			errorElements = errorEstimator.determineCoarseningError(meshUnstructuredRefined_k,meshUnstructuredRefined_k_1,errorEstimationMv_[currentIter_-i],"backwards",markingStrategy_,theta_); 

			errorEstimationMv_[currentIter_-1-i]= errorElements;

		}

		// Setting error of: Mesh_(k-m+1) with the previous error ->downscaling errors
		if(m>1)
			meshUnstructuredRefined_k_m_1 = meshUnstructuredRefined_k_1; 
		else
			meshUnstructuredRefined_k_m_1 =Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainsP1_[currentIter_]->getMesh() , true); 

		// Error of Level l is at l-1
		for(int i=0; i< n; i++){
			iterC = k-m+i;
			if(i==0){
				errorElements = errorEstimator.determineCoarseningError(meshUnstructuredRefined_k_m_1,meshUnstructuredP1[iterC],errorEstimationMv_[iterC+1],"backwards",markingStrategy_,theta_); 
			}
			else{
				errorElements = errorEstimator.determineCoarseningError(meshUnstructuredP1[iterC-1],meshUnstructuredP1[iterC],errorEstimationMv_[iterC-1],"forwards",markingStrategy_,theta_); 
			}
			if(iterC > errorEstimationMv_.size())
				errorEstimationMv_.push_back(errorElements);
			else 
				errorEstimationMv_[iterC]=errorElements;


			errorEstimator.markElements(errorElements,theta_,markingStrategy_, meshUnstructuredP1[iterC]);

   			refinementFactory.refineMesh(meshUnstructuredP1[iterC],iterC, outputMesh, refinementMode_);

			meshUnstructuredP1[iterC+1] = outputMesh;
			outputMesh.reset(new MeshUnstr_Type(domainP1->getComm(),  inputMeshP1_->volumeID_));
		}
		
		outputMesh = meshUnstructuredP1[iterC+1];
	}

	else if( currentIter_ < maxIter_ ){			
		// Estimating the error with the Discretizations Mesh.
		errorElementsMv_ = errorEstimator.estimateError(inputMeshP12_, inputMeshP1_, solution, rhsFunc_, domainP12->getFEType());

		errorEstimationMv_.push_back(errorElementsMv_);

		errorEstimator.markElements(errorElementsMv_,theta_,markingStrategy_, inputMeshP1_);

   		refinementFactory.refineMesh(inputMeshP1_,currentIter_, outputMesh, refinementMode_);
	}


	// Export distribution of elements among processors
	MultiVectorPtr_Type procNumTmp = Teuchos::rcp( new MultiVector_Type(domainP12->getElementMap() , 1 ) );

	procNumTmp->putScalar(comm_->getRank());
	MultiVectorConstPtr_Type vecDecompositionConst = procNumTmp;

	// Error in Nodes	
	MultiVectorConstPtr_Type exactSolution = this->calcExactSolution();

	MultiVectorConstPtr_Type exactSolutionP;
	MultiVectorConstPtr_Type  exportSolutionPMv;
	
	if( calculatePressure_ ){
		exportSolutionPMv = problem->getSolution()->getBlock(1);
		if(exactSolPInput_){
			exactSolutionP = this->calcExactSolutionP();
		}
	}

	calcErrorNorms(exactSolution,solution->getBlock(0), exactSolutionP);

	if(this->exportWithParaview_ && initExporter_==false){
		this->initExporter(  parameterListAll_);
	}
	this->exportSolution( inputMeshP12_, problem->getSolution()->getBlock(0), errorNodesMv_, exactSolution, exportSolutionPMv, exactSolutionP);
	if( currentIter_< maxIter_)
		this->exportError( inputMeshP12_, errorElementsMv_, errorH1ElementsMv_ , difH1EtaElementsMv_ , vecDecompositionConst );


	// Determine all essential values
	maxErrorEl.push_back(errorElementsMv_->getMax());
	maxErrorKn.push_back(errorNodesMv_->getMax());
	numElements.push_back(domainP12_->getElementMap()->getMaxAllGlobalIndex()+1);
	numElementsProc.push_back(domainP12_->getElementsC()->numberElements());
	
	numNodes.push_back(domainP12_->getMapUnique()->getMaxAllGlobalIndex()+1);

	if(currentIter_ == maxIter_){
		writeRefinementInfo();	
		exporterSol_->closeExporter();
	    exporterError_->closeExporter();
	}

    domainRefined->setMesh(outputMesh);
	
	return domainRefined;

}

/*template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::determineCoarsening(){



}*/



/*!
\brief Calculating exact solution if possible with exactSolFunc_

*/

template <class SC, class LO, class GO, class NO>
typename AdaptiveMeshRefinement<SC,LO,GO,NO>:: MultiVectorConstPtr_Type AdaptiveMeshRefinement<SC,LO,GO,NO>::calcExactSolution(){
	
    //if ( !rhsFuncVec_[i].empty() )

	MultiVectorPtr_Type exactSolution = Teuchos::rcp(new MultiVector_Type( solution_->getBlock(0)->getMap() ) ); 
	Teuchos::ArrayRCP<SC> exactSolA = exactSolution->getDataNonConst(0);

	vec2D_dbl_ptr_Type points = inputMeshP12_->getPointsUnique();

	Teuchos::ArrayRCP<SC> exactSol(dofs_);
	for(int i=0; i< points->size(); i++){
		exactSolFunc_(&points->at(i).at(0),&exactSol[0]);
		for(int j=0; j< dofs_ ; j++)
			exactSolA[i*dofs_+j] = exactSol[j];

	}

	MultiVectorConstPtr_Type exactSolConst = exactSolution;

	return exactSolConst;
}

/*!
\brief Calculating exact solution if possible with exactSolPFunc_

*/

template <class SC, class LO, class GO, class NO>
typename AdaptiveMeshRefinement<SC,LO,GO,NO>:: MultiVectorConstPtr_Type AdaptiveMeshRefinement<SC,LO,GO,NO>::calcExactSolutionP(){
	
    //if ( !rhsFuncVec_[i].empty() )

	MultiVectorPtr_Type exactSolution = Teuchos::rcp(new MultiVector_Type( domainP1_->getMapUnique())); 
	Teuchos::ArrayRCP<SC> exactSolA = exactSolution->getDataNonConst(0);

	vec2D_dbl_ptr_Type points = domainP1_->getPointsUnique();

	Teuchos::ArrayRCP<SC> exactSol(dofsP_);
	for(int i=0; i< points->size(); i++){
		exactSolPFunc_(&points->at(i).at(0),&exactSol[0]);
		exactSolA[i] = exactSol[0];
	}

	MultiVectorConstPtr_Type exactSolConst = exactSolution;

	return exactSolConst;
}

/*!
Calculating error norms. If the exact solution is unknown we use approxmated errorNorm and error indicators
@param[in] exact solution if known
@param[in] FE solution
@param[in] error estimation
*/

template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::calcErrorNorms(MultiVectorConstPtr_Type exactSolution, MultiVectorConstPtr_Type solutionP12,MultiVectorConstPtr_Type exactSolutionP){

	// Calculating the error per node
	MultiVectorPtr_Type errorValues = Teuchos::rcp(new MultiVector_Type( solution_->getBlock(0)->getMap() ) ); 
	//this = alpha*A + beta*B + gamma*this
	errorValues->update( 1., exactSolution, -1. ,solutionP12, 0.);

	// Taking abs norm
	MultiVectorConstPtr_Type errorValuesAbs = Teuchos::rcp(new MultiVector_Type(  solution_->getBlock(0)->getMap()) );
	errorValuesAbs = errorValues; 
	errorValues->abs(errorValuesAbs);

	// Absolute Error in nodes
	errorNodesMv_ = errorValues;

	// ---------------------------
	// Calculating H1 Norm
	errorH1.push_back(sqrt(problem_->calculateH1Norm(errorValues)));

	// ---------------------------
	// L2 Norm 
	MultiVectorConstPtr_Type exactSolutionTmp = Teuchos::rcp(new MultiVector_Type( domainP12_ ->getMapUnique() ) ); 
	Teuchos::ArrayRCP<double > exactSolutionTmpA = exactSolutionTmp->getDataNonConst(0);

	MultiVectorConstPtr_Type solutionTmp = Teuchos::rcp(new MultiVector_Type( domainP12_ ->getMapUnique() ) ); 
	Teuchos::ArrayRCP<double > solutionTmpA = solutionTmp->getDataNonConst(0);

	Teuchos::ArrayRCP<double > exactSolutionA = exactSolution->getDataNonConst(0);

	Teuchos::ArrayRCP<double > solutionP12A = solutionP12->getDataNonConst(0);

	double errorL2Tmp=0;
	for(int i=0; i< dofs_ ; i++){

		MultiVectorPtr_Type errorValues = Teuchos::rcp(new MultiVector_Type( domainP12_ ->getMapUnique() ) );
 		for(int j=0; j< solutionTmpA.size(); j++){
			solutionTmpA[j] = solutionP12A[j*dofs_+i];
			exactSolutionTmpA[j] = exactSolutionA[j*dofs_+i];
		}	 

	    //this = alpha*A + beta*B + gamma*this
		errorValues->update( 1., exactSolutionTmp, -1. ,solutionTmp, 0.);

		MultiVectorConstPtr_Type errorValuesAbs = Teuchos::rcp(new MultiVector_Type(  domainP12_ ->getMapUnique()) );
		errorValuesAbs = errorValues; 

		errorValues->abs(errorValuesAbs);

		errorL2Tmp += problem_->calculateL2Norm(errorValues);

	}

	errorL2.push_back(sqrt(errorL2Tmp));

	// -------------------------------
	// Calculating Error bound epsilon
	if(exactSolInput_ == true){
		relError.push_back(sqrt(problem_->calculateH1Norm(errorValues)) / sqrt(problem_->calculateH1Norm(exactSolution)));
	}

	if(exactSolInput_ == true){
		double eta = 0.;
	    Teuchos::ArrayRCP<const double > errorElement = errorElementsMv_->getData(0);
		for(int i=0; i < errorElement.size() ; i++){
			//eta += errorElement[i];
			if(eta < errorElement[i])
				eta=errorElement[i];		
		}
		reduceAll<int, double> (*comm_, REDUCE_MAX, eta, outArg (eta));
		//eta = pow(eta,2);


		eRelError.push_back(sqrt(eta)/ sqrt(problem_->calculateH1Norm(solutionP12)));
	}

	// -------------------------------
	// Calculating H1 Norm elementwise
	// First we need a repeated map that is compatible with a vector solution and error

	errorH1ElementsMv_ =  Teuchos::rcp( new MultiVector_Type( domainP12_ ->getElementMap(), 1 ) );
	errorH1ElementsMv_->putScalar(0.);
	Teuchos::ArrayRCP<SC> errorH1ElementsA = errorH1ElementsMv_->getDataNonConst(0);
	
	difH1EtaElementsMv_ =  Teuchos::rcp( new MultiVector_Type( domainP12_ ->getElementMap(), 1 ) );
	difH1EtaElementsMv_->putScalar(0.);

	if(domainP12_->getElementMap()->getMaxAllGlobalIndex()< 1500){
		ElementsPtr_Type elements = domainP12_->getElementsC();
		MapConstPtr_Type elementMap = domainP12_->getElementMap();
		MapConstPtr_Type mapUnique = domainP12_->getMapUnique();
		MapConstPtr_Type mapRep = domainP12_->getMapRepeated();

		vec_GO_Type repIDsVec(0);
		for(int i=0; i<mapRep->getMaxLocalIndex()+1; i++){
			GO gID = mapRep->getGlobalElement(i);
			for(int d=0; d < dofs_ ; d++)
				repIDsVec.push_back(gID*dofs_+d);
		}
		Teuchos::ArrayView<GO> repIDsVecArray = Teuchos::arrayViewFromVector(repIDsVec);
		// global Ids of Elements' Nodes
		MapPtr_Type mapRepSystem =
				Teuchos::rcp( new Map_Type( elementMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), repIDsVecArray , 0, domainP12_->getComm()) );

		MultiVectorPtr_Type mvValuesError =  Teuchos::rcp( new MultiVector_Type( mapRepSystem, 1 ) );
		Teuchos::ArrayRCP< SC > mvValuesErrorA  = mvValuesError->getDataNonConst(0);	

		MultiVectorPtr_Type mvValuesErrorUnique =  Teuchos::rcp( new MultiVector_Type( solution_->getBlock(0)->getMap(), 1 ) );
		//Teuchos::ArrayRCP< SC > mvValuesErrorA  = mvValuesError->getDataNonConst(0);	


		MultiVectorPtr_Type errorNodesRep =  Teuchos::rcp( new MultiVector_Type( mapRepSystem, 1 ) );
		Teuchos::ArrayRCP< SC > errorNodesRepA  = errorNodesRep->getDataNonConst(0);	
		errorNodesRep->importFromVector(errorNodesMv_,false,"Insert");


		for(int k=0; k< elementMap->getMaxAllGlobalIndex()+1; k++){
			mvValuesError->putScalar(0.);	
			vec_GO_Type notOnMyProc(0);
			vec_dbl_Type notOnMyProcValue(0);
			if(elementMap->getLocalElement(k) != -1){
				vec_int_Type nodeList = elements->getElement(elementMap->getLocalElement(k)).getVectorNodeList();
				for(int j=0; j< nodeList.size(); j++){
					for(int d=0; d < dofs_ ; d++){
						if(mapUnique->getLocalElement(mapRep->getGlobalElement(nodeList[j])) == -1){
							GO gID = mapRep->getGlobalElement(nodeList[j]);
							notOnMyProc.push_back(gID*dofs_+d);
							notOnMyProcValue.push_back(errorNodesRepA[dofs_*nodeList[j]+d]);
						}

						mvValuesErrorA[dofs_*nodeList[j]+d] = errorNodesRepA[dofs_*nodeList[j]+d];
					}
				}
			}
			Teuchos::ArrayView<GO> globalNodeArray = Teuchos::arrayViewFromVector( notOnMyProc);

			// global Ids of Elements' Nodes
			MapPtr_Type mapNodeExport =
				Teuchos::rcp( new Map_Type( elementMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalNodeArray, 0, domainP12_->getComm()) );
					
			MultiVectorPtr_Type notMV  =  Teuchos::rcp( new MultiVector_Type( mapNodeExport, 1 ) );
			Teuchos::ArrayRCP<SC> notMVA = notMV->getDataNonConst(0);
			for(int i=0; i< notMVA.size(); i++)
				notMVA[i] = notOnMyProcValue[i];
			
			mvValuesErrorUnique->importFromVector(mvValuesError,false,"Insert");
			mvValuesErrorUnique->importFromVector(notMV,false,"Insert");
		
			double valueH1 = problem_->calculateH1Norm(mvValuesErrorUnique);
			double valueL2 = 0; // problem_->calculateL2Norm(mvValuesErrorUnique);
			if(elementMap->getLocalElement(k) != -1){
				errorH1ElementsA[elementMap->getLocalElement(k)]= sqrt(valueH1 + valueL2);
			}
		
		}

		MultiVectorConstPtr_Type errorH1 = errorH1ElementsMv_;
		
		MultiVectorConstPtr_Type errorElements = errorElementsMv_;

		difH1EtaElementsMv_->update( 1., errorElements , -1. , errorH1, 0.);
	}
	if( calculatePressure_== true  && exactSolPInput_ == true  ){
		// Calculating the error per node
		MultiVectorPtr_Type errorValuesP = Teuchos::rcp(new MultiVector_Type( domainP1_->getMapUnique() ) ); 

		//this = alpha*A + beta*B + gamma*this
		errorValuesP->update( 1., exactSolutionP, -1. ,solution_->getBlock(1), 0.);

		// Taking abs norm
		MultiVectorConstPtr_Type errorValuesPAbs = Teuchos::rcp(new MultiVector_Type( domainP1_->getMapUnique() ) );
		errorValuesPAbs = errorValuesP; 
		errorValuesP->abs(errorValuesPAbs);

		errorNodesPMv_ = errorValuesP;

		double errorL2Tmp = problem_->calculateL2Norm(errorValuesP,1);

		errorL2P.push_back(errorL2Tmp);
	}


	
}



///
/// ParaView exporter setup
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::initExporter( ParameterListPtr_Type parameterListAll){

	exporterSol_.reset(new ExporterParaViewAMR<SC,LO,GO,NO>());
	exporterError_.reset(new ExporterParaViewAMR<SC,LO,GO,NO>());

	exporterSol_->setup( "RefinementU" , domainP12_->getMesh(),  domainP12_->getFEType(), parameterListAll );

	if(calculatePressure_ ){
		exporterSolP_.reset(new ExporterParaViewAMR<SC,LO,GO,NO>());
		exporterSolP_->setup( "RefinementP" , domainP1_->getMesh(),  domainP1_->getFEType(), parameterListAll );
	}

	exporterError_->setup("ErrorEstimation", domainP1_->getMesh(), "P0",parameterListAll );

	initExporter_=true;

}

///
/// ParaView exporter export of solution on current mesh
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportSolution(MeshUnstrPtr_Type mesh, MultiVectorConstPtr_Type exportSolutionMv, MultiVectorConstPtr_Type errorValues, MultiVectorConstPtr_Type exactSolutionMv,MultiVectorConstPtr_Type exportSolutionPMv,MultiVectorConstPtr_Type exactSolutionPMv){

	string exporterType = "Scalar";
	if(dofs_ >1 )
		exporterType = "Vector";

	if(currentIter_==0){
		
		exporterSol_->addVariable( exportSolutionMv, "u_h", exporterType, dofs_, domainP12_->getMapUnique() );
		exporterSol_->addVariable( exactSolutionMv, "u", exporterType, dofs_, domainP12_->getMapUnique() );
		exporterSol_->addVariable( errorValues, "|u-u_h|", exporterType, dofs_, domainP12_->getMapUnique() );


		if( calculatePressure_ ){
			exporterSolP_->addVariable( exportSolutionPMv, "p_h", "Scalar", dofsP_, domainP1_->getMapUnique() );
			if(exactSolPInput_){
				exporterSolP_->addVariable( exactSolutionPMv, "p", "Scalar", dofsP_, domainP1_->getMapUnique() );
				exporterSolP_->addVariable( errorNodesPMv_, "|p-p_h|", "Scalar", dofsP_, domainP1_->getMapUnique() );
			}
			exporterSolP_->save( (double) currentIter_);
		}
		
	}
	else{
		exporterSol_->reSetup(mesh);
		exporterSol_->updateVariables(exportSolutionMv, "u_h");
		exporterSol_->updateVariables( exactSolutionMv, "u" );
		exporterSol_->updateVariables(errorValues, "|u-u_h|");

		if( calculatePressure_ ){
			exporterSolP_->reSetup(domainP1_->getMesh());
			exporterSolP_->updateVariables( exportSolutionPMv, "p_h");
			if(exactSolPInput_ ){
				exporterSolP_->updateVariables( exactSolutionPMv, "p");
				exporterSolP_->updateVariables( errorNodesPMv_, "|p-p_h|");
			}

		exporterSolP_->save( (double) currentIter_);
		}
	}
			
	exporterSol_->save( (double) currentIter_);

}


///
/// ParaView exporter export of error values and element distribution on current mesh
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportError(MeshUnstrPtr_Type mesh, MultiVectorConstPtr_Type errorElConst, MultiVectorConstPtr_Type errorElConstH1 , MultiVectorConstPtr_Type difH1Eta ,MultiVectorConstPtr_Type vecDecompositionConst ){

	if(currentIter_==0){
		exporterError_->addVariable( errorElConst, "eta_T", "Scalar", 1, domainP1_->getElementMap());
		exporterError_->addVariable( errorElConstH1, "||u-u_h||_H1(T)", "Scalar", 1, domainP1_->getElementMap());
		exporterError_->addVariable( difH1Eta, "eta_T-||u-u_h||_H1(T)", "Scalar", 1, domainP1_->getElementMap());
		exporterError_->addVariable( vecDecompositionConst, "Proc", "Scalar", 1, domainP1_->getElementMap());
	}
	else{
		exporterError_->reSetup(mesh);

		exporterError_->updateVariables(errorElConst,"eta_T");
		exporterError_->updateVariables(errorElConstH1,"||u-u_h||_H1(T)");
		exporterError_->updateVariables(difH1Eta, "eta_T-||u-u_h||_H1(T)");
		exporterError_->updateVariables(vecDecompositionConst,"Proc");
	}

    exporterError_->save((double) currentIter_);


}

/*!
Writing refinement information
*/
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::writeRefinementInfo(){

	vec_GO_Type globalProcs(0);
	for (int i=0; i<= maxRank_; i++)
			globalProcs.push_back(i);

	Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

	vec_GO_Type localProc(0);
	localProc.push_back(comm_->getRank());
	Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

	MapPtr_Type mapGlobalProc =
		Teuchos::rcp( new Map_Type( domainP1_->getElementMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, comm_) );

	// Global IDs of Procs
	MapPtr_Type mapProc =
		Teuchos::rcp( new Map_Type( domainP1_->getElementMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, comm_) );
	
	MultiVectorPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVector_Type( mapProc, 1 ) );

	exportLocalEntry->putScalar( (LO) numElementsProc[currentIter_] );

	MultiVectorPtr_Type elementList= Teuchos::rcp( new MultiVector_Type( mapGlobalProc, 1 ) );
	elementList->putScalar( 0 ); 
	elementList->importFromVector( exportLocalEntry, true, "Insert");

	Teuchos::ArrayRCP<const double > elementProcList = elementList->getData(0);




	if(comm_->getRank() == 0){	
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Summary Mesh Refinement" << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Marking Strategy:	" << markingStrategy_ << endl;
			cout << " Theta:			" << theta_ << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Tolerance:			" << tol_ << endl;
			cout << " Max number of Iterations:	" <<  maxIter_ << endl;
			cout << " Number of Processors:		" << maxRank_ +1 << endl;
			cout << " Number of Refinements:		" << currentIter_ << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Number of elements after Refinement.... " << endl;
			for(int i=1; i<= currentIter_; i++)
				cout <<" "<< i << ":	" << numElements[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Number of Nodes after Refinement.... " << endl;
			for(int i=1; i<=currentIter_ ; i++)
				cout <<" "<< i << ":	" << numNodes[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Errorestimation: max error in Elements according to error Estimator in Refinement.... " << endl;
			for (int i=0; i<currentIter_ ; i++)
				cout <<" "<< i << ":	" << maxErrorEl[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Maximal error in nodes after Refinement. " << endl;
			for (int i=1; i<=currentIter_ ; i++)
				cout <<" "<< i << ":	" << maxErrorKn[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " || u-u_h ||_H1	||	|| u-u_h ||_L2  ||" ;
			if( calculatePressure_== true  && exactSolPInput_ == true  ){
				cout << " 	|| p-p_h||_L2 " << endl;
			}
			else
				cout << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			for (int i=1; i<=currentIter_ ; i++){
				cout <<" "<< i << ":	"  << errorH1[i]<< "	||	" << errorL2[i] ;
				if( calculatePressure_== true  && exactSolPInput_ == true  ){
					cout << "  	||	" <<  errorL2P[i] << endl;
				}
				else
					cout << endl;
			}
			cout << "__________________________________________________________________________________________________________ " << endl;

			cout << " ||u-u_h||_H1 / ||u ||_H1 	||  eta / ||u_h ||_H1	" << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			for (int i=1; i<=currentIter_ ; i++){
				cout <<" "<< i << ":	" << relError[i] << " 		||	" << eRelError[i]  << endl;
			}
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << "Distribution of elements on .. " << endl;
			for(int l=0; l< maxRank_ +1 ; l++)
				cout <<" Processor "<< l << " carries " << elementProcList[l] << " Elements "<< endl; 
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
		}
		
}



}

#endif












