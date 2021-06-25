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
	refinementRestriction_ = parameterListAll->sublist("Mesh Refinement").get("RefinementRestriction","keepRegularity");
	refinement3DDiagonal_ = parameterListAll->sublist("Mesh Refinement").get("3D regular Refinement Diagonal Pick",0);

	writeRefinementTime_ = parameterListAll->sublist("Mesh Refinement").get("Write Refinement Time",true);
	writeMeshQuality_ = parameterListAll->sublist("Mesh Refinement").get("Write Mesh Quality",true);


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

	tol_= parameterListAll->sublist("Mesh Refinement").get("Toleranz",0.001);
	theta_ = parameterListAll->sublist("Mesh Refinement").get("Theta",0.35);
	markingStrategy_ = parameterListAll->sublist("Mesh Refinement").get("RefinementType","Uniform");
	maxIter_ = parameterListAll->sublist("Mesh Refinement").get("MaxIter",3);
	refinementRestriction_ = parameterListAll->sublist("Mesh Refinement").get("Refinement Restriction","keepRegularity");
	refinement3DDiagonal_ = parameterListAll->sublist("Mesh Refinement").get("3D regular Refinement Diagonal Pick",0);

	writeRefinementTime_ = parameterListAll->sublist("Mesh Refinement").get("Write Refinement Time",true);
	writeMeshQuality_ = parameterListAll->sublist("Mesh Refinement").get("Write Mesh Quality",true);

		
	
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
	RefinementFactory<SC,LO,GO,NO> refinementFactory( domainP1->getComm(), inputMeshP1_->volumeID_, refinementRestriction_, refinement3DDiagonal_); 

	// Estimating the error with the Discretizations Mesh.
	int currentLevel =0;
	while(currentLevel < level){		
		errorEstimator.tagArea(inputMeshP1_,area);
		refinementFactory.refineMesh(inputMeshP1_,currentLevel, outputMesh);

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

	// !!!!! Not yes existing function 
	domainRefined->initWithDomain(domainP1);

	inputMeshP12_ = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP12->getMesh() , true);
	inputMeshP12_->FEType_ = domainP12->getFEType();

	this->identifyProblem(solution);

	// Error Estimation object
    ErrorEstimation<SC,LO,GO,NO> errorEstimator (dim_, problemType_ );

	// Refinement Factory object
	RefinementFactory<SC,LO,GO,NO> refinementFactory( domainP1->getComm(), inputMeshP1_->volumeID_, refinementRestriction_, refinement3DDiagonal_); 

	if(dim_ ==3){
		SurfaceElementsPtr_Type surfaceTriangleElements = inputMeshP12_->getSurfaceTriangleElements(); // Surfaces
		if(surfaceTriangleElements.is_null()){
			surfaceTriangleElements.reset(new SurfaceElements()); // Surface
			cout << " Building surfaceTriangleElemenets ... " << endl;
			refinementFactory.buildSurfaceTriangleElements( inputMeshP12_->getElementsC(),inputMeshP12_->getEdgeElements(),surfaceTriangleElements, inputMeshP12_->getEdgeMap(),inputMeshP12_->getElementMap() );
			inputMeshP12_->surfaceTriangleElements_ = surfaceTriangleElements;
			inputMeshP1_->surfaceTriangleElements_ = surfaceTriangleElements;
			cout << " ... done " << endl;
		}
		else if(surfaceTriangleElements->numberElements() ==0){
			cout << " Building surfaceTriangleElemenets " << endl;
			refinementFactory.buildSurfaceTriangleElements( inputMeshP12_->getElementsC(),inputMeshP12_->getEdgeElements(),inputMeshP12_->getSurfaceTriangleElements() , inputMeshP12_->getEdgeMap(),inputMeshP12_->getElementMap() );
			inputMeshP12_->surfaceTriangleElements_ = surfaceTriangleElements;
			inputMeshP1_->surfaceTriangleElements_ = surfaceTriangleElements;
	
			cout << " ... done " << endl;
		}
	}

	// MESH COARSENING
	bool coarsen = false;

	// If coarsen Mesh is false, so consequently we refine the Mesh we go about as folows:
	if(coarsen == false &&  currentIter_ < maxIter_ ){
				
		// Estimating the error with the Discretizations Mesh.
		errorElementsMv_ = errorEstimator.estimateError(inputMeshP12_, inputMeshP1_, solution, rhsFunc_, domainP12->getFEType());

		errorEstimationMv_.push_back(errorElementsMv_);

		errorEstimator.markElements(errorElementsMv_,theta_,markingStrategy_, inputMeshP1_);

   		refinementFactory.refineMesh(inputMeshP1_,currentIter_, outputMesh);


	}
	/*else if(maxIter_ != currentIter_){
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

			errorElements = meshUnstructuredRefined_k_1->determineCoarseningError(meshUnstructuredRefined_k,"backwards"); // (MeshUnstrPtr_Type mesh_k,MeshUnstrPtr_Type mesh_k_m,MultiVectorPtr_Type errorElementMv_k,  string distribution)

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
	}*/

	// Exporting current solution and errorEstimation
	Teuchos::RCP<const MultiVector<SC,LO,GO,NO> >  exportSolutionMv = problem->getSolution()->getBlock(0);
	// Export distribution of elements among processors
	MultiVectorPtr_Type procNumTmp = Teuchos::rcp( new MultiVector_Type(domainP12->getElementMap() , 1 ) );

	procNumTmp->putScalar(comm_->getRank());
	MultiVectorConstPtr_Type vecDecompositionConst = procNumTmp;

	MultiVectorConstPtr_Type errorElConst  = errorElementsMv_ ;

	// Error in Nodes	
	MultiVectorConstPtr_Type exactSolution = this->calcExactSolution();
	calcErrorNorms(exactSolution,solution->getBlock(0));

	MultiVectorConstPtr_Type errorValues = 	errorNodesMv_;; // error of exact vs approx sol

	if(this->exportWithParaview_ && initExporter_==false){
		this->initExporter(  parameterListAll_);
	}
	this->exportSolution( inputMeshP12_, exportSolutionMv, errorValues, exactSolution);
	this->exportError( inputMeshP12_, errorElConst, vecDecompositionConst );
	
	if(currentIter_ == maxIter_){
		writeRefinementInfo();	
		exporterSol_->closeExporter();
	    exporterError_->closeExporter();
	}



	// Determine all essential values
	maxErrorEl.push_back(errorElementsMv_->getMax());
	maxErrorKn.push_back(errorNodesMv_->getMax());
	numElements.push_back(domainP12_->getElementMap()->getMaxAllGlobalIndex()+1);
	numElementsProc.push_back(domainP12_->getElementsC()->numberElements());
	relError.push_back(0);
	eRelError.push_back(0);
	numNodes.push_back(domainP12_->getMapUnique()->getMaxAllGlobalIndex()+1);

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
Calculating error norms. If the exact solution is unknown we use approxmated errorNorm and error indicators
@param[in] exact solution if known
@param[in] FE solution
@param[in] error estimation
*/

template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::calcErrorNorms(MultiVectorConstPtr_Type exactSolution, MultiVectorConstPtr_Type solutionP12){


	MultiVectorPtr_Type errorValues = Teuchos::rcp(new MultiVector_Type( solution_->getBlock(0)->getMap() ) ); 
    //this = alpha*A + beta*B + gamma*this
    errorValues->update( 1., exactSolution, -1. ,solutionP12, 0.);

	MultiVectorConstPtr_Type errorValuesAbs = Teuchos::rcp(new MultiVector_Type(  solution_->getBlock(0)->getMap()) );
	errorValuesAbs = errorValues; 

	errorValues->abs(errorValuesAbs);

	errorNodesMv_ = errorValues;

	errorH1.push_back(sqrt(problem_->calculateH1Norm(errorValues)));

	// L2 Norm is more difficult

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

	//double solElementH1=sqrt(problem_->calculateH1Norm(exactSolution) + problem_->calculateL2Norm(exactSolution));
	//double solhElementH1=sqrt(problem_->calculateH1Norm(solutionP12) + problem_->calculateL2Norm(solutionP12));

}



///
/// ParaView exporter setup
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::initExporter( ParameterListPtr_Type parameterListAll){

	exporterSol_.reset(new ExporterParaViewAMR<SC,LO,GO,NO>());
	exporterError_.reset(new ExporterParaViewAMR<SC,LO,GO,NO>());

	exporterSol_->setup( "Refinement" , domainP12_->getMesh(),  domainP12_->getFEType(), parameterListAll );

	exporterError_->setup("Error_and_Dist", domainP1_->getMesh(), "P0",parameterListAll );

	initExporter_=true;

}

///
/// ParaView exporter export of solution on current mesh
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportSolution(MeshUnstrPtr_Type mesh, MultiVectorConstPtr_Type exportSolutionMv, MultiVectorConstPtr_Type errorValues, MultiVectorConstPtr_Type exactSolutionMv){

	string exporterType = "Scalar";
	if(dofs_ >1 )
		exporterType = "Vector";

	if(currentIter_==0){
		
		exporterSol_->addVariable( exportSolutionMv, "u_h", exporterType, 1, domainP12_->getMapUnique() );
		exporterSol_->addVariable( exactSolutionMv, "u", exporterType, 1, domainP12_->getMapUnique() );
		exporterSol_->addVariable( errorValues, "Error |u-u_h|", exporterType, 1, domainP12_->getMapUnique() );
	}
	else{
		exporterSol_->reSetup(mesh);
		exporterSol_->updateVariables(exportSolutionMv, "u_h");
		exporterSol_->updateVariables( exactSolutionMv, "u" );
		exporterSol_->updateVariables(errorValues, "Error |u-u_h|");
	}
			
	exporterSol_->save( (double) currentIter_);

}


///
/// ParaView exporter export of error values and element distribution on current mesh
///
template <class SC, class LO, class GO, class NO>
void AdaptiveMeshRefinement<SC,LO,GO,NO>::exportError(MeshUnstrPtr_Type mesh, MultiVectorConstPtr_Type errorElConst, MultiVectorConstPtr_Type vecDecompositionConst ){

	if(currentIter_==0){
		exporterError_->addVariable( errorElConst, "ErrorEstimate", "Scalar", 1, domainP1_->getElementMap());
		exporterError_->addVariable( vecDecompositionConst, "Proc", "Scalar", 1, domainP1_->getElementMap());
	}
	else{
		exporterError_->reSetup(mesh);

		exporterError_->updateVariables(errorElConst,"ErrorEstimate");
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

	exportLocalEntry->putScalar( (LO) numElementsProc[currentIter_-1] );

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
			for(int i=1; i< currentIter_; i++)
				cout <<" "<< i << ":	" << numElements[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Number of Nodes after Refinement.... " << endl;
			for(int i=1; i<currentIter_ ; i++)
				cout <<" "<< i << ":	" << numNodes[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Errorestimation: max error in Elements according to error Estimator after Refinement.... " << endl;
			for (int i=1; i<currentIter_ ; i++)
				cout <<" "<< i << ":	" << maxErrorEl[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Maximal error in nodes after Refinement. " << endl;
			for (int i=1; i<currentIter_ ; i++)
				cout <<" "<< i << ":	" << maxErrorKn[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;

			cout << " ||u-u_h||_H1 / ||u ||_H1 	||  eta / ||u_h ||_H1	||	|| u-u_h ||_H1	||	|| u-u_h ||_L2	...." << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			for (int i=1; i<currentIter_ ; i++)
				cout <<" "<< i << ":	" << relError[i] << " 		||	" << eRelError[i] << "  	||	" << errorH1[i]<< "	||	" << errorL2[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
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












