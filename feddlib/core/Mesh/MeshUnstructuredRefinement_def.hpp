#ifndef MESHUNSTRUCTUREDREFINEMENT_def_hpp
#define MESHUNSTRUCTUREDREFINEMENT_def_hpp
#include "MeshUnstructuredRefinement_decl.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector_def.hpp"
/*!
 Definition of MeshUnstructuredRefinement
 
 @brief  MeshUnstructuredRefinement
 @version 1.0

 */

using namespace std;
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::outArg;

namespace FEDD {

template <class SC, class LO, class GO, class NO>
MeshUnstructuredRefinement<SC,LO,GO,NO>::MeshUnstructuredRefinement():
MeshUnstructured<SC,LO,GO,NO>()
{

  
}

template <class SC, class LO, class GO, class NO>
MeshUnstructuredRefinement<SC,LO,GO,NO>::MeshUnstructuredRefinement(CommConstPtr_Type comm, int volumeID):
MeshUnstructured<SC,LO,GO,NO>(comm,volumeID)
{

}

template <class SC, class LO, class GO, class NO>
MeshUnstructuredRefinement<SC,LO,GO,NO>::MeshUnstructuredRefinement(CommConstPtr_Type comm, int volumeID, MeshUnstrPtr_Type meshP1):
MeshUnstructured<SC,LO,GO,NO>(comm,volumeID)
{
	this->dim_ = meshP1->dim_;
	this->FEType_ = meshP1->FEType_;
    this->volumeID_ = meshP1->volumeID_;
	this->rankRange_=meshP1->rankRange_;
    EdgeElementsPtr_Type edgeElements = meshP1->getEdgeElements(); // Edges
    ElementsPtr_Type elements = meshP1->getElementsC(); // Elements
    this->elementMap_ = meshP1->elementMap_;
	this->mapUnique_ = meshP1->mapUnique_;
	this->mapRepeated_=meshP1->mapRepeated_;
	this->edgeMap_ = meshP1->edgeMap_;
    this->edgeElements_.reset(new EdgeElements(*edgeElements)); // Edges
    this->elementsC_.reset(new Elements(*elements));    // Elements set with the elements we already have
    this->pointsRep_.reset(new std::vector<std::vector<double> >(meshP1->pointsRep_->size(),vector<double>(this->dim_,-1.)));
    *this->pointsRep_ = *meshP1->pointsRep_; // Points
	this->pointsUni_.reset(new std::vector<std::vector<double> >( this->mapUnique_->getNodeNumElements(), vector<double>(this->dim_,-1. ) ) );
    *this->pointsUni_ = *meshP1->pointsUni_; 
	this->bcFlagUni_.reset( new std::vector<int> ( this->mapUnique_->getNodeNumElements(), 0 ) );
    *this->bcFlagUni_ = *meshP1->bcFlagUni_; 
	this->bcFlagRep_.reset(new vector<int>(meshP1->bcFlagRep_->size())); // Flags
	*this->bcFlagRep_ = *meshP1->bcFlagRep_;
	this->edgesElementOrder_=meshP1->getEdgeElementOrder();
	this->numElementsGlob_ = meshP1->numElementsGlob_; 
}

template <class SC, class LO, class GO, class NO>
MeshUnstructuredRefinement<SC,LO,GO,NO>::~MeshUnstructuredRefinement(){
}


// Residualbased A-posteriori ErrorEstimation as proposed in Verfuerths' "A Posteriori Error Estimation Techniques for Finite Element Methods"
template <class SC, class LO, class GO, class NO>
vec_dbl_Type MeshUnstructuredRefinement<SC,LO,GO,NO>::errorEstimation(MultiVectorPtrConst_Type valuesSolution, double theta, string strategy){
   		
	auto startError = std::chrono::high_resolution_clock::now();

 	ElementsPtr_Type elements = this->getElementsC();
	Teuchos::ArrayRCP<const SC> valuesLaplace = valuesSolution->getDataNonConst(0);
	int dim = this->dim_;

    EdgeElementsPtr_Type edgeElements = this->getEdgeElements();

	string FEType = this->FEType_;
	MapConstPtr_Type elementMap = this->getElementMap();

	vec2D_dbl_ptr_Type points = this->getPointsRepeated();
	
	int myRank = this->comm_->getRank();


	edgeElements->matchEdgesToElements(this->getElementMap());
	vec2D_GO_Type elementsOfEdgesGlobal = edgeElements->getElementsOfEdgeGlobal();
	vec2D_LO_Type elementsOfEdgesLocal = edgeElements->getElementsOfEdgeLocal();

	// Vector size of elements for the resdual based error estimate
	vec_dbl_Type errorElement(elements->numberElements());
	double maxErrorEl;
	double maxErrorElLoc=0.0;


  if(this->dim_ == 2){ 
	// Jump per Element over edges
	vec_dbl_Type errorEdgesInterior(3); 
	// Edge Numbers of Element
	vec_int_Type edgeNumbers(3);
	// tmp values u_1,u_2 of gradient u of two elements corresponing with edge 
	vec_dbl_Type u1(2),u2(2);
	
	// gradient of u elementwise
	vec2D_dbl_Type u_El(elements->numberElements(),vec_dbl_Type(2));

	vec_dbl_Type u1Tmp(2),u2Tmp(2);

	vec_int_Type kn1(3),kn2(3);
	vec_dbl_Type p1_2(2);
	// necessary entities
	double h_T ; //diameter of Element
	double h_E ; // something with edge
	vec_dbl_Type v_E(2); // normal Vector of edges
	double norm_v_E;

	Teuchos::ArrayRCP<SC> systemVector1;
	Teuchos::ArrayRCP<SC> systemVector2;

	int elTmp1,elTmp2;

    SC detB1;
    SC absDetB1;
    SmallMatrix<SC> B1(dim);
    SmallMatrix<SC> Binv1(dim);  
	int index0,index;
	

	double deriPhi[3][2]= { -1,-1,  1,0,  0,1 };

	double deriPhiT1[3][2];
	double deriPhiT2[3][2]; 

	B1[0][0]=1;
	B1[1][0]=0;
	B1[0][1]=0;
	B1[1][1]=1;

	// We import the values of the solution in order to have a repeated distributed solution on the processors, which we need for the error Estimation
	MultiVectorPtr_Type valuesSolutionRepeated = Teuchos::rcp( new MultiVector_Type( this->getMapRepeated(), 1 ) );
	valuesSolutionRepeated->putScalar( 0);
	valuesSolutionRepeated->importFromVector(valuesSolution , false, "Insert");

	Teuchos::ArrayRCP< SC > valuesSolutionRep  = valuesSolutionRepeated->getDataNonConst(0);

	// First we determine u on each element
	for(int k=0; k< elements->numberElements();k++){
		for (int l=0; l< 3; l++){
				kn1[l]= elements->getElement(k).getNode(l);
					}
		// Transformation Matrices
		 index0 = kn1[0];
		 for (int s=0; s<dim; s++) {
					index = kn1[s+1];
					for (int t=0; t<dim; t++) {
						B1[t][s] = this->getPointsRepeated()->at(index).at(t) - this->getPointsRepeated()->at(index0).at(t);
					}
				}

		detB1 = B1.computeInverse(Binv1);
		detB1 = std::fabs(detB1);

		for(int l=0; l<2; l++){
			for(int p=0;p<3; p++){
				deriPhiT1[p][l]=(deriPhi[p][0]*Binv1[0][l]+deriPhi[p][1]*Binv1[1][l]);
			}
		}
		
		for(int l=0; l<2; l++){
				u_El[k][l]=valuesSolutionRep[kn1[0]]*deriPhiT1[0][l]+valuesSolutionRep[kn1[1]]*deriPhiT1[1][l]+valuesSolutionRep[kn1[2]]*deriPhiT1[2][l];

			}
	}

	// In case we have more than one proc we need to exchange information via the interface. 
	// We design a multivector consisting of u's x and y values, to import and export it easier to other procs only for interface elements

	// Step 1: create elementMap for Interface elements to import and export them

	vec_GO_Type elementImportMap(0);

	for(int i=0; i<edgeElements->numberElements(); i++){
		if(elementsOfEdgesLocal.at(i).size() >1){
			if(elementsOfEdgesLocal.at(i).at(0) == -1){
				elementImportMap.push_back(elementsOfEdgesGlobal.at(i).at(0));
				}
	
			else if(elementsOfEdgesLocal.at(i).at(1) == -1){
				elementImportMap.push_back(elementsOfEdgesGlobal.at(i).at(1));
			}

		}
	}
	sort(elementImportMap.begin(),elementImportMap.end());
	vec_GO_Type::iterator ip = unique( elementImportMap.begin(), elementImportMap.end());
	elementImportMap.resize(distance(elementImportMap.begin(), ip)); 
	Teuchos::ArrayView<GO> globalElementArrayImp = Teuchos::arrayViewFromVector( elementImportMap);

	// global Ids of Elements' Nodes
	MapPtr_Type mapElementImport =
		Teuchos::rcp( new Map_Type( elementMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalElementArrayImp, 0, this->comm_) );


	MultiVectorPtr_Type valuesU_x = Teuchos::rcp( new MultiVector_Type( elementMap, 1 ) );	
	MultiVectorPtr_Type valuesU_y = Teuchos::rcp( new MultiVector_Type( elementMap, 1 ) );	
	Teuchos::ArrayRCP< SC > entriesU_x  = valuesU_x->getDataNonConst(0);
	Teuchos::ArrayRCP< SC > entriesU_y  = valuesU_y->getDataNonConst(0);

	for(int l=0; l<elements->numberElements(); l++){
		entriesU_x[l]= u_El[l][0];
		entriesU_y[l]= u_El[l][1];
	}

	MultiVectorPtr_Type valuesUImp_x = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );
	MultiVectorPtr_Type valuesUImp_y = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );		
	valuesUImp_x->putScalar(0);
	valuesUImp_y->putScalar(0);
	
	valuesUImp_x->importFromVector(valuesU_x,true, "Insert");
	valuesUImp_y->importFromVector(valuesU_y, true, "Insert");

	
	Teuchos::ArrayRCP< SC > entriesUImp_x  = valuesUImp_x->getDataNonConst(0);
	Teuchos::ArrayRCP< SC > entriesUImp_y  = valuesUImp_y->getDataNonConst(0);

	

	// Then we determine the jump over the edges, if the element we need for this is not on our proc, we import the solution u
	for(int k=0; k< elements->numberElements();k++){
		edgeNumbers = edgeElements->getEdgesOfElement(k); // edges of Element k
		for(int i=0;i<3;i++){
			if(elementsOfEdgesGlobal.at(edgeNumbers[i]).size() > 1){  // not a boundary edge
				// Case that both necessary element are on the same Proc
				if(elementsOfEdgesLocal.at(edgeNumbers[i]).at(0) != -1   && elementsOfEdgesLocal.at(edgeNumbers[i]).at(1) != -1){
					elTmp1 = elementsOfEdgesLocal.at(edgeNumbers[i]).at(0);
					elTmp2 = elementsOfEdgesLocal.at(edgeNumbers[i]).at(1);

					u1 = u_El.at(elTmp1);
					u2 = u_El.at(elTmp2);

					}
				// if one of the necessary elements is not on our proc, we need to import the corresponding value of u
				else{
					vec_GO_Type elementInq(2); // the element in question is the one missing on our proc at entry 0
					elTmp1 = elementsOfEdgesGlobal.at(edgeNumbers[i]).at(0);
					elTmp2 = elementsOfEdgesGlobal.at(edgeNumbers[i]).at(1);

					if(elementMap->getLocalElement(elTmp1) !=  -1){
							elementInq[0]=elTmp2;
							elementInq[1]=elTmp1;
		
						}
					else {
							elementInq[0]=elTmp1;
							elementInq[1]=elTmp2;
						}
		
				
				u1[0] = u_El.at(elementMap->getLocalElement(elementInq[1])).at(0);
				u1[1] = u_El.at(elementMap->getLocalElement(elementInq[1])).at(1);

				u2[0] = entriesUImp_x[mapElementImport->getLocalElement(elementInq[0])];
				u2[1] = entriesUImp_y[mapElementImport->getLocalElement(elementInq[0])];
 
				}

			}


			// Normalenvektor bestimmen:
			p1_2[0] = points->at(edgeElements->getElement(edgeNumbers[i]).getNode(0)).at(0) - points->at(edgeElements->getElement(edgeNumbers[i]).getNode(1)).at(0);
			p1_2[1] = points->at(edgeElements->getElement(edgeNumbers[i]).getNode(0)).at(1) - points->at(edgeElements->getElement(edgeNumbers[i]).getNode(1)).at(1);
			v_E.at(0) = p1_2[1];
			v_E.at(1) = -p1_2[0];
			norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2));	    
				
			h_E =  sqrt(pow(p1_2[0],2)+pow(p1_2[1],2));

			errorEdgesInterior[i] =h_E*h_E* pow(( v_E[0]/norm_v_E* (u1[0]-u2[0]) + v_E[1]/norm_v_E*(u1[1]-u2[1])),2);

			u1[0]=0.;
			u1[1]=0.;
			u2[0]=0.;
			u2[1]=0.;
				
		}

		errorElement[k] = sqrt(1./2*(errorEdgesInterior[0]+errorEdgesInterior[1]+errorEdgesInterior[2])); 
		//cout << " Error global ELement k " << elementMap ->getGlobalElement(k) << " Fehler " << errorElement[k] << endl;
		errorEdgesInterior[0]=0.;
		errorEdgesInterior[1]=0.;
		errorEdgesInterior[2]=0.;

		if(maxErrorElLoc < errorElement[k] )
				maxErrorElLoc = errorElement[k];
		
	}
	// Refinement Strategies
	// As no P2 Error Estimation is implemented we use uniform Refinement in that case
	if(FEType == "P2" ){
		strategy = "Uniform";
	}

	// As we decide which element to tag based on the maximum error in the elements globally, we need to communicated this maxErrorEl
    reduceAll<int, double> (*this->comm_, REDUCE_MAX, maxErrorElLoc, outArg (maxErrorEl));
	int flagCount=0;
	if( strategy == "Maximum"){
		for(int k=0; k< elements->numberElements() ; k++){
			if( errorElement[k] > theta * maxErrorEl){
				elements->getElement(k).tagForRefinement();
				flagCount ++;
				}
			}
	}

	else if(strategy == "Doerfler"){
		double thetaSumTmp=0.;
		double thetaSum=0.;
		double muMax=0.;
		double sigmaT=0.;
		vec_bool_Type flagged(elements->numberElements());
		for(int k=0; k< elements->numberElements(); k++){
			thetaSumTmp = thetaSumTmp + pow(errorElement[k],2);
			flagged[k]=false;
		}
    	reduceAll<int, double> (*this->comm_, REDUCE_SUM, thetaSumTmp, outArg (thetaSum));
		while(sigmaT < theta*thetaSum){
			muMax=0.;
			for(int k=0; k< elements->numberElements(); k++){
				if(muMax < errorElement[k] && flagged[k]==false ){
					muMax = errorElement[k];
				}
			}
    		reduceAll<int, double> (*this->comm_, REDUCE_MAX, muMax, outArg (muMax));
			for(int k=0; k< elements->numberElements(); k++){
				if(muMax == errorElement[k] && flagged[k]==false ){
					flagged[k]=true;
					sigmaT = sigmaT + pow(errorElement[k],2);
					}
			}
    	reduceAll<int, double> (*this->comm_, REDUCE_MAX, sigmaT, outArg (sigmaT));
		}
			
	   	for(int k=0; k< elements->numberElements() ; k++){
			if( flagged[k]==true){
				elements->getElement(k).tagForRefinement();
				flagCount++;
				}
			}

	}	
		
	// Uniform Refinement
	else{ 
		 for(int k=0; k< elements->numberElements() ; k++){
				elements->getElement(k).tagForRefinement();
				flagCount++;
			}
	}
    reduceAll<int, int> (*this->comm_, REDUCE_MAX, flagCount , outArg (flagCount));

	auto finishError = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> elapsed_Error = finishError - startError;

	if(this->comm_->getRank() == 0){
		cout << "__________________________________________________________________________________________________________ " << endl;
		cout << " " << endl;
		cout << " The A-posteriori Error Estimation tagged " << flagCount << " Elements for adaptive Refinement with " << strategy << "-Strategy " << endl;
		cout << " The maximal Error of the Elements is " << maxErrorEl << endl;
		cout << " Elapsed Time of Error Estimation is " << elapsed_Error.count()  << " s" << endl;
		cout << "__________________________________________________________________________________________________________ " << endl;
		}
  	}
	else if(this->dim_ == 3){ 
		if(strategy != "Uniform" ){
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
	   		cout << " You selected a Refinement Strategy that is not yet implemented in the 3D Case." << endl;
			cout << " A uniform Refinement will be performed instead. " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			}
		
		for(int k=0; k< elements->numberElements() ; k++){
			elements->getElement(k).tagForRefinement();
			}

   }
	else
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Error Estimation is only available in 2 and 3 Dimensions");

	return errorElement;

}

// Mesh Refinement 
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::refineMesh( MeshUnstrRefPtrArray_Type meshesP1, int iteration, bool checkRestrictions, string restriction){

	MeshUnstrRefPtr_Type meshP1 = meshesP1[iteration];
		
	// not all Edges are assigned a Flag yet -> in Order to save a few steps along the way, we assign all Edges the correct flag now
	// In regular refinement the correct flags are then passed on, so this step is only necessary in Iteration 0
	if(iteration == 0 )
		meshP1->assignEdgeFlags();

	this->dim_ = meshP1->getDimension();
	this->FEType_="P1";
    this->volumeID_ = meshP1->volumeID_;
	this->rankRange_=meshP1->rankRange_;
	// necessary entities
    EdgeElementsPtr_Type edgeElements = meshP1->getEdgeElements(); // Edges
    ElementsPtr_Type elements = meshP1->getElementsC(); // Elements
    vec2D_dbl_ptr_Type points = meshP1->getPointsRepeated(); // Points
    this->elementMap_ = meshP1->elementMap_;
	this->mapUnique_ = meshP1->mapUnique_;
	this->mapRepeated_=meshP1->mapRepeated_;
	this->edgeMap_ = meshP1->edgeMap_;
	// entities for resulting Elements, Points, Edges, Flags
	// - Points, Elements and Flags are added to the existding entities
	// . edges are reset every refinement step
    this->edgeElements_.reset(new EdgeElements()); // Edges
    this->elementsC_.reset(new Elements(*elements));    // Elements
    this->pointsRep_.reset(new std::vector<std::vector<double> >(meshP1->pointsRep_->size(),vector<double>(this->dim_,-1.)));
    *this->pointsRep_ = *meshP1->pointsRep_; // Points
	this->pointsUni_.reset(new std::vector<std::vector<double> >( this->mapUnique_->getNodeNumElements(), vector<double>(this->dim_,-1. ) ) );
    *this->pointsUni_ = *meshP1->pointsUni_; 
	this->bcFlagUni_.reset( new std::vector<int> ( this->mapUnique_->getNodeNumElements(), 0 ) );
    *this->bcFlagUni_ = *meshP1->bcFlagUni_; 

	this->bcFlagRep_.reset(new vector<int>(meshP1->bcFlagRep_->size())); // Flags
	*this->bcFlagRep_ = *meshP1->bcFlagRep_;

	this->edgesElementOrder_=meshP1->getEdgeElementOrder();

	this->numElementsGlob_ = meshP1->numElementsGlob_; 

	

	if (this->dim_ == 2 || this->dim_ == 3){


		if(this->comm_->getRank() == 0){
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Start Iteration " << iteration+1 << " of "<< this->dim_ << "D Mesh Refinement " << endl;
			cout << " Number of Elements:	" << this->elementMap_->getGlobalNumElements() << endl;
			cout << " Number of Nodes:	" << this->mapUnique_->getGlobalNumElements() << endl; 
			cout << " Number of Edges:	" << this->edgeMap_->getGlobalNumElements() << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
		}


		// ------------------------------------------------------------------------------------------------------
		// Part I: Regular Refinement of Mesh
		// We refine the elements that were determined by our error estimation regular
		// ------------------------------------------------------------------------------------------------------	
		auto startPart1 = std::chrono::high_resolution_clock::now();
		const int myRank = this->comm_->getRank();
		// match the Edges to the Elements for being able to tag the edges of Elements in case of refinement -> does not work parallely as it is now; FIX: update as we refine?
		edgeElements->matchEdgesToElements(this->elementMap_);

		// Vector that carry the information which elements belong to an edge in local and global indices
		vec2D_GO_Type elementsOfEdgeGlobal = edgeElements->getElementsOfEdgeGlobal();
		vec2D_LO_Type elementsOfEdgeLocal = edgeElements->getElementsOfEdgeLocal();
		

		// Determine current global Interface IDs
		// As we pass on the interface boolean while refining red,blue,green... we theoretically dont need the whole extend of this function passed the first 
		// refinement, only the globalInterfaceIDs, which can be determined without 'elementsOfEdgeLocal' and 'Global'
		// (If elementsOfEdgeLocal and Global arent working, this can be taken out, yet is also nice for determining if it works correctly)
		vec_GO_Type globalInterfaceIDs;
		if(iteration==0)
			globalInterfaceIDs = edgeElements->determineInterfaceEdges(this->edgeMap_);
		else{
			for(int i=0; i< edgeElements->numberElements(); i++){
				if(edgeElements->getElement(i).isInterfaceElement()){
					globalInterfaceIDs.push_back(this->edgeMap_->getGlobalElement(i));
				}
			}	
			sort(globalInterfaceIDs.begin(), globalInterfaceIDs.end());			

		}

		// counting new Points and Elements:
		int newPoints=0; // total number of new Points
		int newPointsRepeated= 0; // number of new Points that share an other Processor
		int newPointsUnique =0; // number of Points that are uniquely on one Processor

		int newElements=0;	// new Elements on a Processor

		// Counting new Edges
		int newEdgesUnique=0;
		int newEdgesRepeated =0;
		

		// Loop through Elements in order to determine Elements to refine and refine them regular / red 
		int numPoints=0;

		// Depending on dimension we add a certain number of elements in the 2D case it's 3, in the 3D it's 7
		int dimRegRef=0;
		if(this->dim_ == 2)
			dimRegRef = 3;
		if(this->dim_ == 3)
			dimRegRef = 7;

		for(int i=0; i<elements->numberElements();i++){
			if(elements->getElement(i).isTaggedForRefinement()){
				numPoints= this->pointsRep_->size();
				this->refineRegular( meshP1, edgeElements, elements, i);
				newPoints=newPoints + this->pointsRep_->size()-numPoints;
				newElements = newElements +dimRegRef;
			}
		}

		auto finishPart1 = std::chrono::high_resolution_clock::now();
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part II: Communicating the tagged Edges 
		// As it is possible now, that edges were tagged on one processor on the interface between processers
		// we need communicate tagged edges across the interface
		// In order to safe time, we only communicate tagged interface Edges
		// InterfaceEdges can be determined by the vector 'elementsOfEdgesLocal' as the vector carries a -1 for 
		// for each element belonging to an edge that is not on the processor in question
		// ------------------------------------------------------------------------------------------------------
		auto startPart2 = std::chrono::high_resolution_clock::now();
		MapConstPtr_Type edgeMap = this->getEdgeMap();


		// Determine the globalInterfaceIDs of tagged edges
		vec_GO_Type globalInterfaceIDsTagged(0);
		GO indE;
		for(int i=0; i<globalInterfaceIDs.size(); i++){
			indE = edgeMap->getLocalElement(globalInterfaceIDs[i]);
			if(edgeElements->getElement(indE).isTaggedForRefinement()){
				globalInterfaceIDsTagged.push_back(globalInterfaceIDs[i]);
			}	
		}


		// Constructing a map of the global IDs of edges and tagged edges
		Teuchos::ArrayView<GO> globalEdgesInterfaceArray = Teuchos::arrayViewFromVector( globalInterfaceIDs);
		Teuchos::ArrayView<GO> globalEdgesInterfaceTaggedArray = Teuchos::arrayViewFromVector( globalInterfaceIDsTagged);

		MapPtr_Type mapInterfaceEdges =
			Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalEdgesInterfaceArray, 0, this->comm_) );
		MapPtr_Type mapInterfaceEdgesTagged =
			Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalEdgesInterfaceTaggedArray, 0, this->comm_) );

		// Multivector based on interfaceEdges Map with zero entries 
		MultiVectorGOPtr_Type taggedEdgesGlobal = Teuchos::rcp( new MultiVectorGO_Type(mapInterfaceEdges, 1 ) );
		taggedEdgesGlobal->putScalar(0);

		// Multivector based on taggesEdges Map with one as entries
		MultiVectorGOPtr_Type isActiveEdge = Teuchos::rcp( new MultiVectorGO_Type( mapInterfaceEdgesTagged, 1 ) );
		isActiveEdge->putScalar( (LO) 1);

		taggedEdgesGlobal->importFromVector(isActiveEdge, true, "Insert"); // From this we know that -> one signifies a tagged edge -> if we didn't tag it -> refine
		Teuchos::ArrayRCP< const GO >  tags = taggedEdgesGlobal->getData( 0 );


		// Adding Midpoints and tagging the edges that were tagged on other Procs that are on the interface
		LO ind;
		for (int i=0; i<tags.size(); i++) {
			if (tags[i] > 0){
				ind = edgeMap->getLocalElement(globalInterfaceIDs[i]);
				newPointsRepeated ++;
				if(!edgeElements->getElement(ind).isTaggedForRefinement()){
					edgeElements->getElement(ind).tagForRefinement();
					this->addMidpoint(meshP1,edgeElements,ind);
					// Collecting global IDs we didnt already considered
					globalInterfaceIDsTagged.push_back(globalInterfaceIDs[i]);
					newPoints ++;
					}		
			}
		}

		auto finishPart2 = std::chrono::high_resolution_clock::now();				
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part III: Checking the green Tags
		// Before we start refining the elements green, blue or red we can check whether the projected irregular 
		// refinement is fitting for the element
		// -> We only want to refine green, if the longest edge is the one beeing refined
		// -> We only want to refine blue, if the longest edge is among the tagged
		// ------------------------------------------------------------------------------------------------------

		auto startPart3 = std::chrono::high_resolution_clock::now();		
		if(checkRestrictions == true && this->dim_ ==2)	
			this->checkGreenTags(meshP1, elements ,edgeElements, iteration, newPoints, newPointsRepeated, globalInterfaceIDsTagged,mapInterfaceEdges, restriction );

		sort(globalInterfaceIDsTagged.begin(), globalInterfaceIDsTagged.end());
		auto finishPart3 = std::chrono::high_resolution_clock::now();	

		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part IV: Creating and distributing information of interface Edges
		// We now collected all the necessary information concerning the interface Edges
		// Additionally we need to distribute that information to all processors for updating the nodeMap, edgeMap and
		// elementsOfEdgesGlobal and elementsOfEdgesLocal later
		// All Processors need to know:
		// 		- global IDs of tagged interface edges
		// 		- global IDs of untagged interface edges
		// 		- global IDs of interface edges
		// the latter can be extracted from the two prior informations
		// ------------------------------------------------------------------------------------------------------
		
		auto startPart4 = std::chrono::high_resolution_clock::now();	
		// determine global interface IDs of untagged edges 	
		vec_GO_Type globalInterfaceIDsUntagged(0);
		for(int i=0; i<globalInterfaceIDs.size(); i++){
			indE = edgeMap->getLocalElement(globalInterfaceIDs[i]);
			if(!edgeElements->getElement(indE).isTaggedForRefinement()){
				globalInterfaceIDsUntagged.push_back(globalInterfaceIDs[i]);
			}	
		}
		// creating map for global interface IDs of untagged edges
		Teuchos::ArrayView<GO> globalEdgesInterfaceUntaggedArray = Teuchos::arrayViewFromVector( globalInterfaceIDsUntagged);
		MapPtr_Type mapInterfaceEdgesUntagged =
			Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalEdgesInterfaceUntaggedArray, 0, this->comm_) );

		// reset globalEdgesInterfaceTagged map with new information
		globalEdgesInterfaceTaggedArray = Teuchos::arrayViewFromVector( globalInterfaceIDsTagged);
		mapInterfaceEdgesTagged.reset(new Map<LO,GO,NO>( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalEdgesInterfaceTaggedArray, 0, this->comm_) );

		int maxRank = std::get<1>(this->rankRange_);
		// determine unique map
		MapPtr_Type mapInterfaceEdgesTaggedUnique = mapInterfaceEdgesTagged;
		MapPtr_Type mapInterfaceEdgesUntaggedUnique = mapInterfaceEdgesUntagged;

		//mapInterfaceEdgesTagged->print();
	
		if( mapInterfaceEdgesTagged->getGlobalNumElements()>0){
			mapInterfaceEdgesTaggedUnique = mapInterfaceEdgesTagged->buildUniqueMap( this->rankRange_ ); 
		}
		if(mapInterfaceEdgesUntagged->getGlobalNumElements()>0){
			mapInterfaceEdgesUntaggedUnique = mapInterfaceEdgesUntagged->buildUniqueMap( this->rankRange_ ); 
		}

		// Creating map of all procs on each processor, this map is usefull for distributing the information of new nodes, elements etc.
		// each processor will then write its own information in a multivector that is based in the localProc-Map (entry is the own proc number)
		vec_GO_Type globalProcs(0);
		for (int i=0; i<= maxRank; i++)
			globalProcs.push_back(i);

		Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

		vec_GO_Type localProc(0);
		localProc.push_back(this->comm_->getRank());
		Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

		MapPtr_Type mapGlobalProc =
			Teuchos::rcp( new Map_Type( meshP1->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, this->comm_) );

		MapPtr_Type mapProc =
			Teuchos::rcp( new Map_Type( meshP1->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, this->comm_) );
		
		
		// Vectors of global IDs of tagged and untagged edges, that is the same on each processor
		vec_GO_Type edgesInterfaceUntagged(0);
		vec_GO_Type edgesInterfaceTagged(0);
	
		// Function that distributes the local information of tagged and untagged edges to all processors
		this->distributeTaggedAndUntaggedEdges(mapProc, mapGlobalProc,mapInterfaceEdgesTaggedUnique,mapInterfaceEdgesUntaggedUnique, edgesInterfaceUntagged, edgesInterfaceTagged);
		
		// From that information we can construct the global IDs of all Interface Edges
		vec_GO_Type edgesInterface(0);
		sort(edgesInterfaceUntagged.begin(), edgesInterfaceUntagged.end());
		sort(edgesInterfaceTagged.begin(), edgesInterfaceTagged.end());
		edgesInterface.insert( edgesInterface.end(), edgesInterfaceUntagged.begin(), edgesInterfaceUntagged.end() );
		edgesInterface.insert( edgesInterface.end(), edgesInterfaceTagged.begin(), edgesInterfaceTagged.end() );
		sort(edgesInterface.begin(), edgesInterface.end());


		/*for (int i=0; i< edgesInterfaceUntagged.size(); i++)
			cout << " Untagged " << i << " " << edgesInterfaceUntagged[i] << endl;*/

		/*for (int i=0; i< edgesInterfaceTagged.size(); i++)
			cout << " Tagged " << i << " " << edgesInterfaceTagged[i] << endl;*/

		/*for (int i=0; i< edgesInterface.size(); i++)
			cout << " Interface " << i << " " << edgesInterface[i] << endl;*/
		

		auto finishPart4 = std::chrono::high_resolution_clock::now();		
		// ------------------------------------------------------------------------------------------------------


		// ------------------------------------------------------------------------------------------------------		
		// Part V: Communicating the added Points and updating the corresponding Maps
		// generally we keep the nodelist of throughout the refinement steps and only add new points to it
		// consequently we update the maps for those new points depending on nodes on interface and points unqiuely
		// on processors
		// ------------------------------------------------------------------------------------------------------
		auto startPart5 = std::chrono::high_resolution_clock::now();
		// Collecting the number of new Points and repeated new Points in total on each Proc
		int globalCountPoints = 0; 
		int globalCountRepeated = 0;
		reduceAll<int, int> (*this->comm_, REDUCE_SUM, newPoints, outArg (globalCountPoints));
		reduceAll<int, int> (*this->comm_, REDUCE_SUM, newPointsRepeated, outArg (globalCountRepeated));
		// Setting Unique newPoints as value for communication
		MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );

		// Unique Points on each Proc
		newPointsUnique= newPoints - newPointsRepeated;
		exportLocalEntry->putScalar( (LO) newPointsUnique);

		MultiVectorLOPtr_Type isActiveNodeUnique= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		isActiveNodeUnique->putScalar( (LO) 0 ); 
		isActiveNodeUnique->importFromVector( exportLocalEntry, true, "Insert");
		// -> Result: Information on uniquely added Points of each Proc is accessible to all

		Teuchos::ArrayRCP<const LO> nodesUniqueProcList = isActiveNodeUnique-> getData(0);
			
		// Step 1: Add Nodes that are repeated on the processors
		// Generally we allocate the first 'n' Points to Proc 0, the following start with (n+1) on Proc 1 and so on -> ProcOffset

		GO refineOffset = meshP1->mapUnique_->getMaxAllGlobalIndex()+1; // Number of Nodes before Refinement
		int procOffset=0;
		for(int i=0; i< myRank; i++)
			procOffset = procOffset + nodesUniqueProcList[i];

		MapConstPtr_Type mapRepeatedP1 = meshP1->getMapRepeated(); // Maps
		Teuchos::ArrayView<const GO> nodeList = meshP1->getMapRepeated()->getNodeElementList();
		vec_GO_Type vecGlobalIDsOld = Teuchos::createVector( nodeList );
	  
		vec_GO_Type vecGlobalIDs(this->pointsRep_->size(),-1);

		for (int i=0; i<vecGlobalIDsOld.size(); i++){
			vecGlobalIDs[i] = vecGlobalIDsOld[ i] ; // (i + refineOffset+procOffset);
		}

		// Add Nodes that are on the interface -> repeated Points
		int sumPointsUnique = globalCountPoints - globalCountRepeated;

		// loop through the edgesInterfaceTagged vector
		// if the processor carries that tagged interface edge we assign index determined by count, refineOffse and sumPointsUnque

		// !!!! its is important to assign the global ID to the right entry in vecGlobalIDs !!!!
		// that entry can be determined with the edgeMap and midPoint entry, as it is the local ID of the point in question

		int count=0; // Per tagged edge the count goes up, consequently every Proc that shares an edge assigns the same global index
		LO indLO;
		for (int i=0; i<edgesInterfaceTagged.size(); i++){
			if(edgeMap->getLocalElement(edgesInterfaceTagged[i]) != -1 ){
				indLO = edgeMap->getLocalElement(edgesInterfaceTagged[i]);
				vecGlobalIDs[edgeElements->getMidpoint(indLO)] =  refineOffset+sumPointsUnique+count;
			}
			count ++; 	
		}
		// Step 2: Add points that are uniquely on each processor, leftover -1 entries in vecGlobalIDs signal a missing unique entry
		// all left over entries 
		count =0;
	  	for (int i= vecGlobalIDsOld.size(); i < this->pointsRep_->size() ; i++){
			if(vecGlobalIDs[i] == -1){
				vecGlobalIDs[i] = count + refineOffset+procOffset;
				count ++;
			}
		}
		Teuchos::RCP<std::vector<GO> > pointsRepGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDs ) );
		Teuchos::ArrayView<GO> pointsRepGlobMappingArray = Teuchos::arrayViewFromVector( *pointsRepGlobMapping );
		
		this->mapRepeated_.reset(new Map<LO,GO,NO>( meshP1->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), pointsRepGlobMappingArray, 0, this->comm_) );
		this->mapUnique_ = this->mapRepeated_->buildUniqueMap( this->rankRange_ );
		
		//this->mapRepeated_->print();	
	
		// Points and Flags Unique
		this->pointsUni_.reset(new std::vector<std::vector<double> >( this->mapUnique_->getNodeNumElements(), vector<double>(this->dim_,-1. ) ) );
		this->bcFlagUni_.reset( new std::vector<int> ( this->mapUnique_->getNodeNumElements(), 0 ) );
		for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
			GO gid = this->mapUnique_->getGlobalElement( i );
			LO id = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement( i ) );
			this->pointsUni_->at(i) = this->pointsRep_->at(id);
			this->bcFlagUni_->at(i) = this->bcFlagRep_->at(id);
		}
		auto finishPart5 = std::chrono::high_resolution_clock::now();

		// ------------------------------------------------------------------------------------------------------


		// ------------------------------------------------------------------------------------------------------		
		// Part VI: Irregular Refinement
		// Now all Elements, that were supposed to be refined regularly are refined and all neighbouring elements
		// were checked for certain restrictions (see III)
		// Now we have to determine how neighbouring elements are supposed to be refined	
		// -> One Edge of Element 'i' is tagged for refinement -> greenRefinement
		// -> Two Edges of Element 'i' are tagged for refinement -> blueRefinement
		// -> Three Edges of Element 'i' are tagged for refinement -> redRefinement / regular Refinement		
		// ------------------------------------------------------------------------------------------------------
		auto startPart6 = std::chrono::high_resolution_clock::now();
		vec_int_Type tagCounter2(elements->numberElements());
		int edgeNum;
		int edgePerElement = this->dim_+1;
		for(int i=0;i<elements->numberElements() ;i++){
			if(elements->getElement(i).isTaggedForRefinement()==false){ // only looking at the untagged Elements
				for(int j=0;j<edgePerElement;j++){
					edgeNum = edgeElements->getEdgesOfElement(i).at(j);
					if(edgeElements->getElement(edgeNum).isTaggedForRefinement())
						tagCounter2[i] ++;
				}
				if(tagCounter2[i] == 1){
					elements->getElement(i).setFiniteElementRefinementType("green");
					this->refineGreen(meshP1, edgeElements, elements, i);
					newElements ++;
				}
				else if(tagCounter2[i] == 2){
					elements->getElement(i).setFiniteElementRefinementType("blue");
					this->refineBlue(meshP1,edgeElements, elements, i);
					newElements= newElements+2;
					}	
				else if(tagCounter2[i] == 3){
					elements->getElement(i).setFiniteElementRefinementType("red");
					this->refineRed(meshP1, edgeElements, elements, i);
					newElements= newElements+3;
					}
				else {
					// the element in question has not been refined, nor its neighbour -> add edges to edge list
					vec_int_Type edgesOfElement = edgeElements->getEdgesOfElement(i);
					for(int j=0;j<3;j++){
						edgeElements->getElement(edgesOfElement[j]).setFiniteElementRefinementType("unrefined");
						edgeElements->getElement(edgesOfElement[j]).setPredecessorElement(edgeMap->getGlobalElement(edgesOfElement[j]));
						this->edgeElements_->addEdge(edgeElements->getElement(edgesOfElement[j]),i);
						
						if(edgeElements->getElement(edgesOfElement[j]).getFlag() !=0 && edgeElements->getElement(edgesOfElement[j]).getFlag() !=10){
							if ( !this->elementsC_->getElement(i).subElementsInitialized() )
								this->elementsC_->getElement(i).initializeSubElements( this->FEType_, this->dim_ -1) ;
							this->elementsC_->getElement(i).addSubElement(edgeElements->getElement(edgesOfElement[j]));
						}

					}	
					if(this->dim_ == 3){ 
   						TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "3D Algorithm for irregular MeshRefinement is currently not available, please choose uniform Refinement");
					}

				}

			}
		}
	
		// Checking if a tagged element survived
		for(int j=0;j<this->elementsC_->numberElements();j++){
			this->elementsC_->getElement(j).untagForRefinement();
			if(this->elementsC_->getElement(j).isTaggedForRefinement() == true)
				cout<< "tagged element survived " << endl;
		}

		auto finishPart6 = std::chrono::high_resolution_clock::now();
		// ------------------------------------------------------------------------------------------------------
		// Part VII: Updating the Element Map
		// update elementMap by distributing information of locally added elements and add global IDs based on
		// that information
		// the list of elements only extends throughout the refinement steps, so we only neew to allocate the 
		// globalIDs for elements that extend our last element number
		// ------------------------------------------------------------------------------------------------------
		auto startPart7 = std::chrono::high_resolution_clock::now();

		MapConstPtr_Type elementMap = meshP1->getElementMap();
		// information of new elements
		exportLocalEntry->putScalar( (LO) newElements );

		// map equal to the original edgeMap with zero entries
		MultiVectorLOPtr_Type isActiveElement= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		isActiveElement->putScalar( (LO) 0 ); 
		isActiveElement->importFromVector( exportLocalEntry, false, "Insert");

		// number of new elements per Processor
		Teuchos::ArrayRCP< const LO > elementsRank = isActiveElement->getData(0);

		Teuchos::ArrayView<const GO> elementList = elementMap->getNodeElementList();
		std::vector<GO> vecGlobalIDsElements = Teuchos::createVector( elementList );

		// element offset before refinement
		int offsetElements = elementMap->getGlobalNumElements();

		// determine the offset for this processor
		GO procOffsetElements=0;
		for(int i=0; i< myRank; i++)
			procOffsetElements = procOffsetElements + elementsRank[i];

		for (int i=0; i<newElements; i++){
			vecGlobalIDsElements.push_back( i +  offsetElements + procOffsetElements);
		}

		Teuchos::RCP<std::vector<GO> > elementsGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDsElements ) );
		Teuchos::ArrayView<GO> elementsGlobMappingArray = Teuchos::arrayViewFromVector( *elementsGlobMapping);

		this->elementMap_.reset(new Map<LO,GO,NO>(meshP1->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), elementsGlobMappingArray, 0, this->comm_) );

		// determine global number of elements
		this->numElementsGlob_ = this->elementMap_->getMaxAllGlobalIndex()+1; //meshP1->numElementsGlob_; 

		auto finishPart7 = std::chrono::high_resolution_clock::now();
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part VIII: Making the Edges unique, set global IDs and set elements edges
		// we need to make the edge list unique as there are redundant edges
		// furthermore we set up the elementsOfEdgeLocal and elementsOfEdgeGlobal, but the information of other
		// procs is still missing there (this will be finalized in Part X)
		// ------------------------------------------------------------------------------------------------------
		auto startPart8 = std::chrono::high_resolution_clock::now();
		vec2D_GO_Type combinedEdgeElements;
		
		this->edgeElements_->sortUniqueAndSetGlobalIDsParallel(this->elementMap_,combinedEdgeElements);
		auto finishPart8 = std::chrono::high_resolution_clock::now();
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part IX: Updating the EdgeMap 
		// the edges will be newly constructed each refinement step, consequently so will the entire edgeMap
		// the general idea is the same as in previous map updates
		// we first add the unique edges on each processor and determine the global IDs by the usual procOffset
		// determined by the number of unique edges on each proc
		// then we add the globalIDs of untagged edges and tagged interface edges
		// ------------------------------------------------------------------------------------------------------
		auto startPart9 = std::chrono::high_resolution_clock::now();
		auto startPart9a = std::chrono::high_resolution_clock::now();
		// determine the local IDs of new interface Edges an their global node Indices
		vec_GO_Type globalNodeIDsNewInterfaceEdges(0);
		vec2D_int_Type newInterfaceEdgesLocalId = determineNewInterfaceEdgesLocalIds(globalNodeIDsNewInterfaceEdges) ;
		int newInterfaceEdges = newInterfaceEdgesLocalId.size();
		// determine number of new (repeated) interface edges
		newEdgesRepeated= 2*globalInterfaceIDsTagged.size() + globalInterfaceIDsUntagged.size() + newInterfaceEdges  ; // 2* edgesInterfaceTagged.size()+ edgesInterfaceUntagged.size();
	
		// number of new unique edes
		newEdgesUnique = this->edgeElements_->numberElements() - newEdgesRepeated;

		// distribute newEdgesRepeated number 
		exportLocalEntry->putScalar( (LO) newEdgesRepeated );

		MultiVectorLOPtr_Type newEdgesRepeatedProcGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newEdgesRepeatedProcGlobal->putScalar( (LO) 0 ); 
		newEdgesRepeatedProcGlobal->importFromVector( exportLocalEntry, true, "Insert");

		// distribute newEdgesUnique number 
		exportLocalEntry->putScalar( (LO) newEdgesUnique );

		// map equal to the original edgeMap with zero entries
		MultiVectorLOPtr_Type newEdgesUniqueProcGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newEdgesUniqueProcGlobal->putScalar( (LO) 0 ); 
		newEdgesUniqueProcGlobal->importFromVector( exportLocalEntry, true, "Insert");

		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > newEdgesRepeatedList = newEdgesRepeatedProcGlobal->getData(0);
		Teuchos::ArrayRCP< const LO > newEdgesUniqueList = newEdgesUniqueProcGlobal->getData(0);

		GO procOffsetEdgesUnique=0;
		for(int i=0; i< myRank; i++)
			procOffsetEdgesUnique= procOffsetEdgesUnique + newEdgesUniqueList[i];

		int procOffsetEdgesUniqueSum=0;
		for(int i=0; i< maxRank+1; i++)
			procOffsetEdgesUniqueSum= procOffsetEdgesUniqueSum + newEdgesUniqueList[i];

		// global IDs for map
		vec_GO_Type vecGlobalIDsEdges(this->edgeElements_->numberElements());
		// determine the local IDs of the global 'previous' interface IDs
		vec2D_int_Type interfaceEdgesLocalId = determineInterfaceEdgesNewLocalIds(edgeElements, globalInterfaceIDs,edgeMap) ;	

		// Step 1: adding unique global edge IDs
		count=0;
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(!this->edgeElements_->getElement(i).isInterfaceElement()){
				vecGlobalIDsEdges.at(i) = procOffsetEdgesUnique+count;
				count++;
			}
		}
		auto finishPart9a = std::chrono::high_resolution_clock::now();
		auto startPart9b = std::chrono::high_resolution_clock::now();
		// Step 2: adding repeated global edge IDs
		// !!! important to allocate the right IDs to corresponding edges, as tagged edges split in two !!
		this->addRepeatedIDsEdges(vecGlobalIDsEdges, edgesInterfaceUntagged, edgesInterfaceTagged,procOffsetEdgesUniqueSum,interfaceEdgesLocalId);

		Teuchos::RCP<std::vector<GO> > edgesGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDsEdges ) );
		Teuchos::ArrayView<GO> edgesGlobMappingArray = Teuchos::arrayViewFromVector( *edgesGlobMapping);

		this->edgeMap_.reset(new Map<LO,GO,NO>(meshP1->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingArray, 0, this->comm_) );
		//this->edgeMap_->print();
		auto finishPart9b = std::chrono::high_resolution_clock::now();
		// In the 3D Refinement Case we need a third step for updating the edgeMap
		// We might add Edges on the interface, that haven't previously existed 

		// Step 3: adding repeated new global edges that appear in 3D case on interface surfaces
		// In order to correctly set the global IDs we communicate a list of all those new interface edges among the processors
		// With that list a find-operation determines which global ID your own edge in question has
		// (i.e. Processor 0 holds edge [x_0,x_3] and that edge as is at entry 5 in the list, global ID is offset+5, all Procs have the same list and set the same entry for that edge)
		auto startPart9c = std::chrono::high_resolution_clock::now();
		// Communicate the amount of new edges on each Processor 
		exportLocalEntry->putScalar( (LO) newInterfaceEdges );

		MultiVectorLOPtr_Type newInterfaceEdgesProcGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newInterfaceEdgesProcGlobal->putScalar( (LO) 0 ); 
		newInterfaceEdgesProcGlobal->importFromVector( exportLocalEntry, true, "Insert");
		//newInterfaceEdgesProcGlobal->print();

		Teuchos::ArrayRCP<const LO> newInterfaceEdgesGlobalNum = newInterfaceEdgesProcGlobal-> getData(0);

		// Offset
		procOffset=0;
		for(int i=0; i< myRank; i++)
			procOffset = procOffset + newInterfaceEdgesGlobalNum[i];
	
		// Local Map of edges
		vec_GO_Type localNewInterfaceVec(newInterfaceEdges);
		for(int i=0 ; i< newInterfaceEdges; i++)
			localNewInterfaceVec[i]=i + procOffset;

		Teuchos::ArrayView<GO> localNewInterfaceArray = Teuchos::arrayViewFromVector( localNewInterfaceVec);
		MapPtr_Type mapLocalNewInterfaceEdges =
			Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localNewInterfaceArray, 0, this->comm_) );

		// Global Map of Edges
		int sumRepeatedNewInterface=0;
		for(int i=0; i<= maxRank; i++)
			sumRepeatedNewInterface = sumRepeatedNewInterface + newInterfaceEdgesGlobalNum[i];	


		vec_GO_Type globalNewInterfaceVec(sumRepeatedNewInterface);
		for(int i=0; i<sumRepeatedNewInterface; i++)
			globalNewInterfaceVec[i]=i;

		Teuchos::ArrayView<GO> globalNewInterfaceArray = Teuchos::arrayViewFromVector( globalNewInterfaceVec);
		MapPtr_Type mapGlobalNewInterfaceEdges =
			Teuchos::rcp( new Map_Type( edgeMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalNewInterfaceArray, 0, this->comm_) );


		// Export and Import Vectors of new repeated Edges
		MultiVectorLOPtr_Type exportAllGlobalNewInterfaceEdges_1= Teuchos::rcp( new MultiVectorLO_Type( mapLocalNewInterfaceEdges, 1 ) );
		MultiVectorLOPtr_Type exportAllGlobalNewInterfaceEdges_2= Teuchos::rcp( new MultiVectorLO_Type( mapLocalNewInterfaceEdges, 1 ) );

		MultiVectorLOPtr_Type importAllGlobalNewInterfaceEdges_1= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalNewInterfaceEdges, 1 ) );
		MultiVectorLOPtr_Type importAllGlobalNewInterfaceEdges_2= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalNewInterfaceEdges, 1 ) );

		Teuchos::ArrayRCP<LO> exportAllGlobalNewInterfaceEdges_1_entries = exportAllGlobalNewInterfaceEdges_1-> getDataNonConst(0);
		Teuchos::ArrayRCP<LO> exportAllGlobalNewInterfaceEdges_2_entries = exportAllGlobalNewInterfaceEdges_2-> getDataNonConst(0);


		vec2D_int_Type newInterfaceEdgesLocalIdTmp = newInterfaceEdgesLocalId;	
		for(int i=0; i< newInterfaceEdgesLocalId.size() ; i++){
			exportAllGlobalNewInterfaceEdges_1_entries[i] = newInterfaceEdgesLocalIdTmp[i][1]; 
			exportAllGlobalNewInterfaceEdges_2_entries[i] = newInterfaceEdgesLocalIdTmp[i][2];
		}
		
		importAllGlobalNewInterfaceEdges_1->importFromVector(exportAllGlobalNewInterfaceEdges_1,false,"Insert");
		importAllGlobalNewInterfaceEdges_2->importFromVector(exportAllGlobalNewInterfaceEdges_2,false,"Insert");

		Teuchos::ArrayRCP< const LO >  importAllGlobalNewInterfaceEdges_1_entries = importAllGlobalNewInterfaceEdges_1->getData(0); 
		Teuchos::ArrayRCP< const LO >  importAllGlobalNewInterfaceEdges_2_entries = importAllGlobalNewInterfaceEdges_2->getData(0); 

		// Vector that contains all global node IDs of new interface Edges
		vec2D_GO_Type newGlobalInterfaceEdges(sumRepeatedNewInterface);

		for(int i=0; i<sumRepeatedNewInterface; i++){
			 newGlobalInterfaceEdges[i] = {importAllGlobalNewInterfaceEdges_1_entries[i] ,importAllGlobalNewInterfaceEdges_2_entries[i] } ;
			 sort(newGlobalInterfaceEdges[i].begin(),newGlobalInterfaceEdges[i].end());
		}

		sort( newGlobalInterfaceEdges.begin(), newGlobalInterfaceEdges.end() );
		newGlobalInterfaceEdges.erase( unique( newGlobalInterfaceEdges.begin(), newGlobalInterfaceEdges.end() ), newGlobalInterfaceEdges.end() );

		// -> We now have distributed the list of new edges on the interface globally and can determine their ids accordingly
		int offsetNewRep = this->edgeMap_->getMaxAllGlobalIndex()+1;
		LO entry;
		vec_GO_Type edgeTmp;
		for(int i=0 ; i<newInterfaceEdgesLocalId.size(); i++){
			edgeTmp = {newInterfaceEdgesLocalId[i][1],newInterfaceEdgesLocalId[i][2] } ;
			auto it1 = find( newGlobalInterfaceEdges.begin(), newGlobalInterfaceEdges.end() ,edgeTmp);
            entry = distance( newGlobalInterfaceEdges.begin() , it1 );
			vecGlobalIDsEdges[newInterfaceEdgesLocalId[i][0]] = offsetNewRep + entry;
			
		}
		edgesGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDsEdges ) );
		edgesGlobMappingArray = Teuchos::arrayViewFromVector( *edgesGlobMapping);

		this->edgeMap_.reset(new Map<LO,GO,NO>(meshP1->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingArray, 0, this->comm_) );
		//this->edgeMap_->print();
		auto finishPart9c = std::chrono::high_resolution_clock::now();
		auto finishPart9 = std::chrono::high_resolution_clock::now();

		// ------------------------------------------------------------------------------------------------------
		// Part X: Updating elementsOfEdgeGlobal and elementsOfEdgeLocal
		// this step is only necessary if we have more than 1 processor, as the edgeElements function work serially
		// we started the setup before (sortUniqueAndSetGlobalIDsParallel) an now finalize it with the information of other processors
		// the edges on the interface need the global element number of the neighbouring processor
		// ------------------------------------------------------------------------------------------------------
		auto startPart10 = std::chrono::high_resolution_clock::now();

		this->edgeElements_->setElementsEdges( combinedEdgeElements );

		this->edgeElements_->setUpElementsOfEdge( this->elementMap_, this->edgeMap_);
		// the update only works correclty in the 2D case. For P1 Elements it is not necessary to set elementsOfEdgeGlobal and Local correctly
		// For P1 Elements 3D uniform refinement works even without this beeing set correctly. Will be fixed later
		this->updateElementsOfEdgesLocalAndGlobal(interfaceEdgesLocalId, maxRank, edgesInterface, edgeMap);
	

		auto finishPart10 = std::chrono::high_resolution_clock::now();
		// ------------------------------------------------------------------------------------------------------
		std::chrono::duration<double> elapsed_Part1 = finishPart1 - startPart1;
		std::chrono::duration<double> elapsed_Part2 = finishPart2 - startPart2;
		std::chrono::duration<double> elapsed_Part3 = finishPart3 - startPart3;
		std::chrono::duration<double> elapsed_Part4 = finishPart4 - startPart4;
		std::chrono::duration<double> elapsed_Part5 = finishPart5 - startPart5;
		std::chrono::duration<double> elapsed_Part6 = finishPart6 - startPart6;
		std::chrono::duration<double> elapsed_Part7 = finishPart7 - startPart7;
		std::chrono::duration<double> elapsed_Part8 = finishPart8 - startPart8;
		std::chrono::duration<double> elapsed_Part9 = finishPart9 - startPart9;
		std::chrono::duration<double> elapsed_Part9a = finishPart9a - startPart9a;
		std::chrono::duration<double> elapsed_Part9b = finishPart9b - startPart9b;
		std::chrono::duration<double> elapsed_Part9c = finishPart9c - startPart9c;
		std::chrono::duration<double> elapsed_Part10 = finishPart10 - startPart10;

		std::chrono::duration<double> elapsed_Mesh = finishPart10 - startPart1;

		if(this->comm_->getRank() == 0){
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Elapsed Time of different Refinement steps " << endl;
			cout <<" Part 1 - Regular Refinement:							" << elapsed_Part1.count() << " s" << endl;
			cout <<" Part 2 - Communicating tagged interface Edges:					" << elapsed_Part2.count() << " s "<< endl;
			cout <<" Part 3 - Checking Restrictions:						" << elapsed_Part3.count() << " s "<< endl;
			cout <<" Part 4 - Creating and distributing information of interface Edges:		" << elapsed_Part4.count() <<  " s" << endl;
			cout <<" Part 5 - Communicating added Points:						" << elapsed_Part5.count() << " s " << endl;
			cout <<" Part 6 - Irregular Refinement:							" << elapsed_Part6.count() << " s "<< endl;
			cout <<" Part 7 - Updating Element Map:							" << elapsed_Part7.count() << " s" << endl;
			cout <<" Part 8 - Making the Edges unique and set global IDs:				" << elapsed_Part8.count() << " s " << endl;
			cout <<" Part 9 - Updating edgeMap:							" << elapsed_Part9.count() <<  " s" << endl;
			cout <<" ....Part 9a - setting unique IDs 						" << elapsed_Part9a.count() <<  " s" << endl;
			cout <<" ....Part 9b - setting repeated IDs 						" << elapsed_Part9b.count() <<  " s" << endl;
			cout <<" ....Part 9c - setting new repeated ID						" << elapsed_Part9c.count() <<  " s" << endl;
			cout <<" Part 10 -Updating elementsOfEdgesGlobal/Local:					" << elapsed_Part10.count() <<  " s" << endl;
			cout <<" Total elapsed Time of Refinement						" << elapsed_Mesh.count() <<  " s" << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
		}

		if(this->comm_->getRank() == 0){
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " ... finished Iteration " << iteration+1 << " of 2D Mesh Refinement " << endl;
			cout << " Number of new Elements:	 " << this->elementMap_->getGlobalNumElements() - meshP1->elementMap_-> getGlobalNumElements() << endl;
			cout << " Number of new Nodes:		 " << this->mapUnique_->getGlobalNumElements()- meshP1->mapUnique_-> getGlobalNumElements() << endl; 
			cout << " Number of new Edges:	 	 " << this->edgeMap_->getGlobalNumElements()- meshP1->edgeMap_-> getGlobalNumElements() << endl;
			cout << " Elapsed Time of Refinement:	 " << elapsed_Mesh.count()  << " s" << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
		}

		
	}
	else
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Mesh Refinement is only available in 2 and 3 Dimensions");
		

}


// ------------------------------------------------------------------------------------------------------
// 2D Refinement Restrictions
// ------------------------------------------------------------------------------------------------------
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::checkGreenTags(MeshUnstrRefPtr_Type meshP1, ElementsPtr_Type elements ,EdgeElementsPtr_Type edgeElements, int iteration, int& newPoints, int& newPointsCommon, vec_GO_Type& globalInterfaceIDsTagged, MapConstPtr_Type mapInterfaceEdges, string restriction){

	// We determine whether a element that is tagged for green refinement has been refined green in the previous refinement
	// If thats is the case, we choose to refine element blue instead by adding a node to the longest edge
	vec_int_Type tagCounter(elements->numberElements());
    vec2D_dbl_ptr_Type points = meshP1->getPointsRepeated(); // Points
	int edgeNum;

	MapConstPtr_Type edgeMap = meshP1->getEdgeMap();

	// Determine Global Ids of Interface Edges
	vec2D_GO_Type elementsOfEdgesGlobal = edgeElements->getElementsOfEdgeGlobal();
	vec2D_LO_Type elementsOfEdgesLocal = edgeElements->getElementsOfEdgeLocal();
	
	int layer=1;//Indicator for our greenCheckSearch-> as we add Points for blue refinement, new to be green refined elements might accur-> nextGenGreen Elements-> how many Generations do we check?
	int entry;


	for(int i=0;i<elements->numberElements() ;i++){
		for(int j=0;j<3;j++){
			edgeNum = edgeElements->getEdgesOfElement(i).at(j);
			if(edgeElements->getElement(edgeNum).isTaggedForRefinement()){
				tagCounter[i]=tagCounter[i]+1; // We count the tags of element -> one Tag equals green Refinement
			}
		}
	}

	if(restriction != "keepRegularity" && restriction != "checkGreenTags"){
		cout << " !!! The restriction Type you requested is not available, 'keepRegularity' will be performed instead !!! " << endl; 
		restriction = "keepRegularity";
	}

	while(layer >0){
		vec_GO_Type globalEdges(0);
		layer=0;
		if(restriction == "keepRegularity"){
			for(int i=0;i<elements->numberElements() ;i++){
				if(elements->getElement(i).isTaggedForRefinement()==false){ // only looking at the untagged Elements
					if(tagCounter[i]==1){
						entry = this->determineLongestEdge(edgeElements,edgeElements->getEdgesOfElement(i),points); // we determine the edge, we would choose for blue Refinement
						if(!edgeElements->getElement(entry).isTaggedForRefinement()){ // If the longestest edge is already the tagged one, we leave the Element alone
							edgeElements->getElement(entry).tagForRefinement(); // we tag the Element for refinement
							this->addMidpoint(meshP1,edgeElements,entry);	// we add the necessary midpoint
							newPoints ++; // add new Points 
							if(elementsOfEdgesLocal.at(entry).size() > 1){
								if(elementsOfEdgesLocal.at(entry).at(0) != -1)
									tagCounter[elementsOfEdgesLocal.at(entry).at(0)] ++;
								else
									globalEdges.push_back(edgeMap->getGlobalElement(entry));  // and add the index for communication later
								if(elementsOfEdgesLocal.at(entry).at(1) != -1)
									tagCounter[elementsOfEdgesLocal.at(entry).at(1)] ++;
								else
									globalEdges.push_back(edgeMap->getGlobalElement(entry));  // and add the index for communication later
							}
							else
								tagCounter[elementsOfEdgesLocal.at(entry).at(0)] ++;

							tagCounter[i] = -1; // This element now has been refined optimally, so no other checks are necessary
							layer++;
						}
					}
					if(tagCounter[i]==2){
						entry = this->determineLongestEdge(edgeElements,edgeElements->getEdgesOfElement(i),points); // we determine the edge, we would choose for blue Refinement
						if(!edgeElements->getElement(entry).isTaggedForRefinement()){ // If the longestest edge is already the tagged one, we leave the Element alone
							//cout <<" Change to red in k= " << k << endl;
							edgeElements->getElement(entry).tagForRefinement(); // we tag the Element for refinement
							this->addMidpoint(meshP1,edgeElements,entry);	// we add the necessary midpoint
							newPoints ++; // add new Points 
							if(elementsOfEdgesLocal.at(entry).size() > 1){
								if(elementsOfEdgesLocal.at(entry).at(0) != -1)
									tagCounter[elementsOfEdgesLocal.at(entry).at(0)] ++;
								else
									globalEdges.push_back(edgeMap->getGlobalElement(entry));  // and add the index for communication later

								if(elementsOfEdgesLocal.at(entry).at(1) != -1)
									tagCounter[elementsOfEdgesLocal.at(entry).at(1)] ++;
								else
									globalEdges.push_back(edgeMap->getGlobalElement(entry));  // and add the index for communication later

							}
							else
								tagCounter[elementsOfEdgesLocal.at(entry).at(0)] ++;

							tagCounter[i] = -1; // This element now has been refined optimally, so no other checks are necessary
							layer++;
						}
					}
				}
			}
		}

		else if (restriction == "checkGreenTags"){
			for(int i=0;i<elements->numberElements() ;i++){
			
				if(elements->getElement(i).isTaggedForRefinement()==false){ // only looking at the untagged Elements
				
					if(tagCounter[i]==1 && 	elements->getElement(i).getFiniteElementRefinementType() == "green"){
						entry = this->determineLongestEdge(edgeElements,edgeElements->getEdgesOfElement(i),points); // we determine the edge, we would choose for blue Refinement
						if(!edgeElements->getElement(entry).isTaggedForRefinement()){ // If the longestest edge is already the tagged one, we leave the Element alone
							edgeElements->getElement(entry).tagForRefinement(); // we tag the Element for refinement
							this->addMidpoint(meshP1,edgeElements,entry);	// we add the necessary midpoint
							newPoints ++; // add new Points 
							if(elementsOfEdgesLocal.at(entry).size() > 1){
								if(elementsOfEdgesLocal.at(entry).at(0) != -1)
									tagCounter[elementsOfEdgesLocal.at(entry).at(0)] ++;
								else
									globalEdges.push_back(edgeMap->getGlobalElement(entry));  // and add the index for communication later
								if(elementsOfEdgesLocal.at(entry).at(1) != -1)
									tagCounter[elementsOfEdgesLocal.at(entry).at(1)] ++;
								else
									globalEdges.push_back(edgeMap->getGlobalElement(entry));  // and add the index for communication later
							}
							else
								tagCounter[elementsOfEdgesLocal.at(entry).at(0)] ++;

							tagCounter[i] = -1; // This element now has been refined optimally, so no other checks are necessary
							layer++;
						}
					}
					
				}
			}
		}

		reduceAll<int, int> (*this->comm_, REDUCE_SUM, layer, outArg (layer));
		if(this->comm_->getRank() ==0){
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Refinement Tag switched in " << layer << " Elements " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
		}

		// Constructing a map of the global IDs of the tagged Edges	
		Teuchos::ArrayView<GO> globalEdgesArray = Teuchos::arrayViewFromVector( globalEdges);
		MapPtr_Type mapEdgesTagged =
			Teuchos::rcp( new Map_Type( meshP1->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalEdgesArray, 0, this->comm_) );
		// Multivector based on taggesEdgesMap with one as entries
		MultiVectorLOPtr_Type isActiveEdge = Teuchos::rcp( new MultiVectorLO_Type( mapEdgesTagged, 1 ) );
		isActiveEdge->putScalar( (LO) 1);
		//isActiveEdge->print();

		// Creating a Multivector based on the unique edgeMap and repeated edge map with zeros as entries
		MultiVectorLOPtr_Type taggedEdgesLocal = Teuchos::rcp( new MultiVectorLO_Type(mapInterfaceEdges, 1 ) );
		taggedEdgesLocal->putScalar(0);

		taggedEdgesLocal->importFromVector(isActiveEdge, true, "Insert"); // From this we know that -> one signifies a tagged edge -> if we didn't tag it -> refine
	
		//taggedEdgesLocal->print();

		Teuchos::ArrayRCP< const LO >  tags = taggedEdgesLocal->getData( 0 );

		// Adding Midpoints and tagging the edges that were tagged on other Procs an are on the interface
		// Collecting global Ids of tagged Interface Edges
		LO ind;
		GO indG;
		//edgeMap->print();
		for (int i=0; i<tags.size(); i++) {
			if (tags[i] > 0){
				indG = mapInterfaceEdges->getGlobalElement(i);
				ind = edgeMap->getLocalElement(indG);
				globalInterfaceIDsTagged.push_back(indG);
				newPointsCommon ++;
				if(!edgeElements->getElement(ind).isTaggedForRefinement()){
					edgeElements->getElement(ind).tagForRefinement();
					this->addMidpoint(meshP1,edgeElements,ind);
					// globalen Index der Kante -> alle Prozessoren haben diese Kante -> die Knoten die hier hinzugefügt werden sind haben auf allen Proc den gleichen Index
					newPoints ++;

					if(elementsOfEdgesLocal.at(ind).at(0) != -1 ){
						if(tagCounter[elementsOfEdgesLocal.at(ind).at(0)] != -1)
							tagCounter[elementsOfEdgesLocal.at(ind).at(0)] ++;
					}
					if(elementsOfEdgesLocal.at(ind).at(1) != -1){
						if(tagCounter[elementsOfEdgesLocal.at(ind).at(1)] != -1)
							tagCounter[elementsOfEdgesLocal.at(ind).at(1)] ++;
					}
							
				}
			}
	
		}	
	}	
}

// --------------------------------------------------------------------------------
// communicate tagged edges globally
// --------------------------------------------------------------------------------
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::distributeTaggedAndUntaggedEdges(MapConstPtr_Type mapProc, MapConstPtr_Type mapGlobalProc, MapConstPtr_Type mapInterfaceEdgesTaggedUnique, MapConstPtr_Type mapInterfaceEdgesUntaggedUnique, vec_GO_Type& edgesInterfaceUntagged, vec_GO_Type& edgesInterfaceTagged) 
{
		// To Communicate the tagged and untagged edges globally we first determine the number of
		// tagged and untagged interface edges on each Processor and communicate them
		// Then we construct a map based on my procs personal number of those edges and the ones the processor with rank < mine and so on
		// -> Processor 0 holds: (1,...,n) interface edges, Processor 1 holds:(n+1,....,m) ....
		// Then we make this Map unique, as interface edges are on more than one processor at once
		// with 'getGlobalNumElements()' we determine how many interface edges we have and construct a final map on each processor
		// ->(0,...,k) k=getGlobalNumElements()
		// Finally we communicate the Global IDs of the interface edges of each Proc into this vector and have a list of all interface edges
		// (tagged or untagged)
		int maxRank = std::get<1>(this->rankRange_);
		const int myRank = this->comm_->getRank();
		int localNumInterEdgesUnique =  mapInterfaceEdgesTaggedUnique->getNodeNumElements();
		int localNumInterEdgesUntagged = mapInterfaceEdgesUntaggedUnique->getNodeNumElements();
		//---------------------
		// TAGGED EDGES
		// Now we communicate our local number of interface tagged edges
		// -------------------
		MultiVectorLOPtr_Type uniqueTaggedEdgesLocal = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );
		if(maxRank>0)
			uniqueTaggedEdgesLocal->putScalar( (LO) localNumInterEdgesUnique);
		else
			uniqueTaggedEdgesLocal->putScalar( (LO) 0);

		MultiVectorLOPtr_Type uniqueTaggedEdgesGlobal = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		uniqueTaggedEdgesGlobal->putScalar( (LO) 0 ); 

		uniqueTaggedEdgesGlobal->importFromVector( uniqueTaggedEdgesLocal, true, "Insert");

		Teuchos::ArrayRCP<const LO> uniqueTaggedEdgesGlobalNum = uniqueTaggedEdgesGlobal-> getData(0);
				
		GO procOffset=0;
		for(int i=0; i< myRank; i++)
			procOffset = procOffset + uniqueTaggedEdgesGlobalNum[i];

		// Proc 'before me' added 'procOffset' unique tagged Interface Edges, I will add 'mapInterfaceEdgesTaggeUnqiue->getNodeNumElement()' many
		vec_GO_Type localUniqueInterfaceTaggedID(0);
		for(int i=0; i< localNumInterEdgesUnique; i++)
			localUniqueInterfaceTaggedID.push_back(procOffset+i);

		vec_GO_Type globalUniqueInterfaceTaggedID(0);
		for(int i=0; i< mapInterfaceEdgesTaggedUnique->getGlobalNumElements(); i++)
			globalUniqueInterfaceTaggedID.push_back(i);

		// global Map
		Teuchos::ArrayView<GO> globalUniqueInterfaceTaggedIDA = Teuchos::arrayViewFromVector(globalUniqueInterfaceTaggedID);
		MapPtr_Type globalUniqueInterfaceTaggedIDMap =
			Teuchos::rcp( new Map_Type( mapProc->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalUniqueInterfaceTaggedIDA, 0, this->comm_) );

		MultiVectorLOPtr_Type globalUniqueTaggedEdgesGlobalIDs = Teuchos::rcp( new MultiVectorLO_Type(globalUniqueInterfaceTaggedIDMap , 1 ) );
		globalUniqueTaggedEdgesGlobalIDs->putScalar(0);

		// local Map
		Teuchos::ArrayView<GO> localUniqueInterfaceTaggedIDA = Teuchos::arrayViewFromVector(localUniqueInterfaceTaggedID);
		MapPtr_Type localUniqueInterfaceTaggedIDMap =
			Teuchos::rcp( new Map_Type(mapProc->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localUniqueInterfaceTaggedIDA, 0, this->comm_) );

		MultiVectorLOPtr_Type uniqueTaggedEdgesGlobalIDs = Teuchos::rcp( new MultiVectorLO_Type(localUniqueInterfaceTaggedIDMap , 1 ) );

		Teuchos::ArrayRCP<LO> uniqueTaggedEdgesGlobalIDsEntries = uniqueTaggedEdgesGlobalIDs-> getDataNonConst(0);

		for(int i=0; i< localUniqueInterfaceTaggedID.size(); i++)
			uniqueTaggedEdgesGlobalIDsEntries[i] = mapInterfaceEdgesTaggedUnique->getGlobalElement(i);

		globalUniqueTaggedEdgesGlobalIDs->importFromVector(uniqueTaggedEdgesGlobalIDs,false,"Insert");
		//globalUniqueTaggedEdgesGlobalIDs->print();

		Teuchos::ArrayRCP< const LO >  edgesInterfaceTaggedTmp = globalUniqueTaggedEdgesGlobalIDs->getData(0); 

		for(int i=0; i<edgesInterfaceTaggedTmp.size(); i++)
			edgesInterfaceTagged.push_back(edgesInterfaceTaggedTmp[i]);
		

		//globalUniqueTaggedEdgesGlobalIDs->print();
		//---------------------
		// UNTAGGED EDGES
		//---------------------

		MultiVectorLOPtr_Type uniqueUntaggedEdgesLocal = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );
		if(maxRank>0)
			uniqueUntaggedEdgesLocal->putScalar( (LO) localNumInterEdgesUntagged);
		else
			uniqueUntaggedEdgesLocal->putScalar( (LO) 0);

		MultiVectorLOPtr_Type uniqueUntaggedEdgesGlobal = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		uniqueUntaggedEdgesGlobal->putScalar( (LO) 0 ); 

		uniqueUntaggedEdgesGlobal->importFromVector( uniqueUntaggedEdgesLocal, true, "Insert");

		Teuchos::ArrayRCP<const LO> uniqueUntaggedEdgesGlobalNum = uniqueUntaggedEdgesGlobal-> getData(0);
				
		procOffset=0;
		for(int i=0; i< myRank; i++)
			procOffset = procOffset + uniqueUntaggedEdgesGlobalNum[i];

		// Proc 'before me' added procOffset uniqueTaggedInterfaceEdges, i will add 'mapInterfaceEdgesTaggeUnqiue->getNodeNumElement()' many
		vec_GO_Type localUniqueInterfaceUntaggedID(0);
		for(int i=0; i< localNumInterEdgesUntagged ; i++)
			localUniqueInterfaceUntaggedID.push_back(procOffset+i);

		vec_GO_Type globalUniqueInterfaceUntaggedID(0);
		for(int i=0; i<mapInterfaceEdgesUntaggedUnique->getGlobalNumElements(); i++)
			globalUniqueInterfaceUntaggedID.push_back(i);

		Teuchos::ArrayView<GO> globalUniqueInterfaceUntaggedIDA = Teuchos::arrayViewFromVector(globalUniqueInterfaceUntaggedID);
		MapPtr_Type globalUniqueInterfaceUntaggedIDMap =
			Teuchos::rcp( new Map_Type(mapProc->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalUniqueInterfaceUntaggedIDA, 0, this->comm_) );

		MultiVectorLOPtr_Type globalUniqueUntaggedEdgesGlobalIDs = Teuchos::rcp( new MultiVectorLO_Type(globalUniqueInterfaceUntaggedIDMap , 1 ) );
		globalUniqueUntaggedEdgesGlobalIDs->putScalar(0);


		Teuchos::ArrayView<GO> localUniqueInterfaceUntaggedIDA = Teuchos::arrayViewFromVector(localUniqueInterfaceUntaggedID);
		MapPtr_Type localUniqueInterfaceUntaggedIDMap =
			Teuchos::rcp( new Map_Type(mapProc->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localUniqueInterfaceUntaggedIDA, 0, this->comm_) );

		MultiVectorLOPtr_Type uniqueUntaggedEdgesGlobalIDs = Teuchos::rcp( new MultiVectorLO_Type(localUniqueInterfaceUntaggedIDMap , 1 ) );

		Teuchos::ArrayRCP<LO> uniqueUntaggedEdgesGlobalIDsEntries = uniqueUntaggedEdgesGlobalIDs-> getDataNonConst(0);

		for(int i=0; i< localUniqueInterfaceUntaggedID.size(); i++)
			uniqueUntaggedEdgesGlobalIDsEntries[i] = mapInterfaceEdgesUntaggedUnique->getGlobalElement(i);


		globalUniqueUntaggedEdgesGlobalIDs->importFromVector(uniqueUntaggedEdgesGlobalIDs,false,"Insert");

		Teuchos::ArrayRCP< const LO >  edgesInterfaceUntaggedTmp = globalUniqueUntaggedEdgesGlobalIDs->getData(0); 

		for(int i=0; i<edgesInterfaceUntaggedTmp.size(); i++)
			edgesInterfaceUntagged.push_back(edgesInterfaceUntaggedTmp[i]);
		
		//globalUniqueUntaggedEdgesGlobalIDs->print();

	
}
// -----------------------------------------------------------------------------------
// add repeated global IDs of Edges to edgeMap
// -----------------------------------------------------------------------------------
// Functions that adds the repeated global IDs of Edges. We determine the new with 'interfaceEdgesLocalIDs' which tells us the new local IDs of previous interface IDs.
// With that information we can place the right ID to the right edge and in case of a refined edge, the edge with lowest global Node index gets the lower edge index globally
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::addRepeatedIDsEdges(vec_GO_Type& vecGlobalIDsEdges, vec_GO_Type edgesInterfaceUntagged, vec_GO_Type edgesInterfaceTagged,int procOffsetEdgesUniqueSum, vec2D_int_Type interfaceEdgesLocalId)
{

		MapConstPtr_Type edgeMap = this->edgeMap_;
		int countEdge=0;
		for (int i=0; i<edgesInterfaceUntagged.size(); i++){
			if(edgeMap->getLocalElement(edgesInterfaceUntagged[i]) != -1 ){
				for(int j=0; j< interfaceEdgesLocalId.size(); j++){
					if(edgesInterfaceUntagged[i] == interfaceEdgesLocalId.at(j).at(0)){
						vecGlobalIDsEdges.at(interfaceEdgesLocalId.at(j).at(1))= procOffsetEdgesUniqueSum+countEdge;
						}
				} 	
			}
			countEdge ++;
		}	
		GO number1 ,number2;
		GO nodeEdge11,nodeEdge12, nodeEdge21, nodeEdge22;
		LO id1,id2;
		for(int i=0; i< edgesInterfaceTagged.size() ; i++){
			if(edgeMap->getLocalElement(edgesInterfaceTagged[i]) != -1 ){
				for(int j=0; j< interfaceEdgesLocalId.size()-1; j++){
					if(edgesInterfaceTagged[i] == interfaceEdgesLocalId.at(j).at(0)){
						for(int k=j+1; k< interfaceEdgesLocalId.size();k++){		
							if( edgesInterfaceTagged[i] == interfaceEdgesLocalId.at(k).at(0) ){

								// Giving the lower index to the edge with the lowest globalNodeNumber	
								id1=this->edgeElements_->getElement(interfaceEdgesLocalId.at(j).at(1)).getNode(0);
								id2=this->edgeElements_->getElement(interfaceEdgesLocalId.at(j).at(1)).getNode(1);
								nodeEdge11 = this->mapRepeated_->getGlobalElement(id1);
								nodeEdge12 = this->mapRepeated_->getGlobalElement(id2);

								number1= min(nodeEdge11,nodeEdge12);
								
								id1=this->edgeElements_->getElement(interfaceEdgesLocalId.at(k).at(1)).getNode(0);
								id2=this->edgeElements_->getElement(interfaceEdgesLocalId.at(k).at(1)).getNode(1);
								nodeEdge21 = this->mapRepeated_->getGlobalElement(id1);
								nodeEdge22 = this->mapRepeated_->getGlobalElement(id2);

								number2= min(nodeEdge21, nodeEdge22);
				
								if( number1 < number2){
									vecGlobalIDsEdges.at(interfaceEdgesLocalId.at(j).at(1)) = procOffsetEdgesUniqueSum+countEdge;
									vecGlobalIDsEdges.at(interfaceEdgesLocalId.at(k).at(1)) = procOffsetEdgesUniqueSum+countEdge+1;
									}
								else if(number1 > number2){
									vecGlobalIDsEdges.at(interfaceEdgesLocalId.at(k).at(1)) = procOffsetEdgesUniqueSum+countEdge;
									vecGlobalIDsEdges.at(interfaceEdgesLocalId.at(j).at(1))= procOffsetEdgesUniqueSum+countEdge+1;
									}
								else
   									TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "While determining repeated Edge IDs two Edges appear to have the same Node IDs.Please check NodeMap and EdgeMap");
			
							}
						}
					}
				}
			}
			countEdge=countEdge+2;
		}



}


// -----------------------------------------------------------------------------------
// determine new local ID of interfaceEdges 
// -----------------------------------------------------------------------------------
template <class SC, class LO, class GO, class NO>
vec2D_int_Type MeshUnstructuredRefinement<SC,LO,GO,NO>::determineInterfaceEdgesNewLocalIds(EdgeElementsPtr_Type edgeElements, vec_GO_Type edgesInterface,MapConstPtr_Type edgeMap)
{

	vec2D_int_Type interfaceEdgesLocalId(0,vec_int_Type(2));
	vec_int_Type ids(2);

	for(int i=0; i< this->edgeElements_->numberElements(); i++){
		if(this->edgeElements_->getElement(i).isInterfaceElement() == true && this->edgeElements_->getElement(i).getPredecessorElement() !=-1){
		// We check all new edges, that are on the interface and compare them to the recent interfaceEdges 
			// and compare them to the interface edges
			for(int j=0; j<edgesInterface.size();j++){
				// we get the edge the soon to be checked edge might originated from
				if( this->edgeElements_->getElement(i).getPredecessorElement() == edgesInterface[j]){
							ids[0]=edgesInterface[j]; // index of the predecessor global edge
							ids[1]=i; // index of the local edge
							interfaceEdgesLocalId.push_back(ids);
				}			
			}
		}
	}
	sort(interfaceEdgesLocalId.begin(), interfaceEdgesLocalId.end());


return interfaceEdgesLocalId;

}

// -----------------------------------------------------------------------------------
// determine local IDs of new Interface Edges (in 3D new edges appear on the interface)
// -----------------------------------------------------------------------------------
template <class SC, class LO, class GO, class NO>
vec2D_int_Type MeshUnstructuredRefinement<SC,LO,GO,NO>::determineNewInterfaceEdgesLocalIds(vec_GO_Type& globalNodeIDsNewInterfaceEdges) 
{

	vec2D_int_Type newInterfaceEdgesLocalID(0,vec_int_Type(3));
	vec_int_Type ids(3);

	for(int i=0; i< this->edgeElements_->numberElements(); i++){
		if(this->edgeElements_->getElement(i).isInterfaceElement() == true && this->edgeElements_->getElement(i).getPredecessorElement()==-1){
		// We check all new edges, that are on the interface and and don't have a predecessor element 	
				ids[0]=i; // local ID
				ids[1]=this->mapRepeated_->getGlobalElement(this->edgeElements_->getElement(i).getNode(0)); // global node index
				ids[2]=this->mapRepeated_->getGlobalElement(this->edgeElements_->getElement(i).getNode(1));
				sort(ids.begin()+1,ids.end());
				newInterfaceEdgesLocalID.push_back(ids);
				globalNodeIDsNewInterfaceEdges.push_back(ids[1]);
				globalNodeIDsNewInterfaceEdges.push_back(ids[2]);
		}
	}
	// The global IDs List also contains redundant entries, so we make it unqiue
	sort( globalNodeIDsNewInterfaceEdges.begin(), globalNodeIDsNewInterfaceEdges.end() );
	globalNodeIDsNewInterfaceEdges.erase( unique( globalNodeIDsNewInterfaceEdges.begin(), globalNodeIDsNewInterfaceEdges.end() ), globalNodeIDsNewInterfaceEdges.end() );


return newInterfaceEdgesLocalID;

}

// -----------------------------------------------------------------------------------
// Updating ElementsOfEdgesLocal and ElementsOfEdgesGlobal
// -----------------------------------------------------------------------------------
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::updateElementsOfEdgesLocalAndGlobal(vec2D_int_Type interfaceEdgesLocalId, int maxRank, vec_GO_Type edgesInterface, MapConstPtr_Type edgeMap){

	if(maxRank >0){
		vec_GO_Type edgesInterfaceGlobalID(0);
		for(int i=0; i< edgesInterface.size(); i++){
			if(edgeMap->getLocalElement(edgesInterface[i]) != -1 ){
				for(int j=0; j< interfaceEdgesLocalId.size(); j++){
					if(edgesInterface[i] == interfaceEdgesLocalId.at(j).at(0)){
						this->edgeElements_->setElementsOfEdgeLocalEntry(interfaceEdgesLocalId.at(j).at(1),-1);
						edgesInterfaceGlobalID.push_back(this->edgeMap_->getGlobalElement(interfaceEdgesLocalId.at(j).at(1))); // extracting the global IDs of the new interfaceEdges
					}
				}
			}
			
		}

		// communticating elements across interface
		Teuchos::ArrayView<GO> edgesInterfaceGlobalID_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalID);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalID_, 0, this->comm_) );
		//mapGlobalInterface->print();

		// Global IDs of Procs
		// Setting newPoints as to be communicated Values
		MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			interfaceElementsEntries[i] = this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).at(0);
		}

		//interfaceElements->print();

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );

		MultiVectorLOPtr_Type isInterfaceElement_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_imp->putScalar( (LO) 0 ); 
		isInterfaceElement_imp->importFromVector( interfaceElements, false, "Insert");
		//isInterfaceElement_imp->print();

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( interfaceElements, false, "Insert");
		//isInterfaceElement_exp->print();

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_imp, false, "Insert");
		//isInterfaceElement2_imp->print();

		isInterfaceElement2_imp->exportFromVector(isInterfaceElement_exp, false, "Insert");
		// isInterfaceElement2_imp->print();

		interfaceElementsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			this->edgeElements_->setElementsOfEdgeGlobalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),interfaceElementsEntries[i]);
			}

	}

}

// -----------------------------------------------------------------------------------
// adding a Midpoint on an edge
// -----------------------------------------------------------------------------------
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::addMidpoint(MeshUnstrRefPtr_Type meshP1,EdgeElementsPtr_Type edgeElements, int i){

	int midPointInd; // indices of midpoints of edges of soon to be refined element
	int dim = this->dim_;
	vec2D_dbl_Type points2(0, vec_dbl_Type( dim )); // new Points -> depends on already existing point from previously refined neighbouring elements
	vec_dbl_Type point(dim); // midpoint on every edge, that will be added to points, will be added to pointsRep_

	vec_dbl_Type P1(dim),P2(dim); // points we extract from pointsRep_ in order to determine midPoints

	LO p1ID =edgeElements->getElement(i).getNode(0);
	LO p2ID =edgeElements->getElement(i).getNode(1);
	P1 = this->pointsRep_->at(p1ID);
	P2 = this->pointsRep_->at(p2ID);

	for (int d=0; d<dim; d++){
		point[d]= ( (P1)[d] + (P2)[d] ) / 2.;
	}   
	points2.push_back(point); 

	// New Flags:
	this->bcFlagRep_->push_back(edgeElements->getElement(i).getFlag());
		
	// Mittelpunkte der Kanten setzen
	edgeElements->setMidpoint(i, this->pointsRep_->size());

	// We have to keep in mind, that the new added points need to be maped to the associated points
	this->pointsRep_->insert( this->pointsRep_->end(), points2.begin(), points2.end() );
}

// Function to determine longest edge in triangle
template <class SC, class LO, class GO, class NO>
int MeshUnstructuredRefinement<SC,LO,GO,NO>::determineLongestEdge( EdgeElementsPtr_Type edgeElements, vec_int_Type edgeVec, vec2D_dbl_ptr_Type points){
	// We have to determine which edge is longer, as we use the opposite node of the longer edge for the element construction
	vec_dbl_Type length(3);
	vec_dbl_Type P1(2),P2(2);
	double maxLength=0.0;
	int maxEntry=0;
	LO p1ID,p2ID;
	for(int i=0;i<3;i++){
		p1ID =edgeElements->getElement(edgeVec[i]).getNode(0);
		p2ID =edgeElements->getElement(edgeVec[i]).getNode(1);
		P1 = points->at(p1ID);
		P2 = points->at(p2ID);
		length[i] = sqrt(pow(P1[0]-P2[0],2)+pow(P1[1]-P2[1],2));
		if(length[i] > maxLength){
			maxLength = length[i];
			maxEntry= i;
		}
	}
	return edgeVec[maxEntry];
	
}
// Blue Refinement
// refining element according to blue refinement scheme - connecting nodes of shorter edge with midpoint of longer tagged edge and connect that with opposite corner
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::refineBlue( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){
	if(this->dim_ == 2){
		// The necessary Point was already added to the nodelist
		// now we have to figure out, which node it is -> check the tagged edge for midpoint
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement);

		// midpoint index
		vec_int_Type midPointInd(0);
		vec_int_Type taggedEdge(0);
		vec_int_Type untaggedEdge(0);
		int oppositeNodeIndL, oppositeNodeIndS; // we need to determine the nodeIndex of the nodes opposite to the refined edges

		for(int i=0; i<3; i++)	{
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				midPointInd.push_back(edgeElements->getMidpoint(edgeNumbers[i]));		
				taggedEdge.push_back(edgeNumbers[i]);
				}
			else
				untaggedEdge.push_back(edgeNumbers[i]);
			}
		   		
		// We have to determine which edge is longer, as we use the opposite node of the longer edge for the element construction
		double length1, length2;
		int edgeIndexL=0; // index of the longer tagged edge in 'taggedEdge'
		int edgeIndexS=1; // index of the shorter tagged edge
		vec_dbl_Type P1(2),P2(2);
		LO p1ID =edgeElements->getElement(taggedEdge[0]).getNode(0);
		LO p2ID =edgeElements->getElement(taggedEdge[0]).getNode(1);
		P1 = this->pointsRep_->at(p1ID);
		P2 = this->pointsRep_->at(p2ID);
		length1 = sqrt(pow(P1[0]-P2[0],2)+pow(P1[1]-P2[1],2));

		p1ID =edgeElements->getElement(taggedEdge[1]).getNode(0);
		p2ID =edgeElements->getElement(taggedEdge[1]).getNode(1);
		P1 = this->pointsRep_->at(p1ID);
		P2 = this->pointsRep_->at(p2ID);
		length2 = sqrt(pow(P1[0]-P2[0],2)+pow(P1[1]-P2[1],2));

		if(length1 <= length2){
			edgeIndexL=1;
			edgeIndexS=0;		
		}

		// Determine opposite node of the longer Edge
		for(int i=0;i<3; i++){
			if(edgeNumbers[i] != taggedEdge[edgeIndexL]){
				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(taggedEdge[edgeIndexL]).getNode(0))
					oppositeNodeIndL = edgeElements->getElement(edgeNumbers[i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(taggedEdge[edgeIndexL]).getNode(1))
					oppositeNodeIndL = edgeElements->getElement(edgeNumbers[i]).getNode(0);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(taggedEdge[edgeIndexL]).getNode(1))
					oppositeNodeIndL = edgeElements->getElement(edgeNumbers[i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(taggedEdge[edgeIndexL]).getNode(0))
					oppositeNodeIndL = edgeElements->getElement(edgeNumbers[i]).getNode(0);

				i=2;
				}
		   }

        // Determine opposite node of the shorter Edge
		for(int i=0;i<3; i++){
			if(edgeNumbers[i] != taggedEdge[edgeIndexS]){
				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(taggedEdge[edgeIndexS]).getNode(0))
					oppositeNodeIndS = edgeElements->getElement(edgeNumbers[i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(taggedEdge[edgeIndexS]).getNode(1))
					oppositeNodeIndS = edgeElements->getElement(edgeNumbers[i]).getNode(0);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(taggedEdge[edgeIndexS]).getNode(1))
					oppositeNodeIndS = edgeElements->getElement(edgeNumbers[i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(taggedEdge[edgeIndexS]).getNode(0))
					oppositeNodeIndS = edgeElements->getElement(edgeNumbers[i]).getNode(0);

				i=2;
				}
		   }

		// Mutal Node of the two tagged edges
		int mutalNode;
		if(edgeElements->getElement(taggedEdge[0]).getNode(0) == edgeElements->getElement(taggedEdge[1]).getNode(0))
			mutalNode= edgeElements->getElement(taggedEdge[1]).getNode(0);

		if(edgeElements->getElement(taggedEdge[0]).getNode(1) == edgeElements->getElement(taggedEdge[1]).getNode(1))
			mutalNode= edgeElements->getElement(taggedEdge[1]).getNode(1);

		if(edgeElements->getElement(taggedEdge[0]).getNode(0) == edgeElements->getElement(taggedEdge[1]).getNode(1))
			mutalNode= edgeElements->getElement(taggedEdge[1]).getNode(1);

		if(edgeElements->getElement(taggedEdge[0]).getNode(1) == edgeElements->getElement(taggedEdge[1]).getNode(0))
			mutalNode= edgeElements->getElement(taggedEdge[1]).getNode(0);
		


        vec2D_int_Type newElements(3, vec_int_Type( 0 )); // vector for the new elements
		vec2D_int_Type newEdges(9,vec_int_Type(0)); // vector for the new edges
		vec_int_Type edgeFlags(9); // vector for the new flags
        vec2D_LO_Type markedPoints(0);

		vec_bool_Type isInterfaceEdge(9);
		vec_GO_Type predecessorElement(9);;
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =10)

		// Element 1
		(newElements)[0]={oppositeNodeIndL,midPointInd[edgeIndexL],midPointInd[edgeIndexS]};

		(newEdges)[0] = {oppositeNodeIndL ,midPointInd[edgeIndexL]}; 
		(newEdges)[1] = {oppositeNodeIndL,midPointInd[edgeIndexS]}; 
		(newEdges)[2] = {midPointInd[edgeIndexL] ,midPointInd[edgeIndexS]}; 

		edgeFlags[0]=10;
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[edgeIndexS]);
		edgeFlags[2]=10;

		isInterfaceEdge[0] = false;
		isInterfaceEdge[1] = edgeElements->getElement(taggedEdge[edgeIndexS]).isInterfaceElement();
		isInterfaceEdge[2] = false;

		predecessorElement[0] = -1;
		predecessorElement[1] = this->edgeMap_->getGlobalElement(taggedEdge[edgeIndexS]);
		predecessorElement[2] = -1;


		// Element 2
		newElements[1]={oppositeNodeIndL,midPointInd[edgeIndexL],oppositeNodeIndS};

		(newEdges)[3] = {oppositeNodeIndL ,midPointInd[edgeIndexL]}; 
		(newEdges)[4] = {oppositeNodeIndL ,oppositeNodeIndS}; 
		(newEdges)[5] = {midPointInd[edgeIndexL] ,oppositeNodeIndS}; 


		edgeFlags[3]=10;
		edgeFlags[4]=edgeElements->getElement(untaggedEdge[0]).getFlag();
		edgeFlags[5]=this->bcFlagRep_->at(midPointInd[edgeIndexL]);

		isInterfaceEdge[3] = false;
		isInterfaceEdge[4] = edgeElements->getElement(untaggedEdge[0]).isInterfaceElement();
		isInterfaceEdge[5] = edgeElements->getElement(taggedEdge[edgeIndexL]).isInterfaceElement();;

		predecessorElement[3] = -1;
		predecessorElement[4] = this->edgeMap_->getGlobalElement(untaggedEdge[0]);
		predecessorElement[5] = this->edgeMap_->getGlobalElement(taggedEdge[edgeIndexL]);


		// Element 3
		(newElements)[2]={midPointInd[edgeIndexL] , midPointInd[edgeIndexS]  ,mutalNode};

		(newEdges)[6] = {midPointInd[edgeIndexL] ,midPointInd[edgeIndexS]}; 
		(newEdges)[7] = {midPointInd[edgeIndexL] ,mutalNode}; 
		(newEdges)[8] = {midPointInd[edgeIndexS] ,mutalNode}; 

		edgeFlags[6]=10;
		edgeFlags[7]=this->bcFlagRep_->at(midPointInd[edgeIndexL]);
		edgeFlags[8]=this->bcFlagRep_->at(midPointInd[edgeIndexS]);

		isInterfaceEdge[6] = false;
		isInterfaceEdge[7] = edgeElements->getElement(taggedEdge[edgeIndexL]).isInterfaceElement();
		isInterfaceEdge[8] = edgeElements->getElement(taggedEdge[edgeIndexS]).isInterfaceElement();;

		predecessorElement[6] = -1;
		predecessorElement[7] = this->edgeMap_->getGlobalElement(taggedEdge[edgeIndexL]);
		predecessorElement[8] = this->edgeMap_->getGlobalElement(taggedEdge[edgeIndexS]);


        int offsetElements = this->elementsC_->numberElements();
		int offsetEdges = this->edgeElements_->numberElements();

		for( int i=0;i<3; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),10);
			feNew.setFiniteElementRefinementType("blue");
			if(i<2)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		// Kanten hinzufügen
		for( int i=0;i<9; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<6){
				this->edgeElements_->addEdge(feNew,i/3+offsetElements);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=10){
					if ( !this->elementsC_->getElement(i/3+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/3+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/3+offsetElements).addSubElement(feNew);

					}	
				}
			else{
				this->edgeElements_->addEdge(feNew,indexElement);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=10){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	}

				}			
		}	

    }

	else if(this->dim_ == 3){ 
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The irregular Refinement Strategy you requested is only applicable to a 2-dimensional Mesh.");
			}
    


}
// Green Refinement 
// refining the element according to green scheme - connecting node on refined edge with the opposite node (2D)
template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::refineGreen( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){

	if(this->dim_ == 2){
		// The necessary Point was already added to the nodelist
		// now we have to figure out, which node it is -> check the tagged edge for midpoint
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement);

		// midpoint index
		int midPointInd;
		// Id of tagged Edge (id in edgenumbers)
		int taggedEdgeId;
		// Actual Id of edge (id in Edgelist)
		int taggedEdge;
		int oppositeNodeInd; // nodeInd of node opposite to tagged edge
		vec_int_Type sortedEdges(0); // Sorting the edges in edgeNumers as follows: {taggedEdge, untaggedEdge1, untaggedEdge2}

		// Extracting midPoint Ind
		for(int i=0; i<3; i++)	{
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				midPointInd = edgeElements->getMidpoint(edgeNumbers[i]);		
				taggedEdgeId = i;
				taggedEdge = edgeNumbers[i];
		        sortedEdges.push_back(edgeNumbers[i]); // tagged edge first to sortedEdges
				}
			}

		// Determine opposite node of tagged edge
		for(int i=0;i<3; i++){
			if(edgeNumbers[i] != taggedEdge){
				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(taggedEdge).getNode(0))
					oppositeNodeInd = edgeElements->getElement(edgeNumbers[i]).getNode(1);

				else if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(taggedEdge).getNode(1))
					oppositeNodeInd = edgeElements->getElement(edgeNumbers[i]).getNode(0);

				else if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(taggedEdge).getNode(1))
					oppositeNodeInd = edgeElements->getElement(edgeNumbers[i]).getNode(1);

				else if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(taggedEdge).getNode(0))
					oppositeNodeInd = edgeElements->getElement(edgeNumbers[i]).getNode(0);

				i=2;
				}
		   }

		for(int i=0; i<3 ; i++){
			if(i != taggedEdgeId )
				sortedEdges.push_back(edgeNumbers[i]);
		}

		
		// adding Elements and Edges
		vec2D_int_Type newElements(2, vec_int_Type( 0 )); // vector for the new elements
		vec2D_int_Type newEdges(6,vec_int_Type(0)); // vector for the new edges
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =10)
		vec_int_Type edgeFlags(6); // vector for the new flags
        vec2D_LO_Type markedPoints(0);

		vec_bool_Type isInterfaceEdge(6);
		vec_GO_Type predecessorElement(6);

		// Element 1
		(newElements)[0]={midPointInd, edgeElements->getElement(sortedEdges[1]).getNode(0),  edgeElements->getElement(sortedEdges[1]).getNode(1)};

		(newEdges)[0] = {midPointInd, edgeElements->getElement(sortedEdges[1]).getNode(0)}; 
		(newEdges)[1] = {midPointInd, edgeElements->getElement(sortedEdges[1]).getNode(1)}; 
		(newEdges)[2] = {edgeElements->getElement(sortedEdges[1]).getNode(0) ,edgeElements->getElement(sortedEdges[1]).getNode(1) }; 

		if(oppositeNodeInd == edgeElements->getElement(sortedEdges[1]).getNode(0)){
			edgeFlags[0]=10;
			isInterfaceEdge[0] = false;
			predecessorElement[0] = -1;
			}
		else {
			edgeFlags[0]=this->bcFlagRep_->at(midPointInd);
			isInterfaceEdge[0] = edgeElements->getElement(taggedEdge).isInterfaceElement();
			predecessorElement[0] = this->edgeMap_->getGlobalElement(taggedEdge);
		}
	
	    if(oppositeNodeInd == edgeElements->getElement(sortedEdges[1]).getNode(1)){
			edgeFlags[1]=10;
			isInterfaceEdge[1] = false;
			predecessorElement[1] = -1;
		}
		else {
			edgeFlags[1]=this->bcFlagRep_->at(midPointInd);
			isInterfaceEdge[1] = edgeElements->getElement(taggedEdge).isInterfaceElement();
			predecessorElement[1] = this->edgeMap_->getGlobalElement(taggedEdge);
		}
	
		edgeFlags[2]=edgeElements->getElement(sortedEdges[1]).getFlag();

		isInterfaceEdge[2] = edgeElements->getElement(sortedEdges[1]).isInterfaceElement();

		predecessorElement[2] = this->edgeMap_->getGlobalElement(sortedEdges[1]);

		// Element 2
		(newElements)[1]={midPointInd, edgeElements->getElement(sortedEdges[2]).getNode(0),  edgeElements->getElement(sortedEdges[2]).getNode(1)};

		(newEdges)[3] = {midPointInd, edgeElements->getElement(sortedEdges[2]).getNode(0)}; 
		(newEdges)[4] = {midPointInd, edgeElements->getElement(sortedEdges[2]).getNode(1)}; 
		(newEdges)[5] = {edgeElements->getElement(sortedEdges[2]).getNode(0) ,edgeElements->getElement(sortedEdges[2]).getNode(1) }; 

		if(oppositeNodeInd == edgeElements->getElement(sortedEdges[2]).getNode(0)){
			edgeFlags[3]=10;
			isInterfaceEdge[3] = false;
			predecessorElement[3] = -1;

		}	
		else {
			edgeFlags[3]=this->bcFlagRep_->at(midPointInd);
			isInterfaceEdge[3] = edgeElements->getElement(taggedEdge).isInterfaceElement();
			predecessorElement[3] = this->edgeMap_->getGlobalElement(taggedEdge);
		}
	
	    if(oppositeNodeInd == edgeElements->getElement(sortedEdges[2]).getNode(1)){
			edgeFlags[4]=10;
			isInterfaceEdge[4] = false;
			predecessorElement[4] =-1;
		}
		else {
			edgeFlags[4]=this->bcFlagRep_->at(midPointInd);
			isInterfaceEdge[4] = edgeElements->getElement(taggedEdge).isInterfaceElement();
			predecessorElement[4] = this->edgeMap_->getGlobalElement(taggedEdge);
		}
	
		edgeFlags[5]=edgeElements->getElement(sortedEdges[2]).getFlag();

		isInterfaceEdge[5] = edgeElements->getElement(sortedEdges[2]).isInterfaceElement();

		predecessorElement[5] = this->edgeMap_->getGlobalElement(sortedEdges[2]);

		// We add one element and switch one element with the element 'indexElement'
        int offsetElements = this->elementsC_->numberElements(); // for determining which edge corresponds to which element
		int offsetEdges = this->edgeElements_->numberElements(); // just for printing

		for( int i=0;i<2; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),10);
			feNew.setFiniteElementRefinementType("green"); // setting green refinement type in order to check in following refinement stages
			if(i==0)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		// add edges
		for( int i=0;i<6; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<3){
				this->edgeElements_->addEdge(feNew,offsetElements);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=10){
					if ( !this->elementsC_->getElement(offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(offsetElements).addSubElement(feNew);
					}
				}
			else{
				this->edgeElements_->addEdge(feNew,indexElement);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=10){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);
					}
				}
			}				
    }

	else if(this->dim_ == 3){ 
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The irregular Refinement Strategy you requested is only applicable to a 2-dimensional Mesh.");
			}
    


}

// Red Refinement 
// refining the element regularly, aka red -> one element is refined into 4 (2D) 

template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::refineRed( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){

	if(this->dim_ == 2){
		
        // The necessary Point was already added to the nodelist
		// now we have to figure out, which node it is -> check the tagged edge for midpoint
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement);

		// midpoint index
		vec_int_Type midPointInd(3);

		for(int i=0; i<3; i++)	
				midPointInd[i] = edgeElements->getMidpoint(edgeNumbers[i]);
				
		// -> Edge 1: midPoint[0] - mutualNode[0] = mutualNode of Edge 1 and 2 
		// -> Edge 2: midPoint[1] - mutualNode[1] = mutualNode of Edge 1 and 3
		// -> Edge 3: midpoint[3] - mutualNode[2] = mutualNode of Edge 2 and 3
		// Mutal Node of two edges
		vec_int_Type mutualNode(3);
		int ind=0;
		for(int i=0;i<2; i++){
			for(int j=1;j+i<3;j++){
				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(edgeNumbers[j+i]).getNode(0))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(0);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(edgeNumbers[j+i]).getNode(1))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(edgeNumbers[j+i]).getNode(1))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(edgeNumbers[j+i]).getNode(0))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(0);

				ind++;
			}
		}

		// Adding Elements and the corresponding Edges
		vec2D_int_Type newElements(4, vec_int_Type( 0 )); // vector for the new elements
		vec2D_int_Type newEdges(12,vec_int_Type(0)); // vector for the new edges
		vec_int_Type edgeFlags(12); // vector for the new flags
		vec_bool_Type isInterfaceEdge(12); // bool vectot for interfaceEdges
		vec_GO_Type predecessorElement(12);
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =10)

		// Element 1
		(newElements)[0]={mutualNode[0],midPointInd[0],midPointInd[1]};

		(newEdges)[0] = {mutualNode[0] ,midPointInd[0]}; 
		(newEdges)[1] = {mutualNode[0] ,midPointInd[1]}; 
		(newEdges)[2] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[0]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[2]=10;

		isInterfaceEdge[0] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[1] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[2] = false;

		predecessorElement[0] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[1] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[2] = -1;

		// Element 2
		newElements[1]={mutualNode[1],midPointInd[0],midPointInd[2]};

		(newEdges)[3] = {mutualNode[1] ,midPointInd[0]}; 
		(newEdges)[4] = {mutualNode[1] ,midPointInd[2]}; 
		(newEdges)[5] = {midPointInd[0] ,midPointInd[2]}; 

		edgeFlags[3]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[4]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[5]=10;

		isInterfaceEdge[3] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[4] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[5] = false;

		predecessorElement[3] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[4] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[5] = -1;

		// Element 3
		(newElements)[2]={mutualNode[2] , midPointInd[1] ,midPointInd[2]};

		(newEdges)[6] = {mutualNode[2] ,midPointInd[1]}; 
		(newEdges)[7] = {mutualNode[2] ,midPointInd[2]}; 
		(newEdges)[8] = {midPointInd[1] ,midPointInd[2]}; 

		edgeFlags[6]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[7]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[8]=10;

		isInterfaceEdge[6] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[7] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[8] = false;

		predecessorElement[6] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[7] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[8] = -1;
		// Element 4
		(newElements)[3]={midPointInd[0],midPointInd[1],midPointInd[2]};	

		(newEdges)[9] = {midPointInd[0] ,midPointInd[1]}; 
		(newEdges)[10] = {midPointInd[1] ,midPointInd[2]}; 
		(newEdges)[11] = {midPointInd[2] ,midPointInd[0]}; 

		edgeFlags[9]=10;
		edgeFlags[10]=10;
		edgeFlags[11]=10;

		isInterfaceEdge[9] = false;
		isInterfaceEdge[10] = false;
		isInterfaceEdge[11] = false;

		predecessorElement[9] = -1;
		predecessorElement[10] = -1;
		predecessorElement[11] = -1;

		int offsetElements = this->elementsC_->numberElements();
		int offsetEdges = this->edgeElements_->numberElements();
		for( int i=0;i<4; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),0);
			feNew.setFiniteElementRefinementType("red");
				
			if(i<3)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		for( int i=0;i<12; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<9){
				this->edgeElements_->addEdge(feNew,i/3+offsetElements);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=10){
					if ( !this->elementsC_->getElement(i/3+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/3+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/3+offsetElements).addSubElement(feNew);
					}
				}
			else
				this->edgeElements_->addEdge(feNew,indexElement);
			
		}
    }

	else if(this->dim_ == 3){ 
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The irregular Refinement Strategy you requested is only applicable to a 2-dimensional Mesh.");
			}
    
}
// Regular Refinement 
// refining the element regularly, aka red -> one element is refined into 4 (2D) 

template <class SC, class LO, class GO, class NO>
void MeshUnstructuredRefinement<SC,LO,GO,NO>::refineRegular( MeshUnstrRefPtr_Type meshP1, EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){
	if(this->dim_ == 2){
        vec_int_Type midPointInd( 3 ); // indices of midpoints of edges of soon to be refined element
      
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		int k=0;
		for(int i=0; i<3; i++)	{
			if(!edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()) // we tag every edge, after we refine an element -> no tag - no refinement on that edge so far
				{   			
				this->addMidpoint(meshP1,edgeElements,edgeNumbers[i]);
				midPointInd[i] = edgeElements->getMidpoint(edgeNumbers[i]); 
			}

			else
				midPointInd[i] = edgeElements->getMidpoint(edgeNumbers[i]);
		}

		// Mutal Node of two edges
		vec_int_Type mutualNode(3);
		int ind=0;
		for(int i=0;i<2; i++){
			for(int j=1;j+i<3;j++){
				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(edgeNumbers[j+i]).getNode(0))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(0);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(edgeNumbers[j+i]).getNode(1))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(0) == edgeElements->getElement(edgeNumbers[j+i]).getNode(1))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(1);

				if(edgeElements->getElement(edgeNumbers[i]).getNode(1) == edgeElements->getElement(edgeNumbers[j+i]).getNode(0))
					mutualNode[ind]= edgeElements->getElement(edgeNumbers[j+i]).getNode(0);

				ind++;

		}}

		// Adding Elements and the corresponding Edges
		vec2D_int_Type newElements(4, vec_int_Type( 0 )); // vector for the new elements
		vec2D_int_Type newEdges(12,vec_int_Type(0)); // vector for the new edges
		vec_int_Type edgeFlags(12); // vector for the new flags

		vec_bool_Type isInterfaceEdge(12); // bool vector for interfaceEdges
		vec_GO_Type predecessorElement(12);
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =10)

		// Element 1
		(newElements)[0]={mutualNode[0],midPointInd[0],midPointInd[1]};

		(newEdges)[0] = {mutualNode[0] ,midPointInd[0]}; 
		(newEdges)[1] = {mutualNode[0] ,midPointInd[1]}; 
		(newEdges)[2] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[0]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[2]=10;

		isInterfaceEdge[0] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[1] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[2] = false;

		predecessorElement[0] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[1] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[2] = -1;


		// Element 2
		newElements[1]={mutualNode[1],midPointInd[0],midPointInd[2]};

		(newEdges)[3] = {mutualNode[1] ,midPointInd[0]}; 
		(newEdges)[4] = {mutualNode[1] ,midPointInd[2]}; 
		(newEdges)[5] = {midPointInd[0] ,midPointInd[2]}; 

		edgeFlags[3]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[4]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[5]=10;

		isInterfaceEdge[3] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[4] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[5] = false;

		predecessorElement[3] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[4] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[5] = -1;

		// Element 3
		(newElements)[2]={mutualNode[2] , midPointInd[1] ,midPointInd[2]};

		(newEdges)[6] = {mutualNode[2] ,midPointInd[1]}; 
		(newEdges)[7] = {mutualNode[2] ,midPointInd[2]}; 
		(newEdges)[8] = {midPointInd[1] ,midPointInd[2]}; 

		edgeFlags[6]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[7]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[8]=10;

		isInterfaceEdge[6] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[7] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[8] = false;

		predecessorElement[6] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[7] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[8] = -1;


		// Element 4
		(newElements)[3]={midPointInd[0],midPointInd[1],midPointInd[2]};	

		(newEdges)[9] = {midPointInd[0] ,midPointInd[1]}; 
		(newEdges)[10] = {midPointInd[1] ,midPointInd[2]}; 
		(newEdges)[11] = {midPointInd[2] ,midPointInd[0]}; 

		edgeFlags[9]=10;
		edgeFlags[10]=10;
		edgeFlags[11]=10;

		isInterfaceEdge[9] = false;
		isInterfaceEdge[10] = false;
		isInterfaceEdge[11] = false;

		predecessorElement[9] = -1;
		predecessorElement[10] = -1;
		predecessorElement[11] = -1;


		int offsetElements = this->elementsC_->numberElements();
		int offsetEdges = this->edgeElements_->numberElements();
		for( int i=0;i<4; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),0);
			feNew.setFiniteElementRefinementType("regular");	
			if(i<3)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		for( int i=0;i<12; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<9){
				this->edgeElements_->addEdge(feNew,i/3+offsetElements);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=10){
				 	if ( !this->elementsC_->getElement(i/3+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/3+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/3+offsetElements).addSubElement(feNew);

					}		
			}
			else
				this->edgeElements_->addEdge(feNew,indexElement);
			
		}
		
		// tagging the Edges for refinement
	   	for (int d=0; d<3; d++){
			// now we tag the edge as 'refined'
			edgeElements->getElement(edgeNumbers[d]).tagForRefinement();	   
		}
	
	}	
	// -------------------------------------------------------------------
	// 3D Regular Refinement Algorithm 
	// -------------------------------------------------------------------
	else if(this->dim_ == 3){
		// The way we refine the Tetrahedron is defined by how we order the nodes of the tetrahedron
		// (For the algorithm see "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955)

        vec_int_Type midPointInd( 6 ); // indices of midpoints of edges of soon to be refined element
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element

		// Extract the four points of tetraeder
		vec_int_Type nodeInd(0);
		for(int i=0; i<6; i++)	{
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
		}
		sort( nodeInd.begin(), nodeInd.end() );
		nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );
		
		// Right now the Nodes are ordered by the local Indices, which differ depending on the number of Processors. Hence the refinements are not
		// 100% equal when using a different number of processors.
		// If we sort the nodes by their values and not their local Indices, we can solve that problem, as these values don't change

		vec2D_dbl_Type points(4,vec_dbl_Type(4));
		points[0] = {this->pointsRep_->at(nodeInd[0]).at(0),this->pointsRep_->at(nodeInd[0]).at(1),this->pointsRep_->at(nodeInd[0]).at(2), (double) nodeInd[0]};
		points[1] = {this->pointsRep_->at(nodeInd[1]).at(0),this->pointsRep_->at(nodeInd[1]).at(1),this->pointsRep_->at(nodeInd[1]).at(2), (double) nodeInd[1]};
		points[2] = {this->pointsRep_->at(nodeInd[2]).at(0),this->pointsRep_->at(nodeInd[2]).at(1),this->pointsRep_->at(nodeInd[2]).at(2), (double) nodeInd[2]};
		points[3] = {this->pointsRep_->at(nodeInd[3]).at(0),this->pointsRep_->at(nodeInd[3]).at(1),this->pointsRep_->at(nodeInd[3]).at(2), (double) nodeInd[3]};

		sort(points.begin(), points.end());

		nodeInd[0] = (int) points[0][3];
		nodeInd[1] = (int) points[1][3];
		nodeInd[2] = (int) points[2][3];
		nodeInd[3] = (int) points[3][3];

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

		// Tri_0 = [x_0,x_1,x_2]
		// Tri_1 = [x_0,x_1,x_3]
		// Tri_2 = [x_0,x_2,x_3]
		// Tri_3 = [x_1,x_2,x_3]


		// We check if one or more of these triangles are part of the boundary surface and determine there flag

		vec2D_int_Type originTriangles(4,vec_int_Type(3));
		originTriangles[0] = {nodeInd[0], nodeInd[1], nodeInd[2] };
		originTriangles[1] = {nodeInd[0], nodeInd[1], nodeInd[3] };
		originTriangles[2] = {nodeInd[0], nodeInd[2], nodeInd[3] };
		originTriangles[3] = {nodeInd[1], nodeInd[2], nodeInd[3] };
		
		
		vec_int_Type originFlag(4,10); // Triangle Flag

		int numberSubElSurf=0;
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);
		int entry; 
		if (this->elementsC_->getElement(indexElement).subElementsInitialized() ){
			numberSubElSurf = this->elementsC_->getElement(indexElement).getSubElements()->numberElements();
			for(int i=0; i< numberSubElSurf ; i++){
				triTmp = this->elementsC_->getElement(indexElement).getSubElements()->getElement(i).getVectorNodeList();
				for(int j=0; j<4 ; j++){
					originTriangleTmp = originTriangles[j];
					sort(originTriangleTmp.begin(),originTriangleTmp.end());
					if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] ) 
						originFlag[j] = this->elementsC_->getElement(indexElement).getSubElements()->getElement(i).getFlag();
					//auto it1 = find( originTriangles.begin(), originTriangles.end() ,triTmp );
               		//entry = distance( originTriangles.begin() , it1 );
				}
			}
		}

		// Furthermore we have to determine whether the triangles are part of the interface between processors, as we need this information to determine if edges
		// that emerge on the triangles are part of the interface
		// A triangle is part of the interface if all of its edges are part of the interface (the information if edges are part of the interface was determined
		// in the beginning of the Mesh Refinement by 'determineInterfaceEdges')

		vec_bool_Type interfaceSurface(4,false);

		if(edgeElements->getElement(edgeNumbers[0]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[1]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[3]).isInterfaceElement())
			interfaceSurface[0] = true;
		if(edgeElements->getElement(edgeNumbers[0]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[2]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[4]).isInterfaceElement())
			interfaceSurface[1] = true;
		if(edgeElements->getElement(edgeNumbers[1]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[2]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[5]).isInterfaceElement())
			interfaceSurface[2] = true;
		if(edgeElements->getElement(edgeNumbers[3]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[4]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[5]).isInterfaceElement())
			interfaceSurface[3] = true;

		// Finally we need to determine or extract the indices of the edges midpoints. As in the describe algorithm the midpoints are set as follows:
		// Edge 0 = [x_0,x_1] -> x_01
		// Edge 1 = [x_0,x_2] -> x_02
		// Edge 2 = [x_0,x_3] -> x_03
		// Edge 3 = [x_1,x_2] -> x_12
		// Edge 4 = [x_1,x_3] -> x_13
		// Edge 5 = [x_2,x_3] -> x_23

		for(int i=0; i<6; i++)	{
			if(!edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()) // we tag every edge, after we refine an element -> no tag - no refinement on that edge so far
				{   			
				this->addMidpoint(meshP1,edgeElements,edgeNumbers[i]);
				midPointInd[i] = edgeElements->getMidpoint(edgeNumbers[i]); 
			}
			else
				midPointInd[i] = edgeElements->getMidpoint(edgeNumbers[i]);
	
		}
	
		// Now we construct the new Elements as proposed by Bey's Regular Refinement Algorithm

		vec2D_int_Type newElements(8, vec_int_Type( 0 )); // new elements
		vec2D_int_Type newEdges(48,vec_int_Type(0)); // new edges
		vec2D_int_Type newTriangles(32,vec_int_Type(0)); // new Triangles
		vec_int_Type newTrianglesFlag(32) ; // new Triangle Flags
		vec_int_Type edgeFlags(48); // new Edge flags
		vec_bool_Type isInterfaceEdge(48); // bool vector for interfaceEdge Information
		vec_GO_Type predecessorElement(48); // vector that stores the global IDs of the predecessor of each edge 


		// How are Flags determined?
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =10)
		// If an edges emerges on a triangle, the flag is determined by the triangle flag. Opposite to the 2D case, edges that connect midpoints are not automatically interior edges, but are
		// determined by the triangle/surface they are on


		// Element 1: (x_0,x_01,x_02,x_03)
		(newElements)[0]={nodeInd[0],midPointInd[0],midPointInd[1],midPointInd[2]};

		(newEdges)[0] = {nodeInd[0] ,midPointInd[0]}; 
		(newEdges)[1] = {nodeInd[0] ,midPointInd[1]}; 
		(newEdges)[2] = {nodeInd[0] ,midPointInd[2]}; 
		(newEdges)[3] = {midPointInd[0] ,midPointInd[1]}; 
		(newEdges)[4] = {midPointInd[0] ,midPointInd[2]}; 
		(newEdges)[5] = {midPointInd[1] ,midPointInd[2]}; 

		edgeFlags[0]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[2]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[3]=originFlag[0];
		edgeFlags[4]=originFlag[1];
		edgeFlags[5]=originFlag[2];

		isInterfaceEdge[0] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[1] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[2] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[3] = interfaceSurface[0];
		isInterfaceEdge[4] = interfaceSurface[1];
		isInterfaceEdge[5] = interfaceSurface[2];

		predecessorElement[0] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[1] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[2] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[3] = -1;
		predecessorElement[4] = -1;
		predecessorElement[5] = -1;

		// Subelements of thetrahedron
		newTriangles[0]= {nodeInd[0],midPointInd[0],midPointInd[1]};
		newTriangles[1]= {nodeInd[0],midPointInd[0],midPointInd[2]};
		newTriangles[2]= {nodeInd[0],midPointInd[1],midPointInd[2]};
		newTriangles[3]= {midPointInd[0],midPointInd[1],midPointInd[2]};

		newTrianglesFlag[0]= originFlag[0]; 
		newTrianglesFlag[1]= originFlag[1]; 
		newTrianglesFlag[2]= originFlag[2]; 
		newTrianglesFlag[3]= 10;

		// Element 2: (x_1,x_01,x_12,x_13)
		(newElements)[1]={nodeInd[1],midPointInd[0],midPointInd[3],midPointInd[4]};

		(newEdges)[6] = {nodeInd[1] ,midPointInd[0]}; 
		(newEdges)[7] = {nodeInd[1] ,midPointInd[3]}; 
		(newEdges)[8] = {nodeInd[1] ,midPointInd[4]}; 
		(newEdges)[9] = {midPointInd[0] ,midPointInd[3]}; 
		(newEdges)[10] = {midPointInd[0] ,midPointInd[4]}; 
		(newEdges)[11] = {midPointInd[3] ,midPointInd[4]}; 

		edgeFlags[6]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[7]=this->bcFlagRep_->at(midPointInd[3]);
		edgeFlags[8]=this->bcFlagRep_->at(midPointInd[4]);
		edgeFlags[9]=originFlag[0];
		edgeFlags[10]=originFlag[1];
		edgeFlags[11]=originFlag[3];

		isInterfaceEdge[6] =edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[7] =edgeElements->getElement(edgeNumbers[3]).isInterfaceElement();
		isInterfaceEdge[8] =edgeElements->getElement(edgeNumbers[4]).isInterfaceElement();
		isInterfaceEdge[9] = interfaceSurface[0];
		isInterfaceEdge[10] = interfaceSurface[1];
		isInterfaceEdge[11] = interfaceSurface[3];

		predecessorElement[6] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[7] = this->edgeMap_->getGlobalElement(edgeNumbers[3]);
		predecessorElement[8] = this->edgeMap_->getGlobalElement(edgeNumbers[4]);
		predecessorElement[9] = -1;
		predecessorElement[10] = -1;
		predecessorElement[11] = -1;

		// Subelements of tetrahedron
		newTriangles[4]= {nodeInd[1],midPointInd[0],midPointInd[3]};
		newTriangles[5]= {nodeInd[1],midPointInd[0],midPointInd[4]};
		newTriangles[6]= {nodeInd[1],midPointInd[3],midPointInd[4]};
		newTriangles[7]= {midPointInd[0],midPointInd[3],midPointInd[4]};

		newTrianglesFlag[4]= originFlag[0]; 
		newTrianglesFlag[5]= originFlag[1]; 
		newTrianglesFlag[6]= originFlag[3]; 
		newTrianglesFlag[7]= 10;


		// Element 3: (x_2,x_02,x_12,x_23)
		(newElements)[2]={nodeInd[2],midPointInd[1],midPointInd[3],midPointInd[5]};

		(newEdges)[12] = {nodeInd[2] ,midPointInd[1]}; 
		(newEdges)[13] = {nodeInd[2] ,midPointInd[3]}; 
		(newEdges)[14] = {nodeInd[2] ,midPointInd[5]}; 
		(newEdges)[15] = {midPointInd[1] ,midPointInd[3]}; 
		(newEdges)[16] = {midPointInd[1] ,midPointInd[5]}; 
		(newEdges)[17] = {midPointInd[3] ,midPointInd[5]}; 

		edgeFlags[12]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[13]=this->bcFlagRep_->at(midPointInd[3]);
		edgeFlags[14]=this->bcFlagRep_->at(midPointInd[5]);
		edgeFlags[15]=originFlag[0];
		edgeFlags[16]=originFlag[2];
		edgeFlags[17]=originFlag[3];

		isInterfaceEdge[12] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[13] = edgeElements->getElement(edgeNumbers[3]).isInterfaceElement();
		isInterfaceEdge[14] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[15] = interfaceSurface[0];
		isInterfaceEdge[16] = interfaceSurface[2];
		isInterfaceEdge[17] = interfaceSurface[3];

		predecessorElement[12] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[13] = this->edgeMap_->getGlobalElement(edgeNumbers[3]);
		predecessorElement[14] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[15] = -1;
		predecessorElement[16] = -1;
		predecessorElement[17] = -1;

		// Subelements of thetrahedron
		newTriangles[8]= {nodeInd[2],midPointInd[1],midPointInd[3]};
		newTriangles[9]= {nodeInd[2],midPointInd[1],midPointInd[5]};
		newTriangles[10]= {nodeInd[2],midPointInd[3],midPointInd[5]};
		newTriangles[11]= {midPointInd[1],midPointInd[3],midPointInd[5]};

		newTrianglesFlag[8]= originFlag[0]; 
		newTrianglesFlag[9]= originFlag[2]; 
		newTrianglesFlag[10]= originFlag[3]; 
		newTrianglesFlag[11]= 10;


		// Element 4: (x_3,x_03,x_13,x_23)
		(newElements)[3]={nodeInd[3],midPointInd[2],midPointInd[4],midPointInd[5]};

		(newEdges)[18] = {nodeInd[3] ,midPointInd[2]}; 
		(newEdges)[19] = {nodeInd[3] ,midPointInd[4]}; 
		(newEdges)[20] = {nodeInd[3] ,midPointInd[5]}; 
		(newEdges)[21] = {midPointInd[2] ,midPointInd[4]}; 
		(newEdges)[22] = {midPointInd[2] ,midPointInd[5]}; 
		(newEdges)[23] = {midPointInd[4] ,midPointInd[5]}; 

		edgeFlags[18]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[19]=this->bcFlagRep_->at(midPointInd[4]);
		edgeFlags[20]=this->bcFlagRep_->at(midPointInd[5]);
		edgeFlags[21]=originFlag[1];
		edgeFlags[22]=originFlag[2];
		edgeFlags[23]=originFlag[3];

		isInterfaceEdge[18] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[19] = edgeElements->getElement(edgeNumbers[4]).isInterfaceElement();
		isInterfaceEdge[20] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[21] = interfaceSurface[1];
		isInterfaceEdge[22] = interfaceSurface[2];
		isInterfaceEdge[23] = interfaceSurface[3];

		predecessorElement[18] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[19] = this->edgeMap_->getGlobalElement(edgeNumbers[4]);
		predecessorElement[20] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[21] = -1;
		predecessorElement[22] = -1;
		predecessorElement[23] = -1;


		// Subelements of thetrahedron
		newTriangles[12]= {nodeInd[3],midPointInd[2],midPointInd[4]};
		newTriangles[13]= {nodeInd[3],midPointInd[2],midPointInd[5]};
		newTriangles[14]= {nodeInd[3],midPointInd[4],midPointInd[5]};
		newTriangles[15]= {midPointInd[2],midPointInd[4],midPointInd[5]};

		newTrianglesFlag[12]= originFlag[1]; 
		newTrianglesFlag[13]= originFlag[2]; 
		newTrianglesFlag[14]= originFlag[3]; 
		newTrianglesFlag[15]= 10;

		// The following elements are constructed only with the edges midpoints, hence they have the surface flag and interface characteristic
		// Element 5: (x_01,x_02,x_03,x_12)
		(newElements)[4]={midPointInd[0],midPointInd[1],midPointInd[2],midPointInd[4]};

		(newEdges)[24] = {midPointInd[0] ,midPointInd[1]}; 
		(newEdges)[25] = {midPointInd[0] ,midPointInd[2]}; 
		(newEdges)[26] = {midPointInd[0] ,midPointInd[4]}; 
		(newEdges)[27] = {midPointInd[1] ,midPointInd[2]}; 
		(newEdges)[28] = {midPointInd[1] ,midPointInd[4]}; 
		(newEdges)[29] = {midPointInd[2] ,midPointInd[4]}; 

		edgeFlags[24]=originFlag[0];
		edgeFlags[25]=originFlag[1];
		edgeFlags[26]=originFlag[1];
		edgeFlags[27]=originFlag[2];
		edgeFlags[28]=10;
		edgeFlags[29]=originFlag[1];

		isInterfaceEdge[24] = interfaceSurface[0];
		isInterfaceEdge[25] = interfaceSurface[1];
		isInterfaceEdge[26] = interfaceSurface[1];
		isInterfaceEdge[27] = interfaceSurface[2];
		isInterfaceEdge[28] = false;
		isInterfaceEdge[29] = interfaceSurface[1];

		predecessorElement[24] = -1;
		predecessorElement[25] = -1;
		predecessorElement[26] = -1;
		predecessorElement[27] = -1;
		predecessorElement[28] = -1;
		predecessorElement[29] = -1;

		// Subelements of thetrahedron
		newTriangles[16]= {midPointInd[0],midPointInd[1],midPointInd[2]};
		newTriangles[17]= {midPointInd[0],midPointInd[1],midPointInd[4]};
		newTriangles[18]= {midPointInd[0],midPointInd[2],midPointInd[4]};
		newTriangles[19]= {midPointInd[1],midPointInd[2],midPointInd[4]};

		newTrianglesFlag[16]= 10; 
		newTrianglesFlag[17]= 10; 
		newTrianglesFlag[18]= originFlag[1]; 
		newTrianglesFlag[19]= 10;

		// Element 6: (x_01,x_02,x_12,x_13)
		(newElements)[5]={midPointInd[0],midPointInd[1],midPointInd[3],midPointInd[4]};

		(newEdges)[30] = {midPointInd[0],midPointInd[1]}; 
		(newEdges)[31] = {midPointInd[0] ,midPointInd[3]}; 
		(newEdges)[32] = {midPointInd[0],midPointInd[4]}; 
		(newEdges)[33] = {midPointInd[1] ,midPointInd[3]}; 
		(newEdges)[34] = {midPointInd[1] ,midPointInd[4]}; 
		(newEdges)[35] = {midPointInd[3] ,midPointInd[4]}; 

		edgeFlags[30]=originFlag[0];
		edgeFlags[31]=originFlag[0];
		edgeFlags[32]=originFlag[1];
		edgeFlags[33]=originFlag[0];
		edgeFlags[34]=10;
		edgeFlags[35]=originFlag[3];

		isInterfaceEdge[30] = interfaceSurface[0];
		isInterfaceEdge[31] = interfaceSurface[0];
		isInterfaceEdge[32] = interfaceSurface[1];
		isInterfaceEdge[33] = interfaceSurface[0];
		isInterfaceEdge[34] = false;
		isInterfaceEdge[35] = interfaceSurface[3];

		predecessorElement[30] = -1;
		predecessorElement[31] = -1;
		predecessorElement[32] = -1;
		predecessorElement[33] = -1;
		predecessorElement[34] = -1;
		predecessorElement[35] = -1;


		// Subelements of thetrahedron
		newTriangles[20]= {midPointInd[0],midPointInd[1],midPointInd[3]};
		newTriangles[21]= {midPointInd[0],midPointInd[1],midPointInd[4]};
		newTriangles[22]= {midPointInd[0],midPointInd[3],midPointInd[4]};
		newTriangles[23]= {midPointInd[1],midPointInd[3],midPointInd[4]};

		newTrianglesFlag[20]= originFlag[0]; 
		newTrianglesFlag[21]= 10; 
		newTrianglesFlag[22]= 10; 
		newTrianglesFlag[23]= 10;

		// Element 7: (x_02,x_03,x_13,x_23)
		(newElements)[6]={midPointInd[1],midPointInd[2],midPointInd[4],midPointInd[5]};

		(newEdges)[36] = {midPointInd[1] ,midPointInd[2]}; 
		(newEdges)[37] = {midPointInd[1] ,midPointInd[4]}; 
		(newEdges)[38] = {midPointInd[1] ,midPointInd[5]}; 
		(newEdges)[39] = {midPointInd[2] ,midPointInd[4]}; 
		(newEdges)[40] = {midPointInd[2] ,midPointInd[5]}; 
		(newEdges)[41] = {midPointInd[4] ,midPointInd[5]}; 

		edgeFlags[36]=originFlag[2];
		edgeFlags[37]=10;
		edgeFlags[38]=originFlag[2];
		edgeFlags[39]=originFlag[1];
		edgeFlags[40]=originFlag[2];
		edgeFlags[41]=originFlag[3];

		isInterfaceEdge[36] = interfaceSurface[2];
		isInterfaceEdge[37] = false;
		isInterfaceEdge[38] = interfaceSurface[2];
		isInterfaceEdge[39] = interfaceSurface[1];
		isInterfaceEdge[40] = interfaceSurface[2];
		isInterfaceEdge[41] = interfaceSurface[3];

		predecessorElement[36] = -1;
		predecessorElement[37] = -1;
		predecessorElement[38] = -1;
		predecessorElement[39] = -1;
		predecessorElement[40] = -1;
		predecessorElement[41] = -1;

		// Subelements of thetrahedron
		newTriangles[24]= {midPointInd[1],midPointInd[2],midPointInd[4]};
		newTriangles[25]= {midPointInd[1],midPointInd[2],midPointInd[5]};
		newTriangles[26]= {midPointInd[1],midPointInd[4],midPointInd[5]};
		newTriangles[27]= {midPointInd[2],midPointInd[4],midPointInd[5]};

		newTrianglesFlag[24]= 10; 
		newTrianglesFlag[25]= originFlag[2]; 
		newTrianglesFlag[26]= 10; 
		newTrianglesFlag[27]= 10;

		// Element 8: (x_02,x_12,x_13,x_23)
		(newElements)[7]={midPointInd[1],midPointInd[3],midPointInd[4],midPointInd[5]};

		(newEdges)[42] = {midPointInd[1] ,midPointInd[3]}; 
		(newEdges)[43] = {midPointInd[1] ,midPointInd[4]}; 
		(newEdges)[44] = {midPointInd[1] ,midPointInd[5]}; 
		(newEdges)[45] = {midPointInd[3] ,midPointInd[4]}; 
		(newEdges)[46] = {midPointInd[3] ,midPointInd[5]}; 
		(newEdges)[47] = {midPointInd[4] ,midPointInd[5]}; 

		edgeFlags[42]=originFlag[0];
		edgeFlags[43]=10;
		edgeFlags[44]=originFlag[2];
		edgeFlags[45]=originFlag[3];
		edgeFlags[46]=originFlag[3];
		edgeFlags[47]=originFlag[3];

		isInterfaceEdge[42] = interfaceSurface[0];
		isInterfaceEdge[43] = false; 
		isInterfaceEdge[44] = interfaceSurface[2];
		isInterfaceEdge[45] = interfaceSurface[3];
		isInterfaceEdge[46] = interfaceSurface[3];
		isInterfaceEdge[47] = interfaceSurface[3];

		predecessorElement[42] = -1;
		predecessorElement[43] = -1;
		predecessorElement[44] = -1;
		predecessorElement[45] = -1;
		predecessorElement[46] = -1;
		predecessorElement[47] = -1;

		// Subelements of thetrahedron
		newTriangles[28]= {midPointInd[1],midPointInd[3],midPointInd[4]};
		newTriangles[29]= {midPointInd[1],midPointInd[3],midPointInd[5]};
		newTriangles[30]= {midPointInd[1],midPointInd[4],midPointInd[5]};
		newTriangles[31]= {midPointInd[3],midPointInd[4],midPointInd[5]};

		newTrianglesFlag[28]= 10; 
		newTrianglesFlag[29]= 10; 
		newTrianglesFlag[30]= 10; 
		newTrianglesFlag[31]= originFlag[3];

		
		// Now we add the elements, edges and triangles 

		// Adding Elements
		int offsetElements = this->elementsC_->numberElements();
		int offsetEdges = this->edgeElements_->numberElements();
		for( int i=0;i<8; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );

			FiniteElement feNew(newElements.at(i),10);

			feNew.setFiniteElementRefinementType("regular");	
			if(i<7)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		// Adding the edges (they also have to be added to triangles as subelements, but that is not implemented yet)
		for( int i=0;i<48; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<42){
				this->edgeElements_->addEdge(feNew,i/6+offsetElements);
			}
			else
				this->edgeElements_->addEdge(feNew,indexElement);
		}

		// Adding triangles as subelements, if they arent interior triangles
		for( int i=0;i<32; i++){
			sort( newTriangles.at(i).begin(), newTriangles.at(i).end() );
			FiniteElement feNew(newTriangles[i],newTrianglesFlag[i]);
			if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=10){
				if(i<28){
				 	if ( !this->elementsC_->getElement(i/4+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/4+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/4+offsetElements).addSubElement(feNew);
				}
				else{
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);
				}
			}						
		}
		// tagging the Edges for refinement
	   	for (int d=0; d<6; d++){
			// now we tag the edge as 'refined'
			edgeElements->getElement(edgeNumbers[d]).tagForRefinement();	   
		}
	
	}
	else 
       	TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The regular Refinement Algorithm is only available in 2 and 3 Dimensions" );
  }
}
#endif
