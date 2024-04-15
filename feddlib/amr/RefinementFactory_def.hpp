#ifndef RefinementFactory_def_hpp
#define RefinementFactory_def_hpp

#ifndef MESH_TIMER_START
#define MESH_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Mesh Refinement") + std::string(S))));
#endif

#ifndef MESH_TIMER_STOP
#define MESH_TIMER_STOP(A) A.reset();
#endif

#include "RefinementFactory_decl.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector_def.hpp"
#include <chrono> 
/*!
 Definition of RefinementFactory
 
  \brief  RefinementFactory
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
RefinementFactory<SC,LO,GO,NO>::RefinementFactory():
MeshUnstructured<SC,LO,GO,NO>()
{

  
}


/*!
\brief Initiating RefinementFactory via MeshUnstructured.

@param[in] comm CommPtr.
@param[in] volumeID The flag ID of triangles in 2D or tetrahedra in 3D. Usually 10.

*/

template <class SC, class LO, class GO, class NO>
RefinementFactory<SC,LO,GO,NO>::RefinementFactory(CommConstPtr_Type comm, int volumeID):
MeshUnstructured<SC,LO,GO,NO>(comm,volumeID)
{

}

/*!
\brief Initiating RefinementFactory via MeshUnstructured with additional information for mesh refinement

@param[in] comm CommPtr
@param[in] volumeID The flag ID of triangles in 2D or tetrahedra in 3D. Usually 10.
@param[in] refinementRestriction Restriction for repeated refinement steps.
@param[in] refinement3DDiagonal 3D diagonal pick for regular refinement

*/

template <class SC, class LO, class GO, class NO>
RefinementFactory<SC,LO,GO,NO>::RefinementFactory(CommConstPtr_Type comm, int volumeID, ParameterListPtr_Type parameterListAll):
MeshUnstructured<SC,LO,GO,NO>(comm,volumeID)
{
    this->volumeID_ = volumeID;
	this->dim_ = parameterListAll->sublist("Parameter").get("Dimension",2);
	if(this->dim_ == 2){
		refinementRestriction_ = parameterListAll->sublist("Mesh Refinement").get("Refinement Restriction","Bisection");
		writeRefinementTime_ = parameterListAll->sublist("Mesh Refinement").get("Write Refinement Time",true);
	}
	else{
		refinementRestriction_ = parameterListAll->sublist("Mesh Refinement").get("Refinement Restriction","BeyIrregular");
		refinement3DDiagonal_ = parameterListAll->sublist("Mesh Refinement").get("3D regular Refinement Diagonal Pick",0);
		writeRefinementTime_ = parameterListAll->sublist("Mesh Refinement").get("Write Refinement Time",true);
	}

}

template <class SC, class LO, class GO, class NO>
RefinementFactory<SC,LO,GO,NO>::~RefinementFactory(){

}

/*!

 \brief Main function of RefinementFactory, performs one complete mesh refinement, according to red-green refinement (Verfuerth) or tetrahedral grid refinement (Bey). 

@param[in] meshP1 InputMesh P1.
@param[in] iteration Current Iteration.
@param[in] refinementMode In 2D we can choose between 'Regular' (red-green) or 'Bisection. In 3D only 'Regular' is possible, which is also the default mode.

@param[out] outputMesh Refined mesh.

*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineMesh( MeshUnstrPtr_Type meshP1, int iteration, MeshUnstrPtr_Type outputMesh, string refinementMode){
	MESH_TIMER_START(totalTime," Total Time for Mesh Refinement of this Step ");		


	if(meshP1->FEType_ != "P1" && meshP1->FEType_ != "P2"){ 
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Mesh Refinement only works for Triangular Elements");
	}
	
    currentIter_ = iteration;

	refinementMode_ = refinementMode;
	if(refinementMode_ == "Bisection")
		this->refinementRestriction_ = "Bisection";

	this->dim_ = meshP1->getDimension();
	this->FEType_ =meshP1->FEType_;
    this->volumeID_ = meshP1->volumeID_;
	this->rankRange_=meshP1->rankRange_;
	// necessary entities
    //EdgeElementsPtr_Type edgeElementsTmp = meshP1->getEdgeElements(); // Edges
    SurfaceElementsPtr_Type surfaceTriangleElements = meshP1->getSurfaceTriangleElements(); // Surfaces
    ElementsPtr_Type elementsTmp = meshP1->getElementsC(); // Elements

	EdgeElementsPtr_Type edgeElements =meshP1->getEdgeElements(); //Teuchos::rcp( new EdgeElements(*edgeElementsTmp));
	
	//EdgeElementsPtr_Type edgeElements=Teuchos::rcp( new EdgeElements(*edgeElementsTmp));
	ElementsPtr_Type elements =Teuchos::rcp( new Elements(*elementsTmp));

    vec2D_dbl_ptr_Type points = meshP1->getPointsRepeated(); // Points
    this->elementMap_ = meshP1->elementMap_;
	this->mapUnique_ = meshP1->mapUnique_;
	this->mapRepeated_=meshP1->mapRepeated_;
	this->edgeMap_ = meshP1->edgeMap_;
	// entities for resulting Elements, Points, Edges, Flags
	// - Points, Elements and Flags are added to the existding entities
	// . edges are reset every refinement step
    this->edgeElements_.reset(new EdgeElements()); // Edges
    this->surfaceTriangleElements_.reset(new SurfaceElements()); // Surface
    this->elementsC_.reset(new Elements(*elementsTmp));    // Elements 
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

	// not all Edges are assigned a Flag yet -> in Order to save a few steps along the way, we assign all Edges the correct flag now
	// In regular refinement the correct flags are then passed on, so this step is only necessary in Iteration 0
	if( iteration ==0 ){
		meshP1->assignEdgeFlags();
	}

	// Algorithm
	if (this->dim_ == 2 || this->dim_ == 3){


		if(this->comm_->getRank() == 0){
			cout << " 	-- Mesh Refinement -- " << endl;
			cout << "	__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " 	Start Iteration " << iteration+1 << " of "<< this->dim_ << "D Mesh Refinement " << endl;
			cout << " 	Number of Elements:	" << this->elementMap_->getGlobalNumElements() << endl;
			cout << " 	Number of Nodes:	" << this->mapUnique_->getGlobalNumElements() << endl; 
			cout << " 	Number of Edges:	" << this->edgeMap_->getGlobalNumElements() << endl;
			cout << "	__________________________________________________________________________________________________________ " << endl;
		}



		// ------------------------------------------------------------------------------------------------------
		// Part I: Regular Refinement of Mesh
		// We refine the elements that were determined by our error estimation regular
		// ------------------------------------------------------------------------------------------------------	
		MESH_TIMER_START(preprocessingTimer," Step 0:	 Preprocessing");
		const int myRank = this->comm_->getRank();
		// match the Edges to the Elements for being able to tag the edges of Elements in case of refinement
		edgeElements->matchEdgesToElements(this->elementMap_);

		// Vector that carry the information which elements belong to an edge in local and global indices
		vec2D_GO_Type elementsOfEdgeGlobal = edgeElements->getElementsOfEdgeGlobal();
		vec2D_LO_Type elementsOfEdgeLocal = edgeElements->getElementsOfEdgeLocal();
		

		// Determine current global Interface IDs
		// As we pass on the interface boolean while refining red,blue,green... we theoretically dont need the whole extend of this function passed the first 
		// refinement, only the globalInterfaceIDs, which can be determined without 'elementsOfEdgeLocal' and 'Global'
		// (If elementsOfEdgeLocal and Global arent working, this can be taken out, yet is also nice for determining if it works correctly)
		vec_GO_Type globalInterfaceIDs;
		if(iteration ==0 ){
			globalInterfaceIDs = edgeElements->determineInterfaceEdges(this->edgeMap_);
			for(int i=0; i< elements->numberElements(); i++){
				this->elementsC_->getElement(i).setPredecessorElement(i);
			}

		}
		else {
			for(int i=0; i< edgeElements->numberElements(); i++){
				if(edgeElements->getElement(i).isInterfaceElement()){
					globalInterfaceIDs.push_back(this->edgeMap_->getGlobalElement(i));
				}
			}	
			sort(globalInterfaceIDs.begin(), globalInterfaceIDs.end());			
		}
		globalInterfaceIDs_ = globalInterfaceIDs;

		
		if(this->dim_==3){
			if(surfaceTriangleElements.is_null()){
				surfaceTriangleElements.reset(new SurfaceElements()); // Surface
				this->buildSurfaceTriangleElements(elements,edgeElements, surfaceTriangleElements, this->edgeMap_, this->elementMap_ );
				//this->buildTriangleMap();
			}
			else if(surfaceTriangleElements->numberElements() ==0){
				this->buildSurfaceTriangleElements(elements,edgeElements, surfaceTriangleElements, this->edgeMap_, this->elementMap_ );
			} 
			surfaceTriangleElements->matchSurfacesToElements(this->elementMap_);
		}

		// counting new Points and Elements:
		int newPoints=0; // total number of new Points
		int newPointsRepeated= 0; // number of new Points that share an other Processor

		int newElements=0;	// new Elements on a Processor

		// Counting new Edges
		int newEdgesUnique=0;
		int newEdgesRepeated =0;
		
		// Loop through Elements in order to determine Elements to refine and refine them regular / red 
		int numPoints=0;
		MESH_TIMER_STOP(preprocessingTimer);

		MESH_TIMER_START(regularRefinementTimer," Step 1:	 Tagging Edges for Refinement");
		
		// Depending on dimension we add a certain number of elements in the 2D case it's 3, in the 3D it's 7		
		for(int i=0; i<elements->numberElements();i++){
			if(elements->getElement(i).isTaggedForRefinement()){
				numPoints= this->pointsRep_->size();
				this->bisectEdges( edgeElements, elements, i, surfaceTriangleElements);
				newPoints=newPoints + this->pointsRep_->size()-numPoints;
			}
		}

		MESH_TIMER_STOP(regularRefinementTimer);
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part II: Communicating the tagged Edges 
		// As it is possible now, that edges were tagged on one processor on the interface between processers
		// we need communicate tagged edges across the interface
		// In order to safe time, we only communicate tagged interface Edges
		// InterfaceEdges can be determined by the vector 'elementsOfEdgesLocal' as the vector carries a -1 for 
		// for each element belonging to an edge that is not on the processor in question
		// ------------------------------------------------------------------------------------------------------
		MESH_TIMER_START(commEdgesTimer," Step 2:	 Communicating tagged edges");
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
		//taggedEdgesGlobal->print();

		// Adding Midpoints and tagging the edges that were tagged on other Procs that are on the interface
		LO ind;
		for (int i=0; i<tags.size(); i++) {
			if (tags[i] > 0){
				ind = edgeMap->getLocalElement(globalInterfaceIDs[i]);
				newPointsRepeated ++;
				if(!edgeElements->getElement(ind).isTaggedForRefinement()){
					edgeElements->getElement(ind).tagForRefinement();
					this->addMidpoint(edgeElements,ind);
					// Collecting global IDs we didnt already considered
					globalInterfaceIDsTagged.push_back(globalInterfaceIDs[i]);
					newPoints ++;
					}		
			}
		}
		MESH_TIMER_STOP(commEdgesTimer);
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part III: Checking for Refinement Restrictions
		// Before we start refining the elements according to their refinement rules, we need to check certain
		// refinement restrictions
		// ------------------------------------------------------------------------------------------------------

		MESH_TIMER_START(checkTimer," Step 3:	 Checking Restrictions");		
		this->refinementRestrictions(meshP1, elements ,edgeElements, surfaceTriangleElements, newPoints, newPointsRepeated, globalInterfaceIDsTagged, mapInterfaceEdges, newElements);

		sort(globalInterfaceIDsTagged.begin(), globalInterfaceIDsTagged.end());

		MESH_TIMER_STOP(checkTimer);		


		// ------------------------------------------------------------------------------------------------------
		// Part IV: Communicating the added Points and updating the corresponding Maps
		// generally we keep the nodelist of throughout the refinement steps and only add new points to it
		// consequently we update the maps for those new points depending on nodes on interface and points unqiuely
		// on processors
		// ------------------------------------------------------------------------------------------------------
		MESH_TIMER_START(nodeTimer," Step 4:	 Updating Node Map");		
		// determine global interface IDs of untagged edges 	
		int maxRank = std::get<1>(this->rankRange_);
		// determine unique map

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
		
		this->buildNodeMap(edgeElements, mapGlobalProc, mapProc, newPoints, newPointsRepeated);

		MESH_TIMER_STOP(nodeTimer);		

		// ------------------------------------------------------------------------------------------------------


		// ------------------------------------------------------------------------------------------------------		
		// Part V: refineMesh
		// All edges are tagged and checked for additional restricitions. The mesh can now be refined and regular
		// and irregular refinement rules can be performed	
		// ------------------------------------------------------------------------------------------------------
		MESH_TIMER_START(irregRefTimer," Step 5:	 Irregular Refinement");		
		this->refineMeshRegIreg(elements, edgeElements, newElements,edgeMap, surfaceTriangleElements);
		MESH_TIMER_STOP(irregRefTimer);		

		// ------------------------------------------------------------------------------------------------------
		// Part VI: Updating the Element Map
		// update elementMap by distributing information of locally added elements and add global IDs based on
		// that information
		// the list of elements only extends throughout the refinement steps, so we only neew to allocate the 
		// globalIDs for elements that extend our last element number
		// ------------------------------------------------------------------------------------------------------
		MESH_TIMER_START(elementMapTimer," Step 6:	 Updating Element Map");		
		MapConstPtr_Type elementMap = meshP1->getElementMap();
		// information of new elements
		MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );
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
		this->numElementsGlob_ = this->elementMap_->getMaxAllGlobalIndex()+1;  
		MESH_TIMER_STOP(elementMapTimer);		
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part VII: Making the Edges unique, set global IDs and set elements edges
		// we need to make the edge list unique as there are redundant edges
		// furthermore we set up the elementsOfEdgeLocal and elementsOfEdgeGlobal, but the information of other
		// procs is still missing there (this will be finalized in Part X)
		// ------------------------------------------------------------------------------------------------------

		MESH_TIMER_START(uniqueEdgesTimer," Step 7:	 Making Edges Unique");		
		vec2D_GO_Type combinedEdgeElements;
		this->edgeElements_->sortUniqueAndSetGlobalIDsParallel(this->elementMap_,combinedEdgeElements);
		MESH_TIMER_STOP(uniqueEdgesTimer);		
		//this->comm_->barrier();
		// ------------------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------------------
		// Part VIII: Updating the EdgeMap 
		// The general idea is to indentifiy an edge by their global node indices. 
		// This way two edges on different processors can be indentified and given the same global ID
		// First we determine the global IDs of non interface edges 
		// Then we determine those of interface edges with the procedure discribed above
		// ------------------------------------------------------------------------------------------------------

		MESH_TIMER_START(edgeMapTimer," Step 8:	 Creating EdgeMap");		
		this->buildEdgeMap(mapGlobalProc, mapProc);
		MESH_TIMER_STOP(edgeMapTimer);		

		// ------------------------------------------------------------------------------------------------------
		// Part IX: Updating elementsOfEdgeGlobal and elementsOfEdgeLocal
		// this step is only necessary if we have more than 1 processor, as the edgeElements function work serially
		// we started the setup before (sortUniqueAndSetGlobalIDsParallel) an now finalize it with the information of other processors
		// the edges on the interface need the global element number of the neighbouring processor
		// ------------------------------------------------------------------------------------------------------
		MESH_TIMER_START(elementsOfEdgeTimer," Step 9:	 Updating ElementsOfEdgeLocal/Global");		
		this->edgeElements_->setElementsEdges( combinedEdgeElements );

		this->edgeElements_->setUpElementsOfEdge( this->elementMap_, this->edgeMap_);

		this->updateElementsOfEdgesLocalAndGlobal(maxRank, edgeMap);

		MESH_TIMER_STOP(elementsOfEdgeTimer);	

		MESH_TIMER_START(elementsOfSurfaceTimer," Step 10: Updating ElementsOfSurfaceLocal/Global");
		
		vec2D_GO_Type combinedSurfaceElements;

		this->surfaceTriangleElements_->sortUniqueAndSetGlobalIDsParallel(this->elementMap_,combinedSurfaceElements);
		
		this->surfaceTriangleElements_->setElementsSurface( combinedSurfaceElements );

		this->surfaceTriangleElements_->setUpElementsOfSurface( this->elementMap_, this->edgeMap_, this->edgeElements_);

		MESH_TIMER_STOP(elementsOfSurfaceTimer);	

		vec2D_GO_Type elementsOfSurfaceGlobal = this->surfaceTriangleElements_->getElementsOfSurfaceGlobal();

		vec2D_LO_Type elementsOfSurfaceLocal = this->surfaceTriangleElements_->getElementsOfSurfaceLocal();

	
		MESH_TIMER_STOP(totalTime);		

		
		if(this->comm_->getRank() == 0){
			cout << "	__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " 	... finished Iteration " << iteration+1 << " of " << this->dim_ << "D Mesh Refinement " << endl;
			cout << " 	Number of new Elements:	" << this->elementMap_->getGlobalNumElements() - meshP1->elementMap_-> getGlobalNumElements() << endl;
			cout << " 	Number of new Nodes:	" << this->mapUnique_->getGlobalNumElements()- meshP1->mapUnique_-> getGlobalNumElements() << endl; 
			cout << " 	Number of new Edges:	" << this->edgeMap_->getGlobalNumElements()- meshP1->edgeMap_-> getGlobalNumElements() << endl;
			cout << "	__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
		}

		if(writeRefinementTime_ )
   			Teuchos::TimeMonitor::report(cout,"Mesh Refinement");
		
	// Finally we set all changed mesh enteties for outputMesh

	
	outputMesh->dim_ = this->dim_ ;
	outputMesh->FEType_ = this->FEType_ ;
	outputMesh->rankRange_ =  this->rankRange_;

    outputMesh->elementMap_ = this->elementMap_ ;
	outputMesh->mapUnique_ = this->mapUnique_  ;
	outputMesh->mapRepeated_ = this->mapRepeated_;
	outputMesh->edgeMap_  = this->edgeMap_  ;

	outputMesh->elementsC_ = this->elementsC_;
	outputMesh->edgeElements_ = this->edgeElements_;
	outputMesh->surfaceTriangleElements_ = this->surfaceTriangleElements_;

   	outputMesh->pointsRep_ =  this->pointsRep_  ; 
    outputMesh->pointsUni_ = this->pointsUni_; 

    outputMesh->bcFlagUni_ = this->bcFlagUni_ ; 
	outputMesh->bcFlagRep_ = this->bcFlagRep_ ;


	outputMesh->edgesElementOrder_ = this->edgesElementOrder_;
	outputMesh->numElementsGlob_ = this->numElementsGlob_  ; 

	}


}
/*!

 \brief Checking if surfaces are part of the interface. 
Done by checking if all edges of a triangle are part of the interface and if both elements connected to the surface are on different processors.

@param[in] edgeElements Edges.
@param[in] originFlag Flags of surfaces of indexElement.
@param[in] edgeNumbers Numbers of inserted surfaces edges.
@param[in] indexElement Index of element in question.

*/

template <class SC, class LO, class GO, class NO>
vec_bool_Type RefinementFactory<SC,LO,GO,NO>::checkInterfaceSurface( EdgeElementsPtr_Type edgeElements,vec_int_Type originFlag, vec_int_Type edgeNumbers, int indexElement){

		vec_bool_Type interfaceSurface(4);
		
		if(edgeElements->getElement(edgeNumbers[0]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[1]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[3]).isInterfaceElement() && (originFlag[0] == this->volumeID_))
					interfaceSurface[0] = true;
		if(edgeElements->getElement(edgeNumbers[0]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[2]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[4]).isInterfaceElement() && (originFlag[1] == this->volumeID_ ))
					interfaceSurface[1] = true;
		if(edgeElements->getElement(edgeNumbers[1]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[2]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[5]).isInterfaceElement() && (originFlag[2] == this->volumeID_))
					interfaceSurface[2] = true;
		if(edgeElements->getElement(edgeNumbers[3]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[4]).isInterfaceElement() && edgeElements->getElement(edgeNumbers[5]).isInterfaceElement() && (originFlag[3] == this->volumeID_ ))
					interfaceSurface[3] = true;

		vec2D_LO_Type elementsOfEdgeLocal = edgeElements->getElementsOfEdgeLocal();

		vec2D_LO_Type edgeNumTri(4,vec_LO_Type(3));
		edgeNumTri[0] = {edgeNumbers[0],edgeNumbers[1],edgeNumbers[3]};
		edgeNumTri[1] = {edgeNumbers[0],edgeNumbers[2],edgeNumbers[4]};
		edgeNumTri[2] = {edgeNumbers[1],edgeNumbers[2],edgeNumbers[5]};
		edgeNumTri[3] = {edgeNumbers[3],edgeNumbers[4],edgeNumbers[5]};
		

		// We need to check additionally, if maybe all edges are part of the interface, but not the triangle itself.
		// If we find the element and the neighbouring element that share the triangle on our processor, the triangle is not part of the interface.
		for(int i=0; i<4 ; i++){
			if(interfaceSurface[i] == true){
				vec_LO_Type tmpElements(0);
				for(int j=0 ; j< elementsOfEdgeLocal[edgeNumTri[i][0]].size() ; j++){
					if(elementsOfEdgeLocal[edgeNumTri[i][0]][j] != -1)
						tmpElements.push_back( elementsOfEdgeLocal[edgeNumTri[i][0]][j]);
				}

				for(int j=0 ; j< elementsOfEdgeLocal[edgeNumTri[i][1]].size() ; j++){
					if(elementsOfEdgeLocal[edgeNumTri[i][1]][j] != -1)
						tmpElements.push_back( elementsOfEdgeLocal[edgeNumTri[i][1]][j]);
				}


				sort(tmpElements.begin(),tmpElements.end());
				bool found =false;

				for(int j=0; j< tmpElements.size()-1; j++){
					if((tmpElements[j] == tmpElements[j+1] )&& (tmpElements[j] != indexElement) ) { 
						found = true;
					}
				}
				if(found == true)
					interfaceSurface[i] = false;


			}
		}

		return interfaceSurface;
}

/*!

 \brief Building nodemap after refinement.

@param[in] edgeElements Edges.
@param[in] mapGlobalProc Map of global processor numbers.
@param[in] mapProc Map of local processor numbers.
@param[in] newPoints Number of new points per refinement iteration.
@param[in] newPointsRepeated Number of new repeated points per refinement iteration.


*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::buildNodeMap(EdgeElementsPtr_Type edgeElements, MapConstPtr_Type mapGlobalProc, MapConstPtr_Type mapProc, int newPoints, int newPointsRepeated){
		int maxRank = std::get<1>(this->rankRange_);
		const int myRank = this->comm_->getRank();

		// Collecting the number of new Points and repeated new Points in total on each Proc
		int globalCountPoints = 0; 
		int globalCountRepeated = 0;
		reduceAll<int, int> (*this->comm_, REDUCE_SUM, newPoints, outArg (globalCountPoints));
		reduceAll<int, int> (*this->comm_, REDUCE_SUM, newPointsRepeated, outArg (globalCountRepeated));
		// Setting Unique newPoints as value for communication
		MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );

		// Unique Points on each Proc
		int newPointsUnique= newPoints - newPointsRepeated;
		exportLocalEntry->putScalar( (LO) newPointsUnique);

		MultiVectorLOPtr_Type isActiveNodeUnique= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		isActiveNodeUnique->putScalar( (LO) 0 ); 
		isActiveNodeUnique->importFromVector( exportLocalEntry, true, "Insert");
		// -> Result: Information on uniquely added Points of each Proc is accessible to all

		Teuchos::ArrayRCP<const LO> nodesUniqueProcList = isActiveNodeUnique-> getData(0);
			
		// Step 1: Add Nodes that are repeated on the processors
		// Generally we allocate the first 'n' Points to Proc 0, the following start with (n+1) on Proc 1 and so on -> ProcOffset

		GO refineOffset = this->mapUnique_->getMaxAllGlobalIndex()+1; // Number of Nodes before Refinement
		int procOffset=0;
		for(int i=0; i< myRank; i++)
			procOffset = procOffset + nodesUniqueProcList[i];

		MapConstPtr_Type mapRepeatedP1 = this->getMapRepeated(); // Maps
		Teuchos::ArrayView<const GO> nodeList = mapRepeatedP1->getNodeElementList();
		vec_GO_Type vecGlobalIDsOld = Teuchos::createVector( nodeList );
	  
		vec_GO_Type vecGlobalIDs(this->pointsRep_->size(),-1);

		for (int i=0; i<vecGlobalIDsOld.size(); i++){
			vecGlobalIDs[i] = vecGlobalIDsOld[ i] ; // (i + refineOffset+procOffset);
		}

		// Add Nodes that are on the interface -> repeated Points
		int sumPointsUnique = globalCountPoints - globalCountRepeated;
	

		// First we determine a Map only for the interface Nodes
		// This will reduce the size of the Matrix we build later significantly if only look at the interface edges
		int numEdges= edgeElements->numberElements();
		vec2D_GO_Type inzidenzIndices(0,vec_GO_Type(2)); // Vector that stores global IDs of each edge (in Repeated Sense)
		vec_LO_Type localNodeIndex(0); // stores the local ID of node in question 
		vec_GO_Type id(2);


		int interfaceNum=0;
		for(int i=0; i<numEdges; i++ ){
			if(edgeElements->getElement(i).isInterfaceElement() && edgeElements->getElement(i).isTaggedForRefinement()){
				id[0] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(0)); 
				id[1] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(1));
			 	
				sort(id.begin(),id.end());
				inzidenzIndices.push_back(id);

				localNodeIndex.push_back(edgeElements->getMidpoint(i));

			}
		 }


		// This Matrix is row based, where the row is based on mapInterfaceNodesUnqiue
		// We then add a '1' Entry when two global Node IDs form an edge
		MatrixPtr_Type inzidenzMatrix = Teuchos::rcp( new Matrix_Type(this->mapUnique_, 40 ) );
		Teuchos::Array<GO> index(1);
		Teuchos::Array<GO> col(1);
		Teuchos::Array<SC> value(1, Teuchos::ScalarTraits<SC>::one() );

		for(int i=0; i<inzidenzIndices.size(); i++ ){
	
			index[0] = inzidenzIndices[i][0];
			col[0] = inzidenzIndices[i][1];
			inzidenzMatrix->insertGlobalValues(index[0], col(), value());
		
		 }
   		inzidenzMatrix->fillComplete(); 

		// Now we count the row entries on each processor an set global IDs

		Teuchos::ArrayView<const LO> indices;
		Teuchos::ArrayView<const SC> values;
		vec2D_GO_Type inzidenzIndicesUnique(0,vec_GO_Type(2)); // Vector that stores only both global IDs if the first is part of my unique Interface Nodes
		MapConstPtr_Type colMap = inzidenzMatrix->getMap("col");
		MapConstPtr_Type rowMap = inzidenzMatrix->getMap("row");
		int numRows = rowMap->getNodeNumElements();
		int uniqueEdges =0;
		for(int i=0; i<numRows; i++ ){
			inzidenzMatrix->getLocalRowView(i, indices,values); 
			uniqueEdges = uniqueEdges+indices.size();
			vec_GO_Type edgeTmp(2);
			for(int j=0; j<indices.size(); j++){
				edgeTmp[0] = rowMap->getGlobalElement(i);
				edgeTmp[1] = colMap->getGlobalElement(indices[j]);
				inzidenzIndicesUnique.push_back(edgeTmp);
			}
		}
	
		exportLocalEntry->putScalar( uniqueEdges );
		MultiVectorLOPtr_Type newEdgesInterfaceGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newEdgesInterfaceGlobal->putScalar( (LO) 0 ); 
		newEdgesInterfaceGlobal->importFromVector( exportLocalEntry, true, "Insert");

		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > numUniqueInterface = newEdgesInterfaceGlobal->getData(0);

		int procOffsetEdges=0;
		for(int i=0; i< myRank; i++)
			procOffsetEdges= procOffsetEdges + numUniqueInterface[i];

		int numInterfaceEdges=0;
	
		vec_GO_Type uniqueInterfaceIDsList_(inzidenzIndicesUnique.size());
		for(int i=0; i< uniqueInterfaceIDsList_.size(); i++)
			uniqueInterfaceIDsList_[i] = procOffsetEdges + i;

		MatrixPtr_Type indMatrix = Teuchos::rcp( new Matrix_Type(this->mapUnique_, 40 ) );

		for(int i=0; i<inzidenzIndicesUnique.size(); i++ ){
			index[0] = inzidenzIndicesUnique[i][0];
			col[0] = inzidenzIndicesUnique[i][1];
			Teuchos::Array<SC> value2(1,uniqueInterfaceIDsList_[i]);
			indMatrix->insertGlobalValues(index[0], col(), value2());
		 }
   		indMatrix->fillComplete(); 

		MatrixPtr_Type importMatrix = Teuchos::rcp( new Matrix_Type(this->mapRepeated_, 40 ) );
   		
		importMatrix->importFromVector(indMatrix,false,"Insert");
		importMatrix->fillComplete(); 

		// Determine global indices
		GO edgeID=0;
		colMap = importMatrix->getMap("col");
		rowMap = importMatrix->getMap("row");
	
		LO valueID=0;
		bool found = false;
		GO entry =0;
		for(int i=0; i<inzidenzIndices.size(); i++ ){	
			importMatrix->getLocalRowView(rowMap->getLocalElement(inzidenzIndices[i][0]), indices,values); // Indices and values connected to node i / row i in Matrix
			// Entries in 'indices' represent the local entry in 'colmap
			// with 'getGlobalElement' we know the global Node ID that belongs to the first Node that form an edge
			// vector in with entries only for edges belonging to node i;
			vec2D_GO_Type indicesTmp(indices.size(),vec_GO_Type(2));
			vec_GO_Type indTmp(2);
			for(int j=0; j<indices.size(); j++){
				indTmp[0] = colMap->getGlobalElement(indices[j]);
				indTmp[1] = values[j];
				indicesTmp.push_back(indTmp);	// vector with the indices and values belonging to node i
			}
			found = false;
			for(int k=0; k<indicesTmp.size();k++){
				if(inzidenzIndices[i][1] == indicesTmp[k][0]){
					entry =k;
					k = indicesTmp.size();
					edgeID = indicesTmp[entry][1];
					vecGlobalIDs[localNodeIndex[i]] = refineOffset+sumPointsUnique + edgeID;
					found =true;
				}
			}
			if(found == false)
				cout << " Asking for row " << rowMap->getLocalElement(inzidenzIndices[i][0]) << " for Edge [" << inzidenzIndices[i][0] << ",  " << inzidenzIndices[i][1] << "], on Proc " << myRank << " but no Value found " <<endl;
		 }
	// --------------------------------------------------------------------------------------------------
	// --------------------------------------------------------------------------------------------------


		// Step 2: Add points that are uniquely on each processor, leftover -1 entries in vecGlobalIDs signal a missing unique entry
		// all left over entries 
		int count =0;
	  	for (int i= vecGlobalIDsOld.size(); i < this->pointsRep_->size() ; i++){
			if(vecGlobalIDs[i] == -1){
				vecGlobalIDs[i] = count + refineOffset+procOffset;
				count ++;
			}
		}
		Teuchos::RCP<std::vector<GO> > pointsRepGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDs ) );
		Teuchos::ArrayView<GO> pointsRepGlobMappingArray = Teuchos::arrayViewFromVector( *pointsRepGlobMapping );
		
		this->mapRepeated_.reset(new Map<LO,GO,NO>( this->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), pointsRepGlobMappingArray, 0, this->comm_) );
		this->mapUnique_ = this->mapRepeated_->buildUniqueMap( this->rankRange_ );
	
	
		// Points and Flags Unique
		this->pointsUni_.reset(new std::vector<std::vector<double> >( this->mapUnique_->getNodeNumElements(), vector<double>(this->dim_,-1. ) ) );
		this->bcFlagUni_.reset( new std::vector<int> ( this->mapUnique_->getNodeNumElements(), 0 ) );
		for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
			GO gid = this->mapUnique_->getGlobalElement( i );
			LO id = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement( i ) );
			this->pointsUni_->at(i) = this->pointsRep_->at(id);
			this->bcFlagUni_->at(i) = this->bcFlagRep_->at(id);
		}



}

/*!

 \brief Building surface triangle elements, as they are not originally part of the mesh information provided by mesh partitioner.

@param[in] elements Elements.
@param[in] edgeElements Edges.
@param[in] surfaceTriangleElements Pointer which will be filled with surfaceTriangleElements.
@param[in] edgeMap Global Mapping of edges
@param[in] elementMap Global Mapping of elements.


*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::buildSurfaceTriangleElements(ElementsPtr_Type elements, EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceTriangleElements, MapConstPtr_Type edgeMap, MapConstPtr_Type elementMap ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements.is_null(), std::runtime_error, "Elements not initialized!");
    TEUCHOS_TEST_FOR_EXCEPTION( surfaceTriangleElements.is_null(), std::runtime_error, "Surface Triangle Elements not initialized!");


	vec_LO_Type nodeInd(4);
	vec2D_int_Type newTriangles(4,vec_int_Type(0)); // new Triangles
	
	vec_GO_Type globalInterfaceIDs = edgeElements->determineInterfaceEdges(edgeMap);
	//if(edgeElements->getEdgesOfElement(0) ) here we need some sort of test if the function was already used
	edgeElements->matchEdgesToElements(elementMap);

	for(int T=0; T<elements->numberElements(); T++){
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(T); // indeces of edges belonging to element

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

		// We check if one or more of these triangles are part of the boundary surface and determine there flag

		vec2D_int_Type originTriangles(4,vec_int_Type(3));
		originTriangles[0] = {nodeInd[0],nodeInd[1],nodeInd[2]};
		originTriangles[1] = {nodeInd[0],nodeInd[1],nodeInd[3]};
		originTriangles[2] = {nodeInd[0],nodeInd[2],nodeInd[3]};
		originTriangles[3] = {nodeInd[1],nodeInd[2],nodeInd[3]};
		
		
		vec_int_Type originFlag(4,this->volumeID_); // Triangle Flag
		int numberSubElSurf=0;
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);
		int entry; 
		if (elements->getElement(T).subElementsInitialized() ){
			numberSubElSurf = elements->getElement(T).getSubElements()->numberElements();
			for(int k=0; k< numberSubElSurf ; k++){
				//cout << " Flags of Subelements " << elements->getElement(T).getSubElements()->getElement(k).getFlag() << endl;
				triTmp =elements->getElement(T).getSubElements()->getElement(k).getVectorNodeList();
				for(int j=0; j<4 ; j++){
					originTriangleTmp = originTriangles[j];
					sort(originTriangleTmp.begin(),originTriangleTmp.end());
					sort(triTmp.begin(),triTmp.end());
					if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] &&  triTmp[2] == originTriangleTmp[2] ) 
						originFlag[j] = elements->getElement(T).getSubElements()->getElement(k).getFlag();
				
				}
			}
		}

		// Furthermore we have to determine whether the triangles are part of the interface between processors, as we need this information to determine if edges
		// that emerge on the triangles are part of the interface
		// A triangle is part of the interface if all of its edges are part of the interface (the information if edges are part of the interface was determined
		// in the beginning of the Mesh Refinement by 'determineInterfaceEdges')

		vec_bool_Type interfaceSurface = checkInterfaceSurface(edgeElements,originFlag, edgeNumbers,T);
		
		for(int j=0; j<4; j++){	
			sort( newTriangles.at(j).begin(), newTriangles.at(j).end() );
			FiniteElement feNew(originTriangles[j],originFlag[j]);
			feNew.setInterfaceElement(interfaceSurface[j]);
			surfaceTriangleElements->addSurface(feNew, T);
		}
	}
	vec2D_GO_Type combinedSurfaceElements;
	surfaceTriangleElements->sortUniqueAndSetGlobalIDsParallel(elementMap,combinedSurfaceElements);
	
	surfaceTriangleElements->setElementsSurface( combinedSurfaceElements );

	surfaceTriangleElements->setUpElementsOfSurface(elementMap,edgeMap, edgeElements);	
}




/*!

 \brief Building edgeMap after refinement.

@param[in] mapGlobalProc Map of global processor numbers
@param[in] mapProc Map of local processor number

*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::buildEdgeMap(MapConstPtr_Type mapGlobalProc,MapConstPtr_Type mapProc){

		vec2D_int_Type interfaceEdgesLocalId(1,vec_int_Type(1));
		int maxRank = std::get<1>(this->rankRange_);
		const int myRank = this->comm_->getRank();

		MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );

		// (A) First we determine a Map only for the interface Nodes
		// This will reduce the size of the Matrix we build later significantly if only look at the interface edges
		int numEdges= this->edgeElements_->numberElements();
		vec2D_GO_Type inzidenzIndices(0,vec_GO_Type(2)); // Vector that stores global IDs of each edge (in Repeated Sense)
		vec_LO_Type localEdgeIndex(0); // stores the local ID of edges in question 
		vec_GO_Type id(2);
		int edgesUnique=0;
    	EdgeElementsPtr_Type edgeElements = this->edgeElements_; // Edges

		vec2D_dbl_ptr_Type points = this->pointsRep_;

		int interfaceNum=0;
		for(int i=0; i<numEdges; i++ ){
			if(edgeElements->getElement(i).isInterfaceElement()){

				id[0] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(0)); 
				id[1] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(1));
			 	


				sort(id.begin(),id.end());
				inzidenzIndices.push_back(id);

				localEdgeIndex.push_back(i);
				interfaceNum++;
			}
	
			else{
				edgesUnique++;
			}


		 }


		// This Matrix is row based, where the row is based on mapInterfaceNodesUnqiue
		// We then add a '1' Entry when two global Node IDs form an edge
		MatrixPtr_Type inzidenzMatrix = Teuchos::rcp( new Matrix_Type(this->mapUnique_, 40 ) );
		Teuchos::Array<GO> index(1);
		Teuchos::Array<GO> col(1);
		Teuchos::Array<SC> value(1, Teuchos::ScalarTraits<SC>::one() );

		for(int i=0; i<inzidenzIndices.size(); i++ ){
	
			index[0] = inzidenzIndices[i][0];
			col[0] = inzidenzIndices[i][1];
			inzidenzMatrix->insertGlobalValues(index[0], col(), value());
		
		 }
   		inzidenzMatrix->fillComplete(); //mapInterfaceNodesUnique,mapInterfaceNodesUnique);
		
	
		// ---------------------------------------------------
		// 2 .Set unique edges IDs ---------------------------
		// Setting the IDs of Edges that are uniquely on one
		// Processor
		// ---------------------------------------------------
		exportLocalEntry->putScalar( (LO) edgesUnique );

		MultiVectorLOPtr_Type newEdgesUniqueGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newEdgesUniqueGlobal->putScalar( (LO) 0 ); 
		newEdgesUniqueGlobal->importFromVector( exportLocalEntry, true, "Insert");
		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > newEdgesList = newEdgesUniqueGlobal->getData(0);

		GO procOffsetEdges=0;
		for(int i=0; i< myRank; i++)
			procOffsetEdges= procOffsetEdges + newEdgesList[i];

		// global IDs for map
		vec_GO_Type vecGlobalIDsEdges(this->edgeElements_->numberElements()); 
	
		// Step 1: adding unique global edge IDs
		int count=0;
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(!this->edgeElements_->getElement(i).isInterfaceElement()){
				vecGlobalIDsEdges.at(i) = procOffsetEdges+count;
				count++;
			}
		}	
		
		// Now we add the repeated ids, by first turning interfaceEdgesTag into a map
		// Offset for interface IDS:
		GO offsetInterface=0;
		for(int i=0; i< maxRank+1; i++)
			 offsetInterface=  offsetInterface + newEdgesList[i];
		
		//Now we count the row entries on each processor an set global IDs

		Teuchos::ArrayView<const LO> indices;
		Teuchos::ArrayView<const SC> values;
		vec2D_GO_Type inzidenzIndicesUnique(0,vec_GO_Type(2)); // Vector that stores only both global IDs if the first is part of my unique Interface Nodes
		MapConstPtr_Type colMap = inzidenzMatrix->getMap("col");
		MapConstPtr_Type rowMap = inzidenzMatrix->getMap("row");
		int numRows = rowMap->getNodeNumElements();
		int uniqueEdges =0;
		for(int i=0; i<numRows; i++ ){
			inzidenzMatrix->getLocalRowView(i, indices,values); 
			uniqueEdges = uniqueEdges+indices.size();
			vec_GO_Type edgeTmp(2);
			for(int j=0; j<indices.size(); j++){
				edgeTmp[0] = rowMap->getGlobalElement(i);
				edgeTmp[1] = colMap->getGlobalElement(indices[j]);
				inzidenzIndicesUnique.push_back(edgeTmp);
			}
		}
	
		exportLocalEntry->putScalar( uniqueEdges );
		MultiVectorLOPtr_Type newEdgesInterfaceGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newEdgesInterfaceGlobal->putScalar( (LO) 0 ); 
		newEdgesInterfaceGlobal->importFromVector( exportLocalEntry, true, "Insert");

		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > numUniqueInterface = newEdgesInterfaceGlobal->getData(0);

		procOffsetEdges=0;
		for(int i=0; i< myRank; i++)
			procOffsetEdges= procOffsetEdges + numUniqueInterface[i];

		int numInterfaceEdges=0;
		
		vec_GO_Type uniqueInterfaceIDsList_(inzidenzIndicesUnique.size());
		for(int i=0; i< uniqueInterfaceIDsList_.size(); i++)
			uniqueInterfaceIDsList_[i] = procOffsetEdges + i;

		MatrixPtr_Type indMatrix = Teuchos::rcp( new Matrix_Type(this->mapUnique_, 40 ) );

		for(int i=0; i<inzidenzIndicesUnique.size(); i++ ){
			index[0] = inzidenzIndicesUnique[i][0];
			col[0] = inzidenzIndicesUnique[i][1];
			Teuchos::Array<SC> value2(1,uniqueInterfaceIDsList_[i]);
			indMatrix->insertGlobalValues(index[0], col(), value2());
		 }
   		indMatrix->fillComplete(); 

		MatrixPtr_Type importMatrix = Teuchos::rcp( new Matrix_Type(this->mapRepeated_, 40 ) );
   		
		importMatrix->importFromVector(indMatrix,false,"Insert");
		importMatrix->fillComplete(); 		
		
		// Determine global indices
		GO edgeID=0;
		colMap = importMatrix->getMap("col");
		rowMap = importMatrix->getMap("row");
	
		LO valueID=0;
		bool found = false;
		GO entry =0;
		for(int i=0; i<inzidenzIndices.size(); i++ ){
			
			importMatrix->getLocalRowView(rowMap->getLocalElement(inzidenzIndices[i][0]), indices,values); // Indices and values connected to node i / row i in Matrix
			// Entries in 'indices' represent the local entry in 'colmap
			// with 'getGlobalElement' we know the global Node ID that belongs to the first Node that form an edge
			// vector in with entries only for edges belonging to node i;
			vec2D_GO_Type indicesTmp(indices.size(),vec_GO_Type(2));
			vec_GO_Type indTmp(2);
			for(int j=0; j<indices.size(); j++){
				indTmp[0] = colMap->getGlobalElement(indices[j]);
				indTmp[1] = values[j];
				indicesTmp.push_back(indTmp);	// vector with the indices and values belonging to node i
			}
			//sort(indicesTmp.begin(),indicesTmp.end());
			found = false;
			for(int k=0; k<indicesTmp.size();k++){
				if(inzidenzIndices[i][1] == indicesTmp[k][0]){
					entry =k;
					k = indicesTmp.size();
					edgeID = indicesTmp[entry][1];
					vecGlobalIDsEdges.at(localEdgeIndex[i]) = offsetInterface + edgeID;
					found =true;
				}
			}
			if(found == false)
				cout << " Asking for row " << rowMap->getLocalElement(inzidenzIndices[i][0]) << " for Edge [" << inzidenzIndices[i][0] << ",  " << inzidenzIndices[i][1] << "], on Proc " << myRank << " but no Value found " <<endl;
		 }


		Teuchos::RCP<std::vector<GO>> edgesGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDsEdges ) );
		Teuchos::ArrayView<GO> edgesGlobMappingArray = Teuchos::arrayViewFromVector( *edgesGlobMapping);

		this->edgeMap_.reset(new Map<LO,GO,NO>(this->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingArray, 0, this->comm_) );
		//this->edgeMap_->print();
}


/*!

\brief Refinement Restrictions
\brief In 2D we can add some Restrictions to the Mesh Refinement:
\brief Bisection:	this will keep the regularity of the Mesh by only refining whith a irregular strategy 
			  	when the longest edge is involved. If not we add a node to the longest edge, whereby 
					the irregular refinement strategy is changed.
\brief GreenTags:	this will only check tagged green Elements, if its irregular refinement tag from the previous
					refinement is 'green' and if so not refine it green again but add a node to the longest
					edge and thus refine it blue.
\brief In the 3D Case we simply never refine an element irregularly twice, this strategy is called simply 'Bey'.
If an element is refined regular, its refinement tag changes from eventually 'irregular' to regular. If those elements
should still not be refined irregular we use the strategy 'BeyIrregular'.

\brief Furthermore if there is no fitting irrregular refinement strategy (Type(1)-Type(4) don't fit) we refine regular instead.


@param[in] meshP1 P_1 Mesh.
@param[in] elements Element.
@param[in] edgeElements Edges.
@param[in] iteration Current iteration.
@param[in] newPoints Number of new unique points originating from restrictions.
@param[in] newPointsCommon Number of new repeated points originating from restrictions.
@param[in] globalInterfaceIDsTagged List of global IDs of tagged interface edges.
@param[in] mapInterfaceEdges Map of interface edges.
@param[in] restriction The kind of restriction we want to apply.
@param[in] newElements Number of new elements orginating from restrictions.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refinementRestrictions(MeshUnstrPtr_Type meshP1, ElementsPtr_Type elements ,EdgeElementsPtr_Type edgeElements,SurfaceElementsPtr_Type surfaceTriangleElements,int& newPoints, int& newPointsCommon, vec_GO_Type& globalInterfaceIDsTagged, MapConstPtr_Type mapInterfaceEdges,int& newElements){

	vec2D_dbl_ptr_Type points = meshP1->getPointsRepeated(); // Points
	string restriction = refinementRestriction_;

	if(this->dim_ == 2){
		// We determine whether a element that is tagged for green refinement has been refined green in the previous refinement
		// If thats is the case, we choose to refine element blue instead by adding a node to the longest edge
		vec_int_Type tagCounter(elements->numberElements());

		int edgeNum;

		MapConstPtr_Type edgeMap = meshP1->getEdgeMap();

		// Determine Global Ids of Interface Edges
		vec2D_GO_Type elementsOfEdgesGlobal = edgeElements->getElementsOfEdgeGlobal();
		vec2D_LO_Type elementsOfEdgesLocal = edgeElements->getElementsOfEdgeLocal();
		
		int entry;

		for(int i=0;i<elements->numberElements() ;i++){
			for(int j=0;j<3;j++){
				edgeNum = edgeElements->getEdgesOfElement(i).at(j);
				if(edgeElements->getElement(edgeNum).isTaggedForRefinement()){
					tagCounter[i]=tagCounter[i]+1; // We count the tags of element -> one Tag equals green Refinement
				}
			}
		}

		int layer =1;
		while(layer >0){
			vec_GO_Type globalEdges(0);
			layer=0;
			if(refinementRestriction_ == "Bisection"){
				for(int i=0;i<elements->numberElements() ;i++){
					if(elements->getElement(i).isTaggedForRefinement()==false){ // only looking at the untagged Elements
						if(tagCounter[i]==1){
							entry = this->determineLongestEdge(edgeElements,edgeElements->getEdgesOfElement(i),points); // we determine the edge, we would choose for blue Refinement
							if(!edgeElements->getElement(entry).isTaggedForRefinement()){ // If the longestest edge is already the tagged one, we leave the Element alone
								edgeElements->getElement(entry).tagForRefinement(); // we tag the Element for refinement
								this->addMidpoint(edgeElements,entry);	// we add the necessary midpoint
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
								this->addMidpoint(edgeElements,entry);	// we add the necessary midpoint
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

			else if (refinementRestriction_ == "GreenTags"){
				for(int i=0;i<elements->numberElements() ;i++){
				
					if(elements->getElement(i).isTaggedForRefinement()==false){ // only looking at the untagged Elements
					
						if(tagCounter[i]==1 && 	elements->getElement(i).getFiniteElementRefinementType() == "green"){
							entry = this->determineLongestEdge(edgeElements,edgeElements->getEdgesOfElement(i),points); // we determine the edge, we would choose for blue Refinement
							if(!edgeElements->getElement(entry).isTaggedForRefinement()){ // If the longestest edge is already the tagged one, we leave the Element alone
								edgeElements->getElement(entry).tagForRefinement(); // we tag the Element for refinement
								this->addMidpoint(edgeElements,entry);	// we add the necessary midpoint
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
						this->addMidpoint(edgeElements,ind);
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


	else if(this->dim_== 3){
		restriction = refinementRestriction_;
		int alright = 0;
		MapConstPtr_Type edgeMap = meshP1->getEdgeMap();
		int numPoints=0;
		int layer =0;
		int inputLayer = currentIter_;
		while(alright==0){
			alright=1;
			if(layer == inputLayer &&( restriction == "Bey" || restriction == "BeyIrregular") ){
				restriction = "Bey";			
			}
			else if( layer == inputLayer +1)
				restriction = "None";
			layer++;
			vec_GO_Type untaggedIDs(0);
			for(int i=0;i<elements->numberElements() ;i++){
				vec_int_Type nodeInd(0);
				int nodeTag=0;
				int edgeTag=0;
				vec_GO_Type untaggedIDsTmp(0);
				int entry =0;
				if(elements->getElement(i).isTaggedForRefinement()==false){ // only looking at the untagged Elements
					for(int j=0;j<6;j++){
						int edgeNum = edgeElements->getEdgesOfElement(i).at(j);
						if(edgeElements->getElement(edgeNum).isTaggedForRefinement()){
							nodeInd.push_back(edgeElements->getElement(edgeNum).getNode(0));
							nodeInd.push_back(edgeElements->getElement(edgeNum).getNode(1));
							edgeTag++;
						}
						else{
							if(edgeElements->getElement(edgeNum).isInterfaceElement())
								untaggedIDsTmp.push_back(edgeMap->getGlobalElement(edgeNum));
						}
					}
					sort( nodeInd.begin(), nodeInd.end() );
					nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );
					nodeTag = nodeInd.size();
					if(restriction == "BeyIrregular") {
						if(edgeTag > 0 && elements->getElement(i).getFiniteElementRefinementType( ) == "irregularRegular" ){
							numPoints= this->pointsRep_->size();
							elements->getElement(i).tagForRefinement();
							elements->getElement(i).setFiniteElementRefinementType("irregular");
							this->bisectEdges( edgeElements, elements, i,surfaceTriangleElements, "all");
							newPoints=newPoints + this->pointsRep_->size()-numPoints;
							alright=0;
							for(int j=0; j<untaggedIDsTmp.size(); j++)
								untaggedIDs.push_back(untaggedIDsTmp[j]);
				
						}
						else if(edgeTag > 0 && elements->getElement(i).getFiniteElementRefinementType( ) == "irregular"){
							numPoints= this->pointsRep_->size();
							elements->getElement(i).tagForRefinement();
							this->bisectEdges( edgeElements, elements, i,surfaceTriangleElements, "all");
							newPoints=newPoints + this->pointsRep_->size()-numPoints;
							alright=0;
							for(int j=0; j<untaggedIDsTmp.size(); j++)
								untaggedIDs.push_back(untaggedIDsTmp[j]);
				
						}
				    }
					else if(restriction== "Bey") {
						if(edgeTag > 0 && elements->getElement(i).getFiniteElementRefinementType( ) == "irregular" ){
							numPoints= this->pointsRep_->size();
							elements->getElement(i).tagForRefinement();
							this->bisectEdges( edgeElements, elements, i,surfaceTriangleElements, "all");
							newPoints=newPoints + this->pointsRep_->size()-numPoints;
							alright=0;
							for(int j=0; j<untaggedIDsTmp.size(); j++)
								untaggedIDs.push_back(untaggedIDsTmp[j]);
				
						}
				    }
					if(nodeTag >3 && edgeTag >2){
							numPoints= this->pointsRep_->size();
							elements->getElement(i).tagForRefinement();
							this->bisectEdges( edgeElements, elements, i,surfaceTriangleElements, "all");
							newPoints=newPoints + this->pointsRep_->size()-numPoints;
							alright=0;
							for(int j=0; j<untaggedIDsTmp.size(); j++)
								untaggedIDs.push_back(untaggedIDsTmp[j]);
						}
		
				}
			}

			
			sort(untaggedIDs.begin(),untaggedIDs.end());
			untaggedIDs.erase(unique(untaggedIDs.begin(),untaggedIDs.end()),untaggedIDs.end());
			// Constructing a map of the global IDs of the tagged Edges	
			Teuchos::ArrayView<GO> globalEdgesArray = Teuchos::arrayViewFromVector( untaggedIDs);
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

			Teuchos::ArrayRCP< const LO >  tags = taggedEdgesLocal->getData( 0 );

			// Adding Midpoints and tagging the edges that were tagged on other Procs an are on the interface
			// Collecting global Ids of tagged Interface Edges
			LO ind;
			GO indG;

			for (int i=0; i<tags.size(); i++) {
				if (tags[i] > 0){
					indG = mapInterfaceEdges->getGlobalElement(i);
					ind = edgeMap->getLocalElement(indG);
					globalInterfaceIDsTagged.push_back(indG);
					newPointsCommon ++;
					if(!edgeElements->getElement(ind).isTaggedForRefinement()){
						edgeElements->getElement(ind).tagForRefinement();
						this->addMidpoint(edgeElements,ind);
						// globalen Index der Kante -> alle Prozessoren haben diese Kante -> die Knoten die hier hinzugefügt werden sind haben auf allen Proc den gleichen Index
						newPoints ++;	
					}
				}
			}	
			reduceAll<int, int> (*this->comm_, REDUCE_MIN, alright, outArg (alright));
		}
	}	
}

/*!

 \brief Refinement performed according to the set of rules determined by Bey or Verfürth.

@param[in] elements Elements.
@param[in] edgeElements Edges.
@param[in] newElements Number of new elements originating from refinement.
@param[in] edgeMap Map of global edge ids.
@param[in] surfaceTriangleElements Triangle elements (3D case).

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineMeshRegIreg(ElementsPtr_Type elements, EdgeElementsPtr_Type edgeElements, int& newElements,MapConstPtr_Type edgeMap, SurfaceElementsPtr_Type surfaceTriangleElements) 
{
	if(this->dim_==2){
		int edgeNum=0;
		vec_int_Type tagCounter2(elements->numberElements());
		for(int i=0;i<elements->numberElements() ;i++){
			for(int j=0;j<3;j++){
				edgeNum = edgeElements->getEdgesOfElement(i).at(j);
				if(edgeElements->getElement(edgeNum).isTaggedForRefinement())
					tagCounter2[i] ++;
			}
			if(tagCounter2[i] == 1){
				elements->getElement(i).setFiniteElementRefinementType("green");
				this->refineGreen(edgeElements, elements, i);
				newElements ++;
			}
			else if(tagCounter2[i] == 2){
				elements->getElement(i).setFiniteElementRefinementType("blue");
				this->refineBlue(edgeElements, elements, i);
				newElements= newElements+2;
				}	
			else if(tagCounter2[i] == 3){
				elements->getElement(i).setFiniteElementRefinementType("regular");
				if(refinementMode_ == "Bisection"){
					this->bisectElement3(edgeElements, elements, i);
				}
				else	
					this->refineRegular(edgeElements, elements, i, surfaceTriangleElements);
					newElements= newElements+3;
				}
			else {
				// the element in question has not been refined, nor its neighbour -> add edges to edge list
				vec_int_Type edgesOfElement = edgeElements->getEdgesOfElement(i);
				for(int j=0;j<3;j++){
					edgeElements->getElement(edgesOfElement[j]).setFiniteElementRefinementType("unrefined");
					edgeElements->getElement(edgesOfElement[j]).setPredecessorElement(edgeMap->getGlobalElement(edgesOfElement[j]));
					edgeElements->getElement(edgesOfElement[j]).setInterfaceElement( edgeElements->getElement(edgesOfElement[j]).isInterfaceElement());
					this->edgeElements_->addEdge(edgeElements->getElement(edgesOfElement[j]),i);
					
					/*if(edgeElements->getElement(edgesOfElement[j]).getFlag() !=0 && edgeElements->getElement(edgesOfElement[j]).getFlag() !=10){
						if ( !this->elementsC_->getElement(i).subElementsInitialized() )
							this->elementsC_->getElement(i).initializeSubElements( this->FEType_, this->dim_ -1) ;
						this->elementsC_->getElement(i).addSubElement(edgeElements->getElement(edgesOfElement[j]));	
					}*/

				}	

			}
		}
	}
	// Irregular Refinement algorithm according to Bey
	else if(this->dim_==3){
		int edgeNum=0;
		int nodeTag =0;
		int edgeTag =0;
		for(int i=0;i<elements->numberElements() ;i++){
			vec_int_Type nodeInd(0);
			nodeTag=0;
			edgeTag=0;
			for(int j=0;j<6;j++){
				edgeNum = edgeElements->getEdgesOfElement(i).at(j);
				if(edgeElements->getElement(edgeNum).isTaggedForRefinement()){
					nodeInd.push_back(edgeElements->getElement(edgeNum).getNode(0));
					nodeInd.push_back(edgeElements->getElement(edgeNum).getNode(1));
					edgeTag++;
				}
			}
			sort( nodeInd.begin(), nodeInd.end() );
			nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );
			nodeTag = nodeInd.size();
		
			if(nodeTag == 3 && edgeTag ==3){
				//elements->getElement(i).setFiniteElementRefinementType("Type1");
				//cout << " Requesting Type 1 Refinement on Processor " << this->comm_->getRank()<< endl;
				this->refineType1(edgeElements, elements, i,surfaceTriangleElements);
				newElements = newElements+3;
			}
			else if(nodeTag == 2 && edgeTag == 1){
				//elements->getElement(i).setFiniteElementRefinementType("Type2");
				//cout << " Requesting Type 2 Refinement on Processor " << this->comm_->getRank() << endl;
				this->refineType2(edgeElements, elements, i,surfaceTriangleElements);
				newElements ++;

			}
			else if(nodeTag == 3 && edgeTag == 2){
				//elements->getElement(i).setFiniteElementRefinementType("Type3");
				//cout << " Requesting Type 3 Refinement on Processor " << this->comm_->getRank() << endl;
				this->refineType3(edgeElements, elements, i,surfaceTriangleElements);

				newElements = newElements+2;
			}
			else if(nodeTag == 4 && edgeTag == 2){
				//elements->getElement(i).setFiniteElementRefinementType("Type4");
				//cout << " Requesting Type 4 Refinement on Processor " << this->comm_->getRank() << endl;
				this->refineType4(edgeElements, elements, i,surfaceTriangleElements);
				newElements = newElements+3;
			}
			else if(nodeTag == 4 && edgeTag == 6 ){
				this->refineRegular(edgeElements, elements, i, surfaceTriangleElements);
				newElements = newElements+7;
			}
			else if(nodeTag >3 && edgeTag >2){
   					TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "RefineMesh: Requesting 3D Refinement Type, that doesn't exist.");
			}

			else if(edgeTag ==0){
				vec_int_Type edgesOfElement = edgeElements->getEdgesOfElement(i);
				for(int j=0;j<6;j++){
					edgeElements->getElement(edgesOfElement[j]).setFiniteElementRefinementType("unrefined");
					edgeElements->getElement(edgesOfElement[j]).setPredecessorElement(edgeMap->getGlobalElement(edgesOfElement[j]));
					this->edgeElements_->addEdge(edgeElements->getElement(edgesOfElement[j]),i);

				}

				vec_int_Type surfacesOfElement = surfaceTriangleElements->getSurfacesOfElement(i);
				for(int j=0;j<4;j++){
					this->surfaceTriangleElements_->addSurface(surfaceTriangleElements->getElement(surfacesOfElement[j]),i);
				}
			}
		}
	}
						
	for(int j=0;j<this->edgeElements_->numberElements();j++){
		if(this->edgeElements_->getElement(j).isTaggedForRefinement())
			cout<< "tagged edge element somehow made it " << endl;
	}

}


/*!

 \brief Updating ElementsOfEdgesLocal and ElementsOfEdgesGlobal.

@param[in] maxRank The maximal processor rank.
@param[in] edgeMap Map of global edge ids.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::updateElementsOfEdgesLocalAndGlobal(int maxRank,  MapConstPtr_Type edgeMap){

	if(maxRank >0 && this->dim_ == 2){
		vec_GO_Type edgesInterfaceGlobalID(0);
		LO id=0;
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(this->edgeElements_->getElement(i).isInterfaceElement() ){		
				this->edgeElements_->setElementsOfEdgeLocalEntry(i,-1);
				edgesInterfaceGlobalID.push_back(this->edgeMap_->getGlobalElement(i)); // extracting the global IDs of the new interfaceEdges
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

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );

		MultiVectorLOPtr_Type isInterfaceElement_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_imp->putScalar( (LO) 0 ); 
		isInterfaceElement_imp->importFromVector( interfaceElements, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( interfaceElements, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_imp, false, "Insert");

		isInterfaceElement2_imp->exportFromVector(isInterfaceElement_exp, false, "Insert");

		interfaceElementsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			this->edgeElements_->setElementsOfEdgeGlobalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),interfaceElementsEntries[i]);
			}

	}

	// Contrary to the 2D case in 3D an edge can be part of more than two elements
	// We need to determine how many elements are connected to an edge and import the global IDs from different Processors
	if(maxRank >0 && this->dim_ == 3){

		// First we determine the interface edges, which entries in elementsOfEdgesGlobal/Local have to be completed
		vec_GO_Type edgesInterfaceGlobalID(0);
		vec_int_Type numberElements(0); // represents number of elements the edge is connected to on my Processor
		
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(this->edgeElements_->getElement(i).isInterfaceElement() ){
				edgesInterfaceGlobalID.push_back(this->edgeMap_->getGlobalElement(i)); // extracting the global IDs of the new interfaceEdges
			}	
		}
		sort(edgesInterfaceGlobalID.begin(), edgesInterfaceGlobalID.end());

		for(int i=0; i< edgesInterfaceGlobalID.size(); i++){
			numberElements.push_back(this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).size());
		}
			
		
		// from this we build a map
		Teuchos::ArrayView<GO> edgesInterfaceGlobalID_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalID);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalID_, 0, this->comm_) );

		// As edges can be part of multiple elements on different processors we collect the number of elements connected to the edge in total
		MultiVectorLOPtr_Type numberInterfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > numberInterfaceElementsEntries  = numberInterfaceElements->getDataNonConst(0);

		for(int i=0; i< numberInterfaceElementsEntries.size(); i++)
			numberInterfaceElementsEntries[i] = numberElements[i];

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );
	
		// Element are unique to each processor. This means that the number we collect is the number of elements that are connected to my edge on other processors.
		// With the following communication we add up all the entries for a certain global Edge ID
		// Potential causes of error:
		//  - if an edge is not identified as an interface edge, of course it will not import nor export its interface Information, making itself and others incomplete

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( numberInterfaceElements, false, "Add");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_exp, true, "Insert");

		Teuchos::ArrayRCP< LO > numberInterfaceElementsImportEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		vec_int_Type missingEntries(numberInterfaceElementsEntries.size());
		// With this number we can complete the elementsOfEdgeLocal List with -1 for the elements not on our processor
		for(int i=0; i<numberInterfaceElementsEntries.size() ; i++){
			for(int j=0; j< numberInterfaceElementsImportEntries[i] - numberInterfaceElementsEntries[i];j++){
				this->edgeElements_->setElementsOfEdgeLocalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),-1);
				missingEntries[i] = numberInterfaceElementsImportEntries[i] -numberInterfaceElementsEntries[i];
			}
		}

		// Next we need to identify the global Element IDs of those missing entries and communicate them
		// Hey i got the global Elements ... of edge ... -> exchange
		// Elements are uniquely distributed -> you cannot import an element you already have
		// I need x number of entries -> import all i need, export all i have 
		// Global IDs of Procs
		// Setting newPoints as to be communicated Values

		// Communicating max number of necessary values:
		vec_int_Type::iterator it;
		it = max_element(numberElements.begin(), numberElements.end());
		int myNumberElementsMax = numberElements.at(distance(numberElements.begin(), it)); // accumulate(errorElement.begin(), errorElement.end(),0.0);

		reduceAll<int, int> (*this->comm_, REDUCE_MAX,  myNumberElementsMax , outArg ( myNumberElementsMax));

		MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

		vec2D_int_Type importElements(this->edgeElements_->getElementsOfEdgeGlobal().size(),vec_int_Type( 0));

		// We extended this function to also consider the ranks. Before we would only exchange the information, but in case one processor received information from more than one processor at 
		// the same time some of the information would get lost. Now we only send the information one processor holds to the other at the same time and move through the processor only destributing their
		// edgeOfElementGlobal Information in a circle like order
		for(int k=0; k< maxRank+1 ; k++){
			
			vec_GO_Type edgesInterfaceGlobalIDProc;				
			if(this->comm_->getRank() == k ){
				edgesInterfaceGlobalIDProc = edgesInterfaceGlobalID; // extracting the global IDs of the new interfaceEdges
			}
					

			// from this we build a map
			Teuchos::ArrayView<GO> edgesInterfaceGlobalIDProc_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalIDProc);

			MapPtr_Type mapGlobalInterfaceProcs =
				Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalIDProc_, 0, this->comm_) );

			for(int j=0; j< myNumberElementsMax; j++){
				MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceProcs, 1 ) );
				Teuchos::ArrayRCP< LO > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

				for(int i=0; i< interfaceElementsEntries.size() ; i++){		
					if(numberElements[i] > j && this->comm_->getRank() == k )
						interfaceElementsEntries[i] = this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).at(j);
					else
						interfaceElementsEntries[i] = -1; 
				}

				MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
				isInterfaceElement_exp->putScalar( (LO) -1 ); 
				isInterfaceElement_exp->exportFromVector( interfaceElements, false, "Insert");

				if(this->comm_->getRank() == k && mapGlobalInterfaceUnique->getNodeNumElements() > 0){
					Teuchos::ArrayRCP< LO > interfaceElementsEntries_exp  = isInterfaceElement_exp->getDataNonConst(0);
					for(int i=0; i<  interfaceElementsEntries_exp.size() ; i++){
						LO id = mapGlobalInterface->getLocalElement(mapGlobalInterfaceUnique->getGlobalElement(i));
						interfaceElementsEntries_exp[i] = interfaceElementsEntries[id];
					}
									
				}

			
				MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
				isInterfaceElement2_imp->putScalar( (LO) 0 ); 
				isInterfaceElement2_imp->importFromVector(isInterfaceElement_exp, false, "Insert");

				interfaceElementsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

				for(int i=0; i< interfaceElementsEntries.size() ; i++){
					if(this->comm_->getRank() != k && interfaceElementsEntries[i] != -1)
						importElements[i].push_back( interfaceElementsEntries[i]);
				}
				
			}

		}
		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			sort(importElements[i].begin(),importElements[i].end());
			importElements[i].erase( unique(importElements[i].begin(), importElements[i].end() ), importElements[i].end() );
			if(importElements[i].size() != missingEntries[i])
			   cout << " On Processor " << this->comm_->getRank() << " uneven number for edge imported: " << importElements[i].size() << " missing " << missingEntries[i] << " " << edgesInterfaceGlobalID[i] << endl; // " something went wrong while updating elementsOfEdgesGlobal as the imported entries do not match the supposed number of imported entries. Please check." << endl; 	
		}

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			for(int j=0; j < importElements[i].size(); j++){
				if(importElements[i][j] != -1)
					this->edgeElements_->setElementsOfEdgeGlobalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),importElements[i][j]);
			}
		}

	}

}

/*!

 \brief Adding a Midpoint on an edge.

@param[in] edgeElements Edges.
@param[in] edgeID Edge ids where the midpoints is added.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::addMidpoint(EdgeElementsPtr_Type edgeElements, int edgeID){

	int midPointInd; // indices of midpoints of edges of soon to be refined element
	int dim = this->dim_;
	vec2D_dbl_Type points2(0, vec_dbl_Type( dim )); // new Points -> depends on already existing point from previously refined neighbouring elements
	vec_dbl_Type point(dim); // midpoint on every edge, that will be added to points, will be added to pointsRep_

	vec_dbl_Type P1(dim),P2(dim); // points we extract from pointsRep_ in order to determine midPoints

	LO p1ID =edgeElements->getElement(edgeID).getNode(0);
	LO p2ID =edgeElements->getElement(edgeID).getNode(1);
	P1 = this->pointsRep_->at(p1ID);
	P2 = this->pointsRep_->at(p2ID);

	for (int d=0; d<dim; d++){
		point[d]= ( (P1)[d] + (P2)[d] ) / 2.;
	}   
	points2.push_back(point); 

	// New Flags:	
	this->bcFlagRep_->push_back(edgeElements->getElement(edgeID).getFlag());
		
	// Mittelpunkte der Kanten setzen
	edgeElements->setMidpoint(edgeID, this->pointsRep_->size());

	// We have to keep in mind, that the new added points need to be maped to the associated points
	this->pointsRep_->insert( this->pointsRep_->end(), points2.begin(), points2.end() );
}

/*!

 \brief Eetermine longest edge in triangle.

@param[in] edgeElements Edges.
@param[in] edgeVec Vector with edge ids of triangle.
@param[in] points Points.

@param[out] Local edgeID of the longest edge.

*/


template <class SC, class LO, class GO, class NO>
int RefinementFactory<SC,LO,GO,NO>::determineLongestEdge( EdgeElementsPtr_Type edgeElements, vec_int_Type edgeVec, vec2D_dbl_ptr_Type points){
	// We have to determine which edge is longer, as we use the opposite node of the longer edge for the element construction
	vec_dbl_Type length(edgeVec.size());
	vec_dbl_Type P1(this->dim_),P2(this->dim_);
	double maxLength=0.0;
	int maxEntry=0;
	LO p1ID,p2ID;
	vec2D_dbl_Type nodeInd(0,vec_dbl_Type(this->dim_+2));
	for(int i=0;i<edgeVec.size();i++){
		p1ID =edgeElements->getElement(edgeVec[i]).getNode(0);
		p2ID =edgeElements->getElement(edgeVec[i]).getNode(1);
		P1 = points->at(p1ID);
		P2 = points->at(p2ID);
		double sum=0;
		for(int j=0; j< P1.size();j++)
			sum += pow(P1[j]-P2[j],2);

		length[i] = sqrt(sum);
	
		vec2D_dbl_Type tmpN(0,vec_dbl_Type(this->dim_+1));
		vec_dbl_Type tagged(1,2);
		if(edgeElements->getElement(edgeVec[i]).isTaggedForRefinement() == true)
			tagged[0]=0;
		else if(edgeElements->getElement(edgeVec[i]).isMarkedEdge()==true)
			tagged[0]=1;

		tmpN.push_back(P1);
		tmpN.push_back(P2);
		sort(tmpN.begin(),tmpN.end());
		tmpN[0].push_back(i);
		tmpN[0].push_back(length[i]);
		nodeInd.push_back(tmpN[0]);
		nodeInd[i].insert( nodeInd[i].begin(), tagged.begin(), tagged.end() );
	}
	sort(nodeInd.begin(), nodeInd.end());
	int lInd = nodeInd[0].size()-1;
	for(int i=0; i< nodeInd.size(); i++){
		if(nodeInd[i][lInd] > maxLength){
			maxLength = nodeInd[i][lInd];
			maxEntry= (int) nodeInd[i][lInd-1];
		}
	}
			

	return edgeVec[maxEntry];
	
}

/*!

 \brief 2D blue refinement: refining element according to blue refinement scheme - connecting nodes of shorter edge with midpoint of longer tagged edge and connect that with opposite corner

@param[in] edgeElements Edges
@param[in] elements Elements.
@param[in] indexELement Element in question.

*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineBlue(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){
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
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)

		// Element 1
		(newElements)[0]={oppositeNodeIndL,midPointInd[edgeIndexL],midPointInd[edgeIndexS]};

		(newEdges)[0] = {oppositeNodeIndL ,midPointInd[edgeIndexL]}; 
		(newEdges)[1] = {oppositeNodeIndL,midPointInd[edgeIndexS]}; 
		(newEdges)[2] = {midPointInd[edgeIndexL] ,midPointInd[edgeIndexS]}; 

		edgeFlags[0]=this->volumeID_;
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[edgeIndexS]);
		edgeFlags[2]=this->volumeID_;

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


		edgeFlags[3]=this->volumeID_;
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

		edgeFlags[6]=this->volumeID_;
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
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("blue");
			feNew.setPredecessorElement(indexElement);
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
				if(edgeFlags[i]!=0 && edgeFlags[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(i/3+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/3+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/3+offsetElements).addSubElement(feNew);	

					}	
				}
			else{
				this->edgeElements_->addEdge(feNew,indexElement);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	}

				}			
		}	

    }

	else if(this->dim_ == 3){ 
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "3D Algorithm for irregular MeshRefinement is currently not available, please choose uniform Refinement");
			}
    


}
/*!

 \brief 2D green refinement: refining the element according to green scheme - connecting node on refined edge with the opposite node. 

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineGreen(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){

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
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)
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
			edgeFlags[0]=this->volumeID_;
			isInterfaceEdge[0] = false;
			predecessorElement[0] = -1;
			}
		else {
			edgeFlags[0]=this->bcFlagRep_->at(midPointInd);
			isInterfaceEdge[0] = edgeElements->getElement(taggedEdge).isInterfaceElement();
			predecessorElement[0] = this->edgeMap_->getGlobalElement(taggedEdge);
		}
	
	    if(oppositeNodeInd == edgeElements->getElement(sortedEdges[1]).getNode(1)){
			edgeFlags[1]=this->volumeID_;
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
			edgeFlags[3]=this->volumeID_;
			isInterfaceEdge[3] = false;
			predecessorElement[3] = -1;

		}	
		else {
			edgeFlags[3]=this->bcFlagRep_->at(midPointInd);
			isInterfaceEdge[3] = edgeElements->getElement(taggedEdge).isInterfaceElement();
			predecessorElement[3] = this->edgeMap_->getGlobalElement(taggedEdge);
		}
	
	    if(oppositeNodeInd == edgeElements->getElement(sortedEdges[2]).getNode(1)){
			edgeFlags[4]=this->volumeID_;
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
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("green"); // setting green refinement type in order to check in following refinement stages
			feNew.setPredecessorElement(indexElement);
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
				if(edgeFlags[i]!=0 && edgeFlags[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(offsetElements).addSubElement(feNew);		
					}
				}
			else{
				this->edgeElements_->addEdge(feNew,indexElement);
				if(edgeFlags[i]!=0 && edgeFlags[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	
					}
				}
			}				
    }

	else if(this->dim_ == 3){ 
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "3D Algorithm for irregular MeshRefinement is currently not available, please choose uniform Refinement");
			}
    


}
/*!

 \brief 2D red refinement: refining the element red by connecting all tagged edges midpoints. one element is refined into 4.

@param[in] edgeElements Edges
@param[in] elements Elements.
@param[in] indexELement Element in question.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineRed(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){

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
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)

		// Element 1
		(newElements)[0]={mutualNode[0],midPointInd[0],midPointInd[1]};

		(newEdges)[0] = {mutualNode[0] ,midPointInd[0]}; 
		(newEdges)[1] = {mutualNode[0] ,midPointInd[1]}; 
		(newEdges)[2] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[0]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[2]=this->volumeID_;

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
		edgeFlags[5]=this->volumeID_;

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
		edgeFlags[8]=this->volumeID_;

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

		edgeFlags[9]=this->volumeID_;
		edgeFlags[10]=this->volumeID_;
		edgeFlags[11]=this->volumeID_;

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
			feNew.setPredecessorElement(indexElement);
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
				if(edgeFlags[i]!=0 && edgeFlags[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(i/3+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/3+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/3+offsetElements).addSubElement(feNew);	
					}
				}
			else
				this->edgeElements_->addEdge(feNew,indexElement);
			
		}
    }
	else        	
		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The red irregular Refinement Method you requested is only applicable to a 2 dimensional Mesh.");
		
}

/*!

 \brief 3D Type(4) refinement as defined in  "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.
@param[in] surfaceTriangleElements Triangle elements.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineType4(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements){

// Implementation of Type (4) Refinement Type
// We use this Refinement Type 
	if(this->dim_ == 3){ 

		// The way we refine the Tetrahedron is defined by how we order the nodes of the tetrahedron
		// (For the algorithm see "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955)

        vec_int_Type midPointInd( 0 ); // indices of midpoints of edges of soon to be refined element
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		vec_int_Type edgeNumbersUntagged(0);

		// We sort the node as follows:
		// In this refinement type both tagged edges don't share a common node thus we add the nodes of the first edge, then the second and order them per edge
		vec_int_Type nodeInd1(0);
		vec_int_Type nodeInd2(0);
		bool firstEdge = true;
		for(int i=0; i<6; i++)	{

			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement() && firstEdge ==false){
				nodeInd2.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
				nodeInd2.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
			}

			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement() && firstEdge == true){
				nodeInd1.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
				nodeInd1.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
				firstEdge=false;
			}

		}
		sort( nodeInd1.begin(), nodeInd1.end() );
		sort( nodeInd2.begin(), nodeInd2.end() );

		vec_int_Type nodeInd = {nodeInd1[0],nodeInd1[1],nodeInd2[0],nodeInd2[1]};

		// We don't need to sort the nodes in any other way then we already did, as there is not more than one possibility to construct the subtetrahera in refinement type 3

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
		vec_int_Type surfaceElementsIDs = surfaceTriangleElements->getSurfacesOfElement(indexElement); // surfaces of Element k
		vec2D_int_Type originTriangles(4,vec_int_Type(3));

		originTriangles[0] = {nodeInd[0], nodeInd[1], nodeInd[2] };
		originTriangles[1] = {nodeInd[0], nodeInd[1], nodeInd[3] };
		originTriangles[2] = {nodeInd[0], nodeInd[2], nodeInd[3] };
		originTriangles[3] = {nodeInd[1], nodeInd[2], nodeInd[3] };
	
		vec_int_Type originFlag(4,this->volumeID_); // Triangle Flag

		vec_bool_Type interfaceSurface(4);
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);

		for(int i=0; i< 4 ; i++){
			originTriangleTmp = originTriangles[i];
			sort(originTriangleTmp.begin(),originTriangleTmp.end());
			for(int j=0; j<4 ; j++){
				FiniteElement surfaceTmp = surfaceTriangleElements->getElement(surfaceElementsIDs[j]);
				triTmp = surfaceTmp.getVectorNodeList();
				sort(triTmp.begin(),triTmp.end());
				if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] && triTmp[2] == originTriangleTmp[2]  ) {
					originFlag[i] = surfaceTmp.getFlag();
					interfaceSurface[i] = surfaceTmp.isInterfaceElement();
				}
			}
		}

		// Finally we need to determine or extract the indices of the edges midpoints. As in the describe algorithm the midpoints are set as follows:
		// Edge 0 = [x_0,x_1] -> x_01
		// Edge 1 = [x_0,x_2] -> x_02
		// Edge 2 = [x_0,x_3] -> x_03
		// Edge 3 = [x_1,x_2] -> x_12
		// Edge 4 = [x_1,x_3] -> x_13
		// Edge 5 = [x_2,x_3] -> x_23

		for(int i=0; i<6; i++){
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				midPointInd.push_back(edgeElements->getMidpoint(edgeNumbers[i]));
			}
		}

		
	
		// Now we construct the new Elements as proposed by Bey's Regular Refinement Algorithm

		vec2D_int_Type newElements(4, vec_int_Type( 0 )); // new elements
		vec2D_int_Type newEdges(24,vec_int_Type(0)); // new edges
		vec2D_int_Type newTriangles(16,vec_int_Type(0)); // new Triangles
		vec_int_Type newTrianglesFlag(16) ; // new Triangle Flags
		vec_bool_Type isInterfaceSurface(16); // bool vector for interfaceSurface Information
		vec_int_Type edgeFlags(24); // new Edge flags
		vec_bool_Type isInterfaceEdge(24); // bool vector for interfaceEdge Information
		vec_GO_Type predecessorElement(24); // vector that stores the global IDs of the predecessor of each edge 

		vec2D_LO_Type newTriangleEdgeIDs(4,vec_LO_Type(12));

		// How are Flags determined?
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)
		// If an edges emerges on a triangle, the flag is determined by the triangle flag. Opposite to the 2D case, edges that connect midpoints are not automatically interior edges, but are
		// determined by the triangle/surface they are on


		// Element 1: (x_0,x_2,x_01,x_23)
		(newElements)[0]={nodeInd[0],nodeInd[2], midPointInd[0],midPointInd[1]};

		(newEdges)[0] = {nodeInd[0] ,nodeInd[2]}; 
		(newEdges)[1] = {nodeInd[0] ,midPointInd[0]}; 
		(newEdges)[2] = {nodeInd[0] ,midPointInd[1]}; 
		(newEdges)[3] = {nodeInd[2] ,midPointInd[0]}; 
		(newEdges)[4] = {nodeInd[2] ,midPointInd[1]}; 
		(newEdges)[5] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[0]=edgeElements->getElement(edgeNumbers[1]).getFlag();
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[2]=originFlag[2];
		edgeFlags[3]=originFlag[0];
		edgeFlags[4]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[5]=this->volumeID_;

		isInterfaceEdge[0] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[1] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[2] = interfaceSurface[2];
		isInterfaceEdge[3] = interfaceSurface[0];
		isInterfaceEdge[4] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[5] = false;

		predecessorElement[0] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[1] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[2] = -1;
		predecessorElement[3] = -1;
		predecessorElement[4] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[5] = -1;

		// Subelements of thetrahedron
		newTriangles[0]= {nodeInd[0],nodeInd[2],midPointInd[0]};
		newTriangles[1]= {nodeInd[0],nodeInd[2],midPointInd[1]};
		newTriangles[2]= {nodeInd[0],midPointInd[0],midPointInd[1]};
		newTriangles[3]= {nodeInd[2],midPointInd[0],midPointInd[1]};

		newTrianglesFlag[0]= originFlag[0]; 
		newTrianglesFlag[1]= originFlag[2]; 
		newTrianglesFlag[2]= this->volumeID_; 
		newTrianglesFlag[3]= this->volumeID_;

		isInterfaceSurface[0]= interfaceSurface[0];
		isInterfaceSurface[1]= interfaceSurface[2];
		isInterfaceSurface[2]= false;
		isInterfaceSurface[3]= false;

		newTriangleEdgeIDs[0]={0,1,3,0,2,4,1,2,5,3,4,5};

		// Element 2: (x_1,x_2,x_01,x_23)
		(newElements)[1]={nodeInd[1],nodeInd[2],midPointInd[0],midPointInd[1]};

		(newEdges)[6] = {nodeInd[1] ,nodeInd[2]}; 
		(newEdges)[7] = {nodeInd[1] ,midPointInd[0]}; 
		(newEdges)[8] = {nodeInd[1] ,midPointInd[1]}; 
		(newEdges)[9] = {nodeInd[2] ,midPointInd[0]}; 
		(newEdges)[10] = {nodeInd[2] ,midPointInd[1]}; 
		(newEdges)[11] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[6]=edgeElements->getElement(edgeNumbers[3]).getFlag();
		edgeFlags[7]=edgeElements->getElement(edgeNumbers[0]).getFlag();
		edgeFlags[8]=originFlag[3];
		edgeFlags[9]=originFlag[0];
		edgeFlags[10]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[11]=this->volumeID_;

		isInterfaceEdge[6] = edgeElements->getElement(edgeNumbers[3]).isInterfaceElement();
		isInterfaceEdge[7] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[8] = interfaceSurface[3];
		isInterfaceEdge[9] = interfaceSurface[0];
		isInterfaceEdge[10] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[11] = false;

		predecessorElement[6] = this->edgeMap_->getGlobalElement(edgeNumbers[3]);
		predecessorElement[7] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[8] = -1;
		predecessorElement[9] = -1;
		predecessorElement[10] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[11] = -1;

		// Subelements of tetrahedron
		newTriangles[4]= {nodeInd[1],nodeInd[2],midPointInd[0]};
		newTriangles[5]= {nodeInd[1],nodeInd[2],midPointInd[1]};
		newTriangles[6]= {nodeInd[1],midPointInd[0],midPointInd[1]};
		newTriangles[7]= {nodeInd[2],midPointInd[0],midPointInd[1]};

		newTrianglesFlag[4]= originFlag[0]; 
		newTrianglesFlag[5]= originFlag[3]; 
		newTrianglesFlag[6]= this->volumeID_; 
		newTrianglesFlag[7]= this->volumeID_;

		isInterfaceSurface[4]= interfaceSurface[0];
		isInterfaceSurface[5]= interfaceSurface[3];
		isInterfaceSurface[6]= false;
		isInterfaceSurface[7]= false;

		newTriangleEdgeIDs[1]={6,7,9,6,8,10,7,8,11,9,10,11};

		// Element 3: (x_0,x_3,x_01,x_23) 
		(newElements)[2]={nodeInd[0],nodeInd[3],midPointInd[0],midPointInd[1]};

		(newEdges)[12] = {nodeInd[0] ,nodeInd[3]}; 
		(newEdges)[13] = {nodeInd[0] ,midPointInd[0]}; 
		(newEdges)[14] = {nodeInd[0] ,midPointInd[1]}; 
		(newEdges)[15] = {nodeInd[3] ,midPointInd[0]}; 
		(newEdges)[16] = {nodeInd[3] ,midPointInd[1]}; 
		(newEdges)[17] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[12]=edgeElements->getElement(edgeNumbers[2]).getFlag();
		edgeFlags[13]=edgeElements->getElement(edgeNumbers[0]).getFlag();
		edgeFlags[14]=originFlag[2];
		edgeFlags[15]=originFlag[1];
		edgeFlags[16]=edgeElements->getElement(edgeNumbers[5]).getFlag();;
		edgeFlags[17]=this->volumeID_;

		isInterfaceEdge[12] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[13] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[14] = interfaceSurface[2];
		isInterfaceEdge[15] = interfaceSurface[1];
		isInterfaceEdge[16] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[17] = false;

		predecessorElement[12] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[13] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[14] = -1;
		predecessorElement[15] = -1;
		predecessorElement[16] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[17] = -1;

		// Subelements of tetrahedron
		newTriangles[8]= {nodeInd[0],nodeInd[3],midPointInd[0]};
		newTriangles[9]= {nodeInd[0],nodeInd[3],midPointInd[1]};
		newTriangles[10]= {nodeInd[0],midPointInd[0],midPointInd[1]};
		newTriangles[11]= {nodeInd[3],midPointInd[0],midPointInd[1]};

		newTrianglesFlag[8]= originFlag[1]; 
		newTrianglesFlag[9]= originFlag[2]; 
		newTrianglesFlag[10]= this->volumeID_; 
		newTrianglesFlag[11]= this->volumeID_;

		isInterfaceSurface[8]= interfaceSurface[1];
		isInterfaceSurface[9]= interfaceSurface[2];
		isInterfaceSurface[10]= false;
		isInterfaceSurface[11]= false;

		newTriangleEdgeIDs[2]={8,9,11,8,10,12,9,10,13,11,12,13};


		// Element 4: (x_1,x_3,x_01,x_23)
		(newElements)[3]={nodeInd[1],nodeInd[3],midPointInd[0],midPointInd[1]};

		(newEdges)[18] = {nodeInd[1] ,nodeInd[3]}; 
		(newEdges)[19] = {nodeInd[1] ,midPointInd[0]}; 
		(newEdges)[20] = {nodeInd[1] ,midPointInd[1]}; 
		(newEdges)[21] = {nodeInd[3] ,midPointInd[0]}; 
		(newEdges)[22] = {nodeInd[3] ,midPointInd[1]}; 
		(newEdges)[23] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[18]=edgeElements->getElement(edgeNumbers[4]).getFlag();
		edgeFlags[19]=edgeElements->getElement(edgeNumbers[0]).getFlag();
		edgeFlags[20]=originFlag[3];
		edgeFlags[21]=originFlag[1];
		edgeFlags[22]=edgeElements->getElement(edgeNumbers[5]).getFlag();;
		edgeFlags[23]=this->volumeID_;

		isInterfaceEdge[18] = edgeElements->getElement(edgeNumbers[4]).isInterfaceElement();
		isInterfaceEdge[19] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[20] = interfaceSurface[3];
		isInterfaceEdge[21] = interfaceSurface[1];
		isInterfaceEdge[22] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[23] = false;

		predecessorElement[18] = this->edgeMap_->getGlobalElement(edgeNumbers[4]);
		predecessorElement[19] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[20] = -1;
		predecessorElement[21] = -1;
		predecessorElement[22] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[23] = -1;

		// Subelements of tetrahedron
		newTriangles[12]= {nodeInd[1],nodeInd[3],midPointInd[0]};
		newTriangles[13]= {nodeInd[1],nodeInd[3],midPointInd[1]};
		newTriangles[14]= {nodeInd[1],midPointInd[0],midPointInd[1]};
		newTriangles[15]= {nodeInd[3],midPointInd[0],midPointInd[1]};

		newTrianglesFlag[12]= originFlag[1]; 
		newTrianglesFlag[13]= originFlag[3]; 
		newTrianglesFlag[14]= this->volumeID_; 
		newTrianglesFlag[15]= this->volumeID_;

		isInterfaceSurface[12]= interfaceSurface[1];
		isInterfaceSurface[13]= interfaceSurface[3];
		isInterfaceSurface[14]= false;
		isInterfaceSurface[15]= false;

		newTriangleEdgeIDs[3]={18,19,21,18,20,22,19,20,23,21,22,23};
		
		// Now we add the elements, edges and triangles 

		// Adding Elements
		int offsetElements = this->elementsC_->numberElements(); 
		int offsetEdges = this->edgeElements_->numberElements(); 
		for( int i=0;i<4; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("irregular");	
			feNew.setPredecessorElement(indexElement);
			if(i<3)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		// Adding the edges (they also have to be added to triangles as subelements, but that is not implemented yet)
		for( int i=0;i<24; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<18){
				this->edgeElements_->addEdge(feNew,i/6+offsetElements);
			}
			else
				this->edgeElements_->addEdge(feNew,indexElement);		
		}

		// Adding triangles as subelements, if they arent interior triangles
		int offsetSurface=0;
		for( int i=0;i<16; i++){
			sort( newTriangles.at(i).begin(), newTriangles.at(i).end() );
			FiniteElement feNew(newTriangles[i],newTrianglesFlag[i]);
			feNew.setInterfaceElement(isInterfaceSurface[i]);
			if(i<12){
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
				 	if ( !this->elementsC_->getElement(i/4+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/4+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/4+offsetElements).addSubElement(feNew);					
				}
				this->surfaceTriangleElements_->addSurface(feNew, i/4+offsetElements);
			}
			else{
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	
				}
				this->surfaceTriangleElements_->addSurface(feNew, indexElement);
			}							
		}
		FiniteElement element;
		FiniteElement feEdge;
		for( int i=0;i<4; i++){	
			if(i<3)	
				element = this->elementsC_->getElement(i+offsetElements);
			else
				element = this->elementsC_->getElement(indexElement);
			bool init=false;
			for(int j=0; j<24 ; j++){
				FiniteElement feEdge = this->edgeElements_->getElement(j+offsetEdges);
				if(feEdge.getFlag() != this->volumeID_){
					if(init == true)
						element.addSubElement( feEdge );
					else if ( !element.subElementsInitialized() ){
				        element.initializeSubElements( "P1", 1 ); // only P1 for now                
				        element.addSubElement( feEdge );
				        init= true;
				    }
				    else {
				        ElementsPtr_Type surfaces = element.getSubElements();
				        // We set the edge to the corresponding element(s)
				        surfaces->setToCorrectElement( feEdge );
				    }		
				}
			}
		}
	
	}
	else 
		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The Type 4 irregular Refinement Method you requested is only applicable to a 3 dimensional Mesh. Please reconsider.");

  
}
/*!

 \brief 3D Type(3) refinement as defined in  "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.
@param[in] surfaceTriangleElements Triangle elements.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineType3(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements){


// Implementation of Type (3) Refinement Type
// We use this Refinement Type 
	if(this->dim_ == 3){ 

		// The way we refine the Tetrahedron is defined by how we order the nodes of the tetrahedron
		// (For the algorithm see "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955)
		
        vec_int_Type midPointInd(0); // indices of midpoints of edges of soon to be refined element
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		vec_int_Type edgeNumbersUntagged(0);
		// Extract the three points of tetraeder, that connect the tagged edges
		vec_int_Type nodeInd(0);
		for(int i=0; i<6; i++)	{
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
				nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
			}
			else
				edgeNumbersUntagged.push_back(edgeNumbers[i]);
		}
		sort( nodeInd.begin(), nodeInd.end() );
		vec_int_Type nodeIndTmp = nodeInd;
		nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );
		
		// We determine the common node index, this one will be placed at x_0
		int commonNode=-1;
		int entryCommonNode;
		vec_int_Type leftOverNodes(0);
		for(int i=0; i<3; i++){
			if(nodeIndTmp[i] == nodeIndTmp[i+1]){
				commonNode=nodeIndTmp[i];
				entryCommonNode=i;
			}

		}
		for(int i=0; i<3; i++){
			if(nodeInd[i] != commonNode){
				leftOverNodes.push_back(nodeInd[i]);
			}

		}
			
		// Then we determine the longest of the tagged edges, of which the second node will be placed at x_1 and the leftover node at x_2
		vec_dbl_Type length(2);
		vec_dbl_Type P1(3),P2(3);
		double maxLength=0.0;
		int maxEntry=0;
		int minEntry=0;
		LO p1ID,p2ID;
		for(int i=0;i<2;i++){
			p1ID =commonNode;
			p2ID =leftOverNodes[i];
			P1 = this->pointsRep_->at(p1ID);
			P2 = this->pointsRep_->at(p2ID);
			length[i] = sqrt(pow(P1[0]-P2[0],2)+pow(P1[1]-P2[1],2)+pow(P1[2]-P2[2],2));
			if(length[i] > maxLength){
				maxLength = length[i];
				maxEntry= i;
			}
			else
				minEntry=i;
		}
		nodeInd[0] = commonNode;
		nodeInd[1] = leftOverNodes[minEntry];
		nodeInd[2] = leftOverNodes[maxEntry];

		// We now have the two nodes that connect the tagged edge or rather make up the tagged edge
		// The left over nodes are the ones opposite to the tagged edge
		for(int i=0; i<4; i++){
			if(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(0) == nodeInd[0])
				nodeInd.push_back(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(1));
			else if(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(1) == nodeInd[0])
				nodeInd.push_back(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(0));
		}
		
		// We don't need to sort the nodes in any other way then we already did, as there is not more than one possibility to construct the subtetrahera in refinement type 2

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
		vec_int_Type surfaceElementsIDs = surfaceTriangleElements->getSurfacesOfElement(indexElement); // surfaces of Element k
		vec2D_int_Type originTriangles(4,vec_int_Type(3));

		originTriangles[0] = {nodeInd[0], nodeInd[1], nodeInd[2] };
		originTriangles[1] = {nodeInd[0], nodeInd[1], nodeInd[3] };
		originTriangles[2] = {nodeInd[0], nodeInd[2], nodeInd[3] };
		originTriangles[3] = {nodeInd[1], nodeInd[2], nodeInd[3] };
		
		
	
		vec_int_Type originFlag(4,this->volumeID_); // Triangle Flag

		vec_bool_Type interfaceSurface(4);
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);

		for(int i=0; i< 4 ; i++){
			originTriangleTmp = originTriangles[i];
			sort(originTriangleTmp.begin(),originTriangleTmp.end());
			for(int j=0; j<4 ; j++){
				FiniteElement surfaceTmp = surfaceTriangleElements->getElement(surfaceElementsIDs[j]);
				triTmp = surfaceTmp.getVectorNodeList();
				sort(triTmp.begin(),triTmp.end());
				if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] && triTmp[2] == originTriangleTmp[2]  ) {
					originFlag[i] = surfaceTmp.getFlag();
					interfaceSurface[i] = surfaceTmp.isInterfaceElement();
				}
			}
		}

		// Finally we need to determine or extract the indices of the edges midpoints. As in the describe algorithm the midpoints are set as follows:
		// Edge 0 = [x_0,x_1] -> x_01
		// Edge 1 = [x_0,x_2] -> x_02
		// Edge 2 = [x_0,x_3] -> x_03
		// Edge 3 = [x_1,x_2] -> x_12
		// Edge 4 = [x_1,x_3] -> x_13
		// Edge 5 = [x_2,x_3] -> x_23

		for(int i=0; i<6; i++){
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				midPointInd.push_back(edgeElements->getMidpoint(edgeNumbers[i]));
			}
		}
		
	
		// Now we construct the new Elements as proposed by Bey's Regular Refinement Algorithm

		vec2D_int_Type newElements(3, vec_int_Type( 0 )); // new elements
		vec2D_int_Type newEdges(18,vec_int_Type(0)); // new edges
		vec2D_int_Type newTriangles(12,vec_int_Type(0)); // new Triangles
		vec_int_Type newTrianglesFlag(12) ; // new Triangle Flags
		vec_int_Type isInterfaceSurface(12);
		vec_int_Type edgeFlags(18); // new Edge flags
		vec_bool_Type isInterfaceEdge(18); // bool vector for interfaceEdge Information
		vec_GO_Type predecessorElement(18); // vector that stores the global IDs of the predecessor of each edge 

		vec2D_LO_Type newTriangleEdgeIDs(3,vec_LO_Type(12));

		// How are Flags determined?
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)
		// If an edges emerges on a triangle, the flag is determined by the triangle flag. Opposite to the 2D case, edges that connect midpoints are not automatically interior edges, but are
		// determined by the triangle/surface they are on


		// Element 1: (x_0,x_3,x_01,x_02)
		(newElements)[0]={nodeInd[0],nodeInd[3], midPointInd[0],midPointInd[1]};

		(newEdges)[0] = {nodeInd[0] ,nodeInd[3]}; 
		(newEdges)[1] = {nodeInd[0] ,midPointInd[0]}; 
		(newEdges)[2] = {nodeInd[0] ,midPointInd[1]}; 
		(newEdges)[3] = {nodeInd[3] ,midPointInd[0]}; 
		(newEdges)[4] = {nodeInd[3] ,midPointInd[1]}; 
		(newEdges)[5] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[0]=edgeElements->getElement(edgeNumbers[2]).getFlag();
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[2]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[3]=originFlag[1];
		edgeFlags[4]=originFlag[2];
		edgeFlags[5]=originFlag[0];

		isInterfaceEdge[0] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[1] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[2] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[3] = interfaceSurface[1];
		isInterfaceEdge[4] = interfaceSurface[2];
		isInterfaceEdge[5] = interfaceSurface[0];

		predecessorElement[0] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[1] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[2] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[3] = -1;
		predecessorElement[4] = -1;
		predecessorElement[5] = -1;

		// Subelements of thetrahedron
		newTriangles[0]= {nodeInd[0],nodeInd[3],midPointInd[0]};
		newTriangles[1]= {nodeInd[0],nodeInd[3],midPointInd[1]};
		newTriangles[2]= {nodeInd[0],midPointInd[0],midPointInd[1]};
		newTriangles[3]= {nodeInd[3],midPointInd[0],midPointInd[1]};

		newTrianglesFlag[0]= originFlag[1]; 
		newTrianglesFlag[1]= originFlag[2]; 
		newTrianglesFlag[2]= originFlag[0]; 
		newTrianglesFlag[3]= this->volumeID_;

		isInterfaceSurface[0]= interfaceSurface[1];
		isInterfaceSurface[1]= interfaceSurface[2];
		isInterfaceSurface[2]= interfaceSurface[0];
		isInterfaceSurface[3]= false;

		newTriangleEdgeIDs[0]={0,1,3,0,2,4,1,2,5,3,4,5};

		// Element 2: (x_1,x_2,x_3,x_01)
		(newElements)[1]={nodeInd[1],nodeInd[2],nodeInd[3],midPointInd[0]};

		(newEdges)[6] = {nodeInd[1] ,nodeInd[2]}; 
		(newEdges)[7] = {nodeInd[1] ,nodeInd[3]}; 
		(newEdges)[8] = {nodeInd[1] ,midPointInd[0]}; 
		(newEdges)[9] = {nodeInd[2] ,nodeInd[3]}; 
		(newEdges)[10] = {nodeInd[2] ,midPointInd[0]}; 
		(newEdges)[11] = {nodeInd[3] ,midPointInd[0]}; 

		edgeFlags[6]=edgeElements->getElement(edgeNumbers[3]).getFlag();
		edgeFlags[7]=edgeElements->getElement(edgeNumbers[4]).getFlag();
		edgeFlags[8]=edgeElements->getElement(edgeNumbers[0]).getFlag();
		edgeFlags[9]=edgeElements->getElement(edgeNumbers[5]).getFlag();
		edgeFlags[10]=originFlag[0];
		edgeFlags[11]=originFlag[1];

		isInterfaceEdge[6] = edgeElements->getElement(edgeNumbers[3]).isInterfaceElement();
		isInterfaceEdge[7] = edgeElements->getElement(edgeNumbers[4]).isInterfaceElement();
		isInterfaceEdge[8] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[9] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[10] = interfaceSurface[0];
		isInterfaceEdge[11] = interfaceSurface[1];

		predecessorElement[6] = this->edgeMap_->getGlobalElement(edgeNumbers[3]);
		predecessorElement[7] = this->edgeMap_->getGlobalElement(edgeNumbers[4]);
		predecessorElement[8] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[9] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[10] = -1;
		predecessorElement[11] = -1;

		// Subelements of tetrahedron
		newTriangles[4]= {nodeInd[1],nodeInd[2],nodeInd[3]};
		newTriangles[5]= {nodeInd[1],nodeInd[2],midPointInd[0]};
		newTriangles[6]= {nodeInd[1],nodeInd[3],midPointInd[0]};
		newTriangles[7]= {nodeInd[2],nodeInd[3],midPointInd[0]};

		newTrianglesFlag[4]= originFlag[3]; 
		newTrianglesFlag[5]= originFlag[0]; 
		newTrianglesFlag[6]= originFlag[1]; 
		newTrianglesFlag[7]= this->volumeID_;

		isInterfaceSurface[4]= interfaceSurface[3];
		isInterfaceSurface[5]= interfaceSurface[0];
		isInterfaceSurface[6]= interfaceSurface[1];
		isInterfaceSurface[7]= false;

		newTriangleEdgeIDs[1]={6,7,9,6,8,10,7,8,11,9,10,11};
		// Element 3: (x_2,x_3,x_01,x_02)
		(newElements)[2]={nodeInd[2],nodeInd[3],midPointInd[0],midPointInd[1]};

		(newEdges)[12] = {nodeInd[2] ,nodeInd[3]}; 
		(newEdges)[13] = {nodeInd[2] ,midPointInd[0]}; 
		(newEdges)[14] = {nodeInd[2] ,midPointInd[1]}; 
		(newEdges)[15] = {nodeInd[3] ,midPointInd[0]}; 
		(newEdges)[16] = {nodeInd[3] ,midPointInd[1]}; 
		(newEdges)[17] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[12]=edgeElements->getElement(edgeNumbers[5]).getFlag();
		edgeFlags[13]=originFlag[0];
		edgeFlags[14]=edgeElements->getElement(edgeNumbers[1]).getFlag();
		edgeFlags[15]=originFlag[1];
		edgeFlags[16]=originFlag[2];
		edgeFlags[17]=originFlag[0];

		isInterfaceEdge[12] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[13] = interfaceSurface[0];
		isInterfaceEdge[14] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[15] = interfaceSurface[1];
		isInterfaceEdge[16] = interfaceSurface[2];
		isInterfaceEdge[17] = interfaceSurface[0];

		predecessorElement[12] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[13] = -1;
		predecessorElement[14] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[15] = -1;
		predecessorElement[16] = -1;
		predecessorElement[17] = -1;

		// Subelements of tetrahedron
		newTriangles[8]= {nodeInd[2],nodeInd[3],midPointInd[0]};
		newTriangles[9]= {nodeInd[2],nodeInd[3],midPointInd[1]};
		newTriangles[10]= {nodeInd[2],midPointInd[0],midPointInd[1]};
		newTriangles[11]= {nodeInd[3],midPointInd[0],midPointInd[1]};

		newTrianglesFlag[8]=this->volumeID_; 
		newTrianglesFlag[9]= originFlag[2]; 
		newTrianglesFlag[10]= originFlag[0]; 
		newTrianglesFlag[11]= this->volumeID_;

		isInterfaceSurface[8]= false;
		isInterfaceSurface[9]= interfaceSurface[2];
		isInterfaceSurface[10]= interfaceSurface[0];
		isInterfaceSurface[11]= false;

		newTriangleEdgeIDs[2]={8,9,11,8,10,12,9,10,13,11,12,13};

		
		
		// Now we add the elements, edges and triangles 

		// Adding Elements
		int offsetElements = this->elementsC_->numberElements(); 
		int offsetEdges = this->edgeElements_->numberElements(); 
		for( int i=0;i<3; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("irregular");	
			feNew.setPredecessorElement(indexElement);
			if(i<2)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		// Adding the edges (they also have to be added to triangles as subelements, but that is not implemented yet)
		for( int i=0;i<18; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<12){
				this->edgeElements_->addEdge(feNew,i/6+offsetElements);
			}
			else
				this->edgeElements_->addEdge(feNew,indexElement);		
		}

		// Adding triangles as subelements, if they arent interior triangles
		int offsetSurface =0;
		for( int i=0;i<12; i++){
			sort( newTriangles.at(i).begin(), newTriangles.at(i).end() );
			FiniteElement feNew(newTriangles[i],newTrianglesFlag[i]);
			feNew.setInterfaceElement(isInterfaceSurface[i]);
			if(i<8){
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
				 	if ( !this->elementsC_->getElement(i/4+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/4+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/4+offsetElements).addSubElement(feNew);					
				}
				this->surfaceTriangleElements_->addSurface(feNew, i/4+offsetElements);
			}
			else{
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	
				}
				this->surfaceTriangleElements_->addSurface(feNew, indexElement);
			}					
		}
		FiniteElement element;
		FiniteElement feEdge;
		for( int i=0;i<3; i++){	
			if(i<2)	
				element = this->elementsC_->getElement(i+offsetElements);
			else
				element = this->elementsC_->getElement(indexElement);
			bool init=false;
			for(int j=0; j<18 ; j++){
				FiniteElement feEdge = this->edgeElements_->getElement(j+offsetEdges);
				if(feEdge.getFlag() != this->volumeID_){
					if(init == true)
						element.addSubElement( feEdge );
					else if ( !element.subElementsInitialized() ){
				        element.initializeSubElements( "P1", 1 ); // only P1 for now                
				        element.addSubElement( feEdge );
				        init= true;
				    }
				    else {
				        ElementsPtr_Type surfaces = element.getSubElements();
				        // We set the edge to the corresponding element(s)
				        surfaces->setToCorrectElement( feEdge );
				    }		
				}
			}
		}
	
	}
	else 
		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The Type 1 irregular Refinement Method you requested is only applicable to a 3 dimensional Mesh. Please reconsider.");

  
}


/*!

 \brief 3D Type(2) refinement as defined in  "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.
@param[in] surfaceTriangleElements Triangle elements.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineType2(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements){

// Implementation of Type (2) Refinement Type
// We use this Refinement Type 
	if(this->dim_ == 3){ 

		// The way we refine the Tetrahedron is defined by how we order the nodes of the tetrahedron
		// (For the algorithm see "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955)
		// The Type 2 Refinement is similar to a green Refinement in two dimensions, as we connect the midpoint to the opposite points on the same surface
		// The procedure is similar to the regular refinement, we just add less elements

        vec_int_Type midPointInd( 1 ); // indices of midpoints of edges of soon to be refined element
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		vec_int_Type edgeNumbersUntagged(0);
		// Extract the three points of tetraeder, that connect the tagged edges
		vec_int_Type nodeInd(0);
		for(int i=0; i<6; i++)	{
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
				nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
			}
			else
				edgeNumbersUntagged.push_back(edgeNumbers[i]);
		}
		sort( nodeInd.begin(), nodeInd.end() );
		nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );
		
		// We now have the two nodes that connect the tagged edge or rather make up the tagged edge
		// The left over nodes are the ones opposite to the tagged edge
		vec_int_Type nodeIndTmp(0);
		for(int i=0; i<5; i++){
			if(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(0) == nodeInd[0])
				nodeIndTmp.push_back(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(1));
			else if(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(1) == nodeInd[0])
				nodeIndTmp.push_back(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(0));
			else if(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(0) == nodeInd[1])
				nodeIndTmp.push_back(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(1));
			else if(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(1) == nodeInd[1])
				nodeIndTmp.push_back(edgeElements->getElement(edgeNumbersUntagged[i]).getNode(0));
		}

		sort( nodeIndTmp.begin(), nodeIndTmp.end() );
		nodeIndTmp.erase( unique( nodeIndTmp.begin(), nodeIndTmp.end() ), nodeIndTmp.end() );
	
		// Now we add the two remaining points, that dont belong to the tagged edge
		nodeInd.push_back(nodeIndTmp[0]);
		nodeInd.push_back(nodeIndTmp[1]);
		
		// We don't need to sort the nodes in any other way then we already did, as there is not more than one possibility to construct the subtetrahera in refinement type 2

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
		vec_int_Type surfaceElementsIDs = surfaceTriangleElements->getSurfacesOfElement(indexElement); // surfaces of Element k
		vec2D_int_Type originTriangles(4,vec_int_Type(3));

		originTriangles[0] = {nodeInd[0], nodeInd[1], nodeInd[2] };
		originTriangles[1] = {nodeInd[0], nodeInd[1], nodeInd[3] };
		originTriangles[2] = {nodeInd[0], nodeInd[2], nodeInd[3] };
		originTriangles[3] = {nodeInd[1], nodeInd[2], nodeInd[3] };
		
		
	
		vec_int_Type originFlag(4,this->volumeID_); // Triangle Flag

		vec_bool_Type interfaceSurface(4);
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);

		for(int i=0; i< 4 ; i++){
			originTriangleTmp = originTriangles[i];
			sort(originTriangleTmp.begin(),originTriangleTmp.end());
			for(int j=0; j<4 ; j++){
				FiniteElement surfaceTmp = surfaceTriangleElements->getElement(surfaceElementsIDs[j]);
				triTmp = surfaceTmp.getVectorNodeList();
				sort(triTmp.begin(),triTmp.end());
				if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] && triTmp[2] == originTriangleTmp[2]  ) {
					originFlag[i] = surfaceTmp.getFlag();
					interfaceSurface[i] = surfaceTmp.isInterfaceElement();
				}
			}
		}

		// Finally we need to determine or extract the indices of the edges midpoints. As in the describe algorithm the midpoints are set as follows:
		// Edge 0 = [x_0,x_1] -> x_01
		// Edge 1 = [x_0,x_2] -> x_02
		// Edge 2 = [x_0,x_3] -> x_03
		// Edge 3 = [x_1,x_2] -> x_12
		// Edge 4 = [x_1,x_3] -> x_13
		// Edge 5 = [x_2,x_3] -> x_23

		for(int i=0; i<6; i++){
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				midPointInd[0] = edgeElements->getMidpoint(edgeNumbers[i]);
			}
		}
		
	
		// Now we construct the new Elements as proposed by Bey's Regular Refinement Algorithm

		vec2D_int_Type newElements(2, vec_int_Type( 0 )); // new elements
		vec2D_int_Type newEdges(12,vec_int_Type(0)); // new edges
		vec2D_int_Type newTriangles(8,vec_int_Type(0)); // new Triangles
		vec_int_Type newTrianglesFlag(8) ; // new Triangle Flags
		vec_bool_Type isInterfaceSurface(8); // bool vector for interfaceSurface Information
		vec_int_Type edgeFlags(12); // new Edge flags
		vec_bool_Type isInterfaceEdge(12); // bool vector for interfaceEdge Information
		vec_GO_Type predecessorElement(12); // vector that stores the global IDs of the predecessor of each edge 

		vec2D_LO_Type newTriangleEdgeIDs(2,vec_LO_Type(12));

		// How are Flags determined?
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)
		// If an edges emerges on a triangle, the flag is determined by the triangle flag. Opposite to the 2D case, edges that connect midpoints are not automatically interior edges, but are
		// determined by the triangle/surface they are on


		// Element 1: (x_0,x_2,x_3,x_01)
		(newElements)[0]={nodeInd[0],nodeInd[2],nodeInd[3],midPointInd[0]};

		(newEdges)[0] = {nodeInd[0] ,nodeInd[2]}; 
		(newEdges)[1] = {nodeInd[0] ,nodeInd[3]}; 
		(newEdges)[2] = {nodeInd[0] ,midPointInd[0]}; 
		(newEdges)[3] = {nodeInd[2] ,nodeInd[3]}; 
		(newEdges)[4] = {nodeInd[2] ,midPointInd[0]}; 
		(newEdges)[5] = {nodeInd[3] ,midPointInd[0]}; 

		edgeFlags[0]=edgeElements->getElement(edgeNumbers[1]).getFlag();
		edgeFlags[1]=edgeElements->getElement(edgeNumbers[2]).getFlag();
		edgeFlags[2]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[3]=edgeElements->getElement(edgeNumbers[5]).getFlag();
		edgeFlags[4]=originFlag[0];
		edgeFlags[5]=originFlag[1];

		isInterfaceEdge[0] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[1] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[2] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[3] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();;
		isInterfaceEdge[4] = interfaceSurface[0];
		isInterfaceEdge[5] = interfaceSurface[1];

		predecessorElement[0] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[1] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[2] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[3] =  this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[4] = -1;
		predecessorElement[5] = -1;

		// Subelements of thetrahedron
		newTriangles[0]= {nodeInd[0],nodeInd[2],nodeInd[3]};
		newTriangles[1]= {nodeInd[0],nodeInd[2],midPointInd[0]};
		newTriangles[2]= {nodeInd[0],nodeInd[3],midPointInd[0]};
		newTriangles[3]= {nodeInd[2],nodeInd[3],midPointInd[0]};

		newTrianglesFlag[0]= originFlag[2]; 
		newTrianglesFlag[1]= originFlag[0]; 
		newTrianglesFlag[2]= originFlag[1]; 
		newTrianglesFlag[3]= this->volumeID_;

		isInterfaceSurface[0] = interfaceSurface[2];
		isInterfaceSurface[1] = interfaceSurface[0];
		isInterfaceSurface[2] = interfaceSurface[1];
		isInterfaceSurface[3] = false;

		newTriangleEdgeIDs[0]={0,1,3,0,2,4,1,2,5,3,4,5};

		// Element 2: (x_1,x_2,x_3,x_01)
		(newElements)[1]={nodeInd[1],nodeInd[2],nodeInd[3],midPointInd[0]};

		(newEdges)[6] = {nodeInd[1] ,nodeInd[2]}; 
		(newEdges)[7] = {nodeInd[1] ,nodeInd[3]}; 
		(newEdges)[8] = {nodeInd[1] ,midPointInd[0]}; 
		(newEdges)[9] = {nodeInd[2] ,nodeInd[3]}; 
		(newEdges)[10] = {nodeInd[2] ,midPointInd[0]}; 
		(newEdges)[11] = {nodeInd[3] ,midPointInd[0]}; 

		edgeFlags[6]=edgeElements->getElement(edgeNumbers[3]).getFlag();
		edgeFlags[7]=edgeElements->getElement(edgeNumbers[4]).getFlag();
		edgeFlags[8]=edgeElements->getElement(edgeNumbers[0]).getFlag();
		edgeFlags[9]=edgeElements->getElement(edgeNumbers[5]).getFlag();
		edgeFlags[10]=originFlag[0];
		edgeFlags[11]=originFlag[1];

		isInterfaceEdge[6] = edgeElements->getElement(edgeNumbers[3]).isInterfaceElement();
		isInterfaceEdge[7] = edgeElements->getElement(edgeNumbers[4]).isInterfaceElement();
		isInterfaceEdge[8] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[9] = edgeElements->getElement(edgeNumbers[5]).isInterfaceElement();
		isInterfaceEdge[10] = interfaceSurface[0];
		isInterfaceEdge[11] = interfaceSurface[1];

		predecessorElement[6] = this->edgeMap_->getGlobalElement(edgeNumbers[3]);
		predecessorElement[7] = this->edgeMap_->getGlobalElement(edgeNumbers[4]);
		predecessorElement[8] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[9] = this->edgeMap_->getGlobalElement(edgeNumbers[5]);
		predecessorElement[10] = -1;
		predecessorElement[11] = -1;

		// Subelements of tetrahedron
		newTriangles[4]= {nodeInd[1],nodeInd[2],nodeInd[3]};
		newTriangles[5]= {nodeInd[1],nodeInd[2],midPointInd[0]};
		newTriangles[6]= {nodeInd[1],nodeInd[3],midPointInd[0]};
		newTriangles[7]= {nodeInd[2],nodeInd[3],midPointInd[0]};

		newTrianglesFlag[4]= originFlag[3]; 
		newTrianglesFlag[5]= originFlag[0]; 
		newTrianglesFlag[6]= originFlag[1]; 
		newTrianglesFlag[7]= this->volumeID_;

		isInterfaceSurface[4] = interfaceSurface[3];
		isInterfaceSurface[5] = interfaceSurface[0];
		isInterfaceSurface[6] = interfaceSurface[1];
		isInterfaceSurface[7] = false;

		newTriangleEdgeIDs[1]={6,7,9,6,8,10,7,8,11,9,10,11};
		
		// Now we add the elements, edges and triangles 

		// Adding Elements
		int offsetElements = this->elementsC_->numberElements(); 
		int offsetEdges = this->edgeElements_->numberElements(); 
		for( int i=0;i<2; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("irregular");
			if(refinementMode_ == "Bisection")
				feNew.setFiniteElementRefinementType("regular");
			feNew.setPredecessorElement(indexElement);
			if(i<1)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		// Adding the edges (they also have to be added to triangles as subelements, but that is not implemented yet)
		for( int i=0;i<12; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<6){
				this->edgeElements_->addEdge(feNew,i/6+offsetElements);
			}
			else
				this->edgeElements_->addEdge(feNew,indexElement);		
		}

		// Adding triangles as subelements, if they arent interior triangles
		int offsetSurface =0;
		for( int i=0;i<8; i++){
			sort( newTriangles.at(i).begin(), newTriangles.at(i).end() );
			FiniteElement feNew(newTriangles[i],newTrianglesFlag[i]);
			feNew.setInterfaceElement(isInterfaceSurface[i]);
			if(i<4){
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
				 	if ( !this->elementsC_->getElement(i/4+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/4+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/4+offsetElements).addSubElement(feNew);					
				}
				this->surfaceTriangleElements_->addSurface(feNew, i/4+offsetElements);	
			}			
			else{
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	
				}
				this->surfaceTriangleElements_->addSurface(feNew, indexElement);				
			}	

		}
		FiniteElement element;
		FiniteElement feEdge;
		for( int i=0;i<2; i++){	
			if(i<1)	
				element = this->elementsC_->getElement(i+offsetElements);
			else
				element = this->elementsC_->getElement(indexElement);
			bool init=false;
			for(int j=0; j<12 ; j++){
				FiniteElement feEdge = this->edgeElements_->getElement(j+offsetEdges);
				if(feEdge.getFlag() != this->volumeID_){
					if(init == true)
						element.addSubElement( feEdge );
					else if ( !element.subElementsInitialized() ){
				        element.initializeSubElements( "P1", 1 ); // only P1 for now                
				        element.addSubElement( feEdge );
				        init= true;
				    }
				    else {
				        ElementsPtr_Type surfaces = element.getSubElements();
				        // We set the edge to the corresponding element(s)
				        surfaces->setToCorrectElement( feEdge );
				    }		
				}
			}
		}	
	}
	else 
		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The Type 1 irregular Refinement Method you requested is only applicable to a 3 dimensional Mesh. Please reconsider.");

  
}

/*!

 \brief 3D Type(1) refinement as defined in  "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.
@param[in] surfaceTriangleElements Triangle elements.

*/

template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineType1(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements){

// Implementation of Type (1) Refinement Type
// We use this Refinement Type 
	if(this->dim_ == 3){ 

		// The way we refine the Tetrahedron is defined by how we order the nodes of the tetrahedron
		// (For the algorithm see "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955)
		// The Type 1 Refinement is similar to a red Refinement in two dimensions, as we connect the midpoints on one surface and connect them to the remaining point
		// The procedure is similar to the regular refinement, we just add less elements

        vec_int_Type midPointInd(0); // indices of midpoints of edges of soon to be refined element
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		vec_int_Type edgeNumbersUntagged(0);
		// Extract the three points of tetraeder, that connect the tagged edges
		vec_int_Type nodeInd(0);
		int node4=0;
		for(int i=0; i<6; i++)	{
			if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()){
				nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
				nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
			}
			else
				edgeNumbersUntagged.push_back(edgeNumbers[i]);
		}
		sort( nodeInd.begin(), nodeInd.end() );
		nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );
		
		// We now have the three nodes that connect the tagged edges or rather make up the tagged triangle
		// The fourth node that is opposite to this tagged triangle we push back to the nodeInd list
		for(int i=0; i<3; i++){
			if(edgeElements->getElement(edgeNumbersUntagged[0]).getNode(0) == nodeInd[i])
				nodeInd.push_back(edgeElements->getElement(edgeNumbersUntagged[0]).getNode(1));
			if(edgeElements->getElement(edgeNumbersUntagged[0]).getNode(1) == nodeInd[i])
				nodeInd.push_back(edgeElements->getElement(edgeNumbersUntagged[0]).getNode(0));
		}
		
		// We don't need to sort the nodes in any other way then we already did, as there is not more than one possibility to construct the subtetrahera in refinement type 1

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

		vec_int_Type surfaceElementsIDs = surfaceTriangleElements->getSurfacesOfElement(indexElement); // surfaces of Element k
		vec2D_int_Type originTriangles(4,vec_int_Type(3));

		originTriangles[0] = {nodeInd[0], nodeInd[1], nodeInd[2] };
		originTriangles[1] = {nodeInd[0], nodeInd[1], nodeInd[3] };
		originTriangles[2] = {nodeInd[0], nodeInd[2], nodeInd[3] };
		originTriangles[3] = {nodeInd[1], nodeInd[2], nodeInd[3] };
		
		vec_int_Type originFlag(4,this->volumeID_); // Triangle Flag

		vec_bool_Type interfaceSurface(4);
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);

		for(int i=0; i< 4 ; i++){
			originTriangleTmp = originTriangles[i];
			sort(originTriangleTmp.begin(),originTriangleTmp.end());
			for(int j=0; j<4 ; j++){
				FiniteElement surfaceTmp = surfaceTriangleElements->getElement(surfaceElementsIDs[j]);
				triTmp = surfaceTmp.getVectorNodeList();
				sort(triTmp.begin(),triTmp.end());
				if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] && triTmp[2] == originTriangleTmp[2]  ) {
					originFlag[i] = surfaceTmp.getFlag();
					interfaceSurface[i] = surfaceTmp.isInterfaceElement();
				}
			}
		}		
		// Finally we need to determine or extract the indices of the edges midpoints. As in the describe algorithm the midpoints are set as follows:
		// Edge 0 = [x_0,x_1] -> x_01
		// Edge 1 = [x_0,x_2] -> x_02
		// Edge 2 = [x_0,x_3] -> x_03
		// Edge 3 = [x_1,x_2] -> x_12
		// Edge 4 = [x_1,x_3] -> x_13
		// Edge 5 = [x_2,x_3] -> x_23

		for(int i=0; i<6; i++)	{
		  if(edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement())
			midPointInd.push_back(edgeElements->getMidpoint(edgeNumbers[i]));
		}
	
		// Now we construct the new Elements as proposed by Bey's Regular Refinement Algorithm

		vec2D_int_Type newElements(4, vec_int_Type( 0 )); // new elements
		vec2D_int_Type newEdges(24,vec_int_Type(0)); // new edges
		vec2D_int_Type newTriangles(16,vec_int_Type(0)); // new Triangles
		vec_int_Type newTrianglesFlag(16) ; // new Triangle Flags
		vec_bool_Type isInterfaceSurface(16); // bool vector for interfaceSurface Information
		vec_int_Type edgeFlags(24); // new Edge flags
		vec_bool_Type isInterfaceEdge(24); // bool vector for interfaceEdge Information
		vec_GO_Type predecessorElement(24); // vector that stores the global IDs of the predecessor of each edge 

		vec2D_LO_Type newTriangleEdgeIDs(4,vec_LO_Type(12));
		
		// How are Flags determined?
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)
		// If an edges emerges on a triangle, the flag is determined by the triangle flag. Opposite to the 2D case, edges that connect midpoints are not automatically interior edges, but are
		// determined by the triangle/surface they are on


		// Element 1: (x_0,x_01,x_02,x_03)
		(newElements)[0]={nodeInd[0],midPointInd[0],midPointInd[1],nodeInd[3]};

		(newEdges)[0] = {nodeInd[0] ,midPointInd[0]}; 
		(newEdges)[1] = {nodeInd[0] ,midPointInd[1]}; 
		(newEdges)[2] = {nodeInd[0] ,nodeInd[3]}; 
		(newEdges)[3] = {midPointInd[0] ,midPointInd[1]}; 
		(newEdges)[4] = {midPointInd[0] ,nodeInd[3]}; 
		(newEdges)[5] = {midPointInd[1] ,nodeInd[3]}; 

		edgeFlags[0]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[2]=edgeElements->getElement(edgeNumbers[2]).getFlag();
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
		newTriangles[1]= {nodeInd[0],midPointInd[0],nodeInd[3]};
		newTriangles[2]= {nodeInd[0],midPointInd[1],nodeInd[3]};
		newTriangles[3]= {midPointInd[0],midPointInd[1],nodeInd[3]};

		newTrianglesFlag[0]= originFlag[0]; 
		newTrianglesFlag[1]= originFlag[1]; 
		newTrianglesFlag[2]= originFlag[2]; 
		newTrianglesFlag[3]= this->volumeID_;

		isInterfaceSurface[0] = interfaceSurface[0];
		isInterfaceSurface[1] = interfaceSurface[1];
		isInterfaceSurface[2] = interfaceSurface[2];
		isInterfaceSurface[3] = false;

		newTriangleEdgeIDs[0]={0,1,3,0,2,4,1,2,5,3,4,5};

		// Element 2: (x_1,x_01,x_12,x_13)
		(newElements)[1]={nodeInd[1],midPointInd[0],midPointInd[2],nodeInd[3]};

		(newEdges)[6] = {nodeInd[1] ,midPointInd[0]}; 
		(newEdges)[7] = {nodeInd[1] ,midPointInd[2]}; 
		(newEdges)[8] = {nodeInd[1] ,nodeInd[3]}; 
		(newEdges)[9] = {midPointInd[0] ,midPointInd[2]}; 
		(newEdges)[10] = {midPointInd[0] ,nodeInd[3]}; 
		(newEdges)[11] = {midPointInd[2] ,nodeInd[3]}; 

		edgeFlags[6]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[7]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[8]=edgeElements->getElement(edgeNumbers[4]).getFlag();
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
		newTriangles[4]= {nodeInd[1],midPointInd[0],midPointInd[2]};
		newTriangles[5]= {nodeInd[1],midPointInd[0],nodeInd[3]};
		newTriangles[6]= {nodeInd[1],midPointInd[2],nodeInd[3]};
		newTriangles[7]= {midPointInd[0],midPointInd[2],nodeInd[3]};

		newTrianglesFlag[4]= originFlag[0]; 
		newTrianglesFlag[5]= originFlag[1]; 
		newTrianglesFlag[6]= originFlag[3]; 
		newTrianglesFlag[7]= this->volumeID_;

		isInterfaceSurface[4] = interfaceSurface[0];
		isInterfaceSurface[5] = interfaceSurface[1];
		isInterfaceSurface[6] = interfaceSurface[3];
		isInterfaceSurface[7] = false;

		newTriangleEdgeIDs[1]={6,7,9,6,8,10,7,8,11,9,10,11};


		// Element 3: (x_2,x_02,x_12,x_23)
		(newElements)[2]={nodeInd[2],midPointInd[1],midPointInd[2],nodeInd[3]};

		(newEdges)[12] = {nodeInd[2] ,midPointInd[1]}; 
		(newEdges)[13] = {nodeInd[2] ,midPointInd[2]}; 
		(newEdges)[14] = {nodeInd[2] ,nodeInd[3]}; 
		(newEdges)[15] = {midPointInd[1] ,midPointInd[2]}; 
		(newEdges)[16] = {midPointInd[1] ,nodeInd[3]}; 
		(newEdges)[17] = {midPointInd[2] ,nodeInd[3]}; 

		edgeFlags[12]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[13]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[14]=edgeElements->getElement(edgeNumbers[5]).getFlag();
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
		newTriangles[8]= {nodeInd[2],midPointInd[1],midPointInd[2]};
		newTriangles[9]= {nodeInd[2],midPointInd[1],nodeInd[3]};
		newTriangles[10]= {nodeInd[2],midPointInd[2],nodeInd[3]};
		newTriangles[11]= {midPointInd[1],midPointInd[2],nodeInd[3]};

		newTrianglesFlag[8]= originFlag[0]; 
		newTrianglesFlag[9]= originFlag[2]; 
		newTrianglesFlag[10]= originFlag[3]; 
		newTrianglesFlag[11]= this->volumeID_;

		isInterfaceSurface[8] = interfaceSurface[0];
		isInterfaceSurface[9] = interfaceSurface[2];
		isInterfaceSurface[10] = interfaceSurface[3];
		isInterfaceSurface[11] = false;

		newTriangleEdgeIDs[2]={12,13,15,12,14,16,13,14,17,15,16,17};


		// Element 4: (x_3,x_03,x_13,x_23)
		(newElements)[3]={nodeInd[3],midPointInd[0],midPointInd[1],midPointInd[2]};

		(newEdges)[18] = {nodeInd[3] ,midPointInd[0]}; 
		(newEdges)[19] = {nodeInd[3] ,midPointInd[1]}; 
		(newEdges)[20] = {nodeInd[3] ,midPointInd[2]}; 
		(newEdges)[21] = {midPointInd[0] ,midPointInd[1]}; 
		(newEdges)[22] = {midPointInd[0] ,midPointInd[2]}; 
		(newEdges)[23] = {midPointInd[1] ,midPointInd[2]}; 

		edgeFlags[18]=originFlag[1];
		edgeFlags[19]=originFlag[2];
		edgeFlags[20]=originFlag[3];
		edgeFlags[21]=originFlag[0];
		edgeFlags[22]=originFlag[0];
		edgeFlags[23]=originFlag[0];

		isInterfaceEdge[18] = interfaceSurface[1];
		isInterfaceEdge[19] = interfaceSurface[2];
		isInterfaceEdge[20] = interfaceSurface[3];
		isInterfaceEdge[21] = interfaceSurface[0];
		isInterfaceEdge[22] = interfaceSurface[0];
		isInterfaceEdge[23] = interfaceSurface[0];

		predecessorElement[18] = -1;
		predecessorElement[19] = -1;
		predecessorElement[20] = -1;
		predecessorElement[21] = -1;
		predecessorElement[22] = -1;
		predecessorElement[23] = -1;


		// Subelements of thetrahedron
		newTriangles[12]= {nodeInd[3],midPointInd[0],midPointInd[1]};
		newTriangles[13]= {nodeInd[3],midPointInd[0],midPointInd[2]};
		newTriangles[14]= {nodeInd[3],midPointInd[1],midPointInd[2]};
		newTriangles[15]= {midPointInd[0],midPointInd[1],midPointInd[2]};

		newTrianglesFlag[12]= this->volumeID_; 
		newTrianglesFlag[13]= this->volumeID_; 
		newTrianglesFlag[14]= this->volumeID_; 
		newTrianglesFlag[15]= originFlag[0];

		isInterfaceSurface[12] = false;
		isInterfaceSurface[13] = false;
		isInterfaceSurface[14] = false;
		isInterfaceSurface[15] = interfaceSurface[0];
		

		newTriangleEdgeIDs[3]={18,19,21,18,20,22,19,20,23,21,22,23};
		// Now we add the elements, edges and triangles 

		// Adding Elements
		int offsetElements = this->elementsC_->numberElements(); 
		int offsetEdges = this->edgeElements_->numberElements(); 
		for( int i=0;i<4; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("irregular");	
			feNew.setPredecessorElement(indexElement);
			if(i<3)
				this->elementsC_->addElement(feNew);
			else
				this->elementsC_->switchElement(indexElement,feNew);
		}

		// Adding the edges (they also have to be added to triangles as subelements, but that is not implemented yet)
		for( int i=0;i<24; i++){
			sort( newEdges.at(i).begin(), newEdges.at(i).end() );
			FiniteElement feNew(newEdges.at(i),edgeFlags[i]);
			feNew.setInterfaceElement(isInterfaceEdge[i]);
			feNew.setPredecessorElement(predecessorElement[i]);
			if(i<18){
				this->edgeElements_->addEdge(feNew,i/6+offsetElements);
			}
			else
				this->edgeElements_->addEdge(feNew,indexElement);		
		}

		// Adding triangles as subelements, if they arent interior triangles
		int offsetSurface=0;
		for( int i=0;i<16; i++){
			sort( newTriangles.at(i).begin(), newTriangles.at(i).end() );
			FiniteElement feNew(newTriangles[i],newTrianglesFlag[i]);
			feNew.setInterfaceElement(isInterfaceSurface[i]);
			if(i<12){
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
				 	if ( !this->elementsC_->getElement(i/4+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/4+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/4+offsetElements).addSubElement(feNew);	
				}	
				this->surfaceTriangleElements_->addSurface(feNew, i/4+offsetElements);				
			}
			else{
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	
				}
				this->surfaceTriangleElements_->addSurface(feNew, indexElement);				

			}	
		}
		FiniteElement element;
		FiniteElement feEdge;
		for( int i=0;i<4; i++){	
			if(i<3)	
				element = this->elementsC_->getElement(i+offsetElements);
			else
				element = this->elementsC_->getElement(indexElement);
			bool init=false;
			for(int j=0; j<24 ; j++){
				FiniteElement feEdge = this->edgeElements_->getElement(j+offsetEdges);
				if(feEdge.getFlag() != this->volumeID_){
					if(init == true)
						element.addSubElement( feEdge );
					else if ( !element.subElementsInitialized() ){
				        element.initializeSubElements( "P1", 1 ); // only P1 for now                
				        element.addSubElement( feEdge );
				        init= true;
				    }
				    else {
				        ElementsPtr_Type surfaces = element.getSubElements();
				        // We set the edge to the corresponding element(s)
				        surfaces->setToCorrectElement( feEdge );
				    }		
				}
			}
		}
	
	}
	else 
		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The Type 1 irregular Refinement Method you requested is only applicable to a 3 dimensional Mesh. Please reconsider.");

  
}
    
/*!

 \brief 2D and 3D regular refinement. Chosen by error estimator or otherwise elements are refined regular by connecting edge midpoints.
 \brief 3D regular refinement as defined in  "Tetrahedral Grid Refinement" by J. Bey 'Algorithm Regular Refinement' in Computing, Springer Verlag 1955

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.
@param[in] surfaceTriangleElements Triangle elements.

*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::refineRegular(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements){

	if(this->dim_ == 2){
        vec_int_Type midPointInd( 3 ); // indices of midpoints of edges of soon to be refined element
      
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		int k=0;
		for(int i=0; i<3; i++)	{
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
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)

		// Element 1
		(newElements)[0]={mutualNode[0],midPointInd[0],midPointInd[1]};

		(newEdges)[0] = {mutualNode[0] ,midPointInd[0]}; 
		(newEdges)[1] = {mutualNode[0] ,midPointInd[1]}; 
		(newEdges)[2] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[0]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[2]=this->volumeID_;

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
		edgeFlags[5]=this->volumeID_;

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
		edgeFlags[8]=this->volumeID_;

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

		edgeFlags[9]=this->volumeID_;
		edgeFlags[10]=this->volumeID_;
		edgeFlags[11]=this->volumeID_;

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
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("regular");	
			feNew.setPredecessorElement(indexElement);
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
				if(edgeFlags[i]!=0 && edgeFlags[i]!=this->volumeID_){
				 	if ( !this->elementsC_->getElement(i/3+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/3+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/3+offsetElements).addSubElement(feNew);		

					}		
			}
			else
				this->edgeElements_->addEdge(feNew,indexElement);
			
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
		/*vec_int_Type nodeInd(0);
		for(int i=0; i<6; i++)	{
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(0));
			nodeInd.push_back(edgeElements->getElement(edgeNumbers[i]).getNode(1));
		}
		sort( nodeInd.begin(), nodeInd.end() );
		nodeInd.erase( unique( nodeInd.begin(), nodeInd.end() ), nodeInd.end() );
		*/
		vec_int_Type nodeInd = elements->getElement(indexElement).getVectorNodeList();

		vec2D_dbl_ptr_Type pointsRep = this->pointsRep_;
		// Right now the Nodes are ordered by the local Indices, which differ depending on the number of Processors. Hence the refinements are not
		// 100% equal when using a different number of processors.
		// If we sort the nodes by their values and not their local Indices, we can solve that problem, as these values don't change
		// This can also be done by using the global IDs of nodes

		/*vec2D_dbl_Type points(4,vec_dbl_Type(4));
		points[0] = {pointsRep->at(nodeInd[0]).at(0),pointsRep->at(nodeInd[0]).at(1),pointsRep->at(nodeInd[0]).at(2), (double) nodeInd[0]};
		points[1] = {pointsRep->at(nodeInd[1]).at(0),pointsRep->at(nodeInd[1]).at(1),pointsRep->at(nodeInd[1]).at(2), (double) nodeInd[1]};
		points[2] = {pointsRep->at(nodeInd[2]).at(0),pointsRep->at(nodeInd[2]).at(1),pointsRep->at(nodeInd[2]).at(2), (double) nodeInd[2]};
		points[3] = {pointsRep->at(nodeInd[3]).at(0),pointsRep->at(nodeInd[3]).at(1),pointsRep->at(nodeInd[3]).at(2), (double) nodeInd[3]};

		sort(points.begin(), points.end());

		nodeInd[0] = (int) points[0][3];
		nodeInd[1] = (int) points[1][3];
		nodeInd[2] = (int) points[2][3];
		nodeInd[3] = (int) points[3][3];*/

		// NOTE: The above described way of ensuring 100% equal refinement on a varing number of processors is disabled, as we want to keep natural element properties (Sorting such that determinant of element is always > 0)

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

		// We have 4 Triangles in our Tetraedron
		// If one or more of those Triangles are Part of the domains' boundaries, they are added to the element in question as subelement
		// We extract them in follwing pattern:

		// Tri_0 = [x_0,x_1,x_2]
		// Tri_1 = [x_0,x_1,x_3]
		// Tri_2 = [x_0,x_2,x_3]
		// Tri_3 = [x_1,x_2,x_3]

		// We check if one or more of these triangles are part of the boundary surface and determine their flag
		vec_int_Type surfaceElementsIDs = surfaceTriangleElements->getSurfacesOfElement(indexElement); // surfaces of Element k
		vec2D_int_Type originTriangles(4,vec_int_Type(3));

		originTriangles[0] = {nodeInd[0], nodeInd[1], nodeInd[2] };
		originTriangles[1] = {nodeInd[0], nodeInd[1], nodeInd[3] };
		originTriangles[2] = {nodeInd[0], nodeInd[2], nodeInd[3] };
		originTriangles[3] = {nodeInd[1], nodeInd[2], nodeInd[3] };
		
		
	
		vec_int_Type originFlag(4,this->volumeID_); // Triangle Flag

		vec_bool_Type interfaceSurface(4);
		vec_LO_Type triTmp(3);
		vec_int_Type originTriangleTmp(3);

		for(int i=0; i< 4 ; i++){
			originTriangleTmp = originTriangles[i];
			sort(originTriangleTmp.begin(),originTriangleTmp.end());
			for(int j=0; j<4 ; j++){
				FiniteElement surfaceTmp = surfaceTriangleElements->getElement(surfaceElementsIDs[j]);
				triTmp = surfaceTmp.getVectorNodeList();
				sort(triTmp.begin(),triTmp.end());
				if(triTmp[0] == originTriangleTmp[0] && triTmp[1] == originTriangleTmp[1] && triTmp[2] == originTriangleTmp[2]  ) {
					originFlag[i] = surfaceTmp.getFlag();
					interfaceSurface[i] = surfaceTmp.isInterfaceElement();
				}
			}
		}


		// Furthermore we have to determine whether the triangles are part of the interface between processors, as we need this information to determine if edges
		// that emerge on the triangles are part of the interface
		// A triangle is part of the interface if all of its edges are part of the interface (the information if edges are part of the interface was determined
		// at the beginning of the Mesh Refinement by 'determineInterfaceEdges')



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
				this->addMidpoint(edgeElements,edgeNumbers[i]);
				midPointInd[i] = edgeElements->getMidpoint(edgeNumbers[i]); 
			}
			else{
				midPointInd[i] = edgeElements->getMidpoint(edgeNumbers[i]);
				
			}
	
		}
		
	
		// Now we construct the new Elements as proposed by Bey's Regular Refinement Algorithm

		vec2D_int_Type newElements(8, vec_int_Type( 0 )); // new elements
		vec2D_int_Type newEdges(48,vec_int_Type(0)); // new edges
		vec2D_int_Type newTriangles(32,vec_int_Type(0)); // new Triangles
		vec_int_Type newTrianglesFlag(32) ; // new Triangle Flags
		vec_int_Type isInterfaceSurface(32);
		vec_int_Type edgeFlags(48); // new Edge flags
		vec_bool_Type isInterfaceEdge(48); // bool vector for interfaceEdge Information
		vec_GO_Type predecessorElement(48); // vector that stores the global IDs of the predecessor of each edge 

		vec2D_LO_Type newTriangleEdgeIDs(8,vec_LO_Type(12));

		// How are Flags determined?
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)
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
		newTriangles[0]= {nodeInd[0],midPointInd[1],midPointInd[0]};
		newTriangles[1]= {nodeInd[0],midPointInd[0],midPointInd[2]};
		newTriangles[2]= {nodeInd[0],midPointInd[2],midPointInd[1]};
		newTriangles[3]= {midPointInd[0],midPointInd[1],midPointInd[2]};

		newTrianglesFlag[0]= originFlag[0]; 
		newTrianglesFlag[1]= originFlag[1]; 
		newTrianglesFlag[2]= originFlag[2]; 
		newTrianglesFlag[3]= this->volumeID_;

		isInterfaceSurface[0] = interfaceSurface[0];
		isInterfaceSurface[1] = interfaceSurface[1];
		isInterfaceSurface[2] = interfaceSurface[2];
		isInterfaceSurface[3] = false;

		newTriangleEdgeIDs[0]={0,1,3,0,2,4,1,2,5,3,4,5};
		// Element 2: (x_01,x_1,x_12,x_13)
		(newElements)[1]={midPointInd[0],nodeInd[1],midPointInd[3],midPointInd[4]};

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
		newTriangles[5]= {nodeInd[1],midPointInd[4],midPointInd[0]};
		newTriangles[6]= {nodeInd[1],midPointInd[3],midPointInd[4]};
		newTriangles[7]= {midPointInd[0],midPointInd[3],midPointInd[4]};

		newTrianglesFlag[4]= originFlag[0]; 
		newTrianglesFlag[5]= originFlag[1]; 
		newTrianglesFlag[6]= originFlag[3]; 
		newTrianglesFlag[7]= this->volumeID_;

		isInterfaceSurface[4] = interfaceSurface[0];
		isInterfaceSurface[5] = interfaceSurface[1];
		isInterfaceSurface[6] = interfaceSurface[3];
		isInterfaceSurface[7] = false;

		newTriangleEdgeIDs[1]={6,7,9,6,8,10,7,8,11,9,10,11};
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
		newTriangles[8]= {nodeInd[2],midPointInd[3],midPointInd[1]};
		newTriangles[9]= {nodeInd[2],midPointInd[1],midPointInd[5]};
		newTriangles[10]= {nodeInd[2],midPointInd[5],midPointInd[3]};
		newTriangles[11]= {midPointInd[1],midPointInd[3],midPointInd[5]};

		newTrianglesFlag[8]= originFlag[0]; 
		newTrianglesFlag[9]= originFlag[2]; 
		newTrianglesFlag[10]= originFlag[3]; 
		newTrianglesFlag[11]= this->volumeID_;

		isInterfaceSurface[8] = interfaceSurface[0];
		isInterfaceSurface[9] = interfaceSurface[2];
		isInterfaceSurface[10] = interfaceSurface[3];
		isInterfaceSurface[11] = false;

		newTriangleEdgeIDs[2]={12,13,15,12,14,16,13,14,17,15,16,17};

		// Element 4: (x_3,x_03,x_13,x_23)
		(newElements)[3]={midPointInd[2],midPointInd[4],midPointInd[5],nodeInd[3]};

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
		newTriangles[13]= {nodeInd[3],midPointInd[5],midPointInd[2]};
		newTriangles[14]= {nodeInd[3],midPointInd[4],midPointInd[5]};
		newTriangles[15]= {midPointInd[2],midPointInd[4],midPointInd[5]};

		newTrianglesFlag[12]= originFlag[1]; 
		newTrianglesFlag[13]= originFlag[2]; 
		newTrianglesFlag[14]= originFlag[3]; 
		newTrianglesFlag[15]= this->volumeID_;

		isInterfaceSurface[12] = interfaceSurface[1];
		isInterfaceSurface[13] = interfaceSurface[2];
		isInterfaceSurface[14] = interfaceSurface[3];
		isInterfaceSurface[15] = false;


		newTriangleEdgeIDs[3]={18,19,20,18,20,22,19,21,23,21,22,23};

		// The following elements are constructed only with the edges midpoints, hence they have the surface flag and interface characteristic
		// In order to control a certain mesh quality the diagonal cut, which is usually between midpoint[1] and midPoint[4]
		vec_dbl_Type lengthDia(3);

		int diaInd=0;

		lengthDia[0] = sqrt(pow(pointsRep->at(midPointInd[1]).at(0) - pointsRep->at(midPointInd[4]).at(0),2) + pow(pointsRep->at(midPointInd[1]).at(1) - pointsRep->at(midPointInd[4]).at(1),2) +pow(pointsRep->at(midPointInd[1]).at(2) - pointsRep->at(midPointInd[4]).at(2),2) );

		lengthDia[1] = sqrt(pow(pointsRep->at(midPointInd[0]).at(0) - pointsRep->at(midPointInd[5]).at(0),2) + pow(pointsRep->at(midPointInd[0]).at(1) - pointsRep->at(midPointInd[5]).at(1),2) +pow(pointsRep->at(midPointInd[0]).at(2) - pointsRep->at(midPointInd[5]).at(2),2) );

		lengthDia[2] = sqrt(pow(pointsRep->at(midPointInd[2]).at(0) - pointsRep->at(midPointInd[3]).at(0),2) + pow(pointsRep->at(midPointInd[2]).at(1) - pointsRep->at(midPointInd[3]).at(1),2) +pow(pointsRep->at(midPointInd[2]).at(2) - pointsRep->at(midPointInd[3]).at(2),2) );


		vec2D_dbl_Type dia(3,vec_dbl_Type(2));
		dia[0] = { lengthDia[0],0.};
		dia[1] = { lengthDia[1],1.};
		dia[2] = { lengthDia[2],2.};
		sort(dia.begin(),dia.end());		

		// Diagonal 0 represents the shortest
		// Diagonal 1 represents the second shortest
		// Diagonal 2 represents the longest

		// If the Diagonal is not within that range we allways use diaInd =0. Same Diagonal with consistent indexing of elements
		if(refinement3DDiagonal_>2 || refinement3DDiagonal_ == 0)
			diaInd =0;
		else {
			diaInd = (int) dia[refinement3DDiagonal_][1];
			TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The Refinement feature to pick a specific internal Diagonal in regular refinement is not yet orientation preserving!");

		}


		// Element 5: (x_01,x_02,x_03,x_12)
		if(diaInd ==0){
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
			edgeFlags[28]=this->volumeID_;
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
			newTriangles[18]= {midPointInd[0],midPointInd[4],midPointInd[2]};
			newTriangles[19]= {midPointInd[1],midPointInd[2],midPointInd[4]};


			newTrianglesFlag[16]= this->volumeID_; 
			newTrianglesFlag[17]= this->volumeID_; 
			newTrianglesFlag[18]= originFlag[1]; 
			newTrianglesFlag[19]= this->volumeID_;

			isInterfaceSurface[16] = false;
			isInterfaceSurface[17] = false;
			isInterfaceSurface[18] = interfaceSurface[1];
			isInterfaceSurface[19] = false;
		
			newTriangleEdgeIDs[4]={24,25,27,24,26,28,25,26,29,27,28,29};

			// Element 6: (x_01,x_02,x_12,x_13)
			(newElements)[5]={midPointInd[0],midPointInd[3],midPointInd[1],midPointInd[4]};

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
			edgeFlags[34]=this->volumeID_;
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
			newTrianglesFlag[21]= this->volumeID_; 
			newTrianglesFlag[22]= this->volumeID_; 
			newTrianglesFlag[23]= this->volumeID_;

			isInterfaceSurface[20] = interfaceSurface[0];
			isInterfaceSurface[21] = false;
			isInterfaceSurface[22] = false;
			isInterfaceSurface[23] = false;

			newTriangleEdgeIDs[5]={30,31,33,30,32,34,31,32,35,33,34,35};
			// Element 7: (x_02,x_03,x_13,x_23)
			(newElements)[6]={midPointInd[1],midPointInd[2],midPointInd[4],midPointInd[5]};

			(newEdges)[36] = {midPointInd[1] ,midPointInd[2]}; 
			(newEdges)[37] = {midPointInd[1] ,midPointInd[4]}; 
			(newEdges)[38] = {midPointInd[1] ,midPointInd[5]}; 
			(newEdges)[39] = {midPointInd[2] ,midPointInd[4]}; 
			(newEdges)[40] = {midPointInd[2] ,midPointInd[5]}; 
			(newEdges)[41] = {midPointInd[4] ,midPointInd[5]}; 

			edgeFlags[36]=originFlag[2];
			edgeFlags[37]=this->volumeID_;
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

			newTrianglesFlag[24]= this->volumeID_; 
			newTrianglesFlag[25]= originFlag[2]; 
			newTrianglesFlag[26]= this->volumeID_; 
			newTrianglesFlag[27]= this->volumeID_;

			isInterfaceSurface[24] = false;
			isInterfaceSurface[25] = interfaceSurface[2];
			isInterfaceSurface[26] = false;
			isInterfaceSurface[27] = false;

			newTriangleEdgeIDs[6]={36,37,39,36,38,40,37,38,41,39,40,41};

			// Element 8: (x_02,x_12,x_13,x_23)
			(newElements)[7]={midPointInd[1],midPointInd[3],midPointInd[5],midPointInd[4]};

			(newEdges)[42] = {midPointInd[1] ,midPointInd[3]}; 
			(newEdges)[43] = {midPointInd[1] ,midPointInd[4]}; 
			(newEdges)[44] = {midPointInd[1] ,midPointInd[5]}; 
			(newEdges)[45] = {midPointInd[3] ,midPointInd[4]}; 
			(newEdges)[46] = {midPointInd[3] ,midPointInd[5]}; 
			(newEdges)[47] = {midPointInd[4] ,midPointInd[5]}; 

			edgeFlags[42]=originFlag[0];
			edgeFlags[43]=this->volumeID_;
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
			newTriangles[31]= {midPointInd[3],midPointInd[5],midPointInd[4]};

			newTrianglesFlag[28]= this->volumeID_; 
			newTrianglesFlag[29]= this->volumeID_; 
			newTrianglesFlag[30]= this->volumeID_; 
			newTrianglesFlag[31]= originFlag[3];

			isInterfaceSurface[28] = false;
			isInterfaceSurface[29] = false;
			isInterfaceSurface[30] = false;
			isInterfaceSurface[31] = interfaceSurface[3];

			newTriangleEdgeIDs[7]={42,43,45,42,44,46,43,44,47,45,46,47};
		}

		else if(diaInd ==1){
			// Element 5: 
			(newElements)[4]={midPointInd[0],midPointInd[1],midPointInd[2],midPointInd[5]};

			(newEdges)[24] = {midPointInd[0] ,midPointInd[1]}; 
			(newEdges)[25] = {midPointInd[0] ,midPointInd[2]}; 
			(newEdges)[26] = {midPointInd[0] ,midPointInd[5]}; 
			(newEdges)[27] = {midPointInd[1] ,midPointInd[2]}; 
			(newEdges)[28] = {midPointInd[1] ,midPointInd[5]}; 
			(newEdges)[29] = {midPointInd[2] ,midPointInd[5]}; 

			edgeFlags[24]=originFlag[0];
			edgeFlags[25]=originFlag[1];
			edgeFlags[26]=this->volumeID_;
			edgeFlags[27]=originFlag[2];
			edgeFlags[28]=originFlag[2];
			edgeFlags[29]=originFlag[2];

			isInterfaceEdge[24] = interfaceSurface[0];
			isInterfaceEdge[25] = interfaceSurface[1];
			isInterfaceEdge[26] = false;
			isInterfaceEdge[27] = interfaceSurface[2];
			isInterfaceEdge[28] = interfaceSurface[2];
			isInterfaceEdge[29] = interfaceSurface[2];
		
			
			predecessorElement[24] = -1;
			predecessorElement[25] = -1;
			predecessorElement[26] = -1;
			predecessorElement[27] = -1;
			predecessorElement[28] = -1;
			predecessorElement[29] = -1;

			// Subelements of thetrahedron
			newTriangles[16]= {midPointInd[0],midPointInd[1],midPointInd[2]};
			newTriangles[17]= {midPointInd[0],midPointInd[1],midPointInd[5]};
			newTriangles[18]= {midPointInd[0],midPointInd[2],midPointInd[5]};
			newTriangles[19]= {midPointInd[1],midPointInd[2],midPointInd[5]};


			newTrianglesFlag[16]= this->volumeID_; 
			newTrianglesFlag[17]= this->volumeID_; 
			newTrianglesFlag[18]= this->volumeID_; 
			newTrianglesFlag[19]= originFlag[2];

			isInterfaceSurface[16] = false;
			isInterfaceSurface[17] = false;
			isInterfaceSurface[18] = false;
			isInterfaceSurface[19] = interfaceSurface[2];;
		
			newTriangleEdgeIDs[4]={24,25,27,24,26,28,25,26,29,27,28,29};

			// Element 6: (x_01,x_02,x_12,x_13)
			(newElements)[5]={midPointInd[0],midPointInd[3],midPointInd[4],midPointInd[5]};

			(newEdges)[30] = {midPointInd[0],midPointInd[3]}; 
			(newEdges)[31] = {midPointInd[0] ,midPointInd[4]}; 
			(newEdges)[32] = {midPointInd[0],midPointInd[5]}; 
			(newEdges)[33] = {midPointInd[3] ,midPointInd[4]}; 
			(newEdges)[34] = {midPointInd[3] ,midPointInd[5]}; 
			(newEdges)[35] = {midPointInd[4] ,midPointInd[5]}; 

			edgeFlags[30]=originFlag[0];
			edgeFlags[31]=originFlag[1];
			edgeFlags[32]=this->volumeID_;
			edgeFlags[33]=originFlag[3];
			edgeFlags[34]=originFlag[3];
			edgeFlags[35]=originFlag[3];

			isInterfaceEdge[30] = interfaceSurface[0];
			isInterfaceEdge[31] = interfaceSurface[1];
			isInterfaceEdge[32] = false;
			isInterfaceEdge[33] = interfaceSurface[3];
			isInterfaceEdge[34] = interfaceSurface[3];
			isInterfaceEdge[35] = interfaceSurface[3];

			predecessorElement[30] = -1;
			predecessorElement[31] = -1;
			predecessorElement[32] = -1;
			predecessorElement[33] = -1;
			predecessorElement[34] = -1;
			predecessorElement[35] = -1;


			// Subelements of thetrahedron
			newTriangles[20]= {midPointInd[0],midPointInd[3],midPointInd[4]};
			newTriangles[21]= {midPointInd[0],midPointInd[3],midPointInd[5]};
			newTriangles[22]= {midPointInd[0],midPointInd[4],midPointInd[5]};
			newTriangles[23]= {midPointInd[3],midPointInd[4],midPointInd[5]};

			newTrianglesFlag[20]= this->volumeID_; 
			newTrianglesFlag[21]= this->volumeID_; 
			newTrianglesFlag[22]= this->volumeID_; 
			newTrianglesFlag[23]= originFlag[3];

			isInterfaceSurface[20] = false;
			isInterfaceSurface[21] = false;
			isInterfaceSurface[22] = false;
			isInterfaceSurface[23] = interfaceSurface[3];

			newTriangleEdgeIDs[5]={30,31,33,30,32,34,31,32,35,33,34,35};

			// Element 7: (x_02,x_03,x_13,x_23)
			(newElements)[6]={midPointInd[0],midPointInd[1],midPointInd[3],midPointInd[5]};

			(newEdges)[36] = {midPointInd[0] ,midPointInd[1]}; 
			(newEdges)[37] = {midPointInd[0] ,midPointInd[3]}; 
			(newEdges)[38] = {midPointInd[0] ,midPointInd[5]}; 
			(newEdges)[39] = {midPointInd[1] ,midPointInd[3]}; 
			(newEdges)[40] = {midPointInd[1] ,midPointInd[5]}; 
			(newEdges)[41] = {midPointInd[3] ,midPointInd[5]}; 

			edgeFlags[36]=originFlag[0];
			edgeFlags[37]=originFlag[0];
			edgeFlags[38]=this->volumeID_;
			edgeFlags[39]=originFlag[0];
			edgeFlags[40]=originFlag[2];
			edgeFlags[41]=originFlag[3];

			isInterfaceEdge[36] = interfaceSurface[0];
			isInterfaceEdge[37] = interfaceSurface[0];
			isInterfaceEdge[38] = false;
			isInterfaceEdge[39] = interfaceSurface[0];
			isInterfaceEdge[40] = interfaceSurface[2];
			isInterfaceEdge[41] = interfaceSurface[3];

			predecessorElement[36] = -1;
			predecessorElement[37] = -1;
			predecessorElement[38] = -1;
			predecessorElement[39] = -1;
			predecessorElement[40] = -1;
			predecessorElement[41] = -1;

			// Subelements of thetrahedron
			newTriangles[24]= {midPointInd[0],midPointInd[1],midPointInd[3]};
			newTriangles[25]= {midPointInd[0],midPointInd[1],midPointInd[5]};
			newTriangles[26]= {midPointInd[0],midPointInd[3],midPointInd[5]};
			newTriangles[27]= {midPointInd[1],midPointInd[3],midPointInd[5]};

			newTrianglesFlag[24]= originFlag[0]; 
			newTrianglesFlag[25]= this->volumeID_; 
			newTrianglesFlag[26]= this->volumeID_; 
			newTrianglesFlag[27]= this->volumeID_;

			isInterfaceSurface[24] = interfaceSurface[0];
			isInterfaceSurface[25] = false;
			isInterfaceSurface[26] = false;
			isInterfaceSurface[27] = false;

			newTriangleEdgeIDs[6]={36,37,39,36,38,40,37,38,41,39,40,41};

			// Element 8: (x_02,x_12,x_13,x_23)
			(newElements)[7]={midPointInd[0],midPointInd[4],midPointInd[2],midPointInd[5]};

			(newEdges)[42] = {midPointInd[0] ,midPointInd[2]}; 
			(newEdges)[43] = {midPointInd[0] ,midPointInd[4]}; 
			(newEdges)[44] = {midPointInd[0] ,midPointInd[5]}; 
			(newEdges)[45] = {midPointInd[2] ,midPointInd[4]}; 
			(newEdges)[46] = {midPointInd[2] ,midPointInd[5]}; 
			(newEdges)[47] = {midPointInd[4] ,midPointInd[5]}; 

			edgeFlags[42]=originFlag[1];
			edgeFlags[43]=originFlag[1];
			edgeFlags[44]=this->volumeID_;
			edgeFlags[45]=originFlag[1];
			edgeFlags[46]=originFlag[2];
			edgeFlags[47]=originFlag[3];

			isInterfaceEdge[42] = interfaceSurface[1];
			isInterfaceEdge[43] = interfaceSurface[1]; 
			isInterfaceEdge[44] = false;
			isInterfaceEdge[45] = interfaceSurface[1];
			isInterfaceEdge[46] = interfaceSurface[2];
			isInterfaceEdge[47] = interfaceSurface[3];

			predecessorElement[42] = -1;
			predecessorElement[43] = -1;
			predecessorElement[44] = -1;
			predecessorElement[45] = -1;
			predecessorElement[46] = -1;
			predecessorElement[47] = -1;

			// Subelements of thetrahedron
			newTriangles[28]= {midPointInd[0],midPointInd[2],midPointInd[4]};
			newTriangles[29]= {midPointInd[0],midPointInd[2],midPointInd[5]};
			newTriangles[30]= {midPointInd[0],midPointInd[4],midPointInd[5]};
			newTriangles[31]= {midPointInd[2],midPointInd[4],midPointInd[5]};

			newTrianglesFlag[28]= originFlag[1]; 
			newTrianglesFlag[29]= this->volumeID_; 
			newTrianglesFlag[30]= this->volumeID_; 
			newTrianglesFlag[31]= this->volumeID_;

			isInterfaceSurface[28] = interfaceSurface[1];
			isInterfaceSurface[29] = false;
			isInterfaceSurface[30] = false;
			isInterfaceSurface[31] = false;

			newTriangleEdgeIDs[7]={42,43,45,42,44,46,43,44,47,45,46,47};
		}
		else if(diaInd ==2){
			// Element 5: 
			(newElements)[4]={midPointInd[0],midPointInd[1],midPointInd[2],midPointInd[3]};

			(newEdges)[24] = {midPointInd[0] ,midPointInd[1]}; 
			(newEdges)[25] = {midPointInd[0] ,midPointInd[2]}; 
			(newEdges)[26] = {midPointInd[0] ,midPointInd[3]}; 
			(newEdges)[27] = {midPointInd[1] ,midPointInd[2]}; 
			(newEdges)[28] = {midPointInd[1] ,midPointInd[3]}; 
			(newEdges)[29] = {midPointInd[2] ,midPointInd[3]}; 

			edgeFlags[24]=originFlag[0];
			edgeFlags[25]=originFlag[1];
			edgeFlags[26]=originFlag[0];
			edgeFlags[27]=originFlag[2];
			edgeFlags[28]=originFlag[0];
			edgeFlags[29]=this->volumeID_;

			isInterfaceEdge[24] = interfaceSurface[0];
			isInterfaceEdge[25] = interfaceSurface[1];
			isInterfaceEdge[26] = interfaceSurface[0];
			isInterfaceEdge[27] = interfaceSurface[2];
			isInterfaceEdge[28] = interfaceSurface[0];
			isInterfaceEdge[29] = false;
		
			
			predecessorElement[24] = -1;
			predecessorElement[25] = -1;
			predecessorElement[26] = -1;
			predecessorElement[27] = -1;
			predecessorElement[28] = -1;
			predecessorElement[29] = -1;

			// Subelements of thetrahedron
			newTriangles[16]= {midPointInd[0],midPointInd[1],midPointInd[2]};
			newTriangles[17]= {midPointInd[0],midPointInd[1],midPointInd[3]};
			newTriangles[18]= {midPointInd[0],midPointInd[2],midPointInd[3]};
			newTriangles[19]= {midPointInd[1],midPointInd[2],midPointInd[3]};


			newTrianglesFlag[16]= this->volumeID_; 
			newTrianglesFlag[17]= originFlag[0]; 
			newTrianglesFlag[18]= this->volumeID_; 
			newTrianglesFlag[19]= this->volumeID_;

			isInterfaceSurface[16] = false;
			isInterfaceSurface[17] = interfaceSurface[0];
			isInterfaceSurface[18] = false;
			isInterfaceSurface[19] = false;
		
			newTriangleEdgeIDs[4]={24,25,27,24,26,28,25,26,29,27,28,29};

			// Element 6: (x_01,x_02,x_12,x_13)
			(newElements)[5]={midPointInd[0],midPointInd[2],midPointInd[3],midPointInd[4]};

			(newEdges)[30] = {midPointInd[0],midPointInd[2]}; 
			(newEdges)[31] = {midPointInd[0] ,midPointInd[3]}; 
			(newEdges)[32] = {midPointInd[0],midPointInd[4]}; 
			(newEdges)[33] = {midPointInd[2] ,midPointInd[3]}; 
			(newEdges)[34] = {midPointInd[2] ,midPointInd[4]}; 
			(newEdges)[35] = {midPointInd[3] ,midPointInd[4]}; 

			edgeFlags[30]=originFlag[1];
			edgeFlags[31]=originFlag[0];
			edgeFlags[32]=originFlag[1];
			edgeFlags[33]=this->volumeID_;
			edgeFlags[34]=originFlag[1];
			edgeFlags[35]=originFlag[3];

			isInterfaceEdge[30] = interfaceSurface[1];
			isInterfaceEdge[31] = interfaceSurface[0];
			isInterfaceEdge[32] = interfaceSurface[1];
			isInterfaceEdge[33] = false;
			isInterfaceEdge[34] = interfaceSurface[1];
			isInterfaceEdge[35] = interfaceSurface[3];

			predecessorElement[30] = -1;
			predecessorElement[31] = -1;
			predecessorElement[32] = -1;
			predecessorElement[33] = -1;
			predecessorElement[34] = -1;
			predecessorElement[35] = -1;


			// Subelements of thetrahedron
			newTriangles[20]= {midPointInd[0],midPointInd[2],midPointInd[3]};
			newTriangles[21]= {midPointInd[0],midPointInd[2],midPointInd[4]};
			newTriangles[22]= {midPointInd[0],midPointInd[3],midPointInd[4]};
			newTriangles[23]= {midPointInd[2],midPointInd[3],midPointInd[4]};

			newTrianglesFlag[20]= this->volumeID_;
			newTrianglesFlag[21]= originFlag[1]; 
			newTrianglesFlag[22]= this->volumeID_; 
			newTrianglesFlag[23]= this->volumeID_;

			isInterfaceSurface[20] = false;
			isInterfaceSurface[21] = interfaceSurface[1];
			isInterfaceSurface[22] = false;
			isInterfaceSurface[23] = false;

			newTriangleEdgeIDs[5]={30,31,33,30,32,34,31,32,35,33,34,35};
			// Element 7: (x_02,x_03,x_13,x_23)
			(newElements)[6]={midPointInd[1],midPointInd[2],midPointInd[3],midPointInd[5]};

			(newEdges)[36] = {midPointInd[1] ,midPointInd[2]}; 
			(newEdges)[37] = {midPointInd[1] ,midPointInd[3]}; 
			(newEdges)[38] = {midPointInd[1] ,midPointInd[5]}; 
			(newEdges)[39] = {midPointInd[2] ,midPointInd[3]}; 
			(newEdges)[40] = {midPointInd[2] ,midPointInd[5]}; 
			(newEdges)[41] = {midPointInd[3] ,midPointInd[5]}; 

			edgeFlags[36]=originFlag[2];
			edgeFlags[37]=originFlag[0];
			edgeFlags[38]=originFlag[2];
			edgeFlags[39]=this->volumeID_;
			edgeFlags[40]=originFlag[2];
			edgeFlags[41]=originFlag[3];

			isInterfaceEdge[36] = interfaceSurface[2];
			isInterfaceEdge[37] = interfaceSurface[0];
			isInterfaceEdge[38] = interfaceSurface[2];
			isInterfaceEdge[39] = false;
			isInterfaceEdge[40] = interfaceSurface[2];
			isInterfaceEdge[41] = interfaceSurface[3];

			predecessorElement[36] = -1;
			predecessorElement[37] = -1;
			predecessorElement[38] = -1;
			predecessorElement[39] = -1;
			predecessorElement[40] = -1;
			predecessorElement[41] = -1;

			// Subelements of thetrahedron
			newTriangles[24]= {midPointInd[1],midPointInd[2],midPointInd[3]};
			newTriangles[25]= {midPointInd[1],midPointInd[2],midPointInd[5]};
			newTriangles[26]= {midPointInd[1],midPointInd[3],midPointInd[5]};
			newTriangles[27]= {midPointInd[2],midPointInd[3],midPointInd[5]};

			newTrianglesFlag[24]= this->volumeID_; 
			newTrianglesFlag[25]= originFlag[2]; 
			newTrianglesFlag[26]= this->volumeID_; 
			newTrianglesFlag[27]= this->volumeID_;

			isInterfaceSurface[24] = false;
			isInterfaceSurface[25] = interfaceSurface[2];
			isInterfaceSurface[26] = false;
			isInterfaceSurface[27] = false;

			newTriangleEdgeIDs[6]={36,37,39,36,38,40,37,38,41,39,40,41};

			// Element 8: (x_02,x_12,x_13,x_23)
			(newElements)[7]={midPointInd[2],midPointInd[3],midPointInd[4],midPointInd[5]};

			(newEdges)[42] = {midPointInd[2] ,midPointInd[3]}; 
			(newEdges)[43] = {midPointInd[2] ,midPointInd[4]}; 
			(newEdges)[44] = {midPointInd[2] ,midPointInd[5]}; 
			(newEdges)[45] = {midPointInd[3] ,midPointInd[4]}; 
			(newEdges)[46] = {midPointInd[3] ,midPointInd[5]}; 
			(newEdges)[47] = {midPointInd[4] ,midPointInd[5]}; 

			edgeFlags[42]=this->volumeID_;
			edgeFlags[43]=originFlag[1];
			edgeFlags[44]=originFlag[2];
			edgeFlags[45]=originFlag[3];
			edgeFlags[46]=originFlag[3];
			edgeFlags[47]=originFlag[3];

			isInterfaceEdge[42] = false;
			isInterfaceEdge[43] = interfaceSurface[1]; 
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
			newTriangles[28]= {midPointInd[2],midPointInd[3],midPointInd[4]};
			newTriangles[29]= {midPointInd[2],midPointInd[3],midPointInd[5]};
			newTriangles[30]= {midPointInd[2],midPointInd[4],midPointInd[5]};
			newTriangles[31]= {midPointInd[3],midPointInd[4],midPointInd[5]};

			newTrianglesFlag[28]= this->volumeID_; 
			newTrianglesFlag[29]= this->volumeID_; 
			newTrianglesFlag[30]= this->volumeID_; 
			newTrianglesFlag[31]= originFlag[3];

			isInterfaceSurface[28] = false;
			isInterfaceSurface[29] = false;
			isInterfaceSurface[30] = false;
			isInterfaceSurface[31] = interfaceSurface[3];

			newTriangleEdgeIDs[7]={42,43,45,42,44,46,43,44,47,45,46,47};
		}



		// Now we add the elements, edges and triangles 

		// Adding Elements
		int offsetElements = this->elementsC_->numberElements(); 
		int offsetEdges = this->edgeElements_->numberElements(); 
		for( int i=0;i<8; i++){
			//sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("regular");	
			feNew.setPredecessorElement(indexElement);
			if(elements->getElement(indexElement).getFiniteElementRefinementType() == "irregular" || elements->getElement(indexElement).getFiniteElementRefinementType() == "irregularRegular" )
				feNew.setFiniteElementRefinementType("irregularRegular");	
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
		int offsetSurfaces =this->surfaceTriangleElements_->numberElements();
		int offsetSurface =0;
		for( int i=0;i<32; i++){
			//sort( newTriangles.at(i).begin(), newTriangles.at(i).end() );
			FiniteElement feNew(newTriangles[i],newTrianglesFlag[i]);
			feNew.setInterfaceElement(isInterfaceSurface[i]);
			if(i<28){
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
					 if ( !this->elementsC_->getElement(i/4+offsetElements).subElementsInitialized() )
							this->elementsC_->getElement(i/4+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
						this->elementsC_->getElement(i/4+offsetElements).addSubElement(feNew);	
				}
				this->surfaceTriangleElements_->addSurface(feNew, i/4+offsetElements);	

			}
			else{
				if(newTrianglesFlag[i]!=0 && newTrianglesFlag[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(indexElement).subElementsInitialized() )
						this->elementsC_->getElement(indexElement).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(indexElement).addSubElement(feNew);	
				}
				this->surfaceTriangleElements_->addSurface(feNew, indexElement);				
			}

			/*this->surfaceTriangleElements_->getElement(i+offsetSurfaces).initializeSubElements( this->FEType_, this->dim_ -2) ;
			for(int j=0; j< 3 ; j++){
				this->surfaceTriangleElements_->getElement(i+offsetSurfaces).addSubElement(this->edgeElements_->getElement(offsetEdges + newTriangleEdgeIDs[i/4][j+offsetSurface]));	
			}
			offsetSurface += 3;
			if((i % 4) == 0)
				offsetSurface =0;*/	
		}
		FiniteElement element;
		FiniteElement feEdge;
		for( int i=0;i<8; i++){	
			if(i<7)	{
				element = this->elementsC_->getElement(i+offsetElements);
				}
			else{
				element = this->elementsC_->getElement(indexElement);
				}
			bool init=false;
			for(int j=0; j<48 ; j++){
				FiniteElement feEdge = this->edgeElements_->getElement(j+offsetEdges);
				if(feEdge.getFlag() != this->volumeID_){
					if(init == true)
						element.addSubElement( feEdge );
					else if ( !element.subElementsInitialized() ){
				        element.initializeSubElements( "P1", 1 ); // only P1 for now                
				        element.addSubElement( feEdge );
				        init= true;
				    }
				    else {
				        ElementsPtr_Type surfaces = element.getSubElements();
				        // We set the edge to the corresponding element(s)
				        surfaces->setToCorrectElement( feEdge );
				    }		
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
       	TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The Refinement Algorithm for the Dimension at hand is not available");
  }

/*!

 \brief 2D and 3D that bisects the edges of tagged Elements. Chosen by error estimator or otherwise elements are refined regular by connecting edge midpoints.

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.
@param[in] surfaceTriangleElements Triangle elements.


*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::bisectEdges(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement, SurfaceElementsPtr_Type surfaceTriangleElements, string mode){
	
	if(refinementMode_ == "Bisection"){

		if(this->dim_ == 2){
			vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element

			int entry = this->determineLongestEdge(edgeElements,edgeNumbers,this->pointsRep_); // we determine the edge, we would choose for blue Refinement

			if(!edgeElements->getElement(entry).isTaggedForRefinement()) // we tag every edge, after we refine an element -> no tag - no refinement on that edge so far
			{   			
				this->addMidpoint(edgeElements,entry);
				edgeElements->getElement(entry).tagForRefinement();
			}

		}
		else if(this->dim_==3)
			TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Explicit refinement by bisections is only implemented in 2d.");
	}	
	else{
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		
		for(int i=0; i<edgeNumbers.size(); i++){
			if(!edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()) // we tag every edge, after we refine an element -> no tag - no refinement on that edge so far
			{   			
				this->addMidpoint(edgeElements,edgeNumbers[i]);
				edgeElements->getElement(edgeNumbers[i]).tagForRefinement();
			}
		}
	}	
	if(mode == "all"){
		vec_int_Type edgeNumbers = edgeElements->getEdgesOfElement(indexElement); // indeces of edges belonging to element
		
		for(int i=0; i<edgeNumbers.size(); i++){
			if(!edgeElements->getElement(edgeNumbers[i]).isTaggedForRefinement()) // we tag every edge, after we refine an element -> no tag - no refinement on that edge so far
			{   			
				this->addMidpoint(edgeElements,edgeNumbers[i]);
				edgeElements->getElement(edgeNumbers[i]).tagForRefinement();
			}
		}
	}

}
/*!

 \brief 2D refinement by bisection of tagged Elements with three tagged Edges.

@param[in] edgeElements Edges.
@param[in] elements Elements.
@param[in] indexELement Element in question.
@param[in] surfaceTriangleElements Triangle elements.

*/
template <class SC, class LO, class GO, class NO>
void RefinementFactory<SC,LO,GO,NO>::bisectElement3(EdgeElementsPtr_Type edgeElements, ElementsPtr_Type elements, int indexElement){

	if(this->dim_ == 2){
		
        // The necessary Point was already added to the nodelist
		// now we have to figure out, which node it is -> check the tagged edge for midpoint
		vec_int_Type edgeNumbersTmp = edgeElements->getEdgesOfElement(indexElement);


		int entry = this->determineLongestEdge(edgeElements,edgeNumbersTmp,this->pointsRep_); // we determine the edge, we would choose for blue Refinement

		int midPointEntry=0;
		// midpoint index
		vec_int_Type midPointInd(3);
		vec_int_Type edgeNumbers(3);	
		vec_int_Type mutualNode(3);

				
		// -> longest Edge: Edge 1: midPoint[0] - mutualNode[0] = mutualNode of Edge 1 and 2 
		// -> 				Edge 2: midPoint[1] - mutualNode[1] = mutualNode of Edge 1 and 3
		// -> 				Edge 3: midpoint[2] - mutualNode[2] = mutualNode of Edge 2 and 3

		for(int i=0; i<3; i++){	
			if(edgeNumbersTmp[i] == entry)	{
				midPointInd[0] = edgeElements->getMidpoint(edgeNumbersTmp[i]);
				edgeNumbers[0] = edgeNumbersTmp[i];
				mutualNode[0] = edgeElements->getElement(entry).getNode(0);
				mutualNode[2] = edgeElements->getElement(entry).getNode(1);
			}
		}
		int tmp=1;
		for(int i=0; i<3; i++){	
			if(edgeNumbersTmp[i] != entry)	{
				midPointInd[tmp] = edgeElements->getMidpoint(edgeNumbersTmp[i]);
				edgeNumbers[tmp] = edgeNumbersTmp[i];
				tmp++;
				if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(0) == mutualNode[0])
					mutualNode[1] = edgeElements->getElement(edgeNumbersTmp[i]).getNode(1);
				else if(edgeElements->getElement(edgeNumbersTmp[i]).getNode(1) == mutualNode[0])
					mutualNode[1] = edgeElements->getElement(edgeNumbersTmp[i]).getNode(0);
			}
		}


		// Adding Elements and the corresponding Edges
		vec2D_int_Type newElements(4, vec_int_Type( 0 )); // vector for the new elements
		vec2D_int_Type newEdges(12,vec_int_Type(0)); // vector for the new edges
		vec_int_Type edgeFlags(12); // vector for the new flags
		vec_bool_Type isInterfaceEdge(12); // bool vectot for interfaceEdges
		vec_GO_Type predecessorElement(12);
		// Edgeflags are determined by the midpoints flag or by the fact, that they are inside a triangle, which consequently makes them interior edges (flag =this->volumeID_)

		// Element 1
		(newElements)[0]={mutualNode[0],midPointInd[0],midPointInd[1]};

		(newEdges)[0] = {mutualNode[0] ,midPointInd[0]}; 
		(newEdges)[1] = {mutualNode[0] ,midPointInd[1]}; 
		(newEdges)[2] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[0]=this->bcFlagRep_->at(midPointInd[0]);
		edgeFlags[1]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[2]=this->volumeID_;

		isInterfaceEdge[0] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[1] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[2] = false;

		predecessorElement[0] = this->edgeMap_->getGlobalElement(edgeNumbers[0]);
		predecessorElement[1] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[2] = -1;

		// Element 2
		newElements[1]={mutualNode[1],midPointInd[0],midPointInd[1]};

		(newEdges)[3] = {mutualNode[1] ,midPointInd[0]}; 
		(newEdges)[4] = {mutualNode[1] ,midPointInd[1]}; 
		(newEdges)[5] = {midPointInd[0] ,midPointInd[1]}; 

		edgeFlags[3]=this->volumeID_;
		edgeFlags[4]=this->bcFlagRep_->at(midPointInd[1]);
		edgeFlags[5]=this->volumeID_;

		isInterfaceEdge[3] = false;
		isInterfaceEdge[4] = edgeElements->getElement(edgeNumbers[1]).isInterfaceElement();
		isInterfaceEdge[5] = false;

		predecessorElement[3] = -1;
		predecessorElement[4] = this->edgeMap_->getGlobalElement(edgeNumbers[1]);
		predecessorElement[5] = -1;

		// Element 3
		(newElements)[2]={mutualNode[1] , midPointInd[0] ,midPointInd[2]};

		(newEdges)[6] = {mutualNode[1] ,midPointInd[0]}; 
		(newEdges)[7] = {mutualNode[1] ,midPointInd[2]}; 
		(newEdges)[8] = {midPointInd[0] ,midPointInd[2]}; 

		edgeFlags[6]=this->volumeID_;
		edgeFlags[7]=this->bcFlagRep_->at(midPointInd[2]);
		edgeFlags[8]=this->volumeID_;

		isInterfaceEdge[6] = false;
		isInterfaceEdge[7] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();
		isInterfaceEdge[8] = false;

		predecessorElement[6] = -1;
		predecessorElement[7] = this->edgeMap_->getGlobalElement(edgeNumbers[2]);
		predecessorElement[8] = -1;

		// Element 4
		(newElements)[3]={midPointInd[0],midPointInd[2],mutualNode[2]};	

		(newEdges)[9] = {midPointInd[0] ,midPointInd[2]}; 
		(newEdges)[10] = {midPointInd[0] ,mutualNode[2]}; 
		(newEdges)[11] = {midPointInd[2] ,mutualNode[2]}; 

		edgeFlags[9]=this->volumeID_;
		edgeFlags[10]=this->bcFlagRep_->at(midPointInd[0]);;
		edgeFlags[11]=this->bcFlagRep_->at(midPointInd[2]);;

		isInterfaceEdge[9] = false;
		isInterfaceEdge[10] = edgeElements->getElement(edgeNumbers[0]).isInterfaceElement();
		isInterfaceEdge[11] = edgeElements->getElement(edgeNumbers[2]).isInterfaceElement();

		predecessorElement[9] = -1;
		predecessorElement[10] =  this->edgeMap_->getGlobalElement(edgeNumbers[0]);;
		predecessorElement[11] =  this->edgeMap_->getGlobalElement(edgeNumbers[2]);;

		int offsetElements = this->elementsC_->numberElements(); 
		int offsetEdges = this->edgeElements_->numberElements(); 
		for( int i=0;i<4; i++){
			sort( newElements.at(i).begin(), newElements.at(i).end() );
			FiniteElement feNew(newElements.at(i),this->volumeID_);
			feNew.setFiniteElementRefinementType("regular");
			feNew.setPredecessorElement(indexElement);
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
				if(edgeFlags[i]!=0 && edgeFlags[i]!=this->volumeID_){
					if ( !this->elementsC_->getElement(i/3+offsetElements).subElementsInitialized() )
						this->elementsC_->getElement(i/3+offsetElements).initializeSubElements( this->FEType_, this->dim_ -1) ;
					this->elementsC_->getElement(i/3+offsetElements).addSubElement(feNew);	
					}
				}
			else
				this->edgeElements_->addEdge(feNew,indexElement);
			
		}
    }
	else        	
		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "The Refinement by Bisection Method you requested is only applicable to a 2 dimensional Mesh.");
	
 }

}

#endif