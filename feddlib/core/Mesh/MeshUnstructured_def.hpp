#ifndef MESHUNSTRUCTURED_def_hpp
#define MESHUNSTRUCTURED_def_hpp
#include "MeshUnstructured_decl.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector_def.hpp"

/*!
 Definition of MeshUnstructured
 
 @brief  MeshUnstructured
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;

namespace FEDD {

template <class SC, class LO, class GO, class NO>
MeshUnstructured<SC,LO,GO,NO>::MeshUnstructured():
Mesh<SC,LO,GO,NO>(),
meshInterface_(),
volumeID_(10),
edgeElements_(),
surfaceEdgeElements_(),
meshFileName_("fileName.mesh"),
delimiter_(" ")
//#ifdef FULL_TIMER
//,TotalP1ToP2Time_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Total P1 to P2")),
//BuildRedundantInfoTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Build Redundant Info")),
//SortUniqueRedundantInfoTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Sort Unique redundant")),
//CheckingFlaggedAndDeleteTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior")),
//CheckingFlaggedTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged")),
//DeleteInteriorTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Interior")),
//SetGlobalInterfaceIDTime_ (Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set Global IDs")),
//SetP2InfoTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set P2 Information")),
//GatherAllMarkedTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Gather All")),
//UniqueAndFindMarkedTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Unique and MaxAll")),
//UniqueNewFaceTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Unique new Faces")),
//SumAllMarkedTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Flagged & Delete Interior: Flagged: Sum All")),
//SetP2PointsTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set P2 Information: Set Points")),
//SetP2ElementsTime_(Teuchos::TimeMonitor::getNewCounter("FE: Unstructured Mesh: Set P2 Information: Set Elements"))
//#endif
{
    edgeElements_ = Teuchos::rcp( new EdgeElements_Type() );
    surfaceEdgeElements_ = Teuchos::rcp( new Elements_Type() );
    
}

template <class SC, class LO, class GO, class NO>
MeshUnstructured<SC,LO,GO,NO>::MeshUnstructured(CommConstPtr_Type comm, int volumeID):
Mesh<SC,LO,GO,NO>(comm),
meshInterface_(),
volumeID_(volumeID),
edgeElements_(),
surfaceEdgeElements_(),
meshFileName_("fileName.mesh"),
delimiter_(" ")
{
    edgeElements_ = Teuchos::rcp( new EdgeElements_Type() );
    surfaceEdgeElements_ = Teuchos::rcp( new Elements_Type() );
        
}

template <class SC, class LO, class GO, class NO>
MeshUnstructured<SC,LO,GO,NO>::~MeshUnstructured(){
}

//template <class SC, class LO, class GO, class NO>
//vec2D_int_ptr_Type MeshUnstructured<SC,LO,GO,NO>::getElements(){
//
//
//    this->elementsVec_ = Teuchos::rcp( new vec2D_int_Type( this->elementsC_->numberElements() ) );
//    for (int i=0; i<this->elementsVec_->size(); i++)
//        this->elementsVec_->at(i) = this->elementsC_->getElement(i).getVectorNodeList();
//
//    return this->elementsVec_;
//}


template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::getTriangles(int vertex1ID, int vertex2ID, vec_int_Type &vertices3ID){
    vertices3ID.resize(2);
    if (vertex2ID<vertex1ID) {
        int tmp = vertex2ID;
        vertex2ID = vertex1ID;
        vertex1ID = tmp;
    }
    if (vertex1ID == 0) {
        if (vertex2ID == 1) {
            vertices3ID.at(0) = 2;
            vertices3ID.at(1) = 3;
        }
        else if (vertex2ID == 2) {
            vertices3ID.at(0) = 1;
            vertices3ID.at(1) = 3;
        }
        else if (vertex2ID == 3) {
            vertices3ID.at(0) = 1;
            vertices3ID.at(1) = 2;
        }

    }
    else if(vertex1ID == 1){
        if (vertex2ID == 2) {
            vertices3ID.at(0) = 0;
            vertices3ID.at(1) = 3;
        }
        else if (vertex2ID == 3) {
            vertices3ID.at(0) = 0;
            vertices3ID.at(1) = 2;
        }
    }

    else if(vertex1ID == 2){
        if (vertex2ID == 3) {
            vertices3ID.at(0) = 0;
            vertices3ID.at(1) = 1;
        }
    }
    else {
#ifdef ASSERTS_WARNINGS
        MYASSERT(false,"Could not identify triangle.");
#endif
    }
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::buildP2ofP1MeshEdge( MeshUnstrPtr_Type meshP1 ){
    
    // If flags of line segments should be used over surface flags, this functions must be checked
    
    int rank = this->comm_->getRank();
    this->rankRange_ = meshP1->rankRange_;
    bool verbose( this->comm_->getRank() == 0 );
    this->elementMap_ = meshP1->elementMap_;
	this->edgeMap_ = meshP1->edgeMap_;
    this->dim_ = meshP1->getDimension();
    this->FEType_ = "P2";
    this->numElementsGlob_ = meshP1->numElementsGlob_;
	this->surfaceTriangleElements_ = meshP1->surfaceTriangleElements_;
    
    GO P1Offset = meshP1->mapUnique_->getMaxAllGlobalIndex()+1;
    EdgeElementsPtr_Type edgeElements = meshP1->getEdgeElements();
    ElementsPtr_Type elements = meshP1->getElementsC();
    
    if (verbose)
        cout << "-- --  Start building P2 mesh -- -- " << endl;
    
    if (verbose)
        cout << "-- Building edge mid points with edges and setting P2 elements ... " << flush;

    vec2D_dbl_Type newPoints( edgeElements->numberElements(), vec_dbl_Type( this->dim_ ) );
    vec_int_Type newFlags( edgeElements->numberElements(), -1 );
    int newNodesPerElement = -1;
    if (this->dim_==2)
        newNodesPerElement = 3;
    else if (this->dim_==3)
        newNodesPerElement = 6;
    
    vec2D_LO_Type newElementNodes( elements->numberElements(), vec_LO_Type( newNodesPerElement, -1 ) );

    vec2D_dbl_ptr_Type pointsP1 = meshP1->getPointsRepeated();
    MapConstPtr_Type mapRepeatedP1 = meshP1->getMapRepeated();
    vec2D_LO_Type markedPoints(0);
    // loop over all previously created edges
	vec_int_Type markedTrue(edgeElements->numberElements());

	int numberLocalPoints = pointsP1->size();

    for (int i=0; i<edgeElements->numberElements(); i++) {
        
        LO p1ID = edgeElements->getElement(i).getNode( 0 );
        LO p2ID = edgeElements->getElement(i).getNode( 1 );

        GO id1 = mapRepeatedP1->getGlobalElement( p1ID );
        GO id2 = mapRepeatedP1->getGlobalElement( p2ID );
        
        for (int d=0; d<this->dim_; d++)
            newPoints[i][d] = ( (*pointsP1)[p1ID][d] + (*pointsP1)[p2ID][d] ) / 2.;
      
       	newFlags[i] = determineFlagP2( meshP1, p1ID, p2ID, i, markedPoints );
		                
        const vec_LO_Type elementsOfEdge = edgeElements->getElementsOfEdge( i );
        const vec_GO_Type elementsGlobalOfEdge = edgeElements->getElementsOfEdgeGlobal( i );

		if(newFlags[i] != -1){ // questionable point that were given a flag, but that is not certain yet
       		for (int j=0; j<elementsOfEdge.size(); j++) {
           		if ( elementsOfEdge[j] == -1 ) 
					markedTrue[i] =1;
			}
		}	

                
        vec_GO_Type relevantElementsOfEdge(0);
        for (int j=0; j<elementsOfEdge.size(); j++) {
            if ( elementsOfEdge[j] != OTLO::invalid() )
                relevantElementsOfEdge.push_back( elementsOfEdge[j] );
        }

        // We need to determine the correct location in the P2 element
        vec_int_Type positions( relevantElementsOfEdge.size() );
        if ( relevantElementsOfEdge.size() > 0 )
            determinePositionInElementP2( positions, relevantElementsOfEdge, p1ID, p2ID, meshP1 );
    
        int factor = 0;
        if (this->dim_==2)
            factor=3;
        else if(this->dim_==3)
            factor = 4;
        for (int j=0; j<relevantElementsOfEdge.size(); j++)
            newElementNodes[ relevantElementsOfEdge[j] ][ positions[j]-factor ] =  i ;
        
    }
	// It is possible that a Edge is conencted to two surfaces with different Flags, that are also on different Processors
	// This Leads to two different Flags for the same Edge
	// In order to counter that effect we check the interface edges of which we determined the flag via surfaces and check if they have the same flag and if not choose the lower one
	int maxRank = std::get<1>(this->rankRange_);
	if(maxRank >0){
		vec_GO_Type edgeSwitch(0);
		vec_LO_Type flags(0);
		for(int i=0; i<edgeElements->numberElements(); i++){
			if(markedTrue[i]==1){
				edgeSwitch.push_back(this->edgeMap_->getGlobalElement(i));
				flags.push_back(newFlags[i]);
			}
		}

		// communticating elements across interface
		Teuchos::ArrayView<GO> edgeSwitchArray = Teuchos::arrayViewFromVector( edgeSwitch);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgeSwitchArray, 0, this->comm_) );

		// Global IDs of Procs
		// Setting newPoints as to be communicated Values
		MultiVectorLOPtr_Type edgeFlags = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP< LO > edgeFlagsEntries  = edgeFlags->getDataNonConst(0);

		for(int i=0; i< edgeFlagsEntries.size() ; i++){
			edgeFlagsEntries[i] = flags[i] ;
		}
		
		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface;

		if(mapGlobalInterface->getGlobalNumElements()>0){
			mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );
		}

		MultiVectorLOPtr_Type isInterfaceElement_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_imp->putScalar( (LO) 0 ); 
		isInterfaceElement_imp->importFromVector( edgeFlags, false, "Insert");
		//isInterfaceElement_imp->print();

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( edgeFlags, false, "Insert");
		//isInterfaceElement_exp->print();

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_imp, false, "Insert");
		//isInterfaceElement2_imp->print();

		isInterfaceElement2_imp->exportFromVector(isInterfaceElement_exp, false, "Insert");
		//isInterfaceElement2_imp->print();

		edgeFlagsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		for(int i=0; i<edgeFlagsEntries.size(); i++){
			LO entry = this->edgeMap_->getLocalElement(edgeSwitch[i]);
			if(newFlags[entry] > edgeFlagsEntries[i]){
				newFlags[entry] = edgeFlagsEntries[i];

			}
		}
	}
	// In the next Step we need to determine the missing Flags of the 'MarkedMissing' Edges
	// We create a Map of the entries we have and one of the ones we need
	vec_GO_Type edgesNeeded(0); // For the Flag entries we need
	vec_GO_Type edgesActive(0);
	vec_int_Type flagsTmp(0);
	for(int i=0; i<edgeElements->numberElements(); i++){
		if(newFlags[i]==-1){
			edgesNeeded.push_back(this->edgeMap_->getGlobalElement(i));
		}
		else{
			edgesActive.push_back(this->edgeMap_->getGlobalElement(i));
			flagsTmp.push_back(newFlags[i]);
			
		}
	}

	// communticating elements across interface
	Teuchos::ArrayView<GO> edgesNeededArray = Teuchos::arrayViewFromVector( edgesNeeded);
	Teuchos::ArrayView<GO> edgesActiveArray = Teuchos::arrayViewFromVector( edgesActive);
	
	MapPtr_Type mapEdgesNeeded =
		Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesNeededArray, 0, this->comm_) );

	MapPtr_Type mapEdgesActive =
		Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesActiveArray, 0, this->comm_) );

	MultiVectorLOPtr_Type flagsImport = Teuchos::rcp( new MultiVectorLO_Type( mapEdgesNeeded, 1 ) );
	flagsImport->putScalar(10);

	MultiVectorLOPtr_Type flagsExport = Teuchos::rcp( new MultiVectorLO_Type( mapEdgesActive, 1 ) );
	Teuchos::ArrayRCP< LO > flagExportEntries  = flagsExport->getDataNonConst(0);
	for(int i=0; i< flagExportEntries.size(); i++){
		flagExportEntries[i] = flagsTmp[i];
	}
	
	//flagsExport->print();
	flagsImport->importFromVector(flagsExport, false, "Insert");
	//flagsImport->print();

    Teuchos::ArrayRCP< LO > flagImportEntries  = flagsImport->getDataNonConst(0);
	for(int i=0; i<flagImportEntries.size(); i++){
		LO entry = this->edgeMap_->getLocalElement(edgesNeeded[i]);
		if(newFlags[entry] ==-1){
			newFlags[entry] = flagImportEntries[i];
		}
	}

    
    // ##########################################################
    // Set P2 Elements
    this->elementsC_.reset(new Elements());
    
    LO numberLocalP1Nodes = meshP1->getPointsRepeated()->size();
    
    for (int i=0; i<elements->numberElements(); i++) {
        vec_int_Type feNodeList = elements->getElement( i ).getVectorNodeListNonConst(); // get a copy
        for (int j=0; j<newElementNodes[i].size(); j++)
            feNodeList.push_back( newElementNodes[i][j] + numberLocalP1Nodes );
        FiniteElement feP2( feNodeList );
        this->elementsC_->addElement(feP2);
    }
    
    if (verbose)
        cout << "done --" << endl;
    
    if (verbose)
        cout << "-- Setting global point IDs and building repeated P2 map and repeated points ... " << flush;
    
    this->pointsRep_.reset(new std::vector<std::vector<double> >(meshP1->pointsRep_->size(),vector<double>(this->dim_,-1.)));
    *this->pointsRep_ = *meshP1->pointsRep_;
    this->bcFlagRep_.reset(new vector<int>(meshP1->bcFlagRep_->size()));
    *this->bcFlagRep_ = *meshP1->bcFlagRep_;

    this->pointsRep_->insert( this->pointsRep_->end(), newPoints.begin(), newPoints.end() );
    this->bcFlagRep_->insert( this->bcFlagRep_->end(), newFlags.begin(), newFlags.end());
    
    Teuchos::ArrayView<const GO> nodeList = meshP1->getMapRepeated()->getNodeElementList();
    std::vector<GO> vecGlobalIDs = Teuchos::createVector( nodeList );
    
    for (int i=0; i<edgeElements->numberElements(); i++){
        vecGlobalIDs.push_back( edgeElements->getGlobalID( (LO) i ) + P1Offset );
    }
    Teuchos::RCP<std::vector<GO> > pointsRepGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDs ) );
    Teuchos::ArrayView<GO> pointsRepGlobMappingArray = Teuchos::arrayViewFromVector( *pointsRepGlobMapping );
    
    this->mapRepeated_.reset(new Map<LO,GO,NO>( meshP1->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), pointsRepGlobMappingArray, 0, this->comm_) );
    
    if (verbose)
        cout << "done --" << endl;
    
    if (verbose)
        cout << "-- Building unique P2 map, setting unique points and setting P2 elements ... " << flush;

    this->mapUnique_ = this->mapRepeated_->buildUniqueMap( this->rankRange_ );
    
    this->pointsUni_.reset(new std::vector<std::vector<double> >( this->mapUnique_->getNodeNumElements(), vector<double>(this->dim_,-1. ) ) );
    this->bcFlagUni_.reset( new std::vector<int> ( this->mapUnique_->getNodeNumElements(), 0 ) );
    for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
        GO gid = this->mapUnique_->getGlobalElement( i );

        LO id = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement( i ) );
        this->pointsUni_->at(i) = this->pointsRep_->at(id);
        this->bcFlagUni_->at(i) = this->bcFlagRep_->at(id);

    }
	this->edgeElements_ = edgeElements;

    
    if (verbose)
        cout << "done --" << endl;
    
    if (verbose)
        cout << "-- Building P2 surface elements ... " << flush;
    
    setP2SurfaceElements( meshP1 );
    
    if (verbose)
        cout << "done --" << endl;
    
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::setP2SurfaceElements( MeshUnstrPtr_Type meshP1 ){
    // loop over all elements. Inside we loop over all P1 surface elements and set the P2 surface elements accordingly
    ElementsPtr_Type elementsP1 = meshP1->getElementsC();
    ElementsPtr_Type elementsP2 = this->getElementsC();
    vec2D_int_Type localSurfaceIndices;
    int surfaceElementOrder = 2;
    if (this->dim_==3)
        surfaceElementOrder = 3;

    
    vec2D_int_Type surfacePermutation;
    getLocalSurfaceIndices( surfacePermutation, surfaceElementOrder );
    
    vec2D_int_Type surfacePermutation2D; // only used if in 3D, since we need to additionally set edges
    if (this->dim_==3)
        getLocalSurfaceIndices( surfacePermutation2D, 2 );
    
    for (int i=0; i<elementsP1->numberElements(); i++) {
        
        FiniteElement fe = elementsP1->getElement( i );
        FiniteElement feP2 = elementsP2->getElement( i );
        
        ElementsPtr_Type subEl = fe.getSubElements(); //might be null
        for (int j=0; j<fe.numSubElements(); j++) {
            FiniteElement feSurf = subEl->getElement(j);
            this->setSurfaceP2(feP2, feSurf, surfacePermutation, this->dim_);
            
            // set edges for 3D case and if there are any edges, for 2D the edges are handled above
            ElementsPtr_Type subElSurf = feSurf.getSubElements(); //might be null
            if (!subElSurf.is_null()) {
                ElementsPtr_Type subElP2 = feP2.getSubElements();
                FiniteElement feSubP2 = subElP2->getElement(j);
                for (int e=0; e<feSurf.numSubElements(); e++) {
                    FiniteElement feEdge = subElSurf->getElement(e);
                    this->setSurfaceP2( feSubP2, feEdge, surfacePermutation2D, this->dim_-1 );
                    subElP2->switchElement( j, feSubP2 );
                }
            }
            elementsP2->switchElement( i, feP2 );
        }
        
    }
}


template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::setSurfaceP2( FiniteElement &feP2, const FiniteElement &surfFeP1, const vec2D_int_Type &surfacePermutation, int dim ){
    
    if (dim == 2) {
        for (int j=0; j<surfacePermutation.size(); j++) {
            vec_int_Type tmpSurface(2);
            tmpSurface[0] = feP2.getNode( surfacePermutation.at(j).at(0) );
            tmpSurface[1] = feP2.getNode( surfacePermutation.at(j).at(1) );
            
            sort( tmpSurface.begin(), tmpSurface.end() );
            
            vec_int_Type surfaceNodeP2(0);
            const vec_int_Type surfaceNodeP1 = surfFeP1.getVectorNodeList();
            if ( tmpSurface[0] == surfaceNodeP1[0]  && tmpSurface[1] == surfaceNodeP1[1] ) {
                                
                surfaceNodeP2.push_back( surfaceNodeP1[0] );
                surfaceNodeP2.push_back( surfaceNodeP1[1] );
                if (j == 0)
                    surfaceNodeP2.push_back( feP2.getNode( 3 ) );
                else if( j == 1 )
                    surfaceNodeP2.push_back( feP2.getNode( 5 ) );
                else if( j == 2 )
                    surfaceNodeP2.push_back( feP2.getNode( 4 ) );
                
                int flag = surfFeP1.getFlag();
                FiniteElement feP2Surf( surfaceNodeP2, flag );
                
                if ( !feP2.subElementsInitialized() )
                    feP2.initializeSubElements("P2",dim-1);
                feP2.addSubElement(feP2Surf);
            }
            
        }
    }
    else if (dim == 3){
        for (int j=0; j<surfacePermutation.size(); j++) {
            vec_int_Type tmpSurface(3);
            tmpSurface[0] = feP2.getNode( surfacePermutation.at(j).at(0) );
            tmpSurface[1] = feP2.getNode( surfacePermutation.at(j).at(1) );
            tmpSurface[2] = feP2.getNode( surfacePermutation.at(j).at(2) );
                        
            //sort( tmpSurface.begin(), tmpSurface.end() );
            vec_int_Type index(3, 0);
            for (int i = 0 ; i != index.size() ; i++)
                index[i] = i;
            
            sort(index.begin(), index.end(),
                 [&](const int& a, const int& b) {
                     return  tmpSurface[a] < tmpSurface[b];
                    }
                 );
            
            tmpSurface = sort_from_ref(tmpSurface, index);
            
            vec_int_Type surfaceNodeP2(0);
            const vec_int_Type surfaceNodeP1 = surfFeP1.getVectorNodeList();
            
            if ( tmpSurface[0] == surfaceNodeP1[0]  &&
                    tmpSurface[1] == surfaceNodeP1[1] &&
                    tmpSurface[2] == surfaceNodeP1[2]) {
                surfaceNodeP2.push_back( surfaceNodeP1[0] );
                surfaceNodeP2.push_back( surfaceNodeP1[1] );
                surfaceNodeP2.push_back( surfaceNodeP1[2] );
                vec_int_Type additionalP2IDs(3);
                if (j == 0){
                    additionalP2IDs[0] = feP2.getNode( 4 );
                    additionalP2IDs[1] = feP2.getNode( 5 );
                    additionalP2IDs[2] = feP2.getNode( 6 );
                }
                else if( j == 1 ){
                    additionalP2IDs[0] = feP2.getNode( 4 );
                    additionalP2IDs[1] = feP2.getNode( 7 );
                    additionalP2IDs[2] = feP2.getNode( 8 );
                }
                else if( j == 2 ){
                    additionalP2IDs[0] = feP2.getNode( 5 );
                    additionalP2IDs[1] = feP2.getNode( 8 );
                    additionalP2IDs[2] = feP2.getNode( 9 );
                }
                else if( j == 3 ){
                    additionalP2IDs[0] = feP2.getNode( 6 );
                    additionalP2IDs[1] = feP2.getNode( 7 );
                    additionalP2IDs[2] = feP2.getNode( 9 );
                }
                
                additionalP2IDs = this->reorderP2SurfaceIndices(additionalP2IDs, index);
                surfaceNodeP2.push_back( additionalP2IDs[0] );
                surfaceNodeP2.push_back( additionalP2IDs[1] );
                surfaceNodeP2.push_back( additionalP2IDs[2] );
                
                int flag = surfFeP1.getFlag();
                FiniteElement feP2Surf( surfaceNodeP2, flag );
                if ( !feP2.subElementsInitialized() )
                    feP2.initializeSubElements("P2",dim-1);
                feP2.addSubElement(feP2Surf);
            }
            
        }
    }
    
}


template <class SC, class LO, class GO, class NO>
vec_int_Type MeshUnstructured<SC,LO,GO,NO>::reorderP2SurfaceIndices( vec_int_Type& additionalP2IDs, vec_int_Type& index , bool track){
    vec_int_Type reorderedIDs(3);
    // depending on the sorting of P1 surface nodes we have to adjust the new ordering of P2 edge midpoints for surfaces in 3D
    if (index[0] == 0){
        if(index[1] == 1){
            reorderedIDs[0] = 0; reorderedIDs[1] = 1; reorderedIDs[2] = 2;
        }
        else if(index[1] == 2){
            reorderedIDs[0] = 2; reorderedIDs[1] = 1; reorderedIDs[2] = 0;
        }
    }
    else if (index[0] == 1){
        if(index[1] == 0){
            reorderedIDs[0] = 0; reorderedIDs[1] = 2; reorderedIDs[2] = 1;
        }
        else if(index[1] == 2){
            reorderedIDs[0] = 1; reorderedIDs[1] = 2; reorderedIDs[2] = 0;
        }
    }
    else if (index[0] == 2){
        if(index[1] == 0){
            reorderedIDs[0] = 2; reorderedIDs[1] = 0; reorderedIDs[2] = 1;
        }
        else if(index[1] == 1){
            reorderedIDs[0] = 1; reorderedIDs[1] = 0; reorderedIDs[2] = 2;
        }
    }
    
    return sort_from_ref(additionalP2IDs, reorderedIDs);
}

//template <class SC, class LO, class GO, class NO>
//void MeshUnstructured<SC,LO,GO,NO>::setSurfaceP2( const FiniteElement &elementP2 , const vec_int_Type& surfaceNode, vec_int_Type& surfaceNodeP2, const vec2D_int_Type &surfacePermutation ){
//    
//    int loc, id1, id2, id3;
//    if (this->dim_ == 2) {
//        for (int j=0; j<surfacePermutation.size(); j++) {
//            id1 = elementP2.getNode( surfacePermutation.at(j).at(0) );
//            id2 = elementP2.getNode( surfacePermutation.at(j).at(1) );
//            
//            vec_int_Type tmpSurface(2);
//            if (id1 > id2){
//                tmpSurface[0] = id2;
//                tmpSurface[1] = id1;
//            }
//            else{
//                tmpSurface[0] = id1;
//                tmpSurface[1] = id2;
//            }
//            
//            if ( tmpSurface[0] == surfaceNode[0]  && tmpSurface[1] == surfaceNode[1] ) {
//                surfaceNodeP2.push_back( surfaceNode[0] );
//                surfaceNodeP2.push_back( surfaceNode[1] );
//                if (j == 0)
//                    surfaceNodeP2.push_back( elementP2.getNode( 3 ) );
//                else if( j == 1 )
//                    surfaceNodeP2.push_back( elementP2.getNode( 5 ) );
//                else if( j == 2 )
//                    surfaceNodeP2.push_back( elementP2.getNode( 4 ) );
//            }
//        }
//    }
//    else if (this->dim_ == 3){
//        for (int j=0; j<surfacePermutation.size(); j++) {
//            id1 = elementP2.getNode( surfacePermutation.at(j).at(0) );
//            id2 = elementP2.getNode( surfacePermutation.at(j).at(1) );
//            id3 = elementP2.getNode( surfacePermutation.at(j).at(2) );
//            
//            vec_int_Type tmpSurface = {id1 , id2 , id3};
//            sort( tmpSurface.begin(), tmpSurface.end() );
//            
//            if ( tmpSurface[0] == surfaceNode[0]  &&
//                    tmpSurface[1] == surfaceNode[1] &&
//                    tmpSurface[2] == surfaceNode[2]) {
//                surfaceNodeP2.push_back( surfaceNode[0] );
//                surfaceNodeP2.push_back( surfaceNode[1] );
//                surfaceNodeP2.push_back( surfaceNode[2] );
//                if (j == 0){
//                    surfaceNodeP2.push_back( elementP2.getNode( 4 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 5 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 6 ) );
//                }
//                else if( j == 1 ){
//                    surfaceNodeP2.push_back( elementP2.getNode( 4 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 7 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 8 ) );
//                }
//                else if( j == 2 ){
//                    surfaceNodeP2.push_back( elementP2.getNode( 5 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 8 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 9 ) );
//                }
//                else if( j == 3 ){
//                    surfaceNodeP2.push_back( elementP2.getNode( 6 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 7 ) );
//                    surfaceNodeP2.push_back( elementP2.getNode( 9 ) );
//                }
//            }
//        }
//    }
//    
//}
//    
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::getLocalSurfaceIndices(vec2D_int_Type &localSurfaceIndices , int surfaceElementOrder ){
    
    if ( this->dim_ == 2 ) {
        
        if (surfaceElementOrder == 2) { //P1
            localSurfaceIndices.resize(3,vec_int_Type(3,-1));
            localSurfaceIndices.at(0).at(0) = 0;
            localSurfaceIndices.at(0).at(1) = 1;
            localSurfaceIndices.at(1).at(0) = 0;
            localSurfaceIndices.at(1).at(1) = 2;
            localSurfaceIndices.at(2).at(0) = 1;
            localSurfaceIndices.at(2).at(1) = 2;
        }
        else{
#ifdef ASSERTS_WARNINGS
            MYASSERT(false,"no permutation for this surface yet.");
#endif
        }
    }
    else if ( this->dim_ == 3 ){
        if (surfaceElementOrder == 3) {
            localSurfaceIndices.resize(4,vec_int_Type(3,-1));
            localSurfaceIndices.at(0).at(0) = 0;
            localSurfaceIndices.at(0).at(1) = 1;
            localSurfaceIndices.at(0).at(2) = 2;
            localSurfaceIndices.at(1).at(0) = 0;
            localSurfaceIndices.at(1).at(1) = 1;
            localSurfaceIndices.at(1).at(2) = 3;
            localSurfaceIndices.at(2).at(0) = 1;
            localSurfaceIndices.at(2).at(1) = 2;
            localSurfaceIndices.at(2).at(2) = 3;
            localSurfaceIndices.at(3).at(0) = 0;
            localSurfaceIndices.at(3).at(1) = 2;
            localSurfaceIndices.at(3).at(2) = 3;
        }
        else{
#ifdef ASSERTS_WARNINGS
            MYASSERT(false,"no permutation for this surface yet.");
#endif
        }
    }
}
    
    
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::getEdgeCombinations( vec2D_int_Type& edgeCombinations ){

    if (this->dim_ == 2) {
        edgeCombinations[0][0] = 0; edgeCombinations[0][1] = 1;
        edgeCombinations[1][0] = 0; edgeCombinations[1][1] = 2;
        edgeCombinations[2][0] = 1; edgeCombinations[2][1] = 2;
    }
    else if (this->dim_ == 3) {
        edgeCombinations[0][0] = 0; edgeCombinations[0][1] = 1;
        edgeCombinations[1][0] = 0; edgeCombinations[1][1] = 2;
        edgeCombinations[2][0] = 0; edgeCombinations[2][1] = 3;
        edgeCombinations[3][0] = 1; edgeCombinations[3][1] = 2;
        edgeCombinations[4][0] = 1; edgeCombinations[4][1] = 3;
        edgeCombinations[5][0] = 2; edgeCombinations[5][1] = 3;
    }
    
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::determinePositionInElementP2( vec_int_Type& positions, vec_GO_Type& elementsOfEdge, LO p1ID, LO p2ID, MeshUnstrPtr_Type meshP1 ){
    ElementsPtr_Type elements = meshP1->getElementsC();
    for (int i=0; i<elementsOfEdge.size(); i++) {

        const vec_int_Type nodeList = elements->getElement( elementsOfEdge[i] ).getVectorNodeList();

        auto it1 = find( nodeList.begin(), nodeList.end() , p1ID );
        int localElNum1 = distance( nodeList.begin() , it1 );
        auto it2 = find( nodeList.begin(), nodeList.end() , p2ID );
        int localElNum2 = distance( nodeList.begin() , it2 );
        if (localElNum1 > localElNum2) {
            int tmp = localElNum1;
            localElNum1 = localElNum2;
            localElNum2 = tmp;
        }
        if (this->dim_ == 2) {
            if (localElNum1==0) {
                if (localElNum2==1)
                    positions[i] = 3;
                else if(localElNum2==2)
                    positions[i] = 5;
            }
            else
                positions[i] = 4;
        }
        else if (this->dim_ == 3) {
            if (localElNum1==0) {
                if (localElNum2==1)
                    positions[i] = 4;
                else if(localElNum2==2)
                    positions[i] = 6;
                else if(localElNum2==3)
                    positions[i] = 7;
                
            }
            else if(localElNum1==1){
                if (localElNum2==2)
                    positions[i] = 5;
                else if(localElNum2==3)
                    positions[i] = 8;
            }
            else
                positions[i] = 9;
        }
    }
}

template <class SC, class LO, class GO, class NO>
int MeshUnstructured<SC,LO,GO,NO>::determineFlagP2( FiniteElement& fe, LO p1ID, LO p2ID, vec2D_int_Type& permutation ){

    int newFlag = std::numeric_limits<int>::max();
    vec_int_Type newFlags(0);
    bool foundLineSegment;
    const vec_int_Type nodeList = fe.getVectorNodeList();
    
    vec_int_Type localElementNumbering(2);
    auto it1 = find( nodeList.begin(), nodeList.end() , p1ID );
    localElementNumbering[0] = distance( nodeList.begin() , it1 );
    auto it2 = find( nodeList.begin(), nodeList.end() , p2ID );
    localElementNumbering[1] = distance( nodeList.begin() , it2 );

    fe.findEdgeFlagInSubElements( localElementNumbering, newFlags, false /*we are not in a subElement yet*/, permutation, foundLineSegment );
    
    if (newFlags.size() == 0)
        newFlag = volumeID_;
    
    // We use the lowest flag
    for (int k = 0; k < newFlags.size(); k++) {
        if ( newFlag > newFlags[k] )
            newFlag = newFlags[k];
    }
    
    return newFlag;
}
    
template <class SC, class LO, class GO, class NO>
int MeshUnstructured<SC,LO,GO,NO>::determineFlagP2( MeshUnstrPtr_Type meshP1, LO p1ID, LO p2ID, LO localEdgeID, vec2D_LO_Type& markedPoint ){
    

    ElementsPtr_Type elements = meshP1->getElementsC();

    vec2D_int_Type permutation = elements->getElementEdgePermutation();

    EdgeElementsPtr_Type edgeElements = meshP1->getEdgeElements();

    MapConstPtr_Type elementMap = meshP1->getElementMap();

    int rank = this->comm_->getRank();

    int flag1 = ( *meshP1->getBCFlagRepeated() )[p1ID];

    int flag2 = ( *meshP1->getBCFlagRepeated() )[p2ID];

    int newFlag = std::numeric_limits<int>::max();

    if(flag1 == volumeID_ || flag2 == volumeID_ ) // one node is in an inner node, than the new node is an inner node aswell
        newFlag = volumeID_;
    else{

        // check if node 1 and node 2 are part of the same surface. In this case we can use the flag of the corresponding surface. Otherwise we mark the new point as an interior point and it gets the volumeID_.

        const vec_LO_Type elementsOfEdge = edgeElements->getElementsOfEdge( (int) localEdgeID );

        const vec_GO_Type elementsOfEdgeGlobal = edgeElements->getElementsOfEdgeGlobal( (int) localEdgeID );


        bool markPoint = false;
        bool foundFlag = false;
        vec_int_Type localElementNumbering(2);
        vec_int_Type edge(2);
        bool foundLineSegment = false;
        vec_int_Type newFlags(0);
        for (int i=0; i<elementsOfEdge.size() && !foundLineSegment; i++) {
            if ( elementsOfEdge[i] != OTLO::invalid() ) {
                //In elementsOfEdge we can access all elements which have this (the current) edge.
                
                FiniteElement fe = elements->getElement( elementsOfEdge[i] );
                const vec_int_Type nodeList = fe.getVectorNodeList();

                // we need to determine the numbering of p1ID and p2ID corresponding to the current element
                auto it1 = find( nodeList.begin(), nodeList.end() , p1ID );
                localElementNumbering[0] = distance( nodeList.begin() , it1 );
                auto it2 = find( nodeList.begin(), nodeList.end() , p2ID );
                localElementNumbering[1] = distance( nodeList.begin() , it2 );
                edge[0] = p1ID; edge[1] = p2ID;
                sort( edge.begin(), edge.end() );
                
                fe.findEdgeFlagInSubElements( edge, newFlags, false /*we are not in a subElement yet*/, permutation, foundLineSegment );

                //We need to mark this point since it can still be on the surface and another element holds the corresponding surface with the correct flag.
                if (newFlags.size() == 0 && newFlag > volumeID_)
                    newFlag = volumeID_; //do we need this?

                //If we found a line element, the we choose this flag
                if (foundLineSegment){
                    foundFlag = true;
                    newFlag = newFlags [0];
                }
                else {
                    // We use the lowest flag of all surfaces
                    for (int k = 0; k < newFlags.size(); k++) {
                        foundFlag = true;
                        if (newFlag > newFlags[k] )
                            newFlag = newFlags[k];
                    }
                }
            }
            else
                markPoint = true;
        }
        if (markPoint && !foundFlag) {
            // We have not found a surface flag for this point.
            // No part of the element which owns this edge is part of the surface, but there might be a different element which has this edge but does not own it, but this element might have a part of the surface. We mark this point and determine the correct flag later
            markedPoint.push_back( {p1ID, p2ID, localEdgeID } );
            newFlag = -1;
        }
        
    }
    
    return newFlag;
}



template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::findSurfaces( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localSurfaceNodeList_vec, vec_int_Type& loc_vector, bool critical){
    int tmpID;
    loc_vector.resize(0);
    vec2D_int_Type::iterator it;
    int id1 = elementNodeList.at(numbering.at(0));
    int id2 = elementNodeList.at(numbering.at(1));

    if (this->dim_==2) {
        vec_int_Type searchfor(2,-1);
        // Here, we use that the local surfaces are sorted. See MyMeshConvert::ReadMesh(...)
        if ( id1 < id2 ) {
            searchfor.at(0) = id1;
            searchfor.at(1) = id2;
        }
        else{
            searchfor.at(0) = id2;
            searchfor.at(1) = id1;
        }

        it = find( localSurfaceNodeList_vec.begin(), localSurfaceNodeList_vec.end(), searchfor );
        if (it!=localSurfaceNodeList_vec.end()) {
            loc_vector.push_back( distance(localSurfaceNodeList_vec.begin(),it) );
        }

    }
    else if (this->dim_==3){
        
        vec_int_Type vertices3ID(0);
        //we want to identify the possible third surface component of a triangle. If we use other elements in 3D we might change this function.
        getTriangles(numbering.at(0), numbering.at(1), vertices3ID);
        // we have the two  possible combinations as a vector: vertices3ID
        
        vec_int_Type searchfor(3,-1);
        for (int i=0; i<2; i++) {
            searchfor.at(0) = id1;
            searchfor.at(1) = id2;
            searchfor.at(2) = elementNodeList.at( vertices3ID.at(i) );

            sort( searchfor.begin(), searchfor.end() );
            
            
            it = find( localSurfaceNodeList_vec.begin(), localSurfaceNodeList_vec.end(), searchfor );
            if (it!=localSurfaceNodeList_vec.end()){
                loc_vector.push_back( distance(localSurfaceNodeList_vec.begin(),it) );
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::findEdges( const vec_int_Type& elementNodeList, vec_int_Type numbering,  vec2D_int_Type& localEdgeNodeList_vec, vec_int_Type& loc_vector){
    int tmpID;
    loc_vector.resize(0);
    vec2D_int_Type::iterator it;
    int id1 = elementNodeList.at(numbering.at(0));
    int id2 = elementNodeList.at(numbering.at(1));
    vec_int_Type searchEdge(0);
    if (id2 > id1)
        searchEdge = { id1, id2 };
    else
        searchEdge = { id2, id1 };
    
    it = find( localEdgeNodeList_vec.begin(), localEdgeNodeList_vec.end(), searchEdge );
    if ( it != localEdgeNodeList_vec.end() )
        loc_vector.push_back( distance(localEdgeNodeList_vec.begin(),it) );
}

template <class SC, class LO, class GO, class NO>
typename MeshUnstructured<SC,LO,GO,NO>::MeshInterfacePtr_Type MeshUnstructured<SC,LO,GO,NO>::getMeshInterface(){
    TEUCHOS_TEST_FOR_EXCEPTION( meshInterface_.is_null(), std::runtime_error, "MeshInterface is null.");
    return meshInterface_;
}
    
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::buildMeshInterfaceParallelAndDistance( MeshUnstrPtr_Type mesh, vec_int_Type flag_vec, vec_dbl_ptr_Type &distancesToInterface ){
    this->meshInterface_.reset( new MeshInterface<SC,LO,GO,NO> ( this->getComm() ) );

    //    BuildMeshInterface for different flags
    this->meshInterface_->determineInterfaceParallelAndDistance( this->pointsUni_,  mesh->pointsUni_, this->bcFlagUni_, mesh->bcFlagUni_, flag_vec, this->getMapUnique(), mesh->getMapUnique(), distancesToInterface, this->pointsRep_, this->getDimension() );
    
    mesh->meshInterface_.reset( new MeshInterface_Type ( mesh->getComm() ) );
    mesh->meshInterface_->buildFromOtherInterface( meshInterface_ );
    
    //because we had to communicated all interface information in determineInterfaceParallel(), we can now partition the information again.
    this->partitionInterface();
    mesh->partitionInterface();
}
    
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::partitionInterface(){
    FEDD_TIMER_START(interfacePartitionTimer," : Mesh : Partition Interface");
    if (!this->meshInterface_.is_null())
        this->meshInterface_->partitionMeshInterface( this->mapRepeated_, this->mapUnique_);

}
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::setMeshFileName(string meshFileName, string delimiter){

    meshFileName_ = meshFileName;
    delimiter_ = delimiter;
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readMeshSize(){

    int numElement;
    int orderElement;
    int dim;
    int numNode;
    int numSurface = -1;
    int orderSurface = -1;
    int numEdges = 0;
    int orderEdges = 0;
    bool verbose ( this->comm_->getRank() == 0 );

    if (verbose) {
        cout << "\n";
        cout << "  Read data of mesh " << meshFileName_<< ":\n";
    }

    meshReadSize ( meshFileName_, numNode, dim, numElement, orderElement, numSurface, orderSurface, numEdges, orderEdges );
    
    if (verbose) {
        cout << "\n";
        cout << "\n";
        cout << "  Number of nodes = " << numNode << "\n";
        cout << "  Spatial dimension = " << dim << "\n";
        cout << "  Number of elements = " << numElement << "\n";
        cout << "  Element order = " << orderElement << "\n";
        cout << "  Number of surface elements = " << numSurface << "\n";
        cout << "  Surface element order = " << orderSurface << "\n";
        cout << "  Number of edge elements (for 3D) = " << numEdges << "\n";
        cout << "  Edges element order (for 3D) = " << orderEdges << "\n";
        
        cout << "\n";
        cout << "\n";
        cout << " Starting to read the data. \n";
    }

    
    this->elementOrder_ = orderElement;
    this->surfaceElementOrder_ = orderSurface;
    this->edgesElementOrder_ = orderEdges;
    this->numElements_ = numElement;
    this->numSurfaces_ = numSurface;
    this->numEdges_ = numEdges;
    this->numNodes_ = numNode;
    
    this->numElementsGlob_ = numElement;    
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readMeshEntity(string entityType){
    
    if (entityType == "element")
        this->readElements( );
    else if (entityType == "surface")
        this->readSurfaces( );
    else if (entityType == "line")
        this->readLines( );
    else if (entityType == "node")
        this->readNodes( );
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Unknown entity type.");
}

//template <class SC, class LO, class GO, class NO>
//void MeshUnstructured<SC,LO,GO,NO>::readSurfaces(){
//    bool verbose ( this->comm_->getRank() == 0 );
//    if (verbose)
//        cout << "### Starting to read surface data ... " << flush;
//
//    vec_int_Type    surfaceFlags( numSurfaces_, 0 );
//    vec_int_Type    surfacesCont( numSurfaces_* surfaceElementOrder_, 0 );
//
//    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
//    meshReadData ( meshFileName_, "surface", delimiter_, this->getDimension(), numSurfaces_, surfaceElementOrder_, surfacesCont, surfaceFlags );
//
//    if (verbose){
//        cout << "done." << endl;
//        cout << "### Setting surface data ... " << flush;
//    }
//    ElementsPtr_Type surfaceElementsMesh = this->getSurfaceElements();
//
//    for (int i=0; i<numSurfaces_; i++) {
//        vec_int_Type tmp(surfaceElementOrder_);
//        for (int j=0; j<surfaceElementOrder_; j++)
//            tmp.at(j) = surfacesCont.at( i * surfaceElementOrder_ + j ) - 1;// -1 to have start index 0
//
//        sort( tmp.begin(), tmp.end() ); // we sort here in order to identify the corresponding element faster!
//        FiniteElement feSurface( tmp , surfaceFlags[i] );
//        surfaceElementsMesh->addElement( feSurface );
//    }
//
//
//    if (verbose)
//        cout << "done." << endl;
//}

//template <class SC, class LO, class GO, class NO>
//void MeshUnstructured<SC,LO,GO,NO>::readLines(){
//    bool verbose ( this->comm_->getRank() == 0 );
//    if (verbose)
//        cout << "### Starting to read line data ... " << flush;
//
//    vec_int_Type    edgeFlags( numEdges_, 0 );
//    vec_int_Type    edgesCont( numEdges_* edgesElementOrder_, 0 );
//
//    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
//    meshReadData ( meshFileName_, "line", delimiter_, this->getDimension(), numEdges_, edgesElementOrder_, edgesCont, edgeFlags );
//
//
//    if (verbose){
//        cout << "done." << endl;
//        cout << "### Setting line data ... " << flush;
//    }
//
//    ElementsPtr_Type edgeElementsMesh = this->getSurfaceEdgeElements();
//
//    //Continous edge surface elements to FiniteElement object (only relevant in 3D)
//    for (int i=0; i<numEdges_; i++) {
//        vec_int_Type tmp(edgesElementOrder_);
//        for (int j=0; j<edgesElementOrder_; j++)
//            tmp.at(j) = edgesCont.at( i * edgesElementOrder_ + j ) - 1;// -1 to have start index 0
//
//        sort( tmp.begin(), tmp.end() ); // we sort here in order to identify the corresponding element faster!
//        FiniteElement feEdge( tmp , edgeFlags[i] );
//        edgeElementsMesh->addElement( feEdge );
//    }
//
//    if (verbose)
//        cout << "done." << endl;
//}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readNodes(){
    bool verbose ( this->comm_->getRank() == 0 );
    if (verbose)
        cout << "### Starting to read node data ... " << flush;

    vec_dbl_Type nodes(numNodes_ * this->getDimension(), 0.);
    vec_int_Type nodeFlags(numNodes_,0);
    meshReadData ( meshFileName_, "node", delimiter_, this->getDimension(), numNodes_, 3/*order of nodes is always 3*/, nodes, nodeFlags );
    
    if (verbose){
        cout << "done." << endl;
        cout << "### Setting node data ... " << flush;
    }
    //Here, all points are saved on every proc
    this->pointsRep_.reset(new std::vector<std::vector<double> >(numNodes_,std::vector<double>(this->getDimension(),-1.)));
    this->bcFlagRep_.reset(new std::vector<int> (numNodes_,0));

    FEDD_TIMER_START(pointsTimer," : MeshReader : Set Points not partitioned");

    for (int i=0; i<numNodes_ ; i++) {
        for (int j=0; j<this->getDimension(); j++)
            this->pointsRep_->at(i).at(j) = nodes[this->getDimension()*i+j];
        
        this->bcFlagRep_->at(i) = nodeFlags[i];
    }

    if (verbose)
        cout << "done." << endl;
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readElements(){
    bool verbose ( this->comm_->getRank() == 0 );
    if (verbose)
        cout << "### Starting to read element data ... " << flush;

    vec_int_Type    elementFlags( numElements_, 0 );
    vec_int_Type    elementsCont( numElements_* elementOrder_, 0 );
        
    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
    meshReadData ( meshFileName_, "element", delimiter_, this->getDimension(), numElements_, elementOrder_, elementsCont, elementFlags );

    if (verbose){
        cout << "done." << endl;
        cout << "### Setting element data ... " << flush;
    }
    ElementsPtr_Type elementsMesh = this->getElementsC();
    elementsMesh->setFiniteElementType("P1");
    elementsMesh->setDimension(this->getDimension());
    
	
    int id;
    for (int i=0; i<numElements_; i++) {
        vec_int_Type tmp(elementOrder_);
        for (int j=0; j<elementOrder_; j++){
            id = elementsCont.at(i*elementOrder_ + j) - 1;
            tmp.at(j) = id;
        }
        FiniteElement fe( tmp , elementFlags[i] );
        elementsMesh->addElement( fe );
    }
    
    if (verbose)
        cout << "done." << endl;
}


template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readSurfaces(){
    bool verbose ( this->comm_->getRank() == 0 );
    if (verbose)
        cout << "### Starting to read surface data ... " << flush;
    
    vec_int_Type    surfaceFlags( numSurfaces_, 0 );
    vec_int_Type    surfacesCont( numSurfaces_* surfaceElementOrder_, 0 );
        
    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
    meshReadData ( meshFileName_, "surface", delimiter_, this->getDimension(), numSurfaces_, surfaceElementOrder_, surfacesCont, surfaceFlags );
        
    if (verbose){
        cout << "done." << endl;
        cout << "### Setting surface data ... " << flush;
    }
    ElementsPtr_Type surfaceElementsMesh = this->getSurfaceElements();
    
    for (int i=0; i<numSurfaces_; i++) {
        vec_int_Type tmp(surfaceElementOrder_);
        for (int j=0; j<surfaceElementOrder_; j++)
            tmp.at(j) = surfacesCont.at( i * surfaceElementOrder_ + j ) - 1;// -1 to have start index 0
        
        sort( tmp.begin(), tmp.end() ); // we sort here in order to identify the corresponding element faster!
        FiniteElement feSurface( tmp , surfaceFlags[i] );
        surfaceElementsMesh->addElement( feSurface );
    }
    
    
    if (verbose)
        cout << "done." << endl;
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readLines(){
    bool verbose ( this->comm_->getRank() == 0 );

    if (verbose)
        cout << "### Starting to read line data ... " << flush;

    vec_int_Type    edgeFlags( numEdges_, 0 );
    vec_int_Type    edgesCont( numEdges_* edgesElementOrder_, 0 );
        
    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
    meshReadData ( meshFileName_, "line", delimiter_, this->getDimension(), numEdges_, edgesElementOrder_, edgesCont, edgeFlags );

    
    if (verbose){
        cout << "done." << endl;
        cout << "### Setting line data ... " << flush;
    }
    
    ElementsPtr_Type edgeElementsMesh = this->getSurfaceEdgeElements();
    
    //Continous edge surface elements to FiniteElement object (only relevant in 3D)
    for (int i=0; i<numEdges_; i++) {
        vec_int_Type tmp(edgesElementOrder_);
        for (int j=0; j<edgesElementOrder_; j++)
            tmp.at(j) = edgesCont.at( i * edgesElementOrder_ + j ) - 1;// -1 to have start index 0
        
        sort( tmp.begin(), tmp.end() ); // we sort here in order to identify the corresponding element faster!
        FiniteElement feEdge( tmp , edgeFlags[i] );
        edgeElementsMesh->addElement( feEdge );
    }
    
    if (verbose)
        cout << "done." << endl;
}

}
#endif
