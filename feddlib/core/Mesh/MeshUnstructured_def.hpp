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

using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::outArg;

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
	this->surfaceTriangleElements_ = meshP1->surfaceTriangleElements_; // for later
    
	meshP1->assignEdgeFlags(); // Function that determines the flag for each edge. That way the P2 flags can easily be determined
    GO P1Offset = meshP1->mapUnique_->getMaxAllGlobalIndex()+1;
    EdgeElementsPtr_Type edgeElements = meshP1->getEdgeElements();
    ElementsPtr_Type elements = meshP1->getElementsC();
    
    if (verbose)
        std::cout << "-- --  Start building P2 mesh -- -- " << std::endl;
    
    if (verbose)
        std::cout << "-- Building edge mid points with edges and setting P2 elements ... " << std::flush;

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


    // loop over all previously created edges to construct new P2 points 

    for (int i=0; i<edgeElements->numberElements(); i++) {
        
        LO p1ID = edgeElements->getElement(i).getNode( 0 );
        LO p2ID = edgeElements->getElement(i).getNode( 1 );
       
        for (int d=0; d<this->dim_; d++)
            newPoints[i][d] = ( (*pointsP1)[p1ID][d] + (*pointsP1)[p2ID][d] ) / 2.; // New point on middle of nodes
        

       	newFlags[i] = edgeElements->getElement(i).getFlag(); // New flags according to edge flags
		                
        const vec_LO_Type elementsOfEdge = edgeElements->getElementsOfEdge( i );
        const vec_GO_Type elementsGlobalOfEdge = edgeElements->getElementsOfEdgeGlobal( i );
                
        vec_GO_Type relevantElementsOfEdge(0);
        for (int j=0; j<elementsOfEdge.size(); j++) {
            if ( elementsOfEdge[j] != OTLO::invalid() )
                relevantElementsOfEdge.push_back( elementsOfEdge[j] );
        }

        // We need to determine the correct location in the P2 element as prescribed by a specific pattern
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

    // Set P2 Elements, by resetting 'old' elements and adding them with the new nodes per element
    this->elementsC_.reset(new Elements());
    
    LO numberLocalP1Nodes = meshP1->getPointsRepeated()->size();
    
    for (int i=0; i<elements->numberElements(); i++) {
        vec_int_Type feNodeList = elements->getElement( i ).getVectorNodeListNonConst(); // get a copy
        for (int j=0; j<newElementNodes[i].size(); j++)
            feNodeList.push_back( newElementNodes[i][j] + numberLocalP1Nodes );
        FiniteElement feP2( feNodeList,elements->getElement( i ).getFlag() );
        this->elementsC_->addElement(feP2);
    }
    
    if (verbose)
        std::cout << "done --" << std::endl;
    
    if (verbose)
        std::cout << "-- Setting global point IDs and building repeated P2 map and repeated points ... " << std::flush;
    
	// Now the corresponding maps for nodes and flags are updated.
	// Generally the P2 points are added after the P1 points in the node lists. Thus the new IDs just need to be added at the 
	// end of the maps

    this->pointsRep_.reset(new std::vector<std::vector<double> >(meshP1->pointsRep_->size(),std::vector<double>(this->dim_,-1.)));
    *this->pointsRep_ = *meshP1->pointsRep_;
    this->bcFlagRep_.reset(new std::vector<int>(meshP1->bcFlagRep_->size()));
    *this->bcFlagRep_ = *meshP1->bcFlagRep_;

    this->pointsRep_->insert( this->pointsRep_->end(), newPoints.begin(), newPoints.end() );
    this->bcFlagRep_->insert( this->bcFlagRep_->end(), newFlags.begin(), newFlags.end());
    
    Teuchos::ArrayView<const GO> nodeList = meshP1->getMapRepeated()->getNodeElementList();
    std::vector<GO> vecGlobalIDs = Teuchos::createVector( nodeList );
    
	// To correctly determine the global nodeIDs of the new P2 nodes we set them via the global edge IDs.
	// As can be on more than one processor at a time, they are an easy way to determine the correct 
	// nodeIDs via the edges'globalID.
    for (int i=0; i<edgeElements->numberElements(); i++){
        vecGlobalIDs.push_back( edgeElements->getGlobalID( (LO) i ) + P1Offset );
        edgeElements->setMidpoint( i, edgeElements->getGlobalID( (LO) i ) + P1Offset);
    }
    Teuchos::RCP<std::vector<GO> > pointsRepGlobMapping = Teuchos::rcp( new std::vector<GO>( vecGlobalIDs ) );
    Teuchos::ArrayView<GO> pointsRepGlobMappingArray = Teuchos::arrayViewFromVector( *pointsRepGlobMapping );
    
    this->mapRepeated_.reset(new Map<LO,GO,NO>( Teuchos::OrdinalTraits<GO>::invalid(), pointsRepGlobMappingArray, 0, this->comm_) );
    
    if (verbose)
        std::cout << "done --" << std::endl;
    
    if (verbose)
        std::cout << "-- Building unique P2 map, setting unique points and setting P2 elements ... " << std::flush;

    this->mapUnique_ = this->mapRepeated_->buildUniqueMap( this->rankRange_ );
    
    this->pointsUni_.reset(new std::vector<std::vector<double> >( this->mapUnique_->getNodeNumElements(), std::vector<double>(this->dim_,-1. ) ) );
    this->bcFlagUni_.reset( new std::vector<int> ( this->mapUnique_->getNodeNumElements(), 0 ) );
    for (int i=0; i<this->mapUnique_->getNodeNumElements(); i++) {
        GO gid = this->mapUnique_->getGlobalElement( i );

        LO id = this->mapRepeated_->getLocalElement( this->mapUnique_->getGlobalElement( i ) );
        this->pointsUni_->at(i) = this->pointsRep_->at(id);
        this->bcFlagUni_->at(i) = this->bcFlagRep_->at(id);

    }
    // Setting local P2 Node as Midpoint.
    for (int i=0; i<edgeElements->numberElements(); i++){
        edgeElements->setMidpoint( i, this->mapRepeated_->getLocalElement(edgeElements->getGlobalID( (LO) i ) + P1Offset));
    }

	this->edgeElements_ = edgeElements;

    
    if (verbose)
        std::cout << "done --" << std::endl;
    
    if (verbose)
        std::cout << "-- Building P2 surface elements ... " << std::flush;
    
    setP2SurfaceElements( meshP1 );
    
    if (verbose)
        std::cout << "done --" << std::endl;
    
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
            // In some instances a element has only an edge as subelement. In that case this step is not required
            if(feSurf.getVectorNodeList().size() == this->dim_)
                this->addSurfaceP2Nodes(feP2, feSurf, surfacePermutation, this->dim_);
            
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
void MeshUnstructured<SC,LO,GO,NO>::addSurfaceP2Nodes( FiniteElement &feP2, const FiniteElement &surfFeP1, const vec2D_int_Type &surfacePermutation, int dim ){
    
    vec2D_LO_Type edgeElementList = this->edgeElements_->getElementsNodeList();
    vec_int_Type surfaceNodeP2(0);

    if (dim == 2) {
        const vec_int_Type surfaceNodeP1 = surfFeP1.getVectorNodeList();
        
        vec_int_Type tmpSurface = surfFeP1.getVectorNodeList();
        sort(tmpSurface.begin(),tmpSurface.end());
        
        auto it1 = find( edgeElementList.begin(),edgeElementList.end() , tmpSurface );
        int edgeID = distance( edgeElementList.begin() , it1 );
        
        surfaceNodeP2.push_back( surfaceNodeP1[0] );
        surfaceNodeP2.push_back( surfaceNodeP1[1] );
        
        surfaceNodeP2.push_back(this->edgeElements_->getMidpoint( edgeID ) );  
            
        int flag = surfFeP1.getFlag();
        FiniteElement feP2Surf( surfaceNodeP2, flag );
        
        if ( !feP2.subElementsInitialized() )
            feP2.initializeSubElements("P2",dim-1);
        feP2.addSubElement(feP2Surf);
            
    }
    else if (dim == 3){

        const vec_int_Type surfaceNodeP1 = surfFeP1.getVectorNodeList();
        vec2D_LO_Type edges(3,vec_LO_Type(2)); 
        edges[0] = {surfaceNodeP1[0] , surfaceNodeP1[1]};
        edges[1] = {surfaceNodeP1[1],surfaceNodeP1[2]};
        edges[2] =  {surfaceNodeP1[0],surfaceNodeP1[2]};

 
        surfaceNodeP2.push_back( surfaceNodeP1[0] );
        surfaceNodeP2.push_back( surfaceNodeP1[1] );
        surfaceNodeP2.push_back( surfaceNodeP1[2] );


        bool criticalEdge = false;
        for(int i=0; i<3 ; i++){

            sort(edges[i].begin(),edges[i].end());
        
            auto it1 = find( edgeElementList.begin(),edgeElementList.end() , edges[i] );
            int edgeID = distance( edgeElementList.begin() , it1 );
            if(it1 == edgeElementList.end())
               criticalEdge=true; // Some triangles can be -1 -1 -1 thus no edge can be found. Why does that happen??
            else                 
                surfaceNodeP2.push_back(this->edgeElements_->getMidpoint( edgeID ) );  

        }

        if(criticalEdge==false){

            
            int flag = surfFeP1.getFlag();
            FiniteElement feP2Surf( surfaceNodeP2, flag );
            
            if ( !feP2.subElementsInitialized() )
                feP2.initializeSubElements("P2",dim-1);
            feP2.addSubElement(feP2Surf);
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
            const vec_int_Type surfaceNodeP1Unsorted = surfFeP1.getVectorNodeList();
            vec_int_Type surfaceNodeP1 = surfFeP1.getVectorNodeList();
            //sort( surfaceNodeP1.begin(), surfaceNodeP1.end() );

            vec_int_Type index2(3, 0);
            for (int i = 0 ; i != index2.size() ; i++)
                index2[i] = i;
            
            sort(index2.begin(), index2.end(),
                 [&](const int& a, const int& b) {
                     return  surfaceNodeP1[a] < surfaceNodeP1[b];
                    }
                 );
            
            surfaceNodeP1 = sort_from_ref(surfaceNodeP1, index2);

            if ( tmpSurface[0] == surfaceNodeP1[0]  &&
                    tmpSurface[1] == surfaceNodeP1[1] &&
                    tmpSurface[2] == surfaceNodeP1[2]) {
                surfaceNodeP2.push_back( surfaceNodeP1Unsorted[0] );
                surfaceNodeP2.push_back( surfaceNodeP1Unsorted[1] );
                surfaceNodeP2.push_back( surfaceNodeP1Unsorted[2] );
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

                // Resorting according to original IDs from GMSH
                
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
int MeshUnstructured<SC,LO,GO,NO>::determineFlagP2( LO p1ID, LO p2ID, LO localEdgeID, vec2D_LO_Type& markedPoint ){
    

    ElementsPtr_Type elements = this->getElementsC();

    vec2D_int_Type permutation = elements->getElementEdgePermutation();

    EdgeElementsPtr_Type edgeElements = this->getEdgeElements();

    MapConstPtr_Type elementMap = this->getElementMap();

    int rank = this->comm_->getRank();

    int flag1 = ( *this->getBCFlagRepeated() )[p1ID];

    int flag2 = ( *this->getBCFlagRepeated() )[p2ID];

    int newFlag = std::numeric_limits<int>::max();

    if(flag1 == this->volumeID_ || flag2 == this->volumeID_ ) // one node is in an inner node, than the new node is an inner node aswell
        newFlag = this->volumeID_;
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
                if (newFlags.size() == 0 && newFlag > this->volumeID_)
                    newFlag = this->volumeID_; //do we need this?

                //If we found a line element, then we choose this flag
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
void MeshUnstructured<SC,LO,GO,NO>::setMeshFileName(std::string meshFileName, std::string delimiter){

    meshFileName_ = meshFileName;
    delimiter_ = delimiter;
}

/*!

 \brief Not all edges are marked with a flag in the beginning. In order to set the correct flags to new points we assign the edge flag of the edge they originated from, similar to the function determineEdgeFlagP2, but this function uses the edgeMap. 
	@todo Elaboration of Flag Assignment Process. i.e. Lowest flag is used etc.

*/


template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::assignEdgeFlags(){

	vec2D_LO_Type markedPoints(0);
	EdgeElementsPtr_Type edgeElements = this->getEdgeElements();
	ElementsPtr_Type elements = this->getElementsC();
	MapConstPtr_Type mapRepeatedP1 = this->getMapRepeated();
	vec_int_Type newFlags(edgeElements->numberElements());

	MapConstPtr_Type edgeMap = this->getEdgeMap();

	vec_int_Type markedTrue(edgeElements->numberElements());
	for(int i=0; i< edgeElements->numberElements() ; i++){
		LO p1ID =edgeElements->getElement(i).getNode(0);
		LO p2ID =edgeElements->getElement(i).getNode(1);
		newFlags[i]=this->determineFlagP2(p1ID, p2ID, i , markedPoints );
		if(newFlags[i] != -1){ // questionable point that were given a flag, but that is not certain yet
			vec_LO_Type elementsOfEdge = edgeElements->getElementsOfEdge( (int) i );		
	   		for (int j=0; j<elementsOfEdge.size(); j++) {
	       		if ( elementsOfEdge[j] == -1 ) 
					markedTrue[i] =1;
			}
		}	
	}
	// It is possible that a Edge is conencted to two surfaces with different Flags, that are also on different Processors
	// This Leads to two different Flags for the same Edge
	// In order to counter that effect we check the interface edges of which we determined the flag via surfaces and check if they have the same flag and if not choose the lower one
	int maxRank = std::get<1>(this->rankRange_);
	if(maxRank >0 ){
		vec_GO_Type edgeSwitch(0);
		vec_LO_Type flags(0);
		for(int i=0; i<edgeElements->numberElements(); i++){
			if(markedTrue[i]==1){
				edgeSwitch.push_back(edgeMap->getGlobalElement(i));
				flags.push_back(newFlags[i]);
			}
		}

		// communticating elements across interface
		Teuchos::ArrayView<GO> edgeSwitchArray = Teuchos::arrayViewFromVector( edgeSwitch);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type(  Teuchos::OrdinalTraits<GO>::invalid(), edgeSwitchArray, 0, this->comm_) );

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

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar( (LO) 0 ); 
		isInterfaceElement_exp->exportFromVector( edgeFlags, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar( (LO) 0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_imp, false, "Insert");

		isInterfaceElement2_imp->exportFromVector(isInterfaceElement_exp, false, "Insert");

		edgeFlagsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		for(int i=0; i<edgeFlagsEntries.size(); i++){
			LO entry = edgeMap->getLocalElement(edgeSwitch[i]);
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
			edgesNeeded.push_back(edgeMap->getGlobalElement(i));
		}
		else{
			edgesActive.push_back(edgeMap->getGlobalElement(i));
			flagsTmp.push_back(newFlags[i]);
			
		}
	}

	// communticating elements across interface
	Teuchos::ArrayView<GO> edgesNeededArray = Teuchos::arrayViewFromVector( edgesNeeded);
	Teuchos::ArrayView<GO> edgesActiveArray = Teuchos::arrayViewFromVector( edgesActive);
	
	MapPtr_Type mapEdgesNeeded =
		Teuchos::rcp( new Map_Type(  Teuchos::OrdinalTraits<GO>::invalid(), edgesNeededArray, 0, this->comm_) );

	MapPtr_Type mapEdgesActive =
		Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), edgesActiveArray, 0, this->comm_) );

	MultiVectorLOPtr_Type flagsImport = Teuchos::rcp( new MultiVectorLO_Type( mapEdgesNeeded, 1 ) );
	flagsImport->putScalar(this->volumeID_);

	MultiVectorLOPtr_Type flagsExport = Teuchos::rcp( new MultiVectorLO_Type( mapEdgesActive, 1 ) );
	Teuchos::ArrayRCP< LO > flagExportEntries  = flagsExport->getDataNonConst(0);
	for(int i=0; i< flagExportEntries.size(); i++){
		flagExportEntries[i] = flagsTmp[i];
	}
	
	flagsImport->importFromVector(flagsExport, false, "Insert");

    Teuchos::ArrayRCP< LO > flagImportEntries  = flagsImport->getDataNonConst(0);
	for(int i=0; i<flagImportEntries.size(); i++){
		LO entry = edgeMap->getLocalElement(edgesNeeded[i]);
		if(newFlags[entry] ==-1){
			newFlags[entry] = flagImportEntries[i];
		}
	}

	for(int i=0; i< edgeElements->numberElements() ; i++){
		edgeElements->getElement(i).setFlag(newFlags[i]);
	}
    
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::buildEdgeMap(){

		int maxRank = std::get<1>(this->rankRange_);
		vec_GO_Type globalProcs(0);
		for (int i=0; i<= maxRank; i++)
			globalProcs.push_back(i);

		Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

		vec_GO_Type localProc(0);
		localProc.push_back(this->comm_->getRank());
		Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

		MapPtr_Type mapGlobalProc =
			Teuchos::rcp( new Map_Type(  Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, this->comm_) );

		MapPtr_Type mapProc =
			Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, this->comm_) );
		

		vec2D_int_Type interfaceEdgesLocalId(1,vec_int_Type(1));
		const int myRank = this->comm_->getRank();

		MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );

		// First we determine a Map only for the interface Nodes
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
		

		// Set unique edges IDs 
		// Setting the IDs of Edges that are uniquely on one Processor

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
			
		 }


		Teuchos::RCP<std::vector<GO>> edgesGlobMapping = Teuchos::rcp( new std::vector<GO>( vecGlobalIDsEdges ) );
		Teuchos::ArrayView<GO> edgesGlobMappingArray = Teuchos::arrayViewFromVector( *edgesGlobMapping);

		this->edgeMap_.reset(new Map<LO,GO,NO>( Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingArray, 0, this->comm_) );
		//this->edgeMap_->print();
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
        std::cout << "\n";
        std::cout << "  Read data of mesh " << meshFileName_<< ":\n";
    }

    meshReadSize ( meshFileName_, numNode, dim, numElement, orderElement, numSurface, orderSurface, numEdges, orderEdges );
    
    if (verbose) {
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "  Number of nodes = " << numNode << "\n";
        std::cout << "  Spatial dimension = " << dim << "\n";
        std::cout << "  Number of elements = " << numElement << "\n";
        std::cout << "  Element order = " << orderElement << "\n";
        std::cout << "  Number of surface elements = " << numSurface << "\n";
        std::cout << "  Surface element order = " << orderSurface << "\n";
        std::cout << "  Number of edge elements (for 3D) = " << numEdges << "\n";
        std::cout << "  Edges element order (for 3D) = " << orderEdges << "\n";
        
        std::cout << "\n";
        std::cout << "\n";
        std::cout << " Starting to read the data. \n";
    }

    
    this->elementOrder_ = orderElement;
    this->surfaceElementOrder_ = orderSurface;
    this->edgesElementOrder_ = orderEdges;
    this->numSurfaces_ = numSurface;
    this->numEdges_ = numEdges;
    this->numNodes_ = numNode;
    
    this->numElementsGlob_ = numElement;    
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readMeshEntity(std::string entityType){
    
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
//        std::cout << "### Starting to read surface data ... " << std::flush;
//
//    vec_int_Type    surfaceFlags( numSurfaces_, 0 );
//    vec_int_Type    surfacesCont( numSurfaces_* surfaceElementOrder_, 0 );
//
//    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
//    meshReadData ( meshFileName_, "surface", delimiter_, this->getDimension(), numSurfaces_, surfaceElementOrder_, surfacesCont, surfaceFlags );
//
//    if (verbose){
//        std::cout << "done." << std::endl;
//        std::cout << "### Setting surface data ... " << std::flush;
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
//        std::cout << "done." << std::endl;
//}

//template <class SC, class LO, class GO, class NO>
//void MeshUnstructured<SC,LO,GO,NO>::readLines(){
//    bool verbose ( this->comm_->getRank() == 0 );
//    if (verbose)
//        std::cout << "### Starting to read line data ... " << std::flush;
//
//    vec_int_Type    edgeFlags( numEdges_, 0 );
//    vec_int_Type    edgesCont( numEdges_* edgesElementOrder_, 0 );
//
//    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
//    meshReadData ( meshFileName_, "line", delimiter_, this->getDimension(), numEdges_, edgesElementOrder_, edgesCont, edgeFlags );
//
//
//    if (verbose){
//        std::cout << "done." << std::endl;
//        std::cout << "### Setting line data ... " << std::flush;
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
//        std::cout << "done." << std::endl;
//}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readNodes(){
    bool verbose ( this->comm_->getRank() == 0 );
    if (verbose)
        std::cout << "### Starting to read node data ... " << std::flush;

    vec_dbl_Type nodes(numNodes_ * this->getDimension(), 0.);
    vec_int_Type nodeFlags(numNodes_,0);
    meshReadData ( meshFileName_, "node", delimiter_, this->getDimension(), numNodes_, 3/*order of nodes is always 3*/, nodes, nodeFlags );
    
    if (verbose){
        std::cout << "done." << std::endl;
        std::cout << "### Setting node data ... " << std::flush;
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
        std::cout << "done." << std::endl;
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readElements(){
    bool verbose ( this->comm_->getRank() == 0 );
    if (verbose)
        std::cout << "### Starting to read element data ... " << std::flush;

    vec_int_Type    elementFlags( this->numElementsGlob_, 0 );
    vec_int_Type    elementsCont( this->numElementsGlob_* this->elementOrder_, 0 );
        
    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
    meshReadData ( meshFileName_, "element", delimiter_, this->getDimension(), this->numElementsGlob_, this->elementOrder_, elementsCont, elementFlags );

    if (verbose){
        std::cout << "done." << std::endl;
        std::cout << "### Setting element data ... " << std::flush;
    }
    ElementsPtr_Type elementsMesh = this->getElementsC();
    elementsMesh->setFiniteElementType("P1");
    elementsMesh->setDimension(this->getDimension());
    
	
    int id;
    for (int i=0; i<this->numElementsGlob_; i++) {
        vec_int_Type tmp(this->elementOrder_);
        for (int j=0; j<this->elementOrder_; j++){
            id = elementsCont.at(i*this->elementOrder_ + j) - 1;
            tmp.at(j) = id;
        }
        FiniteElement fe( tmp , elementFlags[i] );
        elementsMesh->addElement( fe );
    }
    
    if (verbose)
        std::cout << "done." << std::endl;
}


template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readSurfaces(){
    bool verbose ( this->comm_->getRank() == 0 );
    if (verbose)
        std::cout << "### Starting to read surface data ... " << std::flush;
    
    vec_int_Type    surfaceFlags( numSurfaces_, 0 );
    vec_int_Type    surfacesCont( numSurfaces_* this->surfaceElementOrder_, 0 );
        
    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
    meshReadData ( meshFileName_, "surface", delimiter_, this->getDimension(), numSurfaces_, this->surfaceElementOrder_, surfacesCont, surfaceFlags );
        
    if (verbose){
        std::cout << "done." << std::endl;
        std::cout << "### Setting surface data ... " << std::flush;
    }
    ElementsPtr_Type surfaceElementsMesh = this->getSurfaceElements();
    
    for (int i=0; i<numSurfaces_; i++) {
        vec_int_Type tmp(this->surfaceElementOrder_);
        for (int j=0; j<this->surfaceElementOrder_; j++)
            tmp.at(j) = surfacesCont.at( i * this->surfaceElementOrder_ + j ) - 1;// -1 to have start index 0
        
        //cout << " Reading surfaces " << tmp[0] << " " << tmp[1] << " " << tmp[2] << std::endl;
        //sort( tmp.begin(), tmp.end() ); // we sort here in order to identify the corresponding element faster! !!! We refrain from that so the numbering of surface element nodes remains the consisten with respct to normals.
        FiniteElement feSurface( tmp , surfaceFlags[i] );
        surfaceElementsMesh->addElement( feSurface );
    }
    
    
    if (verbose)
        std::cout << "done." << std::endl;
}

template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::readLines(){
    bool verbose ( this->comm_->getRank() == 0 );

    if (verbose)
        std::cout << "### Starting to read line data ... " << std::flush;

    vec_int_Type    edgeFlags( numEdges_, 0 );
    vec_int_Type    edgesCont( numEdges_* this->edgesElementOrder_, 0 );
        
    // We need the edges of surface elements to use special surface flags, which are only set on 1D line segments. We need to save these edges to determine the correct flag of P2 elements which might be build later
    meshReadData ( meshFileName_, "line", delimiter_, this->getDimension(), numEdges_, this->edgesElementOrder_, edgesCont, edgeFlags );

    
    if (verbose){
        std::cout << "done." << std::endl;
        std::cout << "### Setting line data ... " << std::flush;
    }
    
    ElementsPtr_Type edgeElementsMesh = this->getSurfaceEdgeElements();
    
    //Continous edge surface elements to FiniteElement object (only relevant in 3D)
    for (int i=0; i<numEdges_; i++) {
        vec_int_Type tmp(this->edgesElementOrder_);
        for (int j=0; j<this->edgesElementOrder_; j++)
            tmp.at(j) = edgesCont.at( i * this->edgesElementOrder_ + j ) - 1;// -1 to have start index 0
        
        sort( tmp.begin(), tmp.end() ); // we sort here in order to identify the corresponding element faster!
        FiniteElement feEdge( tmp , edgeFlags[i] );
        edgeElementsMesh->addElement( feEdge );
    }
    
    if (verbose)
        std::cout << "done." << std::endl;
}
template <class SC, class LO, class GO, class NO>
void MeshUnstructured<SC,LO,GO,NO>::exportMesh(MapConstPtr_Type mapUnique, MapConstPtr_Type mapRep, bool exportEdges, bool exportSurface, std::string meshName){

    
    std::ofstream myFile;
    if(this->comm_->getRank() ==0)
        myFile.open (meshName);

    bool verbose = (this->comm_->getRank() == 0);
    if(verbose){
        std::cout << " --------------------------------------" << std::endl;
        std::cout << " ------------ Exporting Mesh ----------" << std::endl;
        std::cout << " --------------------------------------" << std::endl;

    }
    // ################ Vertices #################
    int numberNodes = mapUnique->getGlobalNumElements();


    // ##########################################
    // Import missing nodes
    vec_GO_Type globalImportIDsNodes(0);
    for(int i = 0; i < this->mapUnique_->getGlobalNumElements(); i++)
    {
        if(this->mapUnique_->getLocalElement(i) == -1){
            globalImportIDsNodes.push_back(i);
        }
    }
	Teuchos::ArrayView<GO> globalNodesArrayImp = Teuchos::arrayViewFromVector( globalImportIDsNodes);
    // Map of global IDs with missing Elements
	MapPtr_Type mapNodesImport =
		Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), globalNodesArrayImp, 0, this->getComm()) );

    MultiVectorPtr_Type nodesImp = Teuchos::rcp( new MultiVector_Type( mapNodesImport, 1 ) );	
	Teuchos::ArrayRCP< SC > entriesNodesImp  = nodesImp->getDataNonConst(0);

    MultiVectorPtr_Type nodesExport = Teuchos::rcp( new MultiVector_Type( this->mapUnique_, 1 ) );	
	Teuchos::ArrayRCP< SC > entriesNodesExport  = nodesExport->getDataNonConst(0);

    vec2D_dbl_Type missingNodes(entriesNodesImp.size(),vec_dbl_Type(this->dim_+1));

    for(int j= 0; j < this->dim_; j++)
    {
        for(int i=0; i< entriesNodesExport.size(); i++)
            entriesNodesExport[i] =  this->getPointsUnique()->at(i).at(j);

        nodesImp->importFromVector(nodesExport, false, "Insert");
        for(int i=0; i< entriesNodesImp.size(); i++)
            missingNodes[i][j] =  entriesNodesImp[i];


    }

    // Lastly import flags
    for(int i=0; i< entriesNodesExport.size(); i++)
        entriesNodesExport[i] =  this->bcFlagUni_->at(i);

    nodesImp->importFromVector(nodesExport, false, "Insert");
    for(int i=0; i< entriesNodesImp.size(); i++)
        missingNodes[i][this->dim_] =  entriesNodesImp[i];

    // #############################################

    if(this->comm_->getRank() == 0){
        if(verbose)
            std::cout << " ----- Write Nodes .....";

        myFile << "MeshVersionFormatted 2" << std::endl;
        myFile << "Dimension" << " " << this->dim_ << std::endl;
        myFile << std::endl;
        myFile << "Vertices";
        myFile << std::endl;
        myFile << numberNodes;
        myFile << std::endl;
        for(int i = 0; i < mapUnique->getGlobalNumElements(); i++)
        {

            if(mapUnique->getLocalElement(i) != -1){
                LO id = mapUnique->getLocalElement(i);
                for(int j= 0; j < this->getPointsUnique()->at(id).size(); j++)
                {
                    myFile << this->getPointsUnique()->at(id).at(j);
                    myFile << " ";
                }
                if(this->dim_ == 2){
                    myFile << 0.0;
                    myFile << " ";
                }
                myFile << this->bcFlagUni_->at(id);
                myFile << std::endl;
            }
            else{
                LO id = mapNodesImport->getLocalElement(i);
                for(int j= 0; j < this->dim_; j++)
                {
                    myFile << missingNodes[id][j];
                    myFile << " ";
                }
                if(this->dim_ == 2){
                    myFile << 0.0;
                    myFile << " ";
                }
                myFile << missingNodes[id][this->dim_]; 
                myFile << std::endl;
            }



        }
        myFile << std::endl;
        if(verbose)
            std::cout << ".... done ----- " << '\n';

    }

    // ################ Edges #################
    if(exportEdges){
        if(verbose)
            std::cout << " ----- Write Edges .....";
        int dofsEdges=2;
        if(this->FEType_ == "P2")
            dofsEdges=3;
        // Assumption in 2D no edge is subelement of two elements, as they are only subelement if they are part of the surface
        if(this->dim_ == 2)
        {
            vec2D_GO_Type edgesSubelements(0,vec_GO_Type(dofsEdges+1)); // nodes + 
            ElementsPtr_Type elements = this->getElementsC();
            for(int T=0; T< elements->numberElements() ; T++){
                FiniteElement fe = elements->getElement( T );
                ElementsPtr_Type subEl = fe.getSubElements(); // might be null
                for (int surface=0; surface<fe.numSubElements(); surface++) {
                    FiniteElement feSub = subEl->getElement( surface  );
                    vec_GO_Type edgeTmp;
                    edgeTmp = {feSub.getVectorNodeList().at(0),feSub.getVectorNodeList().at(1),feSub.getFlag()};
                    if(this->FEType_=="P2")
                        edgeTmp = {feSub.getVectorNodeList().at(0),feSub.getVectorNodeList().at(1),feSub.getVectorNodeList().at(2),feSub.getFlag()};

                    edgesSubelements.push_back(edgeTmp);
                }
            
            }
            // Local number of subelement edges:
            int numSubEl = edgesSubelements.size();
            
		    int maxRank = std::get<1>(this->rankRange_);
            vec_GO_Type globalProcs(0);
            for (int i=0; i<= maxRank; i++)
                globalProcs.push_back(i);

            Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

            vec_GO_Type localProc(0);
            localProc.push_back(this->comm_->getRank());
            Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

            MapPtr_Type mapGlobalProc =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, this->comm_) );

            MapPtr_Type mapProc =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, this->comm_) );
            
            MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );
            exportLocalEntry->putScalar( (LO) numSubEl );

            // map equal to the original edgeMap with zero entries
            MultiVectorLOPtr_Type numEdgesProc= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
            numEdgesProc->putScalar( (LO) 0 ); 
            numEdgesProc->importFromVector( exportLocalEntry, false, "Insert");

            // number of new elements per Processor
            Teuchos::ArrayRCP< const LO > edgesRanks = numEdgesProc->getData(0);
		    vec_GO_Type vecGlobalIDsElements(0);

            GO procOffsetElements=0;
            int myRank = this->comm_->getRank();
            for(int i=0; i< myRank; i++)
                procOffsetElements = procOffsetElements + edgesRanks[i];

            for (int i=0; i<numSubEl; i++){
                vecGlobalIDsElements.push_back( i +  procOffsetElements);
            }

            Teuchos::RCP<std::vector<GO> > edgesGlobMapping = Teuchos::rcp( new std::vector<GO>( vecGlobalIDsElements ) );
            Teuchos::ArrayView<GO> edgesGlobMappingArray = Teuchos::arrayViewFromVector( *edgesGlobMapping);
            MapPtr_Type mapEdgesExport =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingArray, 0, this->getComm()) );

            vec_GO_Type vecGlobalIDsEdgesImport(0);
            int maxIndex = mapEdgesExport->getMaxAllGlobalIndex();
            if(this->comm_->getRank() == 0){
                for(int i=numSubEl; i<maxIndex+1 ; i++)
                    vecGlobalIDsEdgesImport.push_back(i);
            }
            Teuchos::RCP<std::vector<GO> > edgesGlobMappingImport = Teuchos::rcp( new std::vector<GO>( vecGlobalIDsEdgesImport ) );
            Teuchos::ArrayView<GO> edgesGlobMappingImportArray = Teuchos::arrayViewFromVector( *edgesGlobMappingImport);
            MapPtr_Type mapEdgesImport =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingImportArray, 0, this->getComm()) );

            int numberEdges = maxIndex+1;

            MultiVectorPtr_Type edgesImp = Teuchos::rcp( new MultiVector_Type( mapEdgesImport, 1 ) );	
            Teuchos::ArrayRCP< SC > entriesEdgesImp  = edgesImp->getDataNonConst(0);

            MultiVectorPtr_Type edgesExport = Teuchos::rcp( new MultiVector_Type( mapEdgesExport , 1 ) );	
            Teuchos::ArrayRCP< SC > entriesEdgesExport  = edgesExport->getDataNonConst(0);


            vec2D_dbl_Type missingEdges(entriesEdgesImp.size(),vec_dbl_Type(dofsEdges+1));

            for(int j= 0; j < dofsEdges+1; j++)
            {
                for(int i=0; i< entriesEdgesExport.size(); i++){
                    if(j<dofsEdges)
                        entriesEdgesExport[i] =  mapRep->getGlobalElement(edgesSubelements[i][j]);
                    else
                        entriesEdgesExport[i] =  edgesSubelements[i][j];
                }
                

                edgesImp->importFromVector(edgesExport, false, "Insert");
                for(int i=0; i< entriesEdgesImp.size(); i++)
                    missingEdges[i][j] =  entriesEdgesImp[i];

            }
            if(this->comm_->getRank() == 0){

                myFile << "Edges";
                myFile << std::endl;
                myFile << numberEdges;
                myFile << std::endl;
                for(int i = 0; i < mapEdgesExport->getGlobalNumElements(); i++)
                {
                    if(mapEdgesExport->getLocalElement(i) != -1){
                        LO id = this->edgeMap_->getLocalElement(i);

                        for(int j= 0; j < dofsEdges; j++)
                        {
                            myFile << mapRep->getGlobalElement(edgesSubelements[i][j])+1;
                            myFile << " ";
                        }
                        
                        myFile << edgesSubelements[i][dofsEdges]; // Flag
                        myFile << std::endl;
                        
                        
                    }
                    else{
                        LO id = mapEdgesImport->getLocalElement(i);
                        for(int j= 0; j < dofsEdges; j++)
                        {
                            myFile << missingEdges[id][j]+1;
                            myFile << " ";
                        }
                        myFile << missingEdges[id][dofsEdges]; // Flag
                        myFile << std::endl;
                        
                    }

                }
                myFile << std::endl;
            }

            if(verbose)
                std::cout << ".... done ---- " << '\n';

        }
        else if(this->dim_ == 3){         

            if(this->FEType_ != "P2") // In case of P2 this function was already called!
                this->assignEdgeFlags();
                
            vec_GO_Type globalImportIDsEdges(0);

            MapConstPtr_Type edgeMapUnique =  this->edgeMap_->buildUniqueMap( this->rankRange_ );

            vec2D_int_Type edgesUnique(edgeMapUnique->getNodeNumElements(),vec_int_Type(2,0));
		    for (int i=0; i < edgeMapUnique->getNodeNumElements(); i++) {
                GO gid = edgeMapUnique->getGlobalElement( i );
                LO id = this->edgeMap_->getLocalElement( gid );
                edgesUnique[i] = this->edgeElements_->getElement(id).getVectorNodeListNonConst();

            }
            if(this->comm_->getRank() == 0){
                for(int i = 0; i < edgeMapUnique->getGlobalNumElements(); i++)
                {
                    if(edgeMapUnique->getLocalElement(i) == -1){
                        globalImportIDsEdges.push_back(i);
                    }
                }
            }
            int numberEdges = edgeMapUnique->getGlobalNumElements();

            Teuchos::ArrayView<GO> globalEdgesArrayImp = Teuchos::arrayViewFromVector( globalImportIDsEdges);
            // Map of global IDs with missing Elements
            MapPtr_Type mapEdgesImport =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), globalEdgesArrayImp, 0, this->getComm()) );

            MultiVectorPtr_Type edgesImp = Teuchos::rcp( new MultiVector_Type( mapEdgesImport, 1 ) );	
            Teuchos::ArrayRCP< SC > entriesEdgesImp  = edgesImp->getDataNonConst(0);

            MultiVectorPtr_Type edgesExport = Teuchos::rcp( new MultiVector_Type( edgeMapUnique, 1 ) );	
            Teuchos::ArrayRCP< SC > entriesEdgesExport  = edgesExport->getDataNonConst(0);

            vec2D_dbl_Type missingEdges(entriesEdgesImp.size(),vec_dbl_Type(dofsEdges+1));

            for(int j= 0; j < 2; j++)
            {
                for(int i=0; i< entriesEdgesExport.size(); i++)
                    entriesEdgesExport[i] =  mapRep->getGlobalElement(edgesUnique[i][j])+1;

                edgesImp->importFromVector(edgesExport, false, "Insert");
                for(int i=0; i< entriesEdgesImp.size(); i++)
                    missingEdges[i][j] =  entriesEdgesImp[i];

            }
            if(this->FEType_ == "P2"){
                for(int i=0; i< entriesEdgesExport.size(); i++){
                    GO gid = edgeMapUnique->getGlobalElement( i );
                    LO id = this->edgeMap_->getLocalElement( gid );
                    entriesEdgesExport[i] =  mapRep->getGlobalElement(this->edgeElements_->getMidpoint(id))+1;
                }

                edgesImp->importFromVector(edgesExport, false, "Insert");
                for(int i=0; i< entriesEdgesImp.size(); i++)
                    missingEdges[i][2] =  entriesEdgesImp[i];

            }
            // Lastly import flags
            for(int i=0; i< entriesEdgesExport.size(); i++){
                GO gid = edgeMapUnique->getGlobalElement( i );
                LO id = this->edgeMap_->getLocalElement( gid);
                entriesEdgesExport[i] =  this->edgeElements_->getElement(id).getFlag();
            }

            edgesImp->importFromVector(edgesExport, false, "Insert");
            for(int i=0; i< entriesEdgesImp.size(); i++)
                missingEdges[i][dofsEdges] =  entriesEdgesImp[i];

             // ########## Writing #######################
            if(this->comm_->getRank() == 0){


                vec2D_int_Type edges(0,vec_int_Type(0));
                for(int i = 0; i < edgeMapUnique->getGlobalNumElements(); i++)
                {
                    vec_int_Type entry(0);
                    if(edgeMapUnique->getLocalElement(i) != -1){
                        LO id = this->edgeMap_->getLocalElement(i);
                        if(this->edgeElements_->getElement(id).getFlag() < this->volumeID_)
                        {
                            for(int j= 0; j < 2; j++)
                            {
                                entry.push_back(mapRep->getGlobalElement(this->edgeElements_->getElement(id).getVectorNodeList().at(j))+1);
                            }
                            if(this->FEType_ == "P2")
                                entry.push_back(mapRep->getGlobalElement(this->edgeElements_->getMidpoint(id))+1);

                            entry.push_back(this->edgeElements_->getElement(id).getFlag());
    
                            edges.push_back(entry);

                        }
                        
                    }
                    else{
                        LO id = mapEdgesImport->getLocalElement(i);
                        if(missingEdges[id][dofsEdges]< this->volumeID_){

                            for(int j= 0; j < dofsEdges; j++)
                            {
                                entry.push_back(missingEdges[id][j]);
                            }
                           entry.push_back(missingEdges[id][dofsEdges]); 

                           edges.push_back(entry);

                        }
                    }
                }


                myFile << "Edges";
                myFile << std::endl;
                myFile << edges.size();
                myFile << std::endl;

                for(int i = 0; i < edges.size(); i++)
                {
                    
                    for(int j= 0; j < 2; j++)
                    {
                        myFile << edges[i][j]; //mapRep->getGlobalElement(this->edgeElements_->getElement(id).getVectorNodeList().at(j))+1;
                        myFile << " ";
                    }
                    if(this->FEType_ == "P2")
                        myFile << edges[i][2]; //mapRep->getGlobalElement(this->edgeElements_->getMidpoint(id))+1 << " ";

                    myFile << edges[i][edges[i].size()-1]; // this->edgeElements_->getElement(id).getFlag(); last entry
                    myFile << std::endl;                    

                }

            }
            myFile << std::endl;
            if(verbose)
                std::cout << ".... done ----- " << '\n';    
        }
        

    }
    this->comm_->barrier();
      // ################ Surfaces #################
    if(exportSurface && this->dim_ >2){
        if(verbose)
            std::cout << " ----- Write Surfaces ...";
        int dofsSurfaces=3;
        if(this->FEType_ == "P2")
            dofsSurfaces=6;
        // Assumption in 2D no edge is subelement of two elements, as they are only subelement if they are part of the surface
        if(this->dim_ == 3)
        {
            vec2D_GO_Type surfacesSubelements(0,vec_GO_Type(dofsSurfaces+1)); // nodes + 
            ElementsPtr_Type elements = this->getElementsC();
            for(int T=0; T< elements->numberElements() ; T++){
                FiniteElement fe = elements->getElement( T );
                for (int surface=0; surface<fe.numSubElements(); surface++) {
                    ElementsPtr_Type subEl = fe.getSubElements(); // might be null
                    if(subEl->getDimension() == this->dim_-1 ){
                        FiniteElement feSub = subEl->getElement( surface  );
                        vec_GO_Type surfaceTmp;
                        surfaceTmp = {feSub.getVectorNodeList().at(0),feSub.getVectorNodeList().at(1),feSub.getVectorNodeList().at(2),feSub.getFlag()};
                        if(this->FEType_=="P2")
                            surfaceTmp = {feSub.getVectorNodeList().at(0),feSub.getVectorNodeList().at(1),feSub.getVectorNodeList().at(2),feSub.getVectorNodeList().at(3),feSub.getVectorNodeList().at(4),feSub.getVectorNodeList().at(5),feSub.getFlag()};

                        surfacesSubelements.push_back(surfaceTmp);
                    }
                }
            }
            // Local number of subelement edges:
            int numSubEl = surfacesSubelements.size();
           
		    int maxRank = std::get<1>(this->rankRange_);
            vec_GO_Type globalProcs(0);
            for (int i=0; i<= maxRank; i++)
                globalProcs.push_back(i);

            Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

            vec_GO_Type localProc(0);
            localProc.push_back(this->comm_->getRank());
            Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

            MapPtr_Type mapGlobalProc =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, this->comm_) );

            MapPtr_Type mapProc =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, this->comm_) );
            
            MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );
            exportLocalEntry->putScalar( (LO) numSubEl );
            // map equal to the original edgeMap with zero entries
            MultiVectorLOPtr_Type numSurfacesProc= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
            numSurfacesProc->putScalar( (LO) 0 ); 
            numSurfacesProc->importFromVector( exportLocalEntry, false, "Insert");
            // number of new elements per Processor
            Teuchos::ArrayRCP< const LO > surfacesRanks = numSurfacesProc->getData(0);
		    vec_GO_Type vecGlobalIDsElements(0);

            GO procOffsetElements=0;
            int myRank = this->comm_->getRank();
            for(int i=0; i< myRank; i++)
                procOffsetElements = procOffsetElements + surfacesRanks[i];

            for (int i=0; i<numSubEl; i++){
                vecGlobalIDsElements.push_back( i +  procOffsetElements);
            }

            Teuchos::RCP<std::vector<GO> > surfacesGlobMapping = Teuchos::rcp( new std::vector<GO>( vecGlobalIDsElements ) );
            Teuchos::ArrayView<GO> surfacesGlobMappingArray = Teuchos::arrayViewFromVector( *surfacesGlobMapping);
            MapPtr_Type mapSurfacesExport =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), surfacesGlobMappingArray, 0, this->getComm()) );


            vec_GO_Type vecGlobalIDsSurfacesImport(0);
            int maxIndex = mapSurfacesExport->getMaxAllGlobalIndex();
            if(this->comm_->getRank() == 0){
                for(int i=numSubEl; i<maxIndex+1 ; i++)
                    vecGlobalIDsSurfacesImport.push_back(i);
            }
            Teuchos::RCP<std::vector<GO> > surfacesGlobMappingImport = Teuchos::rcp( new std::vector<GO>( vecGlobalIDsSurfacesImport ) );
            Teuchos::ArrayView<GO> surfacesGlobMappingImportArray = Teuchos::arrayViewFromVector( *surfacesGlobMappingImport);
            MapPtr_Type mapSurfacesImport =
                Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), surfacesGlobMappingImportArray, 0, this->getComm()) );

            int numberSurfaces = maxIndex+1;

            MultiVectorPtr_Type surfacesImp = Teuchos::rcp( new MultiVector_Type( mapSurfacesImport, 1 ) );	
            Teuchos::ArrayRCP< SC > entriesSurfacesImp  = surfacesImp->getDataNonConst(0);

            MultiVectorPtr_Type surfacesExport = Teuchos::rcp( new MultiVector_Type( mapSurfacesExport , 1 ) );	
            Teuchos::ArrayRCP< SC > entriesSurfacesExport  = surfacesExport->getDataNonConst(0);


            vec2D_dbl_Type missingSurfaces(entriesSurfacesImp.size(),vec_dbl_Type(dofsSurfaces+1));

            for(int j= 0; j < dofsSurfaces+1; j++)
            {
                for(int i=0; i< entriesSurfacesExport.size(); i++){
                    if(j<dofsSurfaces)
                        entriesSurfacesExport[i] =  mapRep->getGlobalElement(surfacesSubelements[i][j]);
                    else{
                        entriesSurfacesExport[i] =  surfacesSubelements[i][j]; // FLAG!
                    }
                }
                surfacesImp->importFromVector(surfacesExport, false, "Insert");
                for(int i=0; i< entriesSurfacesImp.size(); i++)
                    missingSurfaces[i][j] =  entriesSurfacesImp[i];

            }
            if(this->comm_->getRank() == 0){

                myFile << "Triangles";
                myFile << std::endl;
                myFile << numberSurfaces;
                myFile << std::endl;
                for(int i = 0; i < mapSurfacesExport->getGlobalNumElements(); i++)
                {
                    if(mapSurfacesExport->getLocalElement(i) != -1){

                        for(int j= 0; j < dofsSurfaces; j++)
                        {
                            myFile << mapRep->getGlobalElement(surfacesSubelements[i][j])+1;
                            myFile << " ";
                        }
                        

                        myFile << surfacesSubelements[i][dofsSurfaces];
                        myFile << std::endl;
                        
                        
                    }
                    else{
                        LO id = mapSurfacesImport->getLocalElement(i);
                        for(int j= 0; j < dofsSurfaces; j++)
                        {
                            myFile << missingSurfaces[id][j]+1;
                            myFile << " ";
                        }
                        myFile << missingSurfaces[id][dofsSurfaces]; 
                        myFile << std::endl;
                        
                    }

                }
                myFile << std::endl;
            }

         
        }
        if(verbose)
            std::cout << "... done -----" << '\n';

    }


    // ################ Elements #################
    this->comm_->barrier();

    if(verbose)
        std::cout << " ----- Write Elements ...";

    // ###########################################
    // 1. Import missing Elements
    vec_GO_Type globalImportIDs(0);
    for(int i = 0; i < this->elementMap_->getGlobalNumElements(); i++)
    {
        if(this->elementMap_->getLocalElement(i) == -1){
            globalImportIDs.push_back(i);
        }
    }
	Teuchos::ArrayView<GO> globalElementArrayImp = Teuchos::arrayViewFromVector( globalImportIDs);
    // Map of global IDs with missing Elements
	MapPtr_Type mapElementImport =
		Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), globalElementArrayImp, 0, this->getComm()) );
   
    MultiVectorPtr_Type idsElement = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );	
	Teuchos::ArrayRCP< SC > entriesElement  = idsElement->getDataNonConst(0);

    MultiVectorPtr_Type idsElementExport = Teuchos::rcp( new MultiVector_Type( this->elementMap_, 1 ) );	
	Teuchos::ArrayRCP< SC > entriesElementExport  = idsElementExport->getDataNonConst(0);

    int dofsElement = this->getElements()->at(0).size();
    vec2D_GO_Type missingElements(entriesElement.size(),vec_GO_Type(dofsElement+1));

    for(int j= 0; j < this->getElements()->at(0).size(); j++)
    {
        for(int i=0; i< entriesElementExport.size(); i++)
            entriesElementExport[i] =  mapRep->getGlobalElement(this->getElements()->at(i).at(j))+1;

        idsElement->importFromVector(idsElementExport, false, "Insert");
        for(int i=0; i< entriesElement.size(); i++)
            missingElements[i][j] =  entriesElement[i];
    }
    for(int i=0; i< entriesElementExport.size(); i++)
        entriesElementExport[i] =  this->getElementsC()->getElement(i).getFlag();

    idsElement->importFromVector(idsElementExport, false, "Insert");
    for(int i=0; i< entriesElement.size(); i++)
        missingElements[i][dofsElement] =  entriesElement[i];
    // --------------------------------------
    if(this->comm_->getRank() == 0){

        if(this->dim_ == 2)
            myFile << "Triangles";
        else if(this->dim_ ==3)
            myFile << "Tetrahedra";

        myFile << std::endl;

        int numberElements = this->elementMap_->getGlobalNumElements();

        myFile << numberElements;
        myFile << std::endl;

        for(int i = 0; i < this->elementMap_->getGlobalNumElements(); i++)
        {
            if(this->elementMap_->getLocalElement(i) != -1){
                LO id = this->elementMap_->getLocalElement(i);

                for(int j= 0; j < this->getElements()->at(id).size(); j++)
                {
                    myFile << mapRep->getGlobalElement(this->getElements()->at(id).at(j))+1;
                    myFile << " ";
                }
                myFile << this->getElementsC()->getElement(id).getFlag(); 
                myFile << std::endl;
                
            }
            else{
                LO id = mapElementImport->getLocalElement(i);
                for(int j= 0; j < dofsElement; j++)
                {
                    myFile << missingElements[id][j];
                    myFile << " ";
                }
                myFile << missingElements[id][dofsElement]; 
                myFile << std::endl;
            }

        }
        if(verbose)
            std::cout << "... done ----- " << '\n';

    }
    if(this->comm_->getRank() ==0)
       myFile.close();
    if(verbose){
        std::cout << " --------------------------------------" << std::endl;
        std::cout << " ------- Finished exporting Mesh ------" << std::endl;
        std::cout << " ------- File Name: " << meshName << " -------" << std::endl;
        if(this->dim_ == 3)
             std::cout << " - Info: In 3D all edges are exported -" << std::endl;

        std::cout << " --------------------------------------" << std::endl;

    }
    this->comm_->barrier();


}

}
#endif
