#include "EdgeElements.hpp"
/*!
 Definition of Elements
 
 @brief  Elements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {
    
    
EdgeElements::EdgeElements():
Elements(),
elementsOfEdgeGlobal_(0),
elementsOfEdgeLocal_(0),
edgesOfElements_(0),
midPointsInd_(0)
{

};

EdgeElements::EdgeElements( EdgeElements& elements  ):
Elements( *(Teuchos::rcp_dynamic_cast<Elements_Type> ( Teuchos::rcpFromRef( elements ) ) ) ),
elementsOfEdgeGlobal_( elements.elementsOfEdgeGlobal_ ),
elementsOfEdgeLocal_( elements.elementsOfEdgeLocal_ ),
edgesOfElements_(0),
midPointsInd_(0)
{

};

void EdgeElements::addEdge( FiniteElement& fe, GO globalID ){
    this->addElement( fe );
    elementsOfEdgeGlobal_.push_back( vec_GO_Type( 1, globalID ) ); // Global ID refers to the element, to which we add the edge
};
void EdgeElements::setElementsEdges( vec2D_GO_Type& elementsOfEdges ){
    // We assume that elementsOfEdgeGlobal_ has 1 element per redundant edge.
    // SortUnique was already called for the edges, so they are unique now.
    // Here we want to set all elements for an edge. We have the information in elementsOfEdges, which holds the information from unique to redundant edges.
    // and elementsOfEdgeGlobal_, which holds the global element ID.
    
    vec2D_GO_Type newElementsOfEdgeGlobal( this->numberElements(), vec_GO_Type(0) );
    for (int i=0; i<elementsOfEdges.size(); i++) {
        for (int j=0; j<elementsOfEdges[i].size(); j++)
            newElementsOfEdgeGlobal[i].push_back( elementsOfEdgeGlobal_[ elementsOfEdges[i][j] ][0] );
    }
    
    elementsOfEdgeGlobal_ = newElementsOfEdgeGlobal;
};

void EdgeElements::partitionEdges( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated ){
   
    typedef Teuchos::OrdinalTraits<LO> OTLO;
    FE_vec_ptr_Type elementsTmp = Teuchos::rcp( new FE_vec_Type ( *elements_ ) );
    
    // Here it is assumed that elementsOfEdgeGlobal_ is still the list for redundant edges. The correct partitioned and combine list is only set here.
    // We might want to make shift the setup of combined element information to the function sortUniqueAndSetGlobalIDs() in the future
    vec2D_GO_Type elementsOfEdgeGlobalTmp = elementsOfEdgeGlobal_;
    
    elementsOfEdgeGlobal_.resize( 0 );
    elementsOfEdgeLocal_.resize( 0 );

    this->elements_.reset( new FE_vec_Type ( ) );
    vec_GO_Type globaIDs = *(this->globalIDs_);
    this->globalIDs_.reset( new vec_GO_Type(0) );
    bool setEdge = false;
    for (int i=0; i<elementsTmp->size(); i++) {

        vec_GO_Type elementsOfThisEdgeGlobal(0);
        vec_LO_Type elementsOfThisEdgeLocal(0);
        vec_LO_Type localElementIDs(0);
        for (int j=0; j<elementsOfEdgeGlobalTmp[i].size(); j++) {

            LO idLocal = elementMap->getLocalElement( elementsOfEdgeGlobalTmp[i][j] );
            // we need to determine which ancestor elements are owned by this rank
            // the rank which owns the first element of all ancestor elements gets the active edge i
            elementsOfThisEdgeLocal.push_back( idLocal );
            elementsOfThisEdgeGlobal.push_back( elementsOfEdgeGlobalTmp[i][j] );
            if (idLocal != OTLO::invalid())
                setEdge = true;
        }

        if ( setEdge ) {
            
            elementsOfEdgeGlobal_.push_back( elementsOfThisEdgeGlobal );
            elementsOfEdgeLocal_.push_back( elementsOfThisEdgeLocal );
            
            FiniteElement edge = (*elementsTmp)[i];
            vec_int_Type edgeVec = { nodeMapRepeated->getLocalElement( edge.getNode(0) ),
                                     nodeMapRepeated->getLocalElement( edge.getNode(1) )
                                    };
            sort(edgeVec.begin(),edgeVec.end());
            
            this->globalIDs_->push_back( globaIDs[i] );
            FiniteElement edgeLocal( edgeVec );
            elements_->push_back( edgeLocal );
            setEdge = false;
        }
    }

    
};
    
void EdgeElements::sortUniqueAndSetGlobalIDs( vec2D_GO_Type& combinedElements ){

    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! sortUniqueAndSetGlobalIDs( ) not possible.");
    vec2D_int_Type elementsVec( numberElements(), vec_int_Type( nodesPerElement(), -1 ) );
    for (int i=0; i<numberElements(); i++) {
        for (int j=0; j<nodesPerElement(); j++) {
            elementsVec[i][j] = (*elements_)[i].getNode( j );
        }
    }
    vec_GO_Type elementsOfEdgeGlobalTmp( elementsOfEdgeGlobal_.size() );
    for (int i=0; i<elementsOfEdgeGlobalTmp.size(); i++)
        elementsOfEdgeGlobalTmp[i] = elementsOfEdgeGlobal_[i][0];
    
    //we also need to call sort but not unique on the elements belonging to the redundant edges
    makeUniqueWithCombines( elements_, combinedElements, elementsOfEdgeGlobal_ );
    
    globalIDs_->resize(0);
    for (int i=0; i<numberElements(); i++)
        globalIDs_->push_back(i);

}
    
void EdgeElements::makeUniqueWithCombines( FE_vec_ptr_Type& elements, vec2D_GO_Type& combinedElements, vec2D_GO_Type& globaIDs )
{
    // We assume that each inner vector of globalIDs has only one element which means that globalIDs can be represented by a vec_GO_Type. However, we want to use vec2D_GO_Type here because we would have to copy the values otherwise in the method sortUniqueAndSetGlobalIDs(), where this method is used
    {
        std::vector<int> index(elements->size());
        for (int i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        std::sort(index.begin(), index.end(),
                  [&](const int& a, const int& b) {
                      return  (*elements)[a].getVectorNodeList() < (*elements)[b].getVectorNodeList();
                  }
                  );
        elements = sort_from_ref( elements, index );
        globaIDs = sort_from_ref( globaIDs, index );
    }
    {
        std::vector<int> index(elements->size());
        for (int i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        combinedElements.resize( elements->size() );
        
        auto it = uniqueWithCombines( elements->begin(), elements->end(), combinedElements );

        elements->resize( distance( elements->begin(), it ) );

        combinedElements.resize( elements->size() );
    }
};

EdgeElements::FE_vec_ptr_Type EdgeElements::sort_from_ref( FE_vec_ptr_Type& elements,
                                                           std::vector<int> const& reference ) {
    FE_vec_ptr_Type ret = Teuchos::rcp(new FE_vec_Type(elements->size()));
    int const size = elements->size();
    for (long long i = 0; i < size; ++i)
        (*ret)[i] = (*elements)[reference[i]];
    
    return ret;
};

// We assume that vec2D_GO_Type in has only one inner element per outer entry
vec2D_GO_Type EdgeElements::sort_from_ref( vec2D_GO_Type const& in,
                             std::vector<int> const& reference ) {
    vec2D_GO_Type ret(in.size(), vec_GO_Type(1));
    
    int const size = in.size();
    for (int i = 0; i < size; ++i)
        ret[i][0] = in[reference[i]][0];
    
    return ret;
};

const vec_LO_Type& EdgeElements::getElementsOfEdge( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfEdgeLocal_.size()-1 < i, std::logic_error, "No local elements for this edge." );

    return elementsOfEdgeLocal_[i];
};
    
const vec_GO_Type& EdgeElements::getElementsOfEdgeGlobal( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfEdgeGlobal_.size()-1 < i, std::logic_error, "No global elements for this edge." );
    return elementsOfEdgeGlobal_[i];
};


// returns the edges of Element i
const vec_int_Type EdgeElements::getEdgesOfElement( int i ){
	TEUCHOS_TEST_FOR_EXCEPTION( edgesOfElements_.size()-1 < i, std::logic_error, "No edges for this Element, something went wrong :(" );
	return edgesOfElements_.at(i);
};

// function that matches the edges to the right elements, local indexing
void EdgeElements::matchEdgesToElements(MapConstPtr_Type elementMap){

	int numberOfElements=elementMap->getMaxLocalIndex() +1;
	
    vec2D_int_Type newEdgesOfElements(numberOfElements , vec_int_Type(0) );
	int j1;

	for(int i=0; i< numberElements() ; i ++){
		for(int j=0;j< elementsOfEdgeLocal_.at(i).size(); j++){
			if(elementsOfEdgeLocal_.at(i).at(j) != -1){
				j1 = elementsOfEdgeLocal_.at(i).at(j);
				newEdgesOfElements.at(j1).push_back(i);
			}
		}
		
	}
    edgesOfElements_ = newEdgesOfElements;

};

// Sets local Midpoint index 'nodeIndex' for edge 'elementIndex'
void EdgeElements::setMidpoint( int elementIndex, int nodeIndex ){
	if(midPointsInd_.size()<1)
		midPointsInd_.resize( numberElements(),-1 );

	midPointsInd_.at(elementIndex)=nodeIndex;   
};

// Returns local index of edge Midpoint
const int EdgeElements::getMidpoint( int elementIndex){
	
	if(midPointsInd_.size()<1)
		return -1;
	else
		return midPointsInd_.at(elementIndex); 
};

// Functions to set elementsOfEdgeLocal und elementsOfEdgeGlobal entries remotely 
// (used by meshUnstructured Refinement to complete vector updates)
void EdgeElements::setElementsOfEdgeLocalEntry(int index, int entry){
	elementsOfEdgeLocal_[index].push_back(entry);
}

void EdgeElements::setElementsOfEdgeGlobalEntry(int index, int entry){
	elementsOfEdgeGlobal_[index].push_back(entry);

}

// Parallel Functions for building edgelists, elementsOfEdgeGlobal and elementsOfEdgesLocal
// The difference is, that we dont build the elementsVec and use the elementMap to set the right global IDs to elementsOfEdgeGlobal
void EdgeElements::sortUniqueAndSetGlobalIDsParallel(MapConstPtr_Type elementMap, vec2D_GO_Type& combinedElements ){

	LO ind;
    for (int i=0; i<elementsOfEdgeGlobal_.size(); i++){
		ind =  elementMap->getGlobalElement( elementsOfEdgeGlobal_[i][0]);
        elementsOfEdgeGlobal_[i][0] = ind ;
	}

	makeUniqueWithCombines( elements_, combinedElements, elementsOfEdgeGlobal_ );

    globalIDs_->resize(0);
    for (int i=0; i<numberElements(); i++)
        globalIDs_->push_back(i);
	

	// We now have elemenetsOfEdgeGlobal_ but without the entries for elements on other procs. As of now elementsOfEdgeLocal_ is the same just missing the -1 entries, that we add later
	elementsOfEdgeLocal_.resize(elementsOfEdgeGlobal_.size());
	for(int i=0; i < elementsOfEdgeGlobal_.size() ;i++)
	{
		elementsOfEdgeLocal_[i].push_back(elementsOfEdgeGlobal_[i][0]);
		if(elementsOfEdgeGlobal_[i].size()>1)
			elementsOfEdgeLocal_[i].push_back(elementsOfEdgeGlobal_[i][1]);
	}
}

// Function that works similar to 'partitionEdges', but can be used when edges are already partitioned among processors
// IMPORTANT: this alone does not update elementsOfEdgeLocal and elementsOfEdgeGlobal completetly, it only sets the vectors up with information
// that is available on the own processor
// In order to completly set those vectors up, one has to communicate across the interface. This is done in meshUnstructuredRefinement with
// 'updateElementsOfEdgesLocalAndGlobal'
void EdgeElements::setUpElementsOfEdge( MapConstPtr_Type elementMap, MapConstPtr_Type edgeMap){

    // Here it is assumed that elementsOfEdgeGlobal_ is still the list for redundant edges. The correct partitioned and combine list is only set here.
    // We might want to make shift the setup of combined element information to the function sortUniqueAndSetGlobalIDs() in the future

    vec2D_GO_Type elementsOfEdgeGlobalTmp = elementsOfEdgeGlobal_;

    elementsOfEdgeGlobal_.resize( 0 );
    elementsOfEdgeLocal_.resize( 0 );

    vec_GO_Type globaIDs = *(this->globalIDs_);
    this->globalIDs_.reset( new vec_GO_Type(0) );

    for (int i=0; i<this->numberElements(); i++) {

        vec_GO_Type elementsOfThisEdgeGlobal(0);
        vec_LO_Type elementsOfThisEdgeLocal(0);
        vec_LO_Type localElementIDs(0);

        for (int j=0; j<elementsOfEdgeGlobalTmp[i].size(); j++) {

            LO idLocal = elementMap->getLocalElement( elementsOfEdgeGlobalTmp[i][j] );

            // we need to determine which ancestor elements are owned by this rank
            // the rank which owns the first element of all ancestor elements gets the active edge i
            elementsOfThisEdgeLocal.push_back( idLocal );
            elementsOfThisEdgeGlobal.push_back( elementsOfEdgeGlobalTmp[i][j] );

        }
        elementsOfEdgeGlobal_.push_back( elementsOfThisEdgeGlobal );
        elementsOfEdgeLocal_.push_back( elementsOfThisEdgeLocal );
        this->globalIDs_->push_back( edgeMap->getGlobalElement((LO)globaIDs[i]) );
    }
};
// Function that determines the global IDs of interface edges from elementsOfEdgesLocal
// In the seriell case thit would be an empty vector, as there are no edges between processors
vec_GO_Type EdgeElements::determineInterfaceEdges(MapConstPtr_Type edgeMap){
    vec_GO_Type globalInterfaceIDs(0);

    for(int i=0; i<this->numberElements(); i++){
        for(int j=0;j< elementsOfEdgeLocal_.at(i).size(); j++){
            if(elementsOfEdgeLocal_.at(i).at(j) == -1){
                globalInterfaceIDs.push_back(edgeMap->getGlobalElement(i));
                this->getElement(i).setInterfaceElement(true);
                j = elementsOfEdgeLocal_.at(i).size();
            }
        }
    }

    return globalInterfaceIDs;

};


}
