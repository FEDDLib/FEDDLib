//#include "SurfaceElements.hpp"
/*!
 Definition of Elements
 
 @brief  Elements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {
    
    
template<class ForwardIt, class GO>
ForwardIt uniqueWithCombines(ForwardIt first, ForwardIt last, vector<vector<GO> >& combines)
{
    
    if (first == last)
        return last;
    
    ForwardIt firstForDistance = first;
    ForwardIt result = first;
    combines[ distance( firstForDistance, result ) ].push_back( (GO) distance( firstForDistance, first ) );
    while (++first != last) {
        if (!(*result == *first) && ++result != first) {
            *result = std::move(*first);
            // also add the element which is the final unique element (the first in the sorted list)
            combines[ distance( firstForDistance, result ) ].push_back( (GO) distance( firstForDistance, first ) );
        }
        else{
            combines[ distance( firstForDistance, result ) ].push_back( (GO) distance( firstForDistance, first ) );
        }
    }
    return ++result;
}

template <typename T>
vector<T> sort_from_ref(
                        vector<T> const& in,
                        vector<int> const& reference
                        ) {
    vector<T> ret(in.size());
    
    int const size = in.size();
    for (int i = 0; i < size; ++i)
        ret[i] = in[reference[i]];
    
    return ret;
};
    
template <typename T, class GO>
void make_unique( vector<vector<T> >& in, vec2D_GO_Type& combinedElements, vector<GO>& globaIDs )
{
    {
        vector<int> index(in.size(), 0);
        for (int i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        sort(index.begin(), index.end(),
             [&](const int& a, const int& b) {
                 return  in[a] < in[b];
             }
             );
        in = sort_from_ref( in, index );
        globaIDs = sort_from_ref( globaIDs, index );
    }
    {
        vector<int> index(in.size(), 0);
        for (int i = 0 ; i != index.size() ; i++)
            index[i] = i;
        
        combinedElements.resize( in.size() );
        
        auto it = uniqueWithCombines( in.begin(), in.end(), combinedElements );
        
        in.resize( distance( in.begin(), it ) );
        combinedElements.resize( in.size() );
    }
};

SurfaceElements::SurfaceElements():
Elements(),
elementsOfEdgeGlobal_(0),
elementsOfEdgeLocal_(0)
{
//    elements_.reset(new FE_vec_Type());
//    globalIDs_.reset( new vec_GO_Type( 0 ) );
};

SurfaceElements::SurfaceElements( SurfaceElements& elements  ):
Elements( *(Teuchos::rcp_dynamic_cast<Elements_Type> ( Teuchos::rcpFromRef( elements ) ) ) ),
elementsOfEdgeGlobal_( elements.elementsOfEdgeGlobal_ ),
elementsOfEdgeLocal_( elements.elementsOfEdgeLocal_ )
{

};

void SurfaceElements::addEdge( FiniteElement& fe, vec_GO_Type& elementsOfEdge  ){
    this->addElement( fe );
    elementsOfEdgeGlobal_.push_back( elementsOfEdge );
};

void SurfaceElements::addEdge( FiniteElement& fe, GO globalID ){
    this->addElement( fe );
    elementsOfEdgeGlobal_.push_back( vec_GO_Type( 1, globalID ) );
};
void SurfaceElements::setElementsEdges( vec2D_GO_Type& elementsOfEdges ){
    // We assume that elementsOfEdgeGlobal_ has 1 element per redundant edge.
    // SortUnique was already called for the edges, so they are unique now.
    // Here we want to set all elements for an edge. We have the information in elementsOfEdges, which holds the information of the from unique to redundant edges.
    // and elementsOfEdgeGlobal_, which holds the global element ID.
    
    vec2D_GO_Type elementsOfRedundantEdgeGlobal =  elementsOfEdgeGlobal_;

    vec2D_GO_Type newElementsOfEdgeGlobal( this->numberElements(), vec_GO_Type(0) );
    
    for (int i=0; i<newElementsOfEdgeGlobal.size(); i++) {
        for (int j=0; j<elementsOfEdges[i].size(); j++) {
            newElementsOfEdgeGlobal[i].push_back( elementsOfRedundantEdgeGlobal[ elementsOfEdges[i][j] ][0] );
        }
    }
    
    elementsOfEdgeGlobal_ = newElementsOfEdgeGlobal;
    
};

void SurfaceElements::partitionEdges( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated ){
   
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
            
            this->globalIDs_->push_back( globaIDs[i] );
            FiniteElement edgeLocal( edgeVec );
            elements_->push_back( edgeLocal );
            setEdge = false;
        }
    }
    
};
    
void SurfaceElements::sortUniqueAndSetGlobalIDs( vec2D_GO_Type& combinedElements ){

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
    make_unique( elementsVec, combinedElements, elementsOfEdgeGlobalTmp );
    
    for (int i=0; i<elementsOfEdgeGlobalTmp.size(); i++)
        elementsOfEdgeGlobal_[i][0] = elementsOfEdgeGlobalTmp[i];
    
    elements_.reset( new FE_vec_Type( ) );
    // for the globalIDs it is assumed that all global element information exists on every rank
    for (int i=0; i<elementsVec.size(); i++) {
        FiniteElement fe( elementsVec[i] );
        addElement( fe, (GO) i );
    }
}

const vec_LO_Type& SurfaceElements::getElementsOfEdge( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfEdgeLocal_.size()-1 < i, std::logic_error, "No local elements for this edge." );

    return elementsOfEdgeLocal_[i];
};
    
const vec_LO_Type& SurfaceElements::getElementsOfEdgeGlobal( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfEdgeGlobal_.size()-1 < i, std::logic_error, "No global elements for this edge." );
    return elementsOfEdgeGlobal_[i];
}
    
}