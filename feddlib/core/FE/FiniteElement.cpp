#include "FiniteElement.hpp"

/*!
 Definition of FiniteElement
 
 @brief  FiniteElement
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
FiniteElement::FiniteElement():
flag_(0),
subElements_(),
numSubElements_(0)
{

}

FiniteElement::FiniteElement( vec_int_Type& localNodeList ):
flag_(0),
subElements_(),
numSubElements_(0)
{
    setElement( localNodeList );
}

FiniteElement::FiniteElement( vec_int_Type& localNodeList, int elementFlag ):
flag_(elementFlag),
subElements_(),
numSubElements_(0)
{    
    setElement( localNodeList );
}

bool FiniteElement::operator==(const FiniteElement &other){
    return localNodeIDs_ == other.localNodeIDs_;
}

FiniteElement& FiniteElement::operator=(const FiniteElement& in){
    
    localNodeIDs_ = in.localNodeIDs_;
    flag_ = in.flag_;
    numSubElements_ = in.numSubElements_;
    subElements_ = in.subElements_;
    refinementType_ = in.refinementType_;  // Tag of finite Element
	isInterfaceElement_ = in.isInterfaceElement_;
	taggedForRefinement_ = in.taggedForRefinement_;
	predecessorElement_ = in.predecessorElement_; 

    return *this;
}
    

void FiniteElement::setElement( vec_LO_Type& localNodeIDs ){
    localNodeIDs_ = localNodeIDs;
}

// add elementFlag 
void FiniteElement::setFlag( int elementFlag ){
    flag_ = elementFlag;
}

int FiniteElement::getNode( int i ) const{
    TEUCHOS_TEST_FOR_EXCEPTION( i >= localNodeIDs_.size(), std::runtime_error, "Node out of range!");
    return localNodeIDs_[i];
}

int FiniteElement::numSubElements() {
    return numSubElements_;
}

bool FiniteElement::subElementsInitialized() {
    return !subElements_.is_null();
}

void FiniteElement::initializeSubElements( std::string feType, int dim ) {
    subElements_ = Teuchos::rcp( new Elements_Type( feType, dim ) );
}

void FiniteElement::addSubElement( FiniteElement& fe ){
    subElements_->addElement(fe);
    numSubElements_++;
}

void FiniteElement::setSubElements( ElementsPtr_Type& subElements ){
    subElements_ = subElements;
    numSubElements_ = subElements->numberElements();
}

void FiniteElement::globalToLocalIDs( MapConstPtr_Type map ){
    for (int i=0; i<localNodeIDs_.size(); i++)
        localNodeIDs_[i] = map->getLocalElement(localNodeIDs_[i]);
        
    if ( subElementsInitialized() )
        subElements_->globalToLocalIDs( map );
}


void FiniteElement::addSubElementIfPart( FiniteElement& feSub, const vec2D_int_Type& permutation, std::string& feType, int dim ){
    vec_LO_Type partIDsEl( permutation[0].size() );
    vec_LO_Type ids = feSub.getVectorNodeList();
    for (int i=0; i<permutation.size(); i++) {
        for (int j=0; j<partIDsEl.size(); j++)
            partIDsEl[j] = localNodeIDs_[permutation[i][j]];
        std::sort( partIDsEl.begin(),partIDsEl.end() );
        if (ids == partIDsEl){
            if (subElementsInitialized())
                addSubElement( feSub );
            else{
                initializeSubElements( feType, dim-1 );
                addSubElement( feSub );
            }
        }
    }
}

void FiniteElement::findEdgeFlagInSubElements( const vec_LO_Type& edgeIDs, vec_int_Type& flags, bool isSubElement, const vec2D_int_Type& permutation, bool& foundLineSegment ){

    if ( numSubElements_ > 0 ){
        vec2D_int_Type permutationSub = subElements_->getElementEdgePermutation();
        for (int i=0; i<numSubElements_ && !foundLineSegment; i++) {
            // We should go down to the line segments
            subElements_->getElement(i).findEdgeFlagInSubElements( edgeIDs, flags, true, permutationSub, foundLineSegment );
        }
    }
    if(isSubElement && !foundLineSegment){ //we are in a subelement, i.e., a line or a surface starting from a 3D element
        bool found = findEdgeInElement( edgeIDs, flags, permutation );
        if (found && size() == 2 ) { //we found a line segment
            // we delete all prior flags if they were found on a surface
            int tmpFlag = flags[ flags.size()-1 ];
            flags.resize(1);
            flags[0] = tmpFlag;
            foundLineSegment = true;
        }
    }
}

bool FiniteElement::findEdgeInElement( const vec_LO_Type& edgeIDs, vec_int_Type& flags, const vec2D_int_Type& permutation ){
    
    vec_LO_Type partIDsEl( permutation[0].size() );
    for (int i=0; i<permutation.size(); i++) {
        for (int j=0; j<partIDsEl.size(); j++){
            partIDsEl[j] = localNodeIDs_[permutation[i][j]];
        }

        std::sort( partIDsEl.begin(),partIDsEl.end() );
        if (edgeIDs == partIDsEl){        
            flags.push_back(flag_);
            return true;
        }
    }
    return false;
}

void FiniteElement::print(MapConstPtr_Type mapRepeated){
    std::cout << "FiniteElement IDs : " << std::flush;
    for (int i=0; i<localNodeIDs_.size(); i++) {
        std::cout << localNodeIDs_[i] << "(gID:"<< std::flush;
        if (!mapRepeated.is_null()) {
            std::cout << mapRepeated->getGlobalElement(localNodeIDs_[i]) << ")  " << std::flush;
        } else {
            std::cout << "?)  " << std::flush;
        }
    }
    std::cout << std::endl;
}

}
