#include "Elements.hpp"
/*!
 Definition of Elements
 
 @brief  Elements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {

Elements::Elements():
elements_(),
globalIDs_(),
elementsNodeList_(0),
FEType_("no Information"),
dim_(-1),
feDataInitialized_(false)
{
    elements_.reset(new FE_vec_Type());
    globalIDs_.reset( new vec_GO_Type( 0 ) );
};

Elements::Elements( std::string feType ):
elements_(),
globalIDs_(),
elementsNodeList_(0),
FEType_(feType),
dim_(-1),
feDataInitialized_(false)
{
    elements_.reset(new FE_vec_Type());
    globalIDs_.reset( new vec_GO_Type( 0 ) );
};

Elements::Elements( std::string feType, int dim ):
elements_(),
elementsNodeList_(0),
globalIDs_(),
FEType_(feType),
dim_(dim),
feDataInitialized_(false)
{
    elements_.reset(new FE_vec_Type());
    globalIDs_.reset( new vec_GO_Type( 0 ) );
};

Elements::Elements( Elements& Elements  ):
elements_(),
elementsNodeList_(0),
FEType_(Elements.FEType_),
dim_(Elements.dim_),
feDataInitialized_(false)
{
    elements_.reset( new FE_vec_Type( *Elements.elements_ ) );
    globalIDs_.reset( new vec_GO_Type( *Elements.globalIDs_ ) );
};

    
int Elements::numberElements(){
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! numberElements( ) not possible.");
    return elements_->size();
}

void Elements::addElement( FiniteElement& fe ){
    
    elements_->push_back( fe );
    globalIDs_->push_back( -1 );
}

void Elements::switchElement( int loc, FiniteElement& fe ){
    
    elements_->at( loc ) = fe;
}
    
void Elements::addElement( FiniteElement& fe, GO globalID ){
    
    elements_->push_back( fe );
    globalIDs_->push_back( globalID );
}
    
const FiniteElement& Elements::getElement( int i ) const {
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! GetElement( i ) not possible. Elements not initialized");
    TEUCHOS_TEST_FOR_EXCEPTION( elements_->size() - 1 < i, std::runtime_error, "Elements does not exist! GetElement( i ) not possible.");
    return elements_->at(i);
}

FiniteElement& Elements::getElement( int i ) {
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! GetElement( i ) not possible.");
    TEUCHOS_TEST_FOR_EXCEPTION( elements_->size() - 1 < i, std::runtime_error, "Elements does not exist! GetElement( i ) not possible.");
    return elements_->at(i);
}
    
Elements::GO Elements::getGlobalID( LO i ) const {
        TEUCHOS_TEST_FOR_EXCEPTION( globalIDs_->size()-1<i, std::runtime_error, "No global information for requested element.");
    return ( *globalIDs_ )[i];
}
    
int Elements::nodesPerElement( ){

    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! NodesPerElement( ) not possible.");
    
    if (elements_->size()>0)
        return elements_->at(0).size();
    else
        return -1;

}

    
void Elements::sortUnique( ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! sortUnique( ) not possible.");
    vec2D_int_Type ElementsVec( numberElements(), vec_int_Type( nodesPerElement(), -1 ) );
    for (int i=0; i<numberElements(); i++) {
        for (int j=0; j<nodesPerElement(); j++) {
            ElementsVec[i][j] = (*elements_)[i].getNode( j );
        }
    }
  
    make_unique( ElementsVec );
    elements_.reset( new FE_vec_Type( ) );
    for (int i=0; i<ElementsVec.size(); i++) {
        FiniteElement fe( ElementsVec[i] );
        addElement( fe );
    }
    
}

void Elements::sortUnique( vec2D_GO_Type& combinedElements ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! sortUnique( ) not possible.");
    vec2D_int_Type ElementsVec( numberElements(), vec_int_Type( nodesPerElement(), -1 ) );
    for (int i=0; i<numberElements(); i++) {
        for (int j=0; j<nodesPerElement(); j++) {
            ElementsVec[i][j] = (*elements_)[i].getNode( j );
        }
    }
    
    make_unique( ElementsVec, combinedElements );
    
    elements_.reset( new FE_vec_Type( ) );
    for (int i=0; i<ElementsVec.size(); i++) {
        FiniteElement fe( ElementsVec[i] );
        addElement( fe );
    }
    
}

    
void Elements::sortUniqueAndSetGlobalIDs( ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! sortUniqueAndSetGlobalIDs( ) not possible.");
    vec2D_int_Type ElementsVec( numberElements(), vec_int_Type( nodesPerElement(), -1 ) );
    for (int i=0; i<numberElements(); i++) {
        for (int j=0; j<nodesPerElement(); j++) {
            ElementsVec[i][j] = (*elements_)[i].getNode( j );
        }
    }
    
    make_unique( ElementsVec );
    
    elements_.reset( new FE_vec_Type( ) );
    // for the globalIDs it is assumed that all global element information exists on every rank
    for (int i=0; i<ElementsVec.size(); i++) {
        FiniteElement fe( ElementsVec[i] );
        addElement( fe, (GO) i );
    }

}
void Elements::sortUniqueAndSetGlobalIDs( vec2D_GO_Type& combinedElements ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! sortUniqueAndSetGlobalIDs( ) not possible.");
    vec2D_int_Type ElementsVec( numberElements(), vec_int_Type( nodesPerElement(), -1 ) );
    for (int i=0; i<numberElements(); i++) {
        for (int j=0; j<nodesPerElement(); j++) {
            ElementsVec[i][j] = (*elements_)[i].getNode( j );
        }
    }
    
    make_unique( ElementsVec, combinedElements );
    
    elements_.reset( new FE_vec_Type( ) );
    // for the globalIDs it is assumed that all global element information exists on every rank
    for (int i=0; i<ElementsVec.size(); i++) {
        FiniteElement fe( ElementsVec[i] );
        addElement( fe, (GO) i );
    }
}
    
void Elements::print( ){
    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! print( ) not possible.");
    for (int i=0; i<numberElements(); i++) {
        std::cout << "Element " << i << " : ";
        for (int j=0; j<nodesPerElement(); j++)
            std::cout << (*elements_)[i].getNode(j) << " ";
        std::cout << std::endl;
    }
}

void Elements::globalToLocalIDs( MapConstPtr_Type map ){
    
    for (int i=0; i<elements_->size(); i++)
        (*elements_)[i].globalToLocalIDs( map );
}

void Elements::setToCorrectElement( FiniteElement& feSub ){
    
    vec2D_int_Type permutation = getSubElementPermutation();
    for (int i=0; i<elements_->size(); i++)
        (*elements_)[i].addSubElementIfPart( feSub, permutation, FEType_, dim_ );
}

vec2D_int_Type Elements::getSubElementPermutation(){

    if (dim_ == 2) {
        vec2D_int_Type permutation( 3, vec_int_Type(2,0) );
        permutation[0][1] = 1;
        permutation[1][1] = 2;
        permutation[2][0] = 1;
        permutation[2][1] = 2;
        return permutation;
    }
    else if(dim_ == 3){
        if (FEType_ == "P1") {
            vec2D_int_Type permutation( 4, vec_int_Type(3,0) );
            permutation[0][1] = 1;
            permutation[0][2] = 2;
            
            permutation[1][1] = 1;
            permutation[1][2] = 3;
            
            permutation[3][0] = 1;
            permutation[3][1] = 2;
            permutation[3][2] = 3;
            
            permutation[4][1] = 2;
            permutation[4][2] = 3;
            
            return permutation;
        }
        else{
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Not implemented for this element type");
        }
    }
    vec2D_int_Type dummy;
    return dummy;
}

void Elements::setElementsNodeList(){    
    for(int i=0; i<numberElements(); i++)    
        elementsNodeList_.push_back(getElement(i).getVectorNodeList());

}

vec2D_LO_Type Elements::getElementsNodeList()
{
    if(elementsNodeList_.size() < 1)
        setElementsNodeList();

    return elementsNodeList_;
}

vec2D_int_Type Elements::getElementEdgePermutation(){
    if (dim_==1) {
        vec2D_int_Type permutation( 1, vec_int_Type(2,0) );
        permutation[0][1] = 1;
        return permutation;
    }
    else if (dim_ == 2) {
        if (FEType_ == "P1") {
            vec2D_int_Type permutation( 3, vec_int_Type(2,0) );
            permutation[0][1] = 1;
            permutation[1][1] = 2;
            permutation[2][0] = 1;
            permutation[2][1] = 2;
                        
            return permutation;
        }
        else{
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Not implemented for this element type");
        }
    }
    else if(dim_ == 3){
        if (FEType_ == "P1") {
            vec2D_int_Type permutation( 6, vec_int_Type(2,0) );
            permutation[0][1] = 1;
            
            permutation[1][1] = 2;
            
            permutation[2][1] = 3;
            
            permutation[3][0] = 1;
            permutation[3][1] = 2;
            
            permutation[4][0] = 1;
            permutation[4][1] = 3;
            
            permutation[5][0] = 2;
            permutation[5][1] = 3;
            return permutation;
        }
        else{
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Not implemented for this element type");
        }
    }
    vec2D_int_Type dummy;
    return dummy;
}
 
void Elements::initializeFEData( vec2D_dbl_ptr_Type pointsRep ){
    if (!feDataInitialized_) {
        vecBTinv_.resize( numberElements() );
        vecDetBTinv_.resize( numberElements() );
        SM_SC_Type BT(dim_);
        for (int i=0; i<numberElements(); i++) {
            buildTransformation( getElement(i).getVectorNodeList(), pointsRep, BT );
            vecDetBTinv_[i] = std::fabs( BT.computeInverse( vecBTinv_[i] ) );
        }
        feDataInitialized_ = true;
    }
}

void Elements::buildTransformation( const vec_int_Type& element,
                                    vec2D_dbl_ptr_Type pointsRep,
                                    SM_SC_Type& B ){

    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = element.at(0);
    if (FEType_[0]=='P') {
        for (UN j=0; j<B.size(); j++) {
            index = element.at(j+1);
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = pointsRep->at(index).at(i) - pointsRep->at(index0).at(i);
            }
        }
    }
    else if (FEType_[0]=='Q'){
        TEUCHOS_TEST_FOR_EXCEPTION( B.size()!=3, std::logic_error, "Transformation for quadrilateral Elements only in 3D.");
        std::vector<int> indexVec(3);
        indexVec[0] = element[1]; indexVec[1] = element[3]; indexVec[2] = element[4];
        for (UN j=0; j<B.size(); j++) {
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = ( pointsRep->at( indexVec[j] ).at(i) - pointsRep->at( index0 ).at(i) ) / 2.;
            }
        }
    }
}

const typename Elements::SM_SC_Type& Elements::getBTinv(int i){
    return vecBTinv_[i];
}

const double& Elements::getDetBTinv(int i){
    return vecDetBTinv_[i];
}

}
