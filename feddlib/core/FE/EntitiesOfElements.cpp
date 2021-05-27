#include "EntitiesOfElements.hpp"
/*!
 Definition of EntitiesOfElements
 
 @brief  EntitiesOfElements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */
namespace FEDD {

EntitiesOfElements::EntitiesOfElements():
Elements(),
elementsOfEntitiesGlobal_(0),
elementsOfEntitiesLocal_(0)
{
};

EntitiesOfElements::EntitiesOfElements( EntitiesOfElements& elements  ):
Elements( *(Teuchos::rcp_dynamic_cast<Elements_Type> ( Teuchos::rcpFromRef( elements ) ) ) ),
elementsOfEntitiesGlobal_( elements.elementsOfEntitiesGlobal_ ),
elementsOfEntitiesLocal_( elements.elementsOfEntitiesLocal_ )
{

};

void EntitiesOfElements::addEntity( FiniteElement& fe, vec_GO_Type& elementsOfEntity  ){
    this->addElement( fe );
    elementsOfEntitiesGlobal_.push_back( elementsOfEntity );
};

void EntitiesOfElements::addEntity( FiniteElement& fe, GO globalID ){
    this->addElement( fe );
    elementsOfEntitiesGlobal_.push_back( vec_GO_Type( 1, globalID ) );
};
void EntitiesOfElements::setElementsEntities( vec2D_GO_Type& elementsOfEntities ){
    // We assume that elementsOfEntiesGlobal_ has 1 element per redundant entity.
    // SortUnique was already called for the entities, so they are unique now.
    // Here we want to set all elements for an entity. We have the information in elementsOfEntities, which holds the information from unique to redundant entities.
    // and elementsOfEntitiesGlobal_, which holds the global element ID.
    
    vec2D_GO_Type elementsOfRedundantEntitiesGlobal =  elementsOfEntitiesGlobal_;

    vec2D_GO_Type newElementsOfEntitiesGlobal( this->numberElements(), vec_GO_Type(0) );
    
    for (int i=0; i<newElementsOfEntitiesGlobal.size(); i++) {
        for (int j=0; j<elementsOfEntities[i].size(); j++) {
            newElementsOfEntitiesGlobal[i].push_back( elementsOfRedundantEntitiesGlobal[ elementsOfEntities[i][j] ][0] );
        }
    }
    
    elementsOfEntitiesGlobal_ = newElementsOfEntitiesGlobal;
    
};

void EntitiesOfElements::partitionEntities( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated ){
   
    typedef Teuchos::OrdinalTraits<LO> OTLO;
    FE_vec_ptr_Type elementsTmp = Teuchos::rcp( new FE_vec_Type ( *elements_ ) );
    
    // Here it is assumed that elementsOfEntityGlobal_ is still the list for redundant entites. The correct partitioned and combine list is only set here.
    // We might want to make shift the setup of combined element information to the function sortUniqueAndSetGlobalIDs() in the future
    vec2D_GO_Type elementsOfEntitiesGlobalTmp = elementsOfEntitiesGlobal_;
    
    elementsOfEntitiesGlobal_.resize( 0 );
    elementsOfEntitiesLocal_.resize( 0 );

    this->elements_.reset( new FE_vec_Type ( ) );
    vec_GO_Type globaIDs = *(this->globalIDs_);
    this->globalIDs_.reset( new vec_GO_Type(0) );
    bool setEntity = false;
    for (int i=0; i<elementsTmp->size(); i++) {
        vec_GO_Type elementsOfThisEntitiesGlobal(0);
        vec_LO_Type elementsOfThisEntitiesLocal(0);
        vec_LO_Type localElementIDs(0);
        for (int j=0; j<elementsOfEntitiesGlobalTmp[i].size(); j++) {
            LO idLocal = elementMap->getLocalElement( elementsOfEntitiesGlobalTmp[i][j] );
            // we need to determine which ancestor elements are owned by this rank
            // the rank which owns the first element of all ancestor elements gets the active edge i
            elementsOfThisEntitiesLocal.push_back( idLocal );
            elementsOfThisEntitiesGlobal.push_back( elementsOfEntitiesGlobalTmp[i][j] );
            if (idLocal != OTLO::invalid())
                setEntity = true;
        }

        if ( setEntity ) {
            
            elementsOfEntitiesGlobal_.push_back( elementsOfThisEntitiesGlobal );
            elementsOfEntitiesLocal_.push_back( elementsOfThisEntitiesLocal );
            
            FiniteElement entity = (*elementsTmp)[i];
            vec_int_Type entityVec( entity.size() );
            for (int k=0; k<entityVec.size(); k++)
                entityVec[k] = nodeMapRepeated->getLocalElement( entity.getNode(k) );

            this->globalIDs_->push_back( globaIDs[i] );
            FiniteElement entityLocal( entityVec );
            elements_->push_back( entityLocal );
            setEntity = false;
        }
    }
    
};
    
void EntitiesOfElements::sortUniqueAndSetGlobalIDs( vec2D_GO_Type& combinedElements ){

    TEUCHOS_TEST_FOR_EXCEPTION( elements_.is_null(), std::runtime_error, "Elements not initialized! sortUniqueAndSetGlobalIDs( ) not possible.");
    vec2D_int_Type elementsVec( numberElements(), vec_int_Type( nodesPerElement(), -1 ) );
    for (int i=0; i<numberElements(); i++) {
        for (int j=0; j<nodesPerElement(); j++)
            elementsVec[i][j] = (*elements_)[i].getNode( j );
    }
    vec_GO_Type elementsOfEntitiesGlobalTmp( elementsOfEntitiesGlobal_.size() );
    for (int i=0; i<elementsOfEntitiesGlobalTmp.size(); i++)
        elementsOfEntitiesGlobalTmp[i] = elementsOfEntitiesGlobal_[i][0];
    
    //we also need to call sort but not unique on the elements belonging to the redundant edges
    make_unique( elementsVec, combinedElements, elementsOfEntitiesGlobalTmp );
    
    for (int i=0; i<elementsOfEntitiesGlobalTmp.size(); i++)
        elementsOfEntitiesGlobal_[i][0] = elementsOfEntitiesGlobalTmp[i];
    
    elements_.reset( new FE_vec_Type( ) );
    // for the globalIDs it is assumed that all global element information exists on every rank
    for (int i=0; i<elementsVec.size(); i++) {
        FiniteElement fe( elementsVec[i] );
        addElement( fe, (GO) i );
    }
}

const vec_LO_Type& EntitiesOfElements::getElementsOfEntity( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfEntitiesLocal_.size()-1 < i, std::logic_error, "No local elements for this edge." );

    return elementsOfEntitiesLocal_[i];
};
    
const vec_GO_Type& EntitiesOfElements::getElementsOfEntityGlobal( int i ){
    TEUCHOS_TEST_FOR_EXCEPTION( elementsOfEntitiesGlobal_.size()-1 < i, std::logic_error, "No global elements for this edge." );
    return elementsOfEntitiesGlobal_[i];
}
    
}
