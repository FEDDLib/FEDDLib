#ifndef EntitiesOfElements_hpp
#define EntitiesOfElements_hpp
#include "Elements.hpp"
#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"
#include <Teuchos_OrdinalTraits.hpp>
/*!
 Declaration of EntitiesOfElements
 
 @brief  EntitiesOfElements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

class EntitiesOfElements : public Elements {

  public:
    typedef default_lo LO;
    typedef default_go GO;
    typedef default_no NO;
    typedef Elements Elements_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;
    
    EntitiesOfElements();
    
    EntitiesOfElements( EntitiesOfElements& EntitiesOfElements  );

    void addEntity( FiniteElement& fe, GO globalID  );
    
    void addEntity( FiniteElement& fe, vec_GO_Type& elementsOfEntity  );
    
    void setElementsEntities( vec2D_GO_Type&  elementsOfEntities );
    
    void partitionEntities( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated );
    
    void sortUniqueAndSetGlobalIDs(vec2D_GO_Type &combinedElements);
    
    const vec_LO_Type& getElementsOfEntity( int i );
    
    const vec_GO_Type& getElementsOfEntityGlobal( int i );
    
    vec2D_GO_Type getElementsOfEntityGlobal(){return elementsOfEntitiesGlobal_;};
    
  private:
    
    vec2D_GO_Type elementsOfEntitiesGlobal_;
    vec2D_LO_Type elementsOfEntitiesLocal_;
    
};
}
#endif
