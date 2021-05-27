#ifndef SurfaceElements_hpp
#define SurfaceElements_hpp
#include "Elements.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"
#include <Teuchos_OrdinalTraits.hpp>
/*!
 Declaration of SurfaceElements
 
 @brief  SurfaceElements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

class SurfaceElements : public Elements {

  public:
    typedef default_lo LO;
    typedef default_go GO;
    typedef default_no NO;
    typedef Elements Elements_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;
    
    SurfaceElements();
    
    SurfaceElements( SurfaceElements& SurfaceElements  );

    void addSurface( FiniteElement& fe, GO globalID  );
    
    void addSurface( FiniteElement& fe, vec_GO_Type& elementsOfFace  );
    
    void setElementsEdges( vec2D_GO_Type&  elementsOfEdge );
    
    void partitionSurface( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated );
    
  private:
    
    vec2D_GO_Type elementsOfSurfaceGlobal_;
    vec2D_LO_Type elementsOfSurfaceLocal_;
    
};
}
#endif
