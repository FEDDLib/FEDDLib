#ifndef SurfaceElements_hpp
#define SurfaceElements_hpp
#include "Elements.hpp"
#include "EdgeElements.hpp"
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
    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;
    
    SurfaceElements();
    
    SurfaceElements( SurfaceElements& SurfaceElements  );

    void addSurface( FiniteElement& fe, GO globalID  );
    
    void addSurface( FiniteElement& fe, vec_GO_Type& elementsOfFace  );
    
    void setElementsSurface( vec2D_GO_Type&  elementsOfSurface );
    
    void partitionSurfaces( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated );
    
    void sortUniqueAndSetGlobalIDs(vec2D_GO_Type &combinedElements);
    
    void makeUniqueWithCombines( FE_vec_ptr_Type& elements, vec2D_GO_Type& combinedElements, vec2D_GO_Type& globaIDs );
    
    FE_vec_ptr_Type sort_from_ref( FE_vec_ptr_Type& elements, std::vector<int> const& reference );
    
    vec2D_GO_Type sort_from_ref( vec2D_GO_Type const& in, std::vector<int> const& reference );
    
    const vec_LO_Type& getElementsOfSurfaceLocal( int i );
    
    const vec_GO_Type& getElementsOfSurfaceGlobal( int i );
    
    vec2D_GO_Type getElementsOfSurfaceGlobal(){return elementsOfSurfaceGlobal_;}

    vec2D_LO_Type getElementsOfSurfaceLocal(){return elementsOfSurfaceLocal_;}
    
    const vec_int_Type getSurfacesOfElement( int i ); // returns the Surfaces of Element i ( , , )

    void matchSurfacesToElements(MapConstPtr_Type elementMap); // matches the corresponding Surfaces to the elements 

	void sortUniqueAndSetGlobalIDsParallel(MapConstPtr_Type elementMap, vec2D_GO_Type& combinedElements );

	void setElementsOfSurfaceLocalEntry(int index, int entry);

	void setElementsOfSurfaceGlobalEntry(int index, int entry);

	void setUpElementsOfSurface( MapConstPtr_Type elementMap, MapConstPtr_Type edgeMap, EdgeElementsPtr_Type edgeElements);
	

  private:
    
    vec2D_GO_Type elementsOfSurfaceGlobal_;
    vec2D_LO_Type elementsOfSurfaceLocal_;
	vec2D_LO_Type surfacesOfElements_;
    
};
}
#endif
