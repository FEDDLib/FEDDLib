#ifndef EdgeElements_hpp
#define EdgeElements_hpp
#include "Elements.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"
#include <Teuchos_OrdinalTraits.hpp>
/*!
 Declaration of EdgeElements
 
 @brief  EdgeElements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

class EdgeElements : public Elements {

  public:
    typedef default_lo LO;
    typedef default_go GO;
    typedef default_no NO;
    typedef Elements Elements_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;
    typedef Teuchos::RCP<FE_vec_Type> FE_vec_ptr_Type;
    
    EdgeElements();
    
    EdgeElements( EdgeElements& EdgeElements  );

    void addEdge( FiniteElement& fe, GO globalID  );
        
    void setElementsEdges( vec2D_GO_Type&  elementsOfEdge );
    
    void partitionEdges( MapConstPtr_Type elementMap, MapConstPtr_Type nodeMapRepeated);
    
    void sortUniqueAndSetGlobalIDs(vec2D_GO_Type &combinedElements);
    
    void makeUniqueWithCombines( FE_vec_ptr_Type& elements, vec2D_GO_Type& combinedElements, vec2D_GO_Type& globaIDs );
    
    FE_vec_ptr_Type sort_from_ref( FE_vec_ptr_Type& elements, std::vector<int> const& reference );
    
    vec2D_GO_Type sort_from_ref( vec2D_GO_Type const& in, std::vector<int> const& reference );
    
    const vec_LO_Type& getElementsOfEdge( int i );
    
    const vec_GO_Type& getElementsOfEdgeGlobal( int i );
    
    vec2D_GO_Type getElementsOfEdgeGlobal(){return elementsOfEdgeGlobal_;}

    vec2D_LO_Type getElementsOfEdgeLocal(){return elementsOfEdgeLocal_;}
    
    const vec_int_Type getEdgesOfElement( int i ); // returns the edges of Element i ( , , )

    void matchEdgesToElements( MapConstPtr_Type elementMap ); // matches the corresponding edges to the elements

    void setMidpoint( LO edgeId, LO localID );
	
    int getMidpoint( LO edgeId );

    LO findEdgeId ( LO nodeId1, LO nodeID2);

    void sortUniqueAndSetGlobalIDsParallel( MapConstPtr_Type elementMap, vec2D_GO_Type& combinedElements );

    void setElementsOfEdgeLocalEntry(int index, int entry);

    void setElementsOfEdgeGlobalEntry(int index, int entry);

    void setUpElementsOfEdge( MapConstPtr_Type elementMap, MapConstPtr_Type edgeMap );

    vec_GO_Type determineInterfaceEdges( MapConstPtr_Type edgeMap );
	
		
    
  private:
    
    vec2D_GO_Type elementsOfEdgeGlobal_;
    vec2D_LO_Type elementsOfEdgeLocal_;
    vec2D_int_Type edgesOfElements_; // Edges of Element i are being stored in row i of edgesOfElements_
	  vec_LO_Type midPointsLocalID_; // Returns the local index of the Node that is the midpoint of edge 'i'
    
};
}
#endif
