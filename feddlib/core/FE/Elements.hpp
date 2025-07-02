#ifndef Elements_hpp
#define Elements_hpp
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"
#include "FiniteElement.hpp"
#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_RCPBoostSharedPtrConversions.hpp>

/*!
 Declaration of Elements
 
 @brief  Elements
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
class FiniteElement;
class Elements {
  
public:
    
    typedef std::vector<FiniteElement> FE_vec_Type;
    typedef Teuchos::RCP<FE_vec_Type> FE_vec_ptr_Type;
    typedef default_lo LO;
    typedef default_go GO;
    typedef default_sc SC;
    typedef default_no NO;
    
    typedef Teuchos::RCP<const Map<LO,GO,NO> > MapConstPtr_Type;
    
    typedef SmallMatrix<SC> SM_SC_Type;
    typedef std::vector<SM_SC_Type> vecSM_SC_Type;
    
    Elements();
    
    Elements( std::string feType );

    Elements( std::string feType, int dim );
    
    Elements( Elements& Elements  );
    
    int numberElements();
    
    void addElement( FiniteElement& fe  );
    
    void addElement( FiniteElement& fe, GO globalID );

    void switchElement( int loc, FiniteElement& fe  );
    
    const FiniteElement& getElement( int i ) const;

    FiniteElement& getElement( int i );

    GO getGlobalID( LO i ) const;
    
    void setFiniteElementType( std::string feType ){ FEType_ = feType; };

    std::string getFiniteElementType( ){ return FEType_; };
    
    void setDimension( int dim ){ dim_ = dim; };
    
    int getDimension() {return dim_;};
    
    int nodesPerElement();
    
    void sortUnique();

    void sortUnique(vec2D_GO_Type &combinedElements);
    
    void sortUniqueAndSetGlobalIDs();
    
    void sortUniqueAndSetGlobalIDs(vec2D_GO_Type &combinedElements);
    
    void print();
    
    void globalToLocalIDs( MapConstPtr_Type map );
    
    /*! We set the FiniteElement of a lower order to the correct FiniteElement, i.e, an edge is set to the correct triangle(s) */
    void setToCorrectElement( FiniteElement& feSub );
        
    vec2D_int_Type getElementEdgePermutation();
    
    vec2D_int_Type getSubElementPermutation();
    
    vec2D_LO_Type getElementsNodeList();

    void setElementsNodeList();
    /*! Initialized BT^{-1} and det(BT^{-1}) for FE element T */
    void initializeFEData( vec2D_dbl_ptr_Type pointsRep );
    
    void buildTransformation(   const vec_int_Type& element,
                                vec2D_dbl_ptr_Type pointsRep,
                                SM_SC_Type& B );
    
    const SM_SC_Type& getBTinv(int i);
    
    const double& getDetBTinv(int i);
    
    FE_vec_ptr_Type elements_;
    vec_GO_ptr_Type globalIDs_;
    std::string FEType_;
    int dim_;
    
    vecSM_SC_Type vecBTinv_;
    vec_dbl_Type vecDetBTinv_;
    bool feDataInitialized_;

    vec2D_LO_Type elementsNodeList_;
      
};
}
#endif
