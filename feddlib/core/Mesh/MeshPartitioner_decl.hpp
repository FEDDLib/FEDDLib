#ifndef MeshPartitioner_decl_hpp
#define MeshPartitioner_decl_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/FE/Domain.hpp"

#define FEDD_HAVE_METIS
#define FEDD_HAVE_PARMETIS

#ifdef FEDD_HAVE_METIS
#include "metis.h"
#endif
#ifdef FEDD_HAVE_PARMETIS
#include "parmetis.h"
#endif

/*!
 Defintion of MeshPartitioner
 
 @brief  MeshPartitioner
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshPartitioner {
    
public:
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef typename Domain_Type::DomainPtr_Type DomainPtr_Type;
    typedef typename Domain_Type::CommConstPtr_Type CommConstPtr_Type;
    typedef std::vector<DomainPtr_Type> DomainPtrArray_Type;
    typedef typename Domain_Type::MeshUnstr_Type MeshUnstr_Type;
    typedef typename Domain_Type::MeshUnstrPtr_Type MeshUnstrPtr_Type;
    typedef typename Domain_Type::MapConstPtr_Type MapConstPtr_Type;
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;    
    typedef EdgeElements EdgeElements_Type;
    typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;

    typedef std::vector<idx_t> vec_idx_Type; //Metis
    
    MeshPartitioner();

    MeshPartitioner( DomainPtrArray_Type domains, ParameterListPtr_Type pL, std::string feType, int dimension );
    
    ~MeshPartitioner();
    
    void readAndPartition();
        
    /*! Only used in 3D*/
    void setEdgesToSurfaces(int meshNumber);
    
    void setSurfacesToElements(int meshNumber);
    
    void partitionMesh( MeshUnstrPtr_Type& mesh, int meshNumber );

    void buildEdgeListParallel( MeshUnstrPtr_Type mesh, ElementsPtr_Type elementsGlobal );
    
    void buildEdgeList( MeshUnstrPtr_Type mesh, ElementsPtr_Type& elementsGlobal );

    void setLocalEdgeIndices(vec2D_int_Type &localEdgeIndices );

    void determineRanks();
    
    void determineRanksFromNumberRanks(vec_int_Type& ranks);
    
    void determineRanksFromFractions(vec_int_Type& fractions);
    
    void makeContinuousElements(ElementsPtr_Type elements, vec_idx_Type& eind_vec, vec_idx_Type& eptr_vec );
    
    void findAndSetSurfaces( vec2D_int_Type& surfElements_vec, vec_int_Type& surfElementsFlag_vec, FiniteElement& element, vec2D_int_Type& permutation, MapConstPtr_Type mapRepeated);

    void findAndSetSurfacesPartitioned( vec2D_int_Type& surfElements_vec, vec_int_Type& surfElementsFlag_vec, FiniteElement& element, vec2D_int_Type& permutation, vec_GO_Type& linearSurfacePartitionOffset, int globalElID);

    void setLocalSurfaceEdgeIndices( vec2D_int_Type &localSurfaceEdgeIndices, int edgesElementOrder );
    
    void findAndSetSurfaceEdges( vec2D_int_Type& edgeElements_vec, vec_int_Type& edgeElementsFlag_vec, FiniteElement& element, vec2D_int_Type& permutation, MapConstPtr_Type mapRepeated);
    
    int searchInSurfaces( vec2D_int_Type& surfaces, vec_int_Type searchSurface);
    
    void setLocalSurfaceIndices(vec2D_int_Type &localSurfaceIndices, int surfaceElementOrder );
    /* ###################################################################### */
private:    
    
    void readAndPartitionMesh( int meshNumber );
    
    ParameterListPtr_Type pList_;
    DomainPtrArray_Type domains_;
    CommConstPtr_Type comm_;
    std::string feType_;
    std::vector< tuple_intint_Type > rankRanges_;
    int dim_;
    };
}

#endif
