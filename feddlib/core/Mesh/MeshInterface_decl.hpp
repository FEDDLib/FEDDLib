#ifndef MESHINTERFACE_decl_hpp
#define MESHINTERFACE_decl_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

/*!
 Declaration of MeshInterface

 @brief  MeshInterface
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD{
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshInterface  {

public:
    typedef MeshInterface<SC,LO,GO,NO> MeshInterface_Type;
    typedef Teuchos::RCP<MeshInterface_Type> MeshInterfacePtr_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;
    typedef typename Map_Type::Comm_Type Comm_Type;
    typedef typename Map_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

    
    typedef std::vector<GO> vec_GO_Type;
    typedef std::vector<vec_GO_Type> vec2D_GO_Type;
    typedef std::vector<vec2D_GO_Type> vec3D_GO_Type;
    typedef Teuchos::RCP<vec3D_GO_Type> vec3D_GO_ptr_Type;

    /* ###################################################################### */
    //
    MeshInterface();
    
    MeshInterface(CommConstPtr_Type comm);
    
    void determineInterface( vec2D_dbl_ptr_Type pointsRepThis,  vec2D_dbl_ptr_Type pointsRepOther, vec_int_ptr_Type flagThis, vec_int_ptr_Type flagOther, vec_int_Type relevant_flag_vec );

    void determineInterfaceParallelAndDistance( vec2D_dbl_ptr_Type pointsUniThis, vec2D_dbl_ptr_Type pointsUniOther, vec_int_ptr_Type flagThis, vec_int_ptr_Type flagOther, vec_int_Type relevant_flag_vec, MapConstPtr_Type mapUniThis, MapConstPtr_Type mapUniOther, vec_dbl_ptr_Type &distancesToInterface, vec2D_dbl_ptr_Type pointsRepThis, int dim );

    void calculateDistancesToInterfaceParallel( vec_dbl_ptr_Type &distancesToInterface, vec2D_dbl_Type &pointThis/*global interface, every proc has same information*/, vec2D_dbl_ptr_Type sourceNodesRep /*partitioned points*/);
    
    int isPartialCouplingFlag(int flag);
    
    void setPartialCoupling(int flag, std::string type);
    
    int getPartialCouplingFlag(int i);
    
    std::string getPartialCouplingType(int i);
    
    int sizePartialCoupling();
    
    void partitionMeshInterface(MapPtr_Type mapRepeated, MapPtr_Type mapUnique);
    
    void buildFromOtherInterface( MeshInterfacePtr_Type otherMeshInterface);

    void print(CommConstPtr_Type comm);

    vec3D_GO_ptr_Type getIndicesGlobalMatched();

    vec3D_GO_ptr_Type getIndicesGlobalMatchedOrigin();

    vec3D_GO_ptr_Type getIndicesGlobalMatchedUnique();

    private:
    // Am Ende repeated partitioniert!
    // TODO: in indicesGlobalMatchedRepated_ umbenennen.
    vec3D_GO_ptr_Type indicesGlobalMatched_; /* One outer vector entry for each interface flag that is used.
                                              If Interface is determined for 2 flags than  2 = indicesGlobalMatched_.size().
                                             After this, 2 row vectors with global IDs for this and otherMesh.
                                            2 = indicesGlobalMatched_.at(0).size(),
                                            NumberOfInterfaceNodesWithFlag = indicesGlobalMatched_.at(0).at(0).size() */

    vec3D_GO_ptr_Type indicesGlobalMatchedOrigin_; // nicht patitioniert, im Gegensatz zu indicesGlobalMatched_

    vec3D_GO_ptr_Type indicesGlobalMatchedUnique_; // im Gegensatz zu indicesGlobalMatched_ unique partitioniert.

    bool isPartitioned_;

    vec_int_Type partialCouplingFlag_;
    vec_string_Type partialCouplingType_;
    
    CommConstPtr_Type comm_;
    
    };

}
#endif
