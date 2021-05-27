#ifndef BlockMap_DECL_hpp
#define BlockMap_DECL_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "MultiVector.hpp"

/*!
 Declaration of BlockMap
 
 @brief  BlockMap
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class LO = default_lo, class GO = default_go, class NO = default_no>
class BlockMap {
    
public:

    typedef BlockMap<LO,GO,NO> BlockMap_Type;
    typedef Teuchos::RCP<BlockMap_Type> BlockMapPtr_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;
    typedef typename Map_Type::MapConstPtr_Type MapConstPtr_Type;
    
    typedef typename Map_Type::Comm_Type Comm_Type;
    typedef typename Map_Type::CommPtr_Type CommPtr_Type;
    typedef typename Map_Type::CommConstPtr_Type CommConstPtr_Type;

    BlockMap();

    BlockMap(UN size);
    
    ~BlockMap();
    
    void resize(UN size);
    
    void addBlock(MapConstPtr_Type map, int i);
    
    void merge();
    
    std::string getUnderlyingLib( ) const;
        
    MapConstPtr_Type getMergedMap();

    CommConstPtr_Type getComm();

    CommPtr_Type getCommNonConst();
    
    MapPtr_Type getBlock(UN i);
    
    MapConstPtr_Type getBlock(UN i) const;
    
    UN size() const { return blockMap_.size(); };
private:
    
    Teuchos::Array<MapPtr_Type> blockMap_;
    MapPtr_Type mergedMap_;
};
}

#endif
