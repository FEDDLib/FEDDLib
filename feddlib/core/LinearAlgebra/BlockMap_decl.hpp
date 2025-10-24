#ifndef BlockMap_DECL_hpp
#define BlockMap_DECL_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "Map.hpp"

/*!
 Declaration of BlockMap
 
 @brief  BlockMap
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
    /*!
    \class BlockMap
    \brief Block Variant of Map class. 
    
    This becomes relevant for block systems and when we attempt to use monolithic preconditioning or solving to determine the merged map of system. 

    \tparam LO The local ordinal type. The is the index type for local indices
    \tparam GO The global ordinal type. The is the index type for global indices
    \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

    This becomes relevant for block systems and when we attempt to use monolithic preconditioning or solving to determine the merged map of system. 

    i.e. Stokes Problem:

    0 | A   B       BlockMap->getBlock(0) of row-block 0
    -----------
    1 | B^T 0       BlockMap->getBlock(1) of row-block 1

    */
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

    /*! @brief Initializing block maps with the size of system (i.e. 2 for Stokes problem). As systems are distributed row-wise, block i corresponds to row-block i.
     @param size of system 
    */
    BlockMap(UN size);
    
    ~BlockMap();
    
    void resize(UN size);
    
    /// @brief Adding map to corresponding block i. Block i is row-wise distrubted according to the underlying map.
    /// @param map local to global indexing of rows
    /// @param i block id
    void addBlock(MapConstPtr_Type map, int i);
    
    /// @brief Merging the map of different blocks together. Relevant for monolithic solving/precondtioning
    void merge();

    void print();

    void info();
    
    /// @brief Getting merged map of block maps
    /// @return mergedMap
    MapConstPtr_Type getMergedMap();

    /// @brief Get underlying communicator
    /// @return comm_
    CommConstPtr_Type getComm();

    CommPtr_Type getCommNonConst();
    
    MapPtr_Type getBlock(UN i);
    
    MapConstPtr_Type getBlock(UN i) const;
    
    UN size() const { return blockMap_.size(); }
private:
    
    Teuchos::Array<MapPtr_Type> blockMap_;
    MapPtr_Type mergedMap_;
};
}

#endif
