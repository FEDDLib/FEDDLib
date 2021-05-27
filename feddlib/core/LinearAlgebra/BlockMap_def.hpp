#ifndef BlockMap_DEF_hpp
#define BlockMap_DEF_hpp
#include "BlockMap_decl.hpp"
/*!
 Definition of BlockMap
 
 @brief  BlockMap
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class LO, class GO, class NO>
BlockMap<LO,GO,NO>::BlockMap():
blockMap_(0),
mergedMap_()
{

}

template < class LO, class GO, class NO>
BlockMap<LO,GO,NO>::BlockMap(UN size):
blockMap_(size),
mergedMap_()
{
    
}

template <class LO, class GO, class NO>
BlockMap<LO,GO,NO>::~BlockMap(){
    
}

template <class LO, class GO, class NO>
void BlockMap<LO,GO,NO>::resize(UN size) {
    blockMap_.resize( size );
    mergedMap_.reset();
}

template <class LO, class GO, class NO>
void BlockMap<LO,GO,NO>::addBlock(MapConstPtr_Type map, int i){
    TEUCHOS_TEST_FOR_EXCEPTION( map.is_null(), std::runtime_error,"Map which you want to add to BlockMap is null.");
    if (i>blockMap_.size()-1)
        blockMap_.resize( blockMap_.size()+1 );
    if (!blockMap_[i].is_null())
        blockMap_[i].reset();

    MapPtr_Type mapNonConst =  Teuchos::rcp_const_cast<Map_Type>(map);
    blockMap_[i] = mapNonConst;
}

    
template <class LO, class GO, class NO>
void BlockMap<LO,GO,NO>::merge( ){
    if ( mergedMap_.is_null() ) {        
        TEUCHOS_TEST_FOR_EXCEPTION( blockMap_.size()==0, std::logic_error,"BlockMap has no maps - we cannot merge.");
        typedef Teuchos::ScalarTraits<GO> GOST;
        GO globalOffset = GOST::zero();
        LO offset = Teuchos::ScalarTraits<LO>::zero();
        Teuchos::Array<GO> globalElementList(0);
        for (UN i=0; i<blockMap_.size(); i++) {
            TEUCHOS_TEST_FOR_EXCEPTION( blockMap_[i].is_null(), std::runtime_error,"Map in BlockMap is null. This should not happen.");
            Teuchos::ArrayView<const GO> blockGlobElementList = blockMap_[i]->getNodeElementList();
            globalElementList.resize( globalElementList.size() + blockGlobElementList.size() );
            for (UN j=0; j<blockGlobElementList.size(); j++) {
                globalElementList[offset] = blockGlobElementList[j] + globalOffset;
                offset++;
            }
            globalOffset += blockMap_[i]->getMaxAllGlobalIndex() + 1;
        }
        
        CommConstPtr_Type comm = blockMap_[0]->getComm();
        typedef Teuchos::OrdinalTraits<GO> GOOT;

        mergedMap_ = Teuchos::rcp( new Map_Type( blockMap_[0]->getUnderlyingLib(), GOOT::invalid(), globalElementList(), GOST::zero(), comm ) );
    }
}

template <class LO, class GO, class NO>
std::string BlockMap<LO,GO,NO>::getUnderlyingLib( ) const{
    TEUCHOS_TEST_FOR_EXCEPTION(blockMap_.size()==0,std::runtime_error,"BlockMap size is 0, there is no underlying Lib.");
    TEUCHOS_TEST_FOR_EXCEPTION(blockMap_[0].is_null(),std::runtime_error,"BlockMap[0] is null.");
    return blockMap_[0]->getUnderlyingLib();
}

template <class LO, class GO, class NO>
typename BlockMap<LO,GO,NO>::MapConstPtr_Type BlockMap<LO,GO,NO>::getMergedMap() {
    if ( mergedMap_.is_null() )
        this->merge();
    return mergedMap_;
}

template <class LO, class GO, class NO>
typename BlockMap<LO,GO,NO>::CommConstPtr_Type BlockMap<LO,GO,NO>::getComm(){
    TEUCHOS_TEST_FOR_EXCEPTION( blockMap_.size()==0, std::logic_error,"BlockMap has no maps - no Comm available.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMap_[0].is_null(), std::logic_error,"BlockMap[0] is null.");
    return blockMap_[0]->getComm();
}

template <class LO, class GO, class NO>
typename BlockMap<LO,GO,NO>::CommPtr_Type BlockMap<LO,GO,NO>::getCommNonConst(){
    TEUCHOS_TEST_FOR_EXCEPTION( blockMap_.size()==0, std::logic_error,"BlockMap has no maps - no Comm available.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMap_[0].is_null(), std::logic_error,"BlockMap[0] is null.");
    return blockMap_[0]->getCommNonConst();
}
    
template <class LO, class GO, class NO>
typename BlockMap<LO,GO,NO>::MapPtr_Type BlockMap<LO,GO,NO>::getBlock(UN i){
    TEUCHOS_TEST_FOR_EXCEPTION( i > blockMap_.size()-1, std::logic_error,"BlockMap entry does not exist.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMap_[i].is_null(), std::logic_error,"Map in BlockMap entry is null.");
    return blockMap_[i];
}

template <class LO, class GO, class NO>
typename BlockMap<LO,GO,NO>::MapConstPtr_Type BlockMap<LO,GO,NO>::getBlock(UN i) const{
    TEUCHOS_TEST_FOR_EXCEPTION( i > blockMap_.size()-1, std::logic_error,"BlockMap entry does not exist.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMap_[i].is_null(), std::logic_error,"Map in BlockMap entry is null.");
    return blockMap_[i];
}


}
#endif
