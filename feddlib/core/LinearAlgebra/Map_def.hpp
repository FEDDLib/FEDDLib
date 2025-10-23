#ifndef MAP_DEF_hpp
#define MAP_DEF_hpp
#include "Map_decl.hpp"

/*!
 Declaration of Map
 
    June 2024: Transition from XPetra to TPetra. In most instances XPetra interface works through a Factory + Build-call. We replaced this with a direct call 
    to the TPetra functions which can be seen in the corresponding XPetra::Factories. 
 */

namespace FEDD {
template < class LO, class GO, class NO>
Map<LO,GO,NO>::Map():
map_()
{

}

template <class LO, class GO, class NO>
Map<LO,GO,NO>::Map( TpetraMapConstPtrConst_Type& tpetraMapPtrIn ):
map_( tpetraMapPtrIn )
{
    
}

template < class LO, class GO, class NO>
Map<LO,GO,NO>::Map( const Map_Type& mapIn ):
map_()
{
    map_ = Teuchos::RCP(new TpetraMap_Type(mapIn.getGlobalNumElements(), mapIn.getNodeElementList(),mapIn.getIndexBase(), mapIn.getComm()));

}
    
template < class LO, class GO, class NO>
Map<LO,GO,NO>::Map(GO numGlobalElements,
                    const Teuchos::ArrayView<const GO> &elementList,
                    GO indexBase,
                    const CommConstPtr_Type &comm):
map_()
{
    map_ =  Teuchos::RCP(new TpetraMap_Type(numGlobalElements, elementList, indexBase, comm));
}

template < class LO, class GO, class NO>
Map<LO,GO,NO>::Map(GO numGlobalElements,
                   LO numLocalElements,
                   GO indexBase,
                   const CommConstPtr_Type &comm):
map_()
{
    map_ = Teuchos::RCP(new TpetraMap_Type(numGlobalElements, numLocalElements, indexBase, comm));

}
    
template < class LO, class GO, class NO>
Map<LO,GO,NO>::~Map(){
    
}

    
template < class LO, class GO, class NO>
GO Map<LO,GO,NO>::getGlobalElement(LO id) const{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    return map_->getGlobalElement(id);
}

template < class LO, class GO, class NO>
LO Map<LO,GO,NO>::getLocalElement(GO id) const{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    return map_->getLocalElement(id);
}

template < class LO, class GO, class NO>
LO Map<LO,GO,NO>::getNodeNumElements() const{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    return map_->getLocalNumElements();
}

template < class LO, class GO, class NO>
GO Map<LO,GO,NO>::getGlobalNumElements() const{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    return map_->getGlobalNumElements();
}

template < class LO, class GO, class NO>
typename Map<LO,GO,NO>::CommConstPtr_Type Map<LO,GO,NO>::getComm() const{
    
    return map_->getComm();
}

template <class LO,class GO,class NO>
GO Map<LO,GO,NO>::getIndexBase() const{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    return map_->getIndexBase();
}

template < class LO, class GO, class NO>
Teuchos::ArrayView< const GO > Map<LO,GO,NO>::getNodeElementList() const{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    return map_->getLocalElementList();
}

template < class LO, class GO, class NO>
typename Map<LO,GO,NO>::MapPtr_Type Map<LO,GO,NO>::buildVecFieldMap(UN numDofs, std::string ordering) const{

    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    TEUCHOS_TEST_FOR_EXCEPTION(ordering.compare("NodeWise"), std::logic_error,"Select a valid ordering: NodeWise");
    Teuchos::ArrayView<const GO> elementList = map_->getLocalElementList();
    Teuchos::Array<GO> elementListField( numDofs *  elementList.size() );
    for (UN i=0; i<elementList.size(); i++) {
        for (UN dof=0; dof<numDofs; dof++)
            elementListField[ numDofs * i + dof ] = numDofs * elementList[i] + dof;
    }
    typedef Teuchos::OrdinalTraits<GO> GOOT;
    return Teuchos::rcp(new Map_Type(GOOT::invalid(), elementListField(), map_->getIndexBase(), map_->getComm() ) );
}


template < class LO, class GO, class NO>
typename Map<LO,GO,NO>::TpetraMapConstPtr_Type Map<LO,GO,NO>::getTpetraMap() const{
    
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"getTpetraMap(): map_ is null.");
    
    return map_;
}

template < class LO, class GO, class NO>
GO Map<LO,GO,NO>::getMaxAllGlobalIndex() const{
    return map_->getMaxAllGlobalIndex();
}
    
template < class LO, class GO, class NO>
LO Map<LO,GO,NO>::getMaxLocalIndex() const{
    return map_->getMaxLocalIndex();
}
    
template < class LO, class GO, class NO>
void Map<LO,GO,NO>::print(Teuchos::EVerbosityLevel verbLevel) const{

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    map_->describe(*out,verbLevel);
}

template < class LO, class GO, class NO>
typename Map<LO,GO,NO>::ThyraVSBConstPtr_Type Map<LO,GO,NO>::getThyraVectorSpaceBase() const{
    
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"getThyraVectorSpaceBase(): map_ is null.");
    //RCP< const Tpetra::Map< LocalOrdinal, GlobalOrdinal, Node > > tpMap = tpetraMap->getTpetra_Map();
    Teuchos::RCP<Thyra::TpetraVectorSpace<default_sc,LO,GO,NO> > thyraTpetraMap = Thyra::tpetraVectorSpace<default_sc,LO,GO,NO>(map_);
    return thyraTpetraMap;
}

template < class LO, class GO, class NO>
typename Map<LO,GO,NO>::CommPtr_Type Map<LO,GO,NO>::getCommNonConst() {
    CommConstPtr_Type commConst = map_->getComm();
    return Teuchos::rcp_const_cast<Comm_Type>(commConst);
}

template <class LO,class GO,class NO>
Teuchos::RCP<Map<LO,GO,NO> > Map<LO,GO,NO>::buildUniqueMap( int numFreeProcs ) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    if (numFreeProcs==0) {
        
        Teuchos::RCP<Tpetra::Vector<GO,LO,GO,NO> > myIndices = Teuchos::RCP( new Tpetra::Vector<GO,LO,GO,NO>(map_));
        myIndices->putScalar(map_->getComm()->getRank()+1);

        TpetraMapPtr_Type linearMap = Teuchos::RCP(new TpetraMap_Type( map_->getMaxAllGlobalIndex()+1, 0, map_->getComm()));
        Teuchos::RCP<Tpetra::Vector<GO,LO,GO,NO> > globalIndices = Teuchos::RCP( new Tpetra::Vector<GO,LO,GO,NO>(linearMap));

        Teuchos::RCP<Tpetra::Import<LO,GO,NO> > importer = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>( map_,linearMap));
        Teuchos::RCP<Tpetra::Import<LO,GO,NO> > importer2 = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>(linearMap, map_));
        globalIndices->doImport(*myIndices,*importer,Tpetra::INSERT);
        myIndices->putScalar(0);
        myIndices->doImport(*globalIndices,*importer2,Tpetra::ADD);
        
        Teuchos::Array<GO> uniqueVector;
        for (unsigned i=0; i<myIndices->getLocalLength(); i++) {
            if (myIndices->getData(0)[i] == map_->getComm()->getRank()+1) {
                uniqueVector.push_back(map_->getGlobalElement(i));
            }
        }
        TpetraMapPtr_Type mapTpetra = Teuchos::RCP(new TpetraMap_Type(-1, uniqueVector(),0,map_->getComm()));
        Teuchos::RCP<Map<LO,GO,NO> > map = Teuchos::rcp( new Map<LO,GO,NO>( mapTpetra ) );
        return  map;
    }
    else{
        Teuchos::RCP<Tpetra::Vector<GO,LO,GO,NO> > myIndices = Teuchos::RCP( new Tpetra::Vector<GO,LO,GO,NO>(map_));
        myIndices->putScalar(map_->getComm()->getRank()+1);
        GO maxGID = map_->getMaxAllGlobalIndex();
        int numAvailableRanks = map_->getComm()->getSize() - numFreeProcs;
        int numElementsForAvailRank = (int) ( ( (maxGID+1) / numAvailableRanks ) + 100.*std::numeric_limits<double>::epsilon() );
        
        int remainingElement = maxGID+1 - numAvailableRanks * numElementsForAvailRank;
        bool hasOneMoreElement = false;
        
        if ( remainingElement > map_->getComm()->getRank() ) {
            numElementsForAvailRank++;
            hasOneMoreElement = true;
        }
        
        if ( map_->getComm()->getRank() + 1 > map_->getComm()->getSize() - numFreeProcs){
            numElementsForAvailRank = 0;
        }
        
        Teuchos::Array<GO> myElements( numElementsForAvailRank );
        GO offset = numElementsForAvailRank * map_->getComm()->getRank();
        if (!hasOneMoreElement) {
            offset += remainingElement;
        }
        for (int i=0; i<myElements.size(); i++) {
            myElements[i] = i + offset;
        }

        TpetraMapPtr_Type linearMapAvailRanks = Teuchos::RCP(new TpetraMap_Type(Teuchos::OrdinalTraits<GO>::invalid(), myElements(), 0, map_->getComm()));
        Teuchos::RCP<Tpetra::Vector<GO,LO,GO,NO> > globalIndices = Teuchos::RCP( new Tpetra::Vector<GO,LO,GO,NO>(linearMapAvailRanks));

        Teuchos::RCP<Tpetra::Import<LO,GO,NO> > importer = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>(map_, linearMapAvailRanks ));
        Teuchos::RCP<Tpetra::Import<LO,GO,NO> > importer2 = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>( linearMapAvailRanks, map_));
        
        globalIndices->doImport(*myIndices,*importer,Tpetra::INSERT);
        
        myIndices->putScalar(0);
        myIndices->doImport(*globalIndices,*importer2,Tpetra::ADD);
        
        Teuchos::Array<GO> uniqueVector;
        for (unsigned i=0; i<myIndices->getLocalLength(); i++) {
            if (myIndices->getData(0)[i] == map_->getComm()->getRank()+1) {
                uniqueVector.push_back(map_->getGlobalElement(i));
            }
        }
        TpetraMapPtr_Type mapTpetra = Teuchos::RCP(new TpetraMap_Type(-1, uniqueVector(),0,map_->getComm()));
        Teuchos::RCP<Map<LO,GO,NO> > map = Teuchos::rcp( new Map<LO,GO,NO>( mapTpetra ) );
        return  map;
    }
    
}

//merge with above function
template <class LO,class GO,class NO>
Teuchos::RCP<Map<LO,GO,NO> > Map<LO,GO,NO>::buildUniqueMap( tuple_intint_Type rankRange ) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(map_.is_null(),std::runtime_error,"map is null.");
    int rank = map_->getComm()->getRank();
    Teuchos::RCP<Tpetra::Vector<GO,LO,GO,NO> > myIndices = Teuchos::RCP( new Tpetra::Vector<GO,LO,GO,NO>(map_));

    myIndices->putScalar(rank + 1);
    GO maxGID = map_->getMaxAllGlobalIndex();

    int numAvailableRanks = std::get<1>(rankRange) - std::get<0>(rankRange) + 1;
    int numElementsForAvailRank = (int) ( ( (maxGID+1) / numAvailableRanks ) + 100.*std::numeric_limits<double>::epsilon() );
    
    int remainingElement = maxGID+1 - numAvailableRanks * numElementsForAvailRank;
    bool hasOneMoreElement = false;
    
    if ( remainingElement > map_->getComm()->getRank() - std::get<0>(rankRange) ) {
        numElementsForAvailRank++;
        hasOneMoreElement = true;
    }
    if ( map_->getComm()->getRank() < std::get<0>(rankRange) || map_->getComm()->getRank() > std::get<1>(rankRange)  )
        numElementsForAvailRank = 0;
    
    
    Teuchos::Array<GO> myElements( numElementsForAvailRank );
    GO offset = numElementsForAvailRank * ( map_->getComm()->getRank() - std::get<0>(rankRange) );
    
    if (!hasOneMoreElement)
        offset += remainingElement;

    for (int i=0; i<myElements.size(); i++)
        myElements[i] = i + offset;
  
    TpetraMapPtr_Type linearMapAvailRanks = Teuchos::RCP(new TpetraMap_Type( Teuchos::OrdinalTraits<GO>::invalid(), myElements(), 0, map_->getComm()));
    Teuchos::RCP<Tpetra::Vector<GO,LO,GO,NO> > globalIndices = Teuchos::RCP( new Tpetra::Vector<GO,LO,GO,NO>(linearMapAvailRanks));

    Teuchos::RCP<Tpetra::Import<LO,GO,NO> > importer = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>( map_,linearMapAvailRanks));
    Teuchos::RCP<Tpetra::Import<LO,GO,NO> > importer2 = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>(linearMapAvailRanks, map_));

    globalIndices->doImport(*myIndices,*importer,Tpetra::INSERT);
    
    myIndices->putScalar(0);
    myIndices->doImport(*globalIndices,*importer2,Tpetra::ADD);
    
    Teuchos::Array<GO> uniqueVector;
    for (unsigned i=0; i<myIndices->getLocalLength(); i++) {
        if (myIndices->getData(0)[i] == map_->getComm()->getRank()+1) {
            uniqueVector.push_back(map_->getGlobalElement(i));
        }
    }

    TpetraMapPtr_Type mapTpetra = Teuchos::RCP(new TpetraMap_Type(-1, uniqueVector(),0,map_->getComm()));
    Teuchos::RCP<Map<LO,GO,NO> > map = Teuchos::rcp( new Map<LO,GO,NO>( mapTpetra ) );
    return  map;

}
   

}
#endif
