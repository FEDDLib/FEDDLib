#ifndef BlockMultiVector_DEF_hpp
#define BlockMultiVector_DEF_hpp
#include "BlockMultiVector_decl.hpp"

/*!
 Defintion of BlockMultiVector
 
 @brief  BlockMultiVector
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


namespace FEDD {
template <class SC, class LO, class GO, class NO>
BlockMultiVector<SC,LO,GO,NO>::BlockMultiVector():
blockMultiVector_(0),
blockMap_(),
mergedMultiVector_(),
mergedMap_()
{
    blockMap_ = Teuchos::rcp( new BlockMap_Type( 0 ) );
}

template <class SC, class LO, class GO, class NO>
BlockMultiVector<SC,LO,GO,NO>::BlockMultiVector(UN size):
blockMultiVector_(size),
blockMap_(),
mergedMultiVector_(),
mergedMap_()
{
    blockMap_ = Teuchos::rcp( new BlockMap_Type( size ) );
}

template <class SC, class LO, class GO, class NO>
BlockMultiVector<SC,LO,GO,NO>::BlockMultiVector( BlockMultiVectorPtr_Type bMVIn):
blockMultiVector_(bMVIn->size()),
blockMap_(),
mergedMultiVector_(),
mergedMap_()
{
    blockMap_ = Teuchos::rcp( new BlockMap_Type( bMVIn->size() ) );
    for (UN i=0; i<bMVIn->size(); i++) {
        //copy MV and then insert
        MultiVectorConstPtr_Type mvTmp = bMVIn->getBlock(i);
        MultiVectorPtr_Type mvCopy = Teuchos::rcp(new MultiVector_Type( mvTmp ) );
        this->addBlock( mvCopy, i );
    }
}

template <class SC, class LO, class GO, class NO>
BlockMultiVector<SC,LO,GO,NO>::BlockMultiVector( std::vector<MapConstPtr_Type>& maps, int numMV ):
blockMultiVector_( maps.size() ),
blockMap_(),
mergedMultiVector_(),
mergedMap_()
{
    blockMap_ = Teuchos::rcp( new BlockMap_Type( maps.size() ) );
    buildFromMaps( maps, numMV );

    
}
    
template <class SC, class LO, class GO, class NO>
BlockMultiVector<SC,LO,GO,NO>::BlockMultiVector( BlockMapConstPtr_Type blockMap, int numMV ):
blockMultiVector_( blockMap->size() ),
blockMap_( ),
mergedMultiVector_(),
mergedMap_()
{
    
    blockMap_ = Teuchos::rcp_const_cast<BlockMap_Type>( blockMap );
    buildFromBlockMap( blockMap, numMV );
}

template <class SC, class LO, class GO, class NO>
BlockMultiVector<SC,LO,GO,NO>::~BlockMultiVector(){

}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::resize(UN size){

    blockMultiVector_.resize( size );
    blockMap_->resize(size);
    mergedMultiVector_.reset();
    mergedMap_.reset();
    localBlockOffsets_.resize( size );
    globalBlockOffsets_.resize( size );

}


template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::buildFromMaps( std::vector<MapConstPtr_Type>& maps, int numMV ){
    
    for (int i=0; i<maps.size(); i++){
        blockMap_->addBlock( maps[i], i );
        MultiVectorPtr_Type mv = Teuchos::rcp( new MultiVector_Type( maps[i], numMV ) );
        this->addBlock( mv, i );
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::buildFromBlockMap( BlockMapConstPtr_Type blockMap, int numMV ){
    
    for (int i=0; i<blockMap->size(); i++){
        MultiVectorPtr_Type mv = Teuchos::rcp( new MultiVector_Type( blockMap->getBlock( i ), numMV ) );
        this->addBlock( mv, i );
    }
}
    
template <class SC, class LO, class GO, class NO>
UN BlockMultiVector<SC,LO,GO,NO>::getNumVectors() const{
    TEUCHOS_TEST_FOR_EXCEPTION(blockMultiVector_.size()==0, std::logic_error,"No MultiVector in BlockMultiVector.");
    TEUCHOS_TEST_FOR_EXCEPTION(blockMultiVector_[0].is_null(), std::runtime_error,"MultiVector in BlockMultiVector is null.");
    return blockMultiVector_[0]->getNumVectors();
}

template <class SC, class LO, class GO, class NO>
int BlockMultiVector<SC,LO,GO,NO>::size() const{
    return blockMultiVector_.size();
}

template <class SC, class LO, class GO, class NO>
typename BlockMultiVector<SC,LO,GO,NO>::MultiVectorConstPtr_Type BlockMultiVector<SC,LO,GO,NO>::getBlock(int i) const{

    TEUCHOS_TEST_FOR_EXCEPTION( (blockMultiVector_.size()-1) < i, std::logic_error,"Block in BlockMultiVector does not exist.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_[0].is_null(), std::runtime_error, "Block in BlockMultiVector is null.");
    return blockMultiVector_[i];
}

template <class SC, class LO, class GO, class NO>
typename BlockMultiVector<SC,LO,GO,NO>::MultiVectorPtr_Type BlockMultiVector<SC,LO,GO,NO>::getBlockNonConst(int i){

    TEUCHOS_TEST_FOR_EXCEPTION( (blockMultiVector_.size()-1) < i, std::logic_error,"Block in BlockMultiVector does not exist.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_[0].is_null(), std::runtime_error, "Block in BlockMultiVector is null.");
    return blockMultiVector_[i];
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::addBlock(const MultiVectorPtr_Type& multiVector, int i){

    TEUCHOS_TEST_FOR_EXCEPTION( multiVector.is_null(), std::runtime_error,"MultiVector which you want to add to BlockMultiVector is null.");
    if ( blockMultiVector_.size()>0 && !blockMultiVector_[0].is_null() )
        TEUCHOS_TEST_FOR_EXCEPTION( multiVector->getNumVectors()!=blockMultiVector_[0]->getNumVectors(), std::logic_error,"MultiVectors for BlockMultiVector have different numbers of vectors.");

    if (i>blockMultiVector_.size()-1)
        blockMultiVector_.resize( blockMultiVector_.size()+1 );
    //Do we need a warning here?
    blockMultiVector_[i] = multiVector;

    blockMap_->addBlock( multiVector->getMap(), i );

}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::merge(){
    if ( mergedMap_.is_null() ) {
        blockMap_->merge();
        mergedMap_ = Teuchos::rcp_const_cast<Map_Type>(blockMap_->getMergedMap());
    }

    
    mergedMultiVector_ = Teuchos::rcp( new MultiVector_Type( mergedMap_, blockMultiVector_[0]->getNumVectors() ) );
    this->determineLocalOffsets();
    this->determineGlobalOffsets();
    
    for (UN i=0; i<blockMultiVector_.size(); i++) {
        if ( !blockMultiVector_[i].is_null() )
            this->mergeBlock( i );
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "MultiVector in BlockMultiVector is null.");
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::mergeBlock(UN block){

    MultiVectorPtr_Type mv = blockMultiVector_[block];
    Teuchos::ArrayRCP<const SC> values;
    GO offset = globalBlockOffsets_[block];
    UN numVectors = mv->getNumVectors();
    for (UN j=0; j<numVectors; j++) {
        values = mv->getData(j);
        for (UN i=0; i<values.size(); i++) {
            mergedMultiVector_->replaceGlobalValue( mv->getMap()->getGlobalElement(i) + offset, j, values[i]);
        }
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::split(){
    TEUCHOS_TEST_FOR_EXCEPTION( mergedMultiVector_.is_null(), std::runtime_error, "MergedMultiVector is null and we cannot use split on it. Use merge() first to generate the MergedMultiVector.");

    for (UN i=0; i<blockMultiVector_.size(); i++) {
        if ( !blockMultiVector_[i].is_null() )
            this->splitBlock( i );
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "MultiVector in BlockMultiVector is null.");
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::splitBlock(UN block){

    MultiVectorPtr_Type mv = blockMultiVector_[block];
    MapConstPtr_Type map = mv->getMap();
    Teuchos::ArrayRCP<const SC> values;
    GO offset = globalBlockOffsets_[block];
    UN numVectors = mv->getNumVectors();
    UN numElements = (UN) map->getNodeNumElements();
    for (UN j=0; j<numVectors; j++) {
        values = mergedMultiVector_->getData(j);
        for (UN i=0; i<numElements; i++) {
            GO globalIndex =  map->getGlobalElement(i) + offset;
            LO index = mergedMultiVector_->getMap()->getLocalElement( globalIndex );
            mv->replaceGlobalValue( mv->getMap()->getGlobalElement(i), j, values[index]);
        }
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::determineLocalOffsets(){

    typedef Teuchos::ScalarTraits<LO> LOST;
    localBlockOffsets_ = Teuchos::ArrayRCP<LO>( blockMultiVector_.size(), LOST::zero() );

    for (UN i=1; i<blockMultiVector_.size() ; i++)
        localBlockOffsets_[i] = localBlockOffsets_[i-1] + blockMultiVector_[i-1]->getMap()->getMaxLocalIndex() + 1;

}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::determineGlobalOffsets(){

    typedef Teuchos::ScalarTraits<GO> GOST;
    globalBlockOffsets_ = Teuchos::ArrayRCP<GO>( blockMultiVector_.size(), GOST::zero() );

    for (UN i=1; i<blockMultiVector_.size() ; i++)
        globalBlockOffsets_[i] = globalBlockOffsets_[i-1] + blockMultiVector_[i-1]->getMap()->getMaxAllGlobalIndex() + 1;

}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::setMergedVector( MultiVectorPtr_Type& mv ){

    mergedMultiVector_ = mv;
    mergedMap_ = mv->getMapNonConst();
    this->determineLocalOffsets();
    this->determineGlobalOffsets();

}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP< Thyra::MultiVectorBase<SC> > BlockMultiVector<SC,LO,GO,NO>::getThyraMultiVector( ) {
    TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_.size() == 0, std::logic_error,"BlockMultiVector size is 0.");
    Teuchos::RCP<Thyra::MultiVectorBase<SC> > thyraMV;
    if (blockMultiVector_.size() == 1){
        TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_[0].is_null(), std::runtime_error, "Block in BlockMultiVector is null.");
        thyraMV = blockMultiVector_[0]->getThyraMultiVector( );
    }
    else {
        this->merge();
        thyraMV = mergedMultiVector_->getThyraMultiVector( );
    }
    return thyraMV;
}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::MultiVectorBase<SC> > BlockMultiVector<SC,LO,GO,NO>::getThyraMultiVectorConst( ) {
    TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_.size() == 0, std::logic_error,"BlockMultiVector size is 0.");
    Teuchos::RCP<const Thyra::MultiVectorBase<SC> > thyraMV;
    if (blockMultiVector_.size() == 1){
        TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_[0].is_null(), std::runtime_error, "Block in BlockMultiVector is null.");
        thyraMV = blockMultiVector_[0]->getThyraMultiVectorConst( );
    }
    else{
        this->merge();
        thyraMV = mergedMultiVector_->getThyraMultiVector( );
    }
    return thyraMV;
}
    
template <class SC, class LO, class GO, class NO>
Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > BlockMultiVector<SC,LO,GO,NO>::getProdThyraMultiVector( ) {
    
    TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_.size() == 0, std::logic_error,"BlockMultiVector size is 0.");
    
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpaces( blockMultiVector_.size() );
    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > multiVecs( blockMultiVector_.size() );
    for (int i=0; i<blockMultiVector_.size(); i++){
        multiVecs[i] = blockMultiVector_[i]->getThyraMultiVector();
        vectorSpaces[i] = multiVecs[i]->range();
    }

    const Teuchos::RCP< Thyra::DefaultProductVectorSpace<SC> > productSpace = Thyra::productVectorSpace<SC>( vectorSpaces );

    return Thyra::defaultProductMultiVector<SC>( productSpace, multiVecs );
}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > BlockMultiVector<SC,LO,GO,NO>::getThyraProdMultiVectorConst( ) const{
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,"We need to implement this.");
    Teuchos::RCP< const Thyra::ProductMultiVectorBase<SC> > thyraProdMV;
    return thyraProdMV;
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::fromThyraMultiVector( Teuchos::RCP< Thyra::MultiVectorBase<SC> > thyraMV) {
    TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_.size() == 0, std::logic_error,"BlockMultiVector size is 0.");
    if (blockMultiVector_.size() == 1){
        TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_[0].is_null(), std::runtime_error, "Block in BlockMultiVector is null.");
        blockMultiVector_[0]->fromThyraMultiVector( thyraMV );
    }
    else {
        if (mergedMultiVector_.is_null())
            this->merge();
        mergedMultiVector_->fromThyraMultiVector( thyraMV );
        this->split();
    }
}

    
template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::fromThyraProdMultiVector( Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraMV) {
    TEUCHOS_TEST_FOR_EXCEPTION( blockMultiVector_.size() == 0, std::logic_error,"BlockMultiVector size is 0.");
    for (int i=0; i<this->size(); i++) {
        this->getBlockNonConst(i)->fromThyraMultiVector( thyraMV->getNonconstMultiVectorBlock( i ) );
    }
}

    
template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms) const {
    typedef Teuchos::ScalarTraits<SC> ST;
    for (int j=0; j<norms.size(); j++)
        norms[j] = ST::zero();

    for (int i=0; i<size(); i++) {
        Teuchos::Array<SC> partialNorm(norms.size());
        blockMultiVector_[i]->norm2(partialNorm());
        for (int j=0; j<partialNorm.size(); j++)
            norms[j] += partialNorm[j]*partialNorm[j];
    }
    for (int j=0; j<norms.size(); j++)
        norms[j] = ST::squareroot( norms[j] );
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::dot(BlockMultiVectorConstPtr_Type a, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &dots)  {

    for (int j=0; j<dots.size(); j++)
        dots[j] = Teuchos::ScalarTraits<SC>::zero();

    Teuchos::Array<SC> partialDots(dots.size());
    for (int i=0; i<size(); i++) {
        blockMultiVector_[i]->dot( a->getBlock(i), partialDots() );
        for (int j=0; j<partialDots.size(); j++)
            dots[j] += partialDots[j];
    }
}
    
template <class SC, class LO, class GO, class NO>
typename BlockMultiVector<SC,LO,GO,NO>::BlockMultiVectorPtr_Type BlockMultiVector<SC,LO,GO,NO>::sumColumns() const{
    
    BlockMultiVectorPtr_Type sumBlockMV = Teuchos::rcp( new BlockMultiVector_Type ( blockMap_, 1 ) );

    for (int i=0; i<this->size(); i++) {
        MultiVectorConstPtr_Type tmpMV = this->getBlock(i);
        MultiVectorPtr_Type sumMV = tmpMV->sumColumns();
        sumBlockMV->addBlock( sumMV, i );
    }
    return sumBlockMV;
}
    
template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::update( const SC& alpha, const BlockMultiVector_Type& A, const SC& beta) {
    //this = alpha*A + beta*this
    TEUCHOS_TEST_FOR_EXCEPTION( size() != A.size(), std::logic_error,"BlockMultiVector sizes are not equal for update.");

    for (int i=0; i<size(); i++) {
        MultiVectorConstPtr_Type tmpA = A[i];
        blockMultiVector_[i]->update( alpha, *tmpA, beta );
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const SC &alpha, BlockMultiVectorConstPtr_Type &A, BlockMultiVectorConstPtr_Type &B, const SC &beta){
//    if (this->getMap()->getCommNonConst()->getRank()==0)
//        std::cout << "### For testing purposes only." << std::endl;

    for (int i=0; i<this->size(); i++){
        MultiVectorConstPtr_Type a = A->getBlock(i);
        MultiVectorConstPtr_Type b = B->getBlock(i);
        this->getBlockNonConst(i)->multiply( transA, transB, alpha, a, b, beta );
    }
}
    
template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::update( const SC& alpha, const BlockMultiVector_Type& A, const SC& beta , const BlockMultiVector_Type& B, const SC& gamma) {
    //this = alpha*A + beta*B + gamma*this
    TEUCHOS_TEST_FOR_EXCEPTION( size() != A.size(), std::logic_error,"BlockMultiVector sizes are not equal for update.");

    for (int i=0; i<size(); i++) {
        MultiVectorConstPtr_Type tmpA = A[i];
        MultiVectorConstPtr_Type tmpB = B[i];
        blockMultiVector_[i]->update( alpha, *tmpA, beta, *tmpB, gamma);
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::putScalar( const SC& alpha ) {
    for (int i=0; i<size(); i++) {
        blockMultiVector_[i]->putScalar( alpha );
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::scale( const SC& alpha ) {
    for (int i=0; i<size(); i++) {
        blockMultiVector_[i]->scale( alpha );
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::writeMM(std::string fN) const{
    for (int i=0; i<size(); i++) {
        if (!blockMultiVector_[i].is_null() ){
            std::string fileName = fN + to_string(i) + ".mm" ;
            blockMultiVector_[i]->writeMM(fileName);
        }
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMultiVector<SC,LO,GO,NO>::print(Teuchos::EVerbosityLevel verbLevel){

    for (UN i=0; i<blockMultiVector_.size(); i++) {
        if ( !blockMultiVector_[i].is_null() ) {
            blockMultiVector_[i]->print( verbLevel );
        }
    }
}

template <class SC, class LO, class GO, class NO>
typename BlockMultiVector<SC,LO,GO,NO>::BlockMapConstPtr_Type BlockMultiVector<SC,LO,GO,NO>::getMap() const{

    return blockMap_;
}

template <class SC, class LO, class GO, class NO>
typename BlockMultiVector<SC,LO,GO,NO>::MultiVectorConstPtr_Type BlockMultiVector<SC,LO,GO,NO>::getMergedVector(){
    if (this->size()>1) {
        this->merge();
        return mergedMultiVector_;
    }
    else
        return this->getBlockNonConst(0);
}
}
#endif
