#ifndef BlockMultiVector_DECL_hpp
#define BlockMultiVector_DECL_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "BlockMap.hpp"
#include "MultiVector.hpp"
#include <Thyra_ProductVectorSpaceBase.hpp>
#include <Thyra_DefaultProductMultiVector_decl.hpp>
/*!
 Declaration of BlockMultiVector
 
 @brief  BlockMultiVector
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC, class LO, class GO, class NO>
class MultiVector;
template <class LO, class GO, class NO>
class BlockMap;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class BlockMultiVector {

public:

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
    
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    typedef typename MultiVector_Type::Map_Type Map_Type;
    typedef typename MultiVector_Type::MapPtr_Type MapPtr_Type;
    typedef typename MultiVector_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef typename MultiVector_Type::Comm_Type Comm_Type;
    typedef typename MultiVector_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef BlockMap<LO,GO,NO> BlockMap_Type;
    typedef Teuchos::RCP<BlockMap_Type> BlockMapPtr_Type;
    typedef Teuchos::RCP<const BlockMap_Type> BlockMapConstPtr_Type;

    BlockMultiVector();

    BlockMultiVector(UN size);

    BlockMultiVector( BlockMultiVectorPtr_Type bMVIn);

    BlockMultiVector( std::vector<MapConstPtr_Type>& maps, int numMV = 1 );

    BlockMultiVector( BlockMapConstPtr_Type blockMap, int numMV = 1 );
    
    ~BlockMultiVector();

    BlockMultiVector_Type& operator= (const BlockMultiVector_Type& rhs) {
        TEUCHOS_TEST_FOR_EXCEPTION( size() != rhs.size(), std::logic_error,"BlockMultiVector sizes are not equal for deep copy.")
        for (int i=0; i<size(); i++) {
            *blockMultiVector_[i] = *rhs[i];
        }
        return *this;
    }

    MultiVectorPtr_Type operator[] (int i) const{
        TEUCHOS_TEST_FOR_EXCEPTION( i > size()-1, std::logic_error,"The requested MultiVector does not exist in BlockMultiVector.")
        return blockMultiVector_[i];
    }
    
    void resize( UN size );
    
    void buildFromMaps( std::vector<MapConstPtr_Type>& maps, int numMV );

    void buildFromBlockMap( BlockMapConstPtr_Type blockMap, int numMV );
    
    void merge();

    void mergeBlock(UN block);

    void split();

    void splitBlock(UN block);

    void determineLocalOffsets();

    void determineGlobalOffsets();

    void setMergedVector( MultiVectorPtr_Type& mv );

    int size() const;

    UN getNumVectors() const;

    void addBlock(const MultiVectorPtr_Type& multiVector, int i);

    MultiVectorConstPtr_Type getBlock(int i) const;
    
    MultiVectorPtr_Type getBlockNonConst(int i);

    Teuchos::RCP< Thyra::MultiVectorBase<SC> > getThyraMultiVector( );

    Teuchos::RCP<const Thyra::MultiVectorBase<SC> > getThyraMultiVectorConst( );

    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > getProdThyraMultiVector( );

    Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > getThyraProdMultiVectorConst( ) const;

    void fromThyraMultiVector( Teuchos::RCP< Thyra::MultiVectorBase<SC> > thyraMV);
    
    void fromThyraProdMultiVector( Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > thyraMV);
    
    void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms) const;

    void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms) const;

    void dot(BlockMultiVectorConstPtr_Type a, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &dots) ;

    void update( const SC& alpha, const BlockMultiVector_Type& A, const SC& beta );

    void update( const SC& alpha, const BlockMultiVector_Type& A, const SC& beta , const BlockMultiVector_Type& B, const SC& gamma);

    BlockMultiVectorPtr_Type sumColumns() const;
    
    // Matrix-matrix multiplication: this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const SC &alpha, BlockMultiVectorConstPtr_Type &A, BlockMultiVectorConstPtr_Type &B, const SC &beta);
    
    void putScalar( const SC& alpha );

    void scale( const SC& alpha );

    void writeMM(std::string fN = "blockMV") const;

    void print(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME);

    BlockMapConstPtr_Type getMap() const;

    MultiVectorConstPtr_Type getMergedVector();

    /// @brief Return the merged vector as non-const vector. This finds application in pressure projection
    /// @return 
    MultiVectorPtr_Type getMergedVectorNonConst();


private:

    Teuchos::Array<MultiVectorPtr_Type> blockMultiVector_;
    BlockMapPtr_Type blockMap_;
    MultiVectorPtr_Type mergedMultiVector_;
    MapPtr_Type mergedMap_;
    Teuchos::ArrayRCP<LO> localBlockOffsets_;
    Teuchos::ArrayRCP<GO> globalBlockOffsets_;
};
}

#endif
