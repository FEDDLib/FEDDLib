#ifndef BLOCKMATRIX_DECL_hpp
#define BLOCKMATRIX_DECL_hpp

#include <Teuchos_Tuple.hpp>
#include <Thyra_BlockedLinearOpBase.hpp>
#include <Thyra_DefaultBlockedLinearOp_decl.hpp>

#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/SmallMatrix.hpp"
#include "BlockMap.hpp"
#include "BlockMultiVector.hpp"
#include "Matrix.hpp"

/*!
 Declaration of BlockMatrix

 @brief  BlockMatrix
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class BlockMatrix {

public:

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
    typedef Teuchos::RCP<const Matrix_Type> MatrixConstPtr_Type;

    typedef typename Matrix_Type::Map_Type Map_Type;
    typedef typename Matrix_Type::MapPtr_Type MapPtr_Type;
    typedef typename Matrix_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef typename Matrix_Type::MultiVector_Type MultiVector_Type;
    typedef typename Matrix_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Matrix_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    typedef BlockMap<LO,GO,NO> BlockMap_Type;
    typedef Teuchos::RCP<BlockMap_Type> BlockMapPtr_Type;

    typedef Teuchos::Tuple<GO, 2> GOTuple_Type;
    typedef Teuchos::Tuple<LO, 2> LOTuple_Type;

    BlockMatrix();

    BlockMatrix(UN size);

    // TODO
    BlockMatrix(BlockMatrixPtr_Type bMatrixIn);

    ~BlockMatrix();

    int size() const;

    void resize(UN size);
    
    MatrixPtr_Type getBlock(int i, int j);

    MatrixConstPtr_Type getBlockConst(int i, int j) const;
    
    bool blockExists(int i, int j) const;

    void addBlock(const MatrixPtr_Type& matrix, int i, int j);

    void merge();

    void mergeBlockNew(UN blockRow, UN blockCol);

    void determineLocalOffsets();

    void determineGlobalOffsets();

    void print(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME);

    void printMerge(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME);

    void writeMM(std::string fN="blockMat") const;

    Teuchos::RCP<const Thyra::LinearOpBase<SC> > getThyraLinOp();

    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC> > getThyraLinBlockOp() const;

    void apply(const BlockMultiVector_Type& X,
               BlockMultiVector_Type& Y ) const;

    void apply(const BlockMultiVector_Type& X,
               BlockMultiVector_Type& Y,
               const SmallMatrix<SC>& coeff) const;

    void addMatrix( const SmallMatrix<SC>& coeffAlpha, const BlockMatrixPtr_Type& matrix, const SmallMatrix<SC>& coeffbeta );

    MatrixPtr_Type getMergedMatrix();

    BlockMapPtr_Type getMap();
protected:

    SmallMatrix<MatrixPtr_Type> blockMatrix_;
    BlockMapPtr_Type blockMap_;
    MatrixPtr_Type mergedMatrix_;
    MapPtr_Type mergedMap_;
    Teuchos::RCP<SmallMatrix<GOTuple_Type> > globalBlockOffsets_;
    Teuchos::RCP<SmallMatrix<LOTuple_Type> > localBlockOffsets_;
};

}
#endif
