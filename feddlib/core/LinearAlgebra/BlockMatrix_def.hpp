#ifndef BLOCKMATRIX_DEF_hpp
#define BLOCKMATRIX_DEF_hpp
#include "BlockMatrix_decl.hpp"
/*!
 Definition of BlockMatrix

 @brief  BlockMatrix
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

template <class SC, class LO, class GO, class NO>
BlockMatrix<SC,LO,GO,NO>::BlockMatrix():
blockMatrix_(0),
blockMap_(),
mergedMatrix_(),
mergedMap_(),
globalBlockOffsets_(),
localBlockOffsets_()
{
    blockMap_ = Teuchos::rcp( new BlockMap_Type( 0 ) );
}

template <class SC, class LO, class GO, class NO>
BlockMatrix<SC,LO,GO,NO>::BlockMatrix(UN size):
blockMatrix_(size),
blockMap_(),
mergedMatrix_(),
mergedMap_(),
globalBlockOffsets_(),
localBlockOffsets_()
{
    blockMap_ = Teuchos::rcp( new BlockMap_Type( size ) );
}


template <class SC, class LO, class GO, class NO>
BlockMatrix<SC,LO,GO,NO>::BlockMatrix(BlockMatrixPtr_Type bMatrixIn):
blockMatrix_(bMatrixIn->size()),
blockMap_(),
mergedMatrix_(),
mergedMap_(),
globalBlockOffsets_(),
localBlockOffsets_()
{
    blockMap_ = Teuchos::rcp( new BlockMap_Type( bMatrixIn->size() ) );

    for(UN i = 0; i < bMatrixIn->size(); i++)
    {
        for(UN j = 0; j < bMatrixIn->size(); j++)
        {
            // copy matrix and then insert
            if(bMatrixIn->blockExists(i,j))
            {
                MatrixPtr_Type matrixTmp = bMatrixIn->getBlock(i,j);
                MatrixPtr_Type matrixCopy = Teuchos::rcp(new Matrix_Type( matrixTmp ) );
                this->addBlock(matrixCopy, i, j);
            }
        }
    }
}


template <class SC, class LO, class GO, class NO>
BlockMatrix<SC,LO,GO,NO>::~BlockMatrix()
{

}

template <class SC, class LO, class GO, class NO>
int BlockMatrix<SC,LO,GO,NO>::size() const{
    return blockMatrix_.size();
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::resize(UN size) {
    blockMatrix_.resize( size );
    blockMap_->resize( size );
}

template <class SC, class LO, class GO, class NO>
typename BlockMatrix<SC,LO,GO,NO>::MatrixPtr_Type BlockMatrix<SC,LO,GO,NO>::getBlock(int i, int j){

    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_[i][j].is_null(), std::runtime_error,"Block in BlockMatrix which you tried to access is null.");
    return blockMatrix_[i][j];
}

template <class SC, class LO, class GO, class NO>
typename BlockMatrix<SC,LO,GO,NO>::MatrixConstPtr_Type BlockMatrix<SC,LO,GO,NO>::getBlockConst(int i, int j) const{

    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_[i][j].is_null(), std::runtime_error,"Block in BlockMatrix which you tried to access is null.");
    return blockMatrix_[i][j];
}

template <class SC, class LO, class GO, class NO>
bool BlockMatrix<SC,LO,GO,NO>::blockExists(int i, int j) const{
    return ( !blockMatrix_[i][j].is_null() );
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::addBlock(const MatrixPtr_Type& matrix, int i, int j){
    UN size = blockMatrix_.size();
    if (i>size-1 || j>size-1)
        blockMatrix_.resize( max(i,j)+1 );

    if ( blockExists(i,j) )
        blockMatrix_[i][j].reset();

    blockMatrix_[i][j] = matrix;


    blockMap_->addBlock( matrix->getMap(), i);
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::merge(){
    if ( mergedMap_.is_null() ) {
        blockMap_->merge();
        mergedMap_ = Teuchos::rcp_const_cast<Map_Type>(blockMap_->getMergedMap());
    }
    LO maxNumEntries = 0;
    
    for (UN i=0; i<blockMatrix_.size(); i++) {
        LO maxNumEntriesBlockRow = 0;
        for (UN j=0; j<blockMatrix_.size(); j++) {
            if ( !blockMatrix_[i][j].is_null() )
                maxNumEntriesBlockRow += blockMatrix_[i][j]->getGlobalMaxNumRowEntries();
        }
        if (maxNumEntriesBlockRow > maxNumEntries)
            maxNumEntries = maxNumEntriesBlockRow;
    }
    
//    if ( mergedMatrix_.is_null() ){
    mergedMatrix_ = Teuchos::rcp( new Matrix_Type( mergedMap_, maxNumEntries ) ); // we should identify the number of entries per row.
    this->determineLocalOffsets();
    this->determineGlobalOffsets();

    for (UN i=0; i<blockMatrix_.size(); i++) {
        for (UN j=0; j<blockMatrix_.size(); j++) {
            if ( !blockMatrix_[i][j].is_null() ){
                mergeBlockNew( i, j );
            }
        }
    }
    mergedMatrix_->fillComplete( mergedMap_, mergedMap_ ); // only for square matrices!
//    }
//    else {
//        mergedMatrix_->resumeFill();
//        for (UN i=0; i<blockMatrix_.size(); i++) {
//            for (UN j=0; j<blockMatrix_.size(); j++) {
//                if ( !blockMatrix_[i][j].is_null() )
//                    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"We need to implement mergeBlockAgain");
//            }
//        }
//    }
////    mergedMatrix_->print(Teuchos::VERB_EXTREME);
//
//    mergedMatrix_->fillComplete( mergedMap_, mergedMap_ ); // only for square matrices!

}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::determineLocalOffsets(){


    typedef Teuchos::ScalarTraits<LO> LOST;

    Teuchos::Array<LO> offsetRows( blockMatrix_.size() - 1, LOST::zero() );
    Teuchos::Array<LO> offsetCols( blockMatrix_.size() - 1, LOST::zero() );
    for (UN row=0; row<blockMatrix_.size() - 1; row++) {
        bool foundOffset = false;
        for (UN col=0; col<blockMatrix_.size() && !foundOffset; col++) {
            if ( !blockMatrix_[row][col].is_null() ){
                LO mli = blockMatrix_[row][col]->getMap("row")->getMaxLocalIndex();
                if (mli>=0)
                    offsetRows[row] = mli + 1;
            }
        }
    }
    for (UN col=0; col<blockMatrix_.size() - 1; col++) {
        bool foundOffset = false;
        for (UN row=0; row<blockMatrix_.size() && !foundOffset; row++) {
            if ( !blockMatrix_[row][col].is_null() ){
                LO mli = blockMatrix_[row][col]->getMap("col")->getMaxLocalIndex();
                if (mli>=0)
                    offsetCols[col] = mli + 1;
            }
        }
    }

    localBlockOffsets_ =
        Teuchos::rcp( new SmallMatrix< LOTuple_Type >( blockMatrix_.size() ) );

    LO offsetRow = LOST::zero();
    for (UN row=0; row<blockMatrix_.size(); row++) {
        LO offsetCol = LOST::zero();
        for (UN col=0; col<blockMatrix_.size(); col++) {
            LOTuple_Type tuple = Teuchos::tuple( offsetRow, offsetCol );
            (*localBlockOffsets_)[row][col] = tuple;
            if (offsetCols.size() > col)
                offsetCol += offsetCols[col];
        }
        if (offsetRows.size() > row)
            offsetRow += offsetRows[row];
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::determineGlobalOffsets(){

    typedef Teuchos::ScalarTraits<GO> GOST;

    Teuchos::Array<GO> offsetRows( blockMatrix_.size() - 1, GOST::zero() );
    Teuchos::Array<GO> offsetCols( blockMatrix_.size() - 1, GOST::zero() );
    for (UN row=0; row<blockMatrix_.size() - 1; row++) {
        bool foundOffset = false;
        for (UN col=0; col<blockMatrix_.size() && !foundOffset; col++) {
            if ( !blockMatrix_[row][col].is_null() )
            {
                offsetRows[row] = blockMatrix_[row][col]->getMap("row")->getMaxAllGlobalIndex() + 1;
                foundOffset = true;
            }
        }
    }
    for (UN col=0; col<blockMatrix_.size() - 1; col++) {
        bool foundOffset = false;
        for (UN row=0; row<blockMatrix_.size() && !foundOffset; row++) {
            if ( !blockMatrix_[row][col].is_null() )
            {
                offsetCols[col] = blockMatrix_[row][col]->getMap("col")->getMaxAllGlobalIndex() + 1;
                foundOffset = true;
            }
        }
    }

    globalBlockOffsets_ =
        Teuchos::rcp( new SmallMatrix< GOTuple_Type >( blockMatrix_.size() ) );

    GO offsetRow = GOST::zero();
    for (UN row=0; row<blockMatrix_.size(); row++) {
        GO offsetCol = GOST::zero();
        for (UN col=0; col<blockMatrix_.size(); col++) {
            GOTuple_Type tuple = Teuchos::tuple( offsetRow, offsetCol );
            (*globalBlockOffsets_)[row][col] = tuple;
            if (offsetCols.size() > col)
                offsetCol += offsetCols[col];
        }
        if (offsetRows.size() > row)
            offsetRow += offsetRows[row];
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::mergeBlockNew(UN blockRow, UN blockCol){

    MatrixPtr_Type matrix = blockMatrix_[blockRow][blockCol];
    if ( matrix->isLocallyIndexed() ) {

        Teuchos::ArrayView<const SC> values;
        Teuchos::ArrayView<const LO> indices;
        MapConstPtr_Type colMap = matrix->getMap("col");
        MapConstPtr_Type rowMap = matrix->getMap("row");

        for (UN i=0; i<matrix->getNodeNumRows(); i++) {
            matrix->getLocalRowView( i, indices, values );
            Teuchos::Array<GO> indicesGlobal( indices.size() );
            for (UN j=0; j<indices.size(); j++) {
                indicesGlobal[j] = colMap->getGlobalElement( indices[j] ) + (*globalBlockOffsets_)[blockRow][blockCol][1];
//                if (indicesGlobal[j]<0) {
//                    std::cout << "blockRow:"<< blockRow << " blockCol:" << blockCol << " "<<i << " "<< j << " indicesGlobal[j]:"<< indicesGlobal[j]<<" indices.size():"<<indices.size() << " indices[j]:"<< indices[j] <<" gBO:" << (*globalBlockOffsets_)[blockRow][blockCol][1] << std::endl;
//                }
            }
            GO offset = (*globalBlockOffsets_)[blockRow][blockCol][0];
            mergedMatrix_->insertGlobalValues( rowMap->getGlobalElement( i ) + offset, indicesGlobal(), values );
        }

    }
//    else if ( blockMatrix_[blockRow][blockCol]->isGloballyIndexed() ) {
//        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"We need to implement mergeBlockNew for isGloballyIndexed()==true.");
//    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Call fillComplete() before you merge the matrix blocks.");

}


template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::print(Teuchos::EVerbosityLevel verbLevel){

    for (UN i=0; i<blockMatrix_.size(); i++) {
        for (UN j=0; j<blockMatrix_.size(); j++) {
            if ( !blockMatrix_[i][j].is_null() ) {
                blockMatrix_[i][j]->print( verbLevel );
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::printMerge(Teuchos::EVerbosityLevel verbLevel){

    if ( !mergedMatrix_.is_null() )
        mergedMatrix_->print( verbLevel );

}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::LinearOpBase<SC> > BlockMatrix<SC,LO,GO,NO>::getThyraLinOp(){

    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() == 0, std::logic_error,"BlockMatrix size is 0.");
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > thyraLinOp;
    if (blockMatrix_.size() == 1){
        TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_[0][0].is_null(), std::runtime_error, "Block in BlockMatrix is null.");
        thyraLinOp = blockMatrix_[0][0]->getThyraLinOp( );
    }
    else{
        this->merge();
        thyraLinOp = mergedMatrix_->getThyraLinOp( );
    }
    return thyraLinOp;
}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP< const Thyra::BlockedLinearOpBase<SC> > BlockMatrix<SC,LO,GO,NO>::getThyraLinBlockOp() const{
    
    Teuchos::RCP< Thyra::DefaultBlockedLinearOp< SC > > dblOp = Thyra::defaultBlockedLinearOp<SC> ();
    dblOp->beginBlockFill( this->size(),  this->size() );
    
    for (int i=0; i<this->size(); i++) {
        for (int j=0; j<this->size(); j++) {
            if ( this->blockExists(i,j) ) {
                dblOp->setBlock( i, j,  blockMatrix_[i][j]->getThyraLinOp( ) );
            }
        }
    }
    
    dblOp->endBlockFill();
    Teuchos::RCP< const  Thyra::BlockedLinearOpBase< SC > > blOpConst = Teuchos::rcp_dynamic_cast<const Thyra::BlockedLinearOpBase< SC > > (dblOp);
//    Teuchos::RCP< const Thyra::BlockedLinearOpBase< SC > > blOpConst = Teuchos::rcp_const_cast<const Thyra::BlockedLinearOpBase< SC > > (blOp);
    return blOpConst;
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::apply(const BlockMultiVector_Type& X,
                                     BlockMultiVector_Type& Y ) const{

    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() == 0, std::logic_error,"BlockMatrix size is 0.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() != X.size(), std::logic_error,"BlockMatrix and BlockMultiVector have different sizes.");
    TEUCHOS_TEST_FOR_EXCEPTION( Y.size() != X.size(), std::logic_error,"BlockMultiVectors have different sizes.");

    for (int i=0; i<Y.size(); i++)
        Y[i]->putScalar( Teuchos::ScalarTraits<SC>::zero() );

    typedef Teuchos::ScalarTraits<SC> ST;

    for (int i=0; i<blockMatrix_.size(); i++) {
        for (int j=0; j<blockMatrix_.size(); j++){
            MultiVectorConstPtr_Type tmpX = X[j];
            if (!blockMatrix_[i][j].is_null())
            {
                blockMatrix_[i][j]->apply( *tmpX, *Y[i], Teuchos::NO_TRANS, ST::one(), ST::one() );
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::apply(const BlockMultiVector_Type& X,
                                     BlockMultiVector_Type& Y,
                                     const SmallMatrix<SC>& coeff) const{

    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() != coeff.size(), std::logic_error,"BlockMatrix size and coefficient size are different.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() == 0, std::logic_error,"BlockMatrix size is 0.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() != X.size(), std::logic_error,"BlockMatrix and BlockMultiVector have different sizes.");
    TEUCHOS_TEST_FOR_EXCEPTION( Y.size() != X.size(), std::logic_error,"BlockMultiVectors have different sizes.");

    for (int i=0; i<Y.size(); i++)
        Y[i]->putScalar( Teuchos::ScalarTraits<SC>::zero() );

    typedef Teuchos::ScalarTraits<SC> ST;

    for (int i=0; i<blockMatrix_.size(); i++) {
        for (int j=0; j<blockMatrix_.size(); j++){
            MultiVectorConstPtr_Type tmpX = X[j];
            if (!blockMatrix_[i][j].is_null())
                blockMatrix_[i][j]->apply( *tmpX, *Y[i], Teuchos::NO_TRANS, coeff[i][j], ST::one() );
        }
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::addMatrix(const SmallMatrix<SC>& coeffAlpha, const BlockMatrixPtr_Type &B, const SmallMatrix<SC>& coeffBeta){
    //B = alpha*A + beta*B.
    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() != coeffAlpha.size(), std::logic_error,"BlockMatrix size and coefficient size are different.");
    TEUCHOS_TEST_FOR_EXCEPTION( coeffAlpha.size() != coeffBeta.size(), std::logic_error,"Coefficients have different sizes.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() == 0, std::logic_error,"BlockMatrix size is 0.");
    TEUCHOS_TEST_FOR_EXCEPTION( blockMatrix_.size() != B->size(), std::logic_error,"BlockMatrices have different sizes.");

    // CAREFUL TO NOT IGNORE ANYTHING
    // B must have the correct blocks!
    for (int i=0; i<blockMatrix_.size(); i++) {
        for (int j=0; j<blockMatrix_.size(); j++){
            if (blockMatrix_[i][j].is_null()){
                if ( coeffBeta[i][j] != Teuchos::ScalarTraits<SC>::one() && B->blockExists(i,j) )
                {
                    B->getBlock(i,j)->scale( coeffBeta[i][j] );
                }
            }
            else
                blockMatrix_[i][j]->addMatrix( coeffAlpha[i][j], B->getBlock(i,j), coeffBeta[i][j] );
        }
    }
}

template <class SC, class LO, class GO, class NO>
void BlockMatrix<SC,LO,GO,NO>::writeMM(std::string fN) const{
    for (int i=0; i<blockMatrix_.size(); i++) {
        for (int j=0; j<blockMatrix_.size(); j++){
            if (!blockMatrix_[i][j].is_null() ){
                std::string fileName = fN + to_string(i) + to_string(j) + ".mm" ;
                blockMatrix_[i][j]->writeMM(fileName);
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
typename BlockMatrix<SC,LO,GO,NO>::MatrixPtr_Type BlockMatrix<SC,LO,GO,NO>::getMergedMatrix() {
    if (this->size()>1) {
        this->merge();
        return mergedMatrix_;
    }
    else
        return this->getBlock(0,0);
}

template <class SC, class LO, class GO, class NO>
typename BlockMatrix<SC,LO,GO,NO>::BlockMapPtr_Type BlockMatrix<SC,LO,GO,NO>::getMap() {

    return blockMap_;
}

}
#endif
