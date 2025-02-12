#ifndef MATRIX_DEF_hpp
#define MATRIX_DEF_hpp
#include "Matrix_decl.hpp"

/*!
 Defintion of Matrix

 @brief  Matrix
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {
//using namespace Teuchos;
template <class SC, class LO, class GO, class NO>
Matrix<SC,LO,GO,NO>::Matrix():
	matrix_()
{

}

template <class SC, class LO, class GO, class NO>
Matrix<SC,LO,GO,NO>::Matrix( TpetraMatrixPtr_Type& tpetraMatPtrIn ):
matrix_( tpetraMatPtrIn )
{

}

template <class SC, class LO, class GO, class NO>
Matrix<SC,LO,GO,NO>::Matrix( MapConstPtr_Type map , LO numEntries):
matrix_(  )
{
    matrix_ = Teuchos::RCP( new TpetraMatrix_Type(map->getTpetraMap(), numEntries)); //, plist) );

}


template <class SC, class LO, class GO, class NO>
Matrix<SC,LO,GO,NO>::Matrix( MatrixPtr_Type matrixIn ):
matrix_( )
{
	matrix_ =Teuchos::RCP( new TpetraMatrix_Type( matrixIn->getMap()->getTpetraMap(), matrixIn->getGlobalMaxNumRowEntries() ));
	if(matrixIn->isLocallyIndexed())
	{
		Teuchos::ArrayView<const SC> values;
        Teuchos::ArrayView<const LO> indices;
        
        MapConstPtr_Type colMap = matrixIn->getMap("col");
        MapConstPtr_Type rowMap = matrixIn->getMap("row");

        for (UN i=0; i<matrixIn->getNodeNumRows(); i++)
		{
            matrixIn->getLocalRowView( i, indices, values );
            Teuchos::Array<GO> indicesGlobal( indices.size() );
            for (UN j=0; j<indices.size(); j++)
			{
                indicesGlobal[j] = colMap->getGlobalElement( indices[j] );

            }
            matrix_->insertGlobalValues( rowMap->getGlobalElement( i ), indicesGlobal(), values );
        }
	}

	matrix_->fillComplete( matrixIn->getMap("domain")->getTpetraMap(), matrixIn->getMap("range")->getTpetraMap() );
}

template <class SC, class LO, class GO, class NO>
Matrix<SC,LO,GO,NO>::~Matrix(){

}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::insertGlobalValues(GO globalRow, const Teuchos::ArrayView< const GO > &cols, const Teuchos::ArrayView< const SC > &vals){

    TEUCHOS_TEST_FOR_EXCEPTION(matrix_.is_null(),std::runtime_error,"");
    matrix_->insertGlobalValues( globalRow, cols, vals );
}

template <class SC, class LO, class GO, class NO>
LO Matrix<SC,LO,GO,NO>::getNodeNumRows() const{
    return matrix_->getLocalNumRows();
}


template <class SC, class LO, class GO, class NO>
typename  Matrix<SC,LO,GO,NO>::MapConstPtr_Type Matrix<SC,LO,GO,NO>::getMap(std::string map_string){

    TEUCHOS_TEST_FOR_EXCEPTION(matrix_.is_null(),std::runtime_error,"RCP<Matrix> is null.");
    TpetraMapConstPtr_Type tpetraMap;
    if (!map_string.compare("row")) {
        tpetraMap = matrix_->getRowMap();
    }
    else if (!map_string.compare("col")) {
        tpetraMap = matrix_->getColMap();
    }
    else if (!map_string.compare("domain")) {
        tpetraMap = matrix_->getDomainMap();
    }
    else if (!map_string.compare("range")) {
        tpetraMap = matrix_->getRangeMap();
    }
    else
        tpetraMap = matrix_->getMap();

    return Teuchos::rcp( new Map_Type(tpetraMap) );
}

template <class SC, class LO, class GO, class NO>
typename  Matrix<SC,LO,GO,NO>::MapConstPtr_Type Matrix<SC,LO,GO,NO>::getMap(std::string map_string) const{

    TEUCHOS_TEST_FOR_EXCEPTION(matrix_.is_null(),std::runtime_error,"RCP<Matrix> is null.");
    TpetraMapConstPtr_Type tpetraMap;
    if (!map_string.compare("row")) {
        tpetraMap = matrix_->getRowMap();
    }
    else if (!map_string.compare("col")) {
        tpetraMap = matrix_->getColMap();
    }
    else if (!map_string.compare("domain")) {
        tpetraMap = matrix_->getDomainMap();
    }
    else if (!map_string.compare("range")) {
        tpetraMap = matrix_->getRangeMap();
    }
    else
        tpetraMap = matrix_->getMap();

    return Teuchos::rcp( new Map_Type(tpetraMap) );
}


template <class SC, class LO, class GO, class NO>
typename  Matrix<SC,LO,GO,NO>::TpetraMapConstPtr_Type Matrix<SC,LO,GO,NO>::getMapTpetra(std::string map_string){

    TEUCHOS_TEST_FOR_EXCEPTION(matrix_.is_null(),std::runtime_error,"RCP<Matrix> is null.");

    if (!map_string.compare("row")) {
        return matrix_->getRowMap();
    }
    else if (!map_string.compare("col")) {
        return matrix_->getColMap();
    }
    else if (!map_string.compare("domain")) {
        return matrix_->getDomainMap();
    }
    else if (!map_string.compare("range")) {
        return matrix_->getRangeMap();
    }
    return matrix_->getMap();
}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::LinearOpBase<SC> > Matrix<SC,LO,GO,NO>::getThyraLinOp() const{
    //Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsWrapMatrix = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*matrix_);
    //return Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyra(crsWrapMatrix.getCrsMatrix());

    Teuchos::RCP<Thyra::LinearOpBase<SC> > thyraOp = Teuchos::null;

    Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> >  tpCrsMat = matrix_; // = xTpCrsMat->getTpetra_CrsMatrixNonConst();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));
    Teuchos::RCP<Tpetra::RowMatrix<SC,LO,GO,NO> > tpRowMat   = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<SC,LO,GO,NO> >(tpCrsMat, true);
    Teuchos::RCP<Tpetra::Operator <SC,LO,GO,NO> > tpOperator = Teuchos::rcp_dynamic_cast<Tpetra::Operator<SC,LO,GO,NO> >(tpRowMat, true);

    thyraOp = Thyra::createLinearOp(tpOperator);

    return thyraOp;

}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > Matrix<SC,LO,GO,NO>::getThyraLinOpNonConst() {
    //Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsWrapMatrix = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*matrix_);
    //return Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> > (Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyra(crsWrapMatrix.getCrsMatrix()) );
    Teuchos::RCP<Thyra::LinearOpBase<SC> > thyraOp = Teuchos::null;

    Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> >  tpCrsMat = matrix_; // = xTpCrsMat->getTpetra_CrsMatrixNonConst();
    TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(tpCrsMat));
    Teuchos::RCP<Tpetra::RowMatrix<SC,LO,GO,NO> > tpRowMat   = Teuchos::rcp_dynamic_cast<Tpetra::RowMatrix<SC,LO,GO,NO> >(tpCrsMat, true);
    Teuchos::RCP<Tpetra::Operator <SC,LO,GO,NO> > tpOperator = Teuchos::rcp_dynamic_cast<Tpetra::Operator<SC,LO,GO,NO> >(tpRowMat, true);

    thyraOp = Thyra::createLinearOp(tpOperator);

    return Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> > (thyraOp); // is this now const or not?? 
}
    
template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::print(Teuchos::EVerbosityLevel verbLevel){

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    matrix_->describe(*out,verbLevel);
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::resumeFill(){
    matrix_->resumeFill();
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::fillComplete(){
    matrix_->fillComplete();
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::fillComplete(MapConstPtr_Type domainMap, MapConstPtr_Type rangeMap){
    matrix_->fillComplete( domainMap->getTpetraMap(), rangeMap->getTpetraMap() );
}

template <class SC, class LO, class GO, class NO>
bool Matrix<SC,LO,GO,NO>::isFillComplete(){
    return matrix_->isFillComplete();
}

template <class SC, class LO, class GO, class NO>
bool Matrix<SC,LO,GO,NO>::isLocallyIndexed(){
    return matrix_->isLocallyIndexed();
}

//typename Matrix<SC,LO,GO,NO>::ThyraLinOpPtr_Type Matrix<SC,LO,GO,NO>::getThyraLinOp(){
//
//    return ;
//}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::getGlobalRowView(GO globalRow, Teuchos::ArrayView< const GO > &indices, Teuchos::ArrayView< const SC > &values) const{
    TEUCHOS_TEST_FOR_EXCEPTION(matrix_->isLocallyIndexed(),std::logic_error,"Underlying matrix is locally indexed and we can not use a global row view. Global row copy is valid here and can be implemented.");
    
    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::global_inds_host_view_type Indices;  //ArrayView< const LO > indices
    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type Values;
    //typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::local_inds_host_view_type indices;
    //typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::values_host_view_type values;

    matrix_->getGlobalRowView(globalRow, Indices, Values);
    indices = Teuchos::ArrayView<const GO> (Indices.data(), Indices.extent(0));
    values = Teuchos::ArrayView<const SC> (reinterpret_cast<const SC*>(Values.data()), Values.extent(0));
  
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::getLocalRowView(LO localRow, Teuchos::ArrayView< const LO > &indices, Teuchos::ArrayView< const SC > &values) const{
    TEUCHOS_TEST_FOR_EXCEPTION(matrix_->isGloballyIndexed(),std::logic_error,"Underlying matrix is globally indexed and we can not use a local row view.");
    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type Indices;  //ArrayView< const LO > indices
    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type Values;
    //typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::local_inds_host_view_type indices;
    //typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::values_host_view_type values;

    matrix_->getLocalRowView(localRow, Indices, Values);
    indices = Teuchos::ArrayView<const LO> (Indices.data(), Indices.extent(0));
    values = Teuchos::ArrayView<const SC> (reinterpret_cast<const SC*>(Values.data()), Values.extent(0));
  
}


template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::replaceGlobalValues(GO globalRow, const Teuchos::ArrayView< const GO > &indices, const Teuchos::ArrayView< const SC > &values){
    matrix_->replaceGlobalValues(globalRow, indices, values);
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::replaceLocalValues(LO localRow, const Teuchos::ArrayView< const LO > &indices, const Teuchos::ArrayView< const SC > &values){
    matrix_->replaceLocalValues(localRow, indices, values);
}

template <class SC, class LO, class GO, class NO>
typename  Matrix<SC,LO,GO,NO>::TpetraMatrixConstPtr_Type Matrix<SC,LO,GO,NO>::getTpetraMatrix() const{
    return matrix_;
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::apply(const MultiVector_Type& X,
                                MultiVector_Type& Y,
                                Teuchos::ETransp mode,
                                SC alpha,
                                SC beta) const{

    TpetraMVPtr_Type yTpetra = Y.getTpetraMultiVectorNonConst();

    matrix_->apply( *X.getTpetraMultiVector(), *yTpetra, mode, alpha, beta );
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::scale(const SC& alpha) {

    matrix_->scale( alpha );

}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::writeMM(std::string fileName) const{
    TEUCHOS_TEST_FOR_EXCEPTION( matrix_.is_null(), std::runtime_error,"Matrix in writeMM is null.");
   // TEUCHOS_TEST_FOR_EXCEPTION( !(matrix_->getMap()->lib()==Xpetra::UseTpetra), std::logic_error,"Only available for Tpetra underlying lib.");
    //typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;
    //typedef Teuchos::RCP<TpetraCrsMatrix> TpetraCrsMatrixPtr;

    //Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*matrix_);
    //Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());

    TpetraMatrixPtr_Type tpetraMat = matrix_;

    Tpetra::MatrixMarket::Writer< TpetraMatrix_Type > tpetraWriter;

    tpetraWriter.writeSparseFile(fileName, tpetraMat, "matrix", "");
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::addMatrix(SC alpha, const MatrixPtr_Type &B, SC beta){
    //B = alpha*A + beta*B.
    if (B->isFillComplete())
        B->resumeFill();
    TEUCHOS_TEST_FOR_EXCEPTION( B->isLocallyIndexed(), std::runtime_error,"Matrix in is locally index but Trilinos Epetra/Tpetra can not add to a matrix at this stage.");
    
    //const Tpetra::CrsMatrix<SC,LO,GO,NO>& tpA = Xpetra::Helpers<SC,LO,GO,NO>::Op2TpetraCrs(A);
    //Tpetra::CrsMatrix<SC,LO,GO,NO>& tpB = Xpetra::Helpers<SC,LO,GO,NO>::Op2NonConstTpetraCrs(B);

    Tpetra::MatrixMatrix::Add(*matrix_, false, alpha, *B->matrix_, beta);

    //Xpetra::MatrixMatrix<SC,LO,GO,NO>::TwoMatrixAdd( *matrix_, false, alpha, *B->matrix_, beta );
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::toMV( MultiVectorPtr_Type& mv ){
    
    MapConstPtr_Type map = this->getMap("row");
    MapConstPtr_Type mapCol = this->getMap("col");
    // number of matrix columns will be the number of multivectors
    int numMV = mapCol->getMaxAllGlobalIndex() + 1;
    
    mv = Teuchos::rcp( new MultiVector_Type ( map, numMV ) );
    Teuchos::ArrayView< const LO > indices;
    Teuchos::ArrayView< const SC > values;
    Teuchos::ArrayRCP< SC > valuesMV;
    GO globalCol = -1;
    for (int i=0; i<this->getNodeNumRows(); i++) {
        this->getLocalRowView( i, indices, values );
        for (int j=0; j<indices.size(); j++) {
            globalCol = mapCol->getGlobalElement( indices[j] );
            valuesMV = mv->getDataNonConst( globalCol );
            valuesMV[ i ] = values[ j ];
        }
    }
}
template <class SC, class LO, class GO, class NO>
LO Matrix<SC,LO,GO,NO>::getGlobalMaxNumRowEntries() const{
    return matrix_->getGlobalMaxNumRowEntries();
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::insertLocalValues(LO localRow, const Teuchos::ArrayView< const LO > &cols, const Teuchos::ArrayView< const SC > &vals){

    TEUCHOS_TEST_FOR_EXCEPTION(matrix_.is_null(),std::runtime_error,"");
    matrix_->insertLocalValues( localRow, cols, vals );
}


template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::importFromVector( MatrixPtr_Type mvIn, bool reuseImport, std::string combineMode, std::string type) {

    //TEUCHOS_TEST_FOR_EXCEPTION( getNodeNumCols() != mvIn->getNodeNumCols(), std::logic_error,"MultiVectors for fillFromVector have different number of vectors.");

    if ( importer_.is_null() || !reuseImport) {
        if (type=="Forward")
            importer_ = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>( mvIn->getMapTpetra(), this->getMapTpetra() ));
        else if(type=="Reverse")
            importer_ = Teuchos::RCP(new Tpetra::Import<LO, GO, NO>( this->getMapTpetra(), mvIn->getMapTpetra() ));
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for import. Choose Forward or Reverse");
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION( !importer_->getSourceMap()->isSameAs( *mvIn->getMap()->getTpetraMap() ), std::logic_error,"Source maps of Importer and Matrix are not the same.");
        TEUCHOS_TEST_FOR_EXCEPTION( !importer_->getTargetMap()->isSameAs( *this->getMap()->getTpetraMap() ), std::logic_error,"Target maps of Importer and Matrix are not the same.");
    }

        
    if (type=="Forward") {
        if ( !combineMode.compare("Insert") )
            matrix_->doImport ( *mvIn->getTpetraMatrix(), *importer_, Tpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            matrix_->doImport ( *mvIn->getTpetraMatrix(), *importer_, Tpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else if(type=="Reverse"){
        if ( !combineMode.compare("Insert") )
            matrix_->doExport ( *mvIn->getTpetraMatrix(), *importer_, Tpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            matrix_->doExport ( *mvIn->getTpetraMatrix(), *importer_, Tpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for import. Choose Forward or Reverse");
}

template <class SC, class LO, class GO, class NO>
void Matrix<SC,LO,GO,NO>::exportFromVector( MatrixPtr_Type mvIn, bool reuseExport, std::string combineMode, std::string type) {

    //TEUCHOS_TEST_FOR_EXCEPTION( getNodeNumCols() != mvIn->getNodeNumCols(), std::logic_error,"MultiVectors for fillFromVector have different number of vectors.");

    if ( exporter_.is_null() || !reuseExport) {
        if (type=="Forward")
            exporter_ = Teuchos::RCP(new Tpetra::Export<LO, GO, NO>( mvIn->getMapTpetra(), this->getMapTpetra() ));
        else if(type=="Reverse")
            exporter_ = Teuchos::RCP(new Tpetra::Export<LO, GO, NO>( this->getMapTpetra(), mvIn->getMapTpetra() ));
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for import. Choose Forward or Reverse");
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION( !exporter_->getSourceMap()->isSameAs( *this->getMap()->getTpetraMap() ), std::logic_error,"Source maps of Exporter and Multivector are not the same.");
        TEUCHOS_TEST_FOR_EXCEPTION( !exporter_->getTargetMap()->isSameAs( *mvIn->getMap()->getTpetraMap() ), std::logic_error,"Target maps of Exporter and Multivector are not the same.");
    }

        
    if (type=="Forward") {
        if ( !combineMode.compare("Insert") )
            matrix_->doExport ( *mvIn->getTpetraMatrix(), *exporter_, Tpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            matrix_->doExport ( *mvIn->getTpetraMatrix(), *exporter_, Tpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else if(type=="Reverse"){
        if ( !combineMode.compare("Insert") )
            matrix_->doImport ( *mvIn->getTpetraMatrix(), *exporter_, Tpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            matrix_->doImport ( *mvIn->getTpetraMatrix(), *exporter_, Tpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for import. Choose Forward or Reverse");
}

}
#endif
