#ifndef MULTIVECTOR_DEF_hpp
#define MULTIVECTOR_DEF_hpp
#include "MultiVector_decl.hpp"

/*!
 Defintion of MultiVector
 
 @brief  MultiVector
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {

extern template class MultiVector<default_sc, default_lo, default_go, default_no>;

template <class SC, class LO, class GO, class NO>
MultiVector<SC,LO,GO,NO>::MultiVector():
multiVector_(),
map_(),
importer_(),
exporter_()
{

}

//template <class SC, class LO, class GO, class NO>
//MultiVector<SC,LO,GO,NO>::MultiVector(const EpetraMat_Type& epetraMatIn):
//    MultiVector_()
//{
//#ifdef FEDD_UseTpetra
//    Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
//#else
//    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
//#endif
//    CommConstPtr_Type teuchosComm = rcp( new Teuchos::MpiComm<int> ( epetraMatIn.Map().Comm() ) );
//    XpetraMultiVectorPtr_Type MultiVector_ = FROSch::ConvertToXpetra( lib, epetraMatIn, teuchosComm );
//}

template <class SC, class LO, class GO, class NO>
MultiVector<SC,LO,GO,NO>::MultiVector( MapConstPtr_Type map, UN nmbVectors ):
multiVector_( ),
map_(map),
importer_(),
exporter_()
{
    multiVector_ = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build( map->getXpetraMap(), nmbVectors );
}

template <class SC, class LO, class GO, class NO>
MultiVector<SC,LO,GO,NO>::MultiVector( XpetraMultiVectorPtr_Type& xpetraMVPtrIn ):
multiVector_( xpetraMVPtrIn ),
map_(),
importer_(),
exporter_()
{
    map_.reset( new Map_Type( xpetraMVPtrIn->getMap() ) );
}

template <class SC, class LO, class GO, class NO>
MultiVector<SC,LO,GO,NO>::MultiVector( MultiVectorConstPtr_Type mvIn ):
multiVector_( ),
map_(),
importer_(),
exporter_()
{
    multiVector_ = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build( mvIn->getMap()->getXpetraMap(), mvIn->getNumVectors() );
    map_.reset( new Map_Type( *mvIn->getMap() ) );
    for (UN j=0; j<this->getNumVectors(); j++) {
        Teuchos::ArrayRCP< const SC > valuesIn = mvIn->getData(j);
        Teuchos::ArrayRCP< SC > valuesThis = this->getDataNonConst(j);
        for (UN i=0; i<valuesThis.size(); i++)//can this be quicker?
            valuesThis[i] = valuesIn[i];
    }

}

template <class SC, class LO, class GO, class NO>
MultiVector<SC,LO,GO,NO>::~MultiVector(){

}
template <class SC, class LO, class GO, class NO>
typename MultiVector<SC,LO,GO,NO>::MapConstPtr_Type MultiVector<SC,LO,GO,NO>::getMap() const{
    return map_;
}

template <class SC, class LO, class GO, class NO>
typename MultiVector<SC,LO,GO,NO>::MapPtr_Type MultiVector<SC,LO,GO,NO>::getMapNonConst() {
    return Teuchos::rcp_const_cast<Map_Type>( map_ );
}
    
template <class SC, class LO, class GO, class NO>
typename MultiVector<SC,LO,GO,NO>::XpetraMapConstPtr_Type MultiVector<SC,LO,GO,NO>::getMapXpetra() const{
    TEUCHOS_TEST_FOR_EXCEPTION(multiVector_.is_null(),std::runtime_error,"RCP<MultiVector> is null.");

    return multiVector_->getMap();
}

// checking whether multivector is null
template <class SC, class LO, class GO, class NO>
bool MultiVector<SC,LO,GO,NO>::is_null() const{
    return multiVector_.is_null();
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::replaceGlobalValue (GO globalRow, UN vectorIndex, const SC &value){
    multiVector_->replaceGlobalValue( globalRow, vectorIndex, value );
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::sumIntoGlobalValue (GO globalRow, UN vectorIndex, const SC &value){
    multiVector_->sumIntoGlobalValue( globalRow, vectorIndex, value );
}

template <class SC, class LO, class GO, class NO>
LO MultiVector<SC,LO,GO,NO>::getLocalLength() const{
    return multiVector_->getLocalLength();
}

template <class SC, class LO, class GO, class NO>
Teuchos::ArrayRCP< const SC > MultiVector<SC,LO,GO,NO>::getData(UN i) const{
    return multiVector_->getData(i);
}

template <class SC, class LO, class GO, class NO>
Teuchos::ArrayRCP< SC > MultiVector<SC,LO,GO,NO>::getDataNonConst(UN i) const{
    return multiVector_->getDataNonConst(i);
}


template <class SC, class LO, class GO, class NO>
UN MultiVector<SC,LO,GO,NO>::getNumVectors() const{
    return multiVector_->getNumVectors();
}

template <class SC, class LO, class GO, class NO>
typename MultiVector<SC,LO,GO,NO>::XpetraMultiVectorConstPtr_Type MultiVector<SC,LO,GO,NO>::getXpetraMultiVector() const{
    return multiVector_;
}

template <class SC, class LO, class GO, class NO>
typename MultiVector<SC,LO,GO,NO>::XpetraMultiVectorPtr_Type MultiVector<SC,LO,GO,NO>::getXpetraMultiVectorNonConst() {
    return multiVector_;
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::print(Teuchos::EVerbosityLevel verbLevel) const{

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    multiVector_->describe(*out,verbLevel);
}

template <class SC, class LO, class GO, class NO>
Teuchos::RCP<Thyra::MultiVectorBase<SC> > MultiVector<SC,LO,GO,NO>::getThyraMultiVector( ) {
    Teuchos::RCP<Thyra::MultiVectorBase<SC> > mv = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<SC> >(Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(multiVector_));
    return mv;
}


template <class SC, class LO, class GO, class NO>
Teuchos::RCP<const Thyra::MultiVectorBase<SC> > MultiVector<SC,LO,GO,NO>::getThyraMultiVectorConst( ) const{

     return Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector( multiVector_ );
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::fromThyraMultiVector( Teuchos::RCP< Thyra::MultiVectorBase<SC> > thyraMV){
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    TEUCHOS_TEST_FOR_EXCEPTION( multiVector_.is_null(), std::runtime_error,"MultiVector is null but we need to know the underlying lib for function fromThyraMultiVector().");

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;

    typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> TOVE_Type;
    typedef Tpetra::MultiVector<SC,LO,GO,NO> TpetraMultiVector_Type;
    typedef Xpetra::TpetraMultiVector<SC,LO,GO,NO> XTpetraMultiVector_Type;

    if (!this->getMap()->getUnderlyingLib().compare("Tpetra") ) {
        Teuchos::RCP<TpetraMultiVector_Type> tMV = TOVE_Type::getTpetraMultiVector( thyraMV );
        Teuchos::RCP<XTpetraMultiVector_Type> xTMV = Teuchos::rcp( new XTpetraMultiVector_Type( tMV ) );
        multiVector_ = xTMV;
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,"we must implement fromThyraMultiVector() for other underlying libs.");
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::norm2(const Teuchos::ArrayView< typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms) const {
    TEUCHOS_TEST_FOR_EXCEPTION( multiVector_.is_null(), std::runtime_error,"MultiVector in norm2 is null.");
    multiVector_->norm2(norms);
}

// Inf Norm of Multivector
template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::normInf(const Teuchos::ArrayView< typename Teuchos::ScalarTraits<SC>::magnitudeType> &normsInf) const {
    TEUCHOS_TEST_FOR_EXCEPTION( multiVector_.is_null(), std::runtime_error,"MultiVector in normInf is null.");
    multiVector_->normInf(normsInf);
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::dot(MultiVectorConstPtr_Type a, const Teuchos::ArrayView< typename Teuchos::ScalarTraits<SC>::magnitudeType> &dots) const {
    TEUCHOS_TEST_FOR_EXCEPTION( multiVector_.is_null(), std::runtime_error,"MultiVector in dot is null.");
    multiVector_->dot( *a->getXpetraMultiVector(), dots );
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::abs(MultiVectorConstPtr_Type a) {
    TEUCHOS_TEST_FOR_EXCEPTION( multiVector_.is_null(), std::runtime_error,"MultiVector in abs is null.");
    multiVector_->abs( *a->getXpetraMultiVector());
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::update( const SC& alpha, const MultiVector_Type& A, const SC& beta) {
    TEUCHOS_TEST_FOR_EXCEPTION( getNumVectors() != A.getNumVectors(), std::logic_error,"MultiVectors for update have different number of vectors.");

    multiVector_->update( alpha, *A.getXpetraMultiVector(), beta );

}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::update( const SC& alpha, const MultiVector_Type& A, const SC& beta , const MultiVector_Type& B, const SC& gamma) {
    TEUCHOS_TEST_FOR_EXCEPTION( getNumVectors() != A.getNumVectors(), std::logic_error,"MultiVectors for update have different number of vectors.");

    multiVector_->update( alpha, *A.getXpetraMultiVector(), beta, *B.getXpetraMultiVector(), gamma );

}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::putScalar( const SC& alpha ){
    multiVector_->putScalar( alpha );
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::scale( const SC& alpha ){
    multiVector_->scale( alpha );
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const SC &alpha, MultiVectorConstPtr_Type &A, MultiVectorConstPtr_Type &B, const SC &beta){
    multiVector_->multiply( transA, transB, alpha, *A->getXpetraMultiVector(), *B->getXpetraMultiVector(), beta );
}
template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const SC &alpha, BlockMultiVectorConstPtr_Type &A, BlockMultiVectorConstPtr_Type &B, const SC &beta){
//    if (this->getMap()->getCommNonConst()->getRank()==0)
//        std::cout << "### For testing purposes only." << std::endl;
    
    for (int i=0; i<A->size(); i++){
        MultiVectorConstPtr_Type a = A->getBlock(i);
        MultiVectorConstPtr_Type b = B->getBlock(i);
        if ( i==0 )
            this->multiply( transA, transB, alpha, a, b, beta );
        else
            this->multiply( transA, transB, alpha, a, b, Teuchos::ScalarTraits<SC>::one() );
    }

}
    
template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::importFromVector( MultiVectorConstPtr_Type mvIn, bool reuseImport, std::string combineMode, std::string type) {
        
    TEUCHOS_TEST_FOR_EXCEPTION( getNumVectors() != mvIn->getNumVectors(), std::logic_error,"MultiVectors for fillFromVector have different number of vectors.");

    if ( importer_.is_null() || !reuseImport) {
        if (type=="Forward")
            importer_ = Xpetra::ImportFactory<LO,GO,NO>::Build ( mvIn->getMapXpetra(), this->getMapXpetra() );
        else if(type=="Reverse")
            importer_ = Xpetra::ImportFactory<LO,GO,NO>::Build ( this->getMapXpetra(), mvIn->getMapXpetra() );
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for import. Choose Forward or Reverse");
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION( !importer_->getSourceMap()->isSameAs( *mvIn->getMap()->getXpetraMap() ), std::logic_error,"Source maps of Importer and Multivector are not the same.");
        TEUCHOS_TEST_FOR_EXCEPTION( !importer_->getTargetMap()->isSameAs( *this->getMap()->getXpetraMap() ), std::logic_error,"Target maps of Importer and Multivector are not the same.");
    }

        
    if (type=="Forward") {
        if ( !combineMode.compare("Insert") )
            multiVector_->doImport ( *mvIn->getXpetraMultiVector(), *importer_, Xpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            multiVector_->doImport ( *mvIn->getXpetraMultiVector(), *importer_, Xpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else if(type=="Reverse"){
        if ( !combineMode.compare("Insert") )
            multiVector_->doExport ( *mvIn->getXpetraMultiVector(), *importer_, Xpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            multiVector_->doExport ( *mvIn->getXpetraMultiVector(), *importer_, Xpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for import. Choose Forward or Reverse");
}

template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::exportFromVector( MultiVectorConstPtr_Type mvIn, bool reuseExport, std::string combineMode, std::string type) {
    TEUCHOS_TEST_FOR_EXCEPTION( getNumVectors() != mvIn->getNumVectors(), std::logic_error,"MultiVectors for exportToVector have different number of vectors.");

    if ( exporter_.is_null() || !reuseExport) {
        if (type=="Forward")
            exporter_ = Xpetra::ExportFactory<LO,GO,NO>::Build ( mvIn->getMapXpetra(), this->getMapXpetra() );
        else if(type=="Reverse")
            exporter_ = Xpetra::ExportFactory<LO,GO,NO>::Build ( this->getMapXpetra(), mvIn->getMapXpetra() );
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for export. Choose Forward or Reverse");
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION( !exporter_->getSourceMap()->isSameAs( *this->getMap()->getXpetraMap() ), std::logic_error,"Source maps of Exporter and Multivector are not the same.");
        TEUCHOS_TEST_FOR_EXCEPTION( !exporter_->getTargetMap()->isSameAs( *mvIn->getMap()->getXpetraMap() ), std::logic_error,"Target maps of Exporter and Multivector are not the same.");
    }
    if (type=="Forward") {
        if ( !combineMode.compare("Insert") )
            multiVector_->doExport ( *mvIn->getXpetraMultiVector(), *exporter_, Xpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            multiVector_->doExport ( *mvIn->getXpetraMultiVector(), *exporter_, Xpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else if (type=="Reverse") {
        if ( !combineMode.compare("Insert") )
            multiVector_->doImport ( *mvIn->getXpetraMultiVector(), *exporter_, Xpetra::INSERT);
        else if ( !combineMode.compare("Add") )
            multiVector_->doImport ( *mvIn->getXpetraMultiVector(), *exporter_, Xpetra::ADD);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown combine mode.");
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,"Unknown type for export. Choose Forward or Reverse");
}
    
template <class SC, class LO, class GO, class NO>
void MultiVector<SC,LO,GO,NO>::writeMM(std::string fileName) const{
    TEUCHOS_TEST_FOR_EXCEPTION( multiVector_.is_null(), std::runtime_error,"MultiVector in writeMM is null.");
    TEUCHOS_TEST_FOR_EXCEPTION( !(multiVector_->getMap()->lib()==Xpetra::UseTpetra), std::logic_error,"Only available for Tpetra underlying lib.");

    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;
    typedef Teuchos::RCP<TpetraCrsMatrix> TpetraCrsMatrixPtr;

    typedef Tpetra::MultiVector<SC,LO,GO,NO> TpetraMultiVector;
    typedef Teuchos::RCP<TpetraMultiVector> TpetraMultiVectorPtr;

    const Xpetra::TpetraMultiVector<SC,LO,GO,NO>& xTpetraMultiVector = dynamic_cast<const Xpetra::TpetraMultiVector<SC,LO,GO,NO> &>(*multiVector_);

    TpetraMultiVectorPtr tpetraMultiVector = xTpetraMultiVector.getTpetra_MultiVector();

    Tpetra::MatrixMarket::Writer< TpetraCrsMatrix > tpetraWriter;

    tpetraWriter.writeDenseFile(fileName, tpetraMultiVector, "multivector", "");
}
    
template <class SC, class LO, class GO, class NO>
typename MultiVector<SC,LO,GO,NO>::MultiVectorConstPtr_Type MultiVector<SC,LO,GO,NO>::getVector( int i ) const{

    XpetraMultiVectorConstPtr_Type xpetraMV = multiVector_->getVector( i );
    XpetraMultiVectorPtr_Type xpetraMVNonConst = Teuchos::rcp_const_cast<XpetraMultiVector_Type>( xpetraMV );
    MultiVectorConstPtr_Type singleMV = Teuchos::rcp( new const MultiVector_Type ( xpetraMVNonConst ) );
    return singleMV;
    
}
    
template <class SC, class LO, class GO, class NO>
typename MultiVector<SC,LO,GO,NO>::MultiVectorPtr_Type MultiVector<SC,LO,GO,NO>::sumColumns() const{
    
    MultiVectorPtr_Type sumMV = Teuchos::rcp( new MultiVector_Type ( map_, 1 ) );
    sumMV->putScalar(0.);
    for (int i=0; i<this->getNumVectors(); i++)
        sumMV->getXpetraMultiVectorNonConst()->getVectorNonConst(0)->update( 1., *multiVector_->getVector(i), 1. );

    return sumMV;
}

template <class SC, class LO, class GO, class NO>
SC MultiVector<SC,LO,GO,NO>::getMax() const{
    TEUCHOS_TEST_FOR_EXCEPTION( this->getNumVectors() > 1, std::runtime_error, "numMultiVector>1: max function not implemented!");
    return multiVector_->getVector(0)->normInf();
}

//typename MultiVector<SC,LO,GO,NO>::ThyraLinOpPtr_Type MultiVector<SC,LO,GO,NO>::getThyraLinOp(){
//
//    return ;
//}

}
#endif
