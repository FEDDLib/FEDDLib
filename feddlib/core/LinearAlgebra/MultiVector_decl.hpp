#ifndef MULTIVECTOR_DECL_hpp
#define MULTIVECTOR_DECL_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "Map.hpp"
#include "BlockMap.hpp"
#include "BlockMultiVector.hpp"
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Thyra_LinearOpBase_decl.hpp>
#include "Xpetra_ThyraUtils.hpp"
#include <Teuchos_VerboseObject.hpp>
#include <MatrixMarket_Tpetra.hpp>

/*!
 Declaration of MultiVector

 @brief  MultiVector
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {
template <class SC, class LO, class GO, class NO>
class BlockMultiVector;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MultiVector {

public:
    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const XpetraMapConstPtr_Type XpetraMapConstPtrConst_Type;

    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
    
    typedef Xpetra::MultiVector<SC,LO,GO,NO> XpetraMultiVector_Type;
    typedef Teuchos::RCP<XpetraMultiVector_Type> XpetraMultiVectorPtr_Type;
    typedef Teuchos::RCP<const XpetraMultiVector_Type> XpetraMultiVectorConstPtr_Type;
    typedef const XpetraMultiVectorConstPtr_Type XpetraMultiVectorConstPtrConst_Type;

    typedef Tpetra::Map<LO,GO,NO> TpetraMap_Type;
    typedef Teuchos::RCP<const TpetraMap_Type> TpetraMapConstPtr_Type;

    typedef Tpetra::MultiVector<SC,LO,GO,NO> TpetraMultiVector_Type;
    typedef Teuchos::RCP<TpetraMultiVector_Type> TpetraMultiVectorPtr_Type;
    typedef Teuchos::RCP<const TpetraMultiVector_Type> TpetraMultiVectorConstPtr_Type;

    typedef Xpetra::Import<LO,GO,NO> XpetraImport_Type;
    typedef Teuchos::RCP<XpetraImport_Type> XpetraImportPtr_Type;

    typedef Xpetra::Export<LO,GO,NO> XpetraExport_Type;
    typedef Teuchos::RCP<XpetraExport_Type> XpetraExportPtr_Type;
    
    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<Comm_Type> CommPtr_Type;    
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    MultiVector();

    MultiVector( MapConstPtr_Type map, UN nmbVectors=1 );

    MultiVector( XpetraMultiVectorPtr_Type& xpetraMVPtrIn );

    MultiVector( MultiVectorConstPtr_Type mvIn );

    ~MultiVector();

    MultiVector_Type& operator= (const MultiVector_Type& rhs) {
        *multiVector_ = *rhs.getXpetraMultiVector();
        return *this;
    }

	bool is_null() const;

    MapConstPtr_Type getMap() const;
    
    MapPtr_Type getMapNonConst();

    XpetraMapConstPtr_Type getMapXpetra() const;

    void replaceGlobalValue (GO globalRow, UN vectorIndex, const SC &value);

    void sumIntoGlobalValue (GO globalRow, UN vectorIndex, const SC &value);

    LO getLocalLength() const;

    Teuchos::ArrayRCP< const SC >  getData(UN i) const;

    Teuchos::ArrayRCP< SC > getDataNonConst(UN i) const;

    UN getNumVectors() const;

    void print(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME) const;

    XpetraMultiVectorConstPtr_Type getXpetraMultiVector() const;

    XpetraMultiVectorPtr_Type getXpetraMultiVectorNonConst();

    Teuchos::RCP< Thyra::MultiVectorBase<SC> > getThyraMultiVector( );

    Teuchos::RCP<const Thyra::MultiVectorBase<SC> > getThyraMultiVectorConst( ) const;

    void fromThyraMultiVector( Teuchos::RCP< Thyra::MultiVectorBase<SC> > thyraMV);

    void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms) const;

    void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms) const;

    void dot(MultiVectorConstPtr_Type a, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &dots) const;

	// Calculate absolute value of Multivector
	void abs(MultiVectorConstPtr_Type a);
    //this = alpha*A + beta*this
    void update( const SC& alpha, const MultiVector_Type& A, const SC& beta );

    //this = alpha*A + beta*B + gamma*this
    void update( const SC& alpha, const MultiVector_Type& A, const SC& beta , const MultiVector_Type& B, const SC& gamma);

    // Matrix-matrix multiplication: this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const SC &alpha, MultiVectorConstPtr_Type &A, MultiVectorConstPtr_Type &B, const SC &beta);

    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const SC &alpha, BlockMultiVectorConstPtr_Type &A, BlockMultiVectorConstPtr_Type &B, const SC &beta);
    
    void putScalar( const SC& alpha );

    void scale( const SC& alpha );

    void importFromVector( MultiVectorConstPtr_Type mvIn, bool reuseImport = false, std::string combineMode = "Insert", std::string type="Forward" );
    
    void exportFromVector( MultiVectorConstPtr_Type mvIn, bool reuseExport = false, std::string combineMode = "Insert", std::string type="Forward" );
    
    void writeMM(std::string fileName="mv.mm") const;
    
    MultiVectorConstPtr_Type getVector( int i ) const;
    
    MultiVectorPtr_Type sumColumns() const;
    
    SC getMax() const;
    
private:

    XpetraMultiVectorPtr_Type multiVector_;
    MapConstPtr_Type map_;
    XpetraImportPtr_Type importer_;
    XpetraExportPtr_Type exporter_;
};
}

#endif
