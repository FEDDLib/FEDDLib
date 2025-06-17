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

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Import.hpp>

/*!
 Declaration of MultiVector

 @brief  MultiVector
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC, class LO, class GO, class NO>
class BlockMultiVector;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MultiVector {

public:
   /*typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const XpetraMapConstPtr_Type XpetraMapConstPtrConst_Type;*/

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::RCP<const BlockMultiVector_Type> BlockMultiVectorConstPtr_Type;
    
   // typedef Xpetra::MultiVector<SC,LO,GO,NO> XpetraMultiVector_Type;
   // typedef Teuchos::RCP<XpetraMultiVector_Type> XpetraMultiVectorPtr_Type;
   // typedef Teuchos::RCP<const XpetraMultiVector_Type> XpetraMultiVectorConstPtr_Type;
   // typedef const XpetraMultiVectorConstPtr_Type XpetraMultiVectorConstPtrConst_Type;


   // typedef Xpetra::Import<LO,GO,NO> XpetraImport_Type;
   // typedef Teuchos::RCP<XpetraImport_Type> XpetraImportPtr_Type;

   // typedef Xpetra::Export<LO,GO,NO> XpetraExport_Type;
   // typedef Teuchos::RCP<XpetraExport_Type> XpetraExportPtr_Type;
    
    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<Comm_Type> CommPtr_Type;    
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    // -------------
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;

    typedef Tpetra::Map<LO,GO,NO> TpetraMap_Type;
    typedef Teuchos::RCP<TpetraMap_Type> TpetraMapPtr_Type;
    typedef Teuchos::RCP<const TpetraMap_Type> TpetraMapConstPtr_Type;
    typedef const TpetraMapConstPtr_Type TpetraMapConstPtrConst_Type;

    typedef Tpetra::MultiVector<SC,LO,GO,NO> TpetraMultiVector_Type;
    typedef Teuchos::RCP<TpetraMultiVector_Type> TpetraMultiVectorPtr_Type;
    typedef Teuchos::RCP<const TpetraMultiVector_Type> TpetraMultiVectorConstPtr_Type;
    typedef const TpetraMultiVectorConstPtr_Type TpetraMultiVectorConstPtrConst_Type;

    typedef Tpetra::Import<LO,GO,NO> TpetraImport_Type;
    typedef Teuchos::RCP<TpetraImport_Type> TpetraImportPtr_Type;

    typedef Tpetra::Export<LO,GO,NO> TpetraExport_Type;
    typedef Teuchos::RCP<TpetraExport_Type> TpetraExportPtr_Type;

 

    MultiVector();

    
    /// @brief Initialize tpetra multivector based on underyling map and number of vectors within. In almost all cases nmbVectors is 1.
    /// @param map parallel local to global indexing of row entries
    /// @param nmbVectors number of vectors in multivector (seen maybe as different columns)
    MultiVector( MapConstPtr_Type map, UN nmbVectors=1 );

    /// @brief Initialize tpetra multivector based on input multivector. Uses underlying map. !! Probably, this is not a deep copy. Both mv have the same pointer.
    /// @param TpetraMVPtrIn tpetra multivector used to build new multivector 
    MultiVector( TpetraMultiVectorPtr_Type& TpetraMVPtrIn );

    /// @brief Initialize tpetra multivector based on input multivector. Uses underlying map and value information to construct new mv
    /// @param mvIn multivector used to build new multivector
    MultiVector( MultiVectorConstPtr_Type mvIn );

    /// @brief Destructor
    ~MultiVector();


    /// @brief This will replace *this contents with the rhs input. Updated to deep copy, as this is necessary so both this and rhs would NOT have the same pointer
    /// @param rhs source for copy
    /// @return destination/result of copy
    MultiVector_Type& operator= (const MultiVector_Type& rhs) {
        //*multiVector_ = *rhs.getTpetraMultiVector(); // old version which worked with xpetra
        FEDDLIB_NOTIFICATION("MultiVector_decl",rhs.getMap()->getComm()->getRank() == 0, " '=' creating a deep copy of input vector into this.");
        Tpetra::deep_copy<SC,LO,GO,NO>(*multiVector_,*rhs.getTpetraMultiVector()); // (destination, source)
        return *this;
    }

    /// @brief checking whether the multiVector exists
	bool is_null() const;

    /// @brief Return underlying map
    /// @return  MapConstPtr_Type
    MapConstPtr_Type getMap() const;
    
    /// @brief Return non constant version of underlying map
    /// @return MapPtr_Type
    MapPtr_Type getMapNonConst();

    /// @brief Return direct tpetra map of underlying map
    /// @return TpetraMapConstPtr_Type
    TpetraMapConstPtr_Type getMapTpetra() const;

    /// @brief Replace global value in mv 
    /// @param globalRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// @param vectorIndex [in] Column index of the entry to modify.
    /// @param value [in] Incoming value to add to the entry.
    void replaceGlobalValue (GO globalRow, UN vectorIndex, const SC &value);

    /// \param lclRow [in] Local row index of the entry to modify.
    ///   Must be a valid local index in this MultiVector's Map on the
    ///   calling process.
    /// \param vectorIndex [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    void replaceLocalValue (LO localRow, UN vectorIndex, const SC &value);

    /// \brief Update (+=) a value in host memory, using global row index.
    ///
    /// Add the given value to the existing value at row \c gblRow (a
    /// global index) and column \c col.  The column index is zero
    /// based.
    /// \param globalRow [in] Global row index of the entry to modify.
    ///   This <i>must</i> be a valid global row index on the calling
    ///   process with respect to the MultiVector's Map.
    /// \param vectorIndex [in] Column index of the entry to modify.
    /// \param value [in] Incoming value to add to the entry.
    void sumIntoGlobalValue (GO globalRow, UN vectorIndex, const SC &value);

    LO getLocalLength() const;

    /// @brief Get data of multivector
    /// @param i 'column' of multivector
    /// @return Array format of const entries (on my processor)
    Teuchos::ArrayRCP< const SC >  getData(UN i) const;

    /// @brief Get data of multivector
    /// @param i 'column' of multivector
    /// @return Array format of entries (on my processor)
    Teuchos::ArrayRCP< SC > getDataNonConst(UN i) const;

    /// @brief Get number of multivector (columns)
    UN getNumVectors() const;

    /// @brief Printing mv. Depending on verbosity level, output increases
    void print(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME) const;

    TpetraMultiVectorConstPtr_Type getTpetraMultiVector() const;

    TpetraMultiVectorPtr_Type getTpetraMultiVectorNonConst();

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
    
    void readMM(std::string fileName) const;
    
    MultiVectorConstPtr_Type getVector( int i ) const;
    
    MultiVectorPtr_Type sumColumns() const;
    
    SC getMax() const;
    
private:

    TpetraMultiVectorPtr_Type multiVector_;
    MapConstPtr_Type map_;
    TpetraImportPtr_Type importer_;
    TpetraExportPtr_Type exporter_;
};
}

#endif
