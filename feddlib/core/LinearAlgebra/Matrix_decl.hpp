#ifndef MATRIX_DECL_hpp
#define MATRIX_DECL_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "Map.hpp"
#include "MultiVector.hpp"
#include "BlockMultiVector.hpp"

#include <Xpetra_MatrixFactory.hpp>
#include "Xpetra_ThyraUtils.hpp"
#include <Teuchos_VerboseObject.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

/*!
 Declaration of Matrix

 @brief  Matrix
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

using namespace std;
namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Matrix {

public:

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtrConst_Type;

    typedef Xpetra::Matrix<SC,LO,GO,NO> XpetraMatrix_Type;
    typedef Teuchos::RCP<XpetraMatrix_Type> XpetraMatrixPtr_Type;
    typedef Teuchos::RCP<const XpetraMatrix_Type> XpetraMatrixConstPtr_Type;
    typedef const Teuchos::RCP<XpetraMatrixConstPtr_Type> XpetraMatrixConstPtrConst_Type;

    typedef Xpetra::MultiVector<SC,LO,GO,NO> XpetraMV_Type;
    typedef Teuchos::RCP<XpetraMV_Type> XpetraMVPtr_Type;

    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
//    typedef Teuchos::RCP<const MultiVector_Type> BlockMultiVectorConstPtr_Type;
    
    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    typedef Xpetra::Import<LO,GO,NO> XpetraImport_Type;
    typedef Teuchos::RCP<XpetraImport_Type> XpetraImportPtr_Type;

    typedef Xpetra::Export<LO,GO,NO> XpetraExport_Type;
    typedef Teuchos::RCP<XpetraExport_Type> XpetraExportPtr_Type;

    Matrix();

    Matrix( XpetraMatrixPtr_Type& xpetraMatPtrIn );

    Matrix( MapConstPtr_Type map , LO numEntries);
//    Matrix(const EpetraMat_Type& epetraMatIn);

    Matrix( MatrixPtr_Type matrixIn );

    ~Matrix();

    Matrix& operator+=(const Matrix& matIn);

    void insertGlobalValues (GO globalRow, const Teuchos::ArrayView< const GO > &cols, const Teuchos::ArrayView< const SC > &vals);

    LO getNodeNumRows() const;

    MapConstPtr_Type getMap(string map_string="");

    MapConstPtr_Type getMap(string map_string="") const;

    XpetraMapConstPtr_Type getMapXpetra(string map_string="");

    Teuchos::RCP<const Thyra::LinearOpBase<SC> > getThyraLinOp() const;

    Teuchos::RCP<Thyra::LinearOpBase<SC> > getThyraLinOpNonConst();

    void print(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME);

    void resumeFill();

    void fillComplete();

    void fillComplete(MapConstPtr_Type domainMap, MapConstPtr_Type rangeMap);

    bool isLocallyIndexed();

    bool isFillComplete();
//    void fillComplete(with maps);

//    ThyraLinOpPtr_Type getThyraLinOp();

    void getGlobalRowView(GO globalRow, Teuchos::ArrayView< const GO > &indices, Teuchos::ArrayView< const SC > &values) const;

    void getLocalRowView(LO localRow, Teuchos::ArrayView< const LO > &indices, Teuchos::ArrayView< const SC > &values) const;

    void replaceGlobalValues(GO globalRow, const Teuchos::ArrayView< const GO > &indices, const Teuchos::ArrayView< const SC > &values);

    void replaceLocalValues(LO localRow, const Teuchos::ArrayView< const LO > &indices, const Teuchos::ArrayView< const SC > &values);

    XpetraMatrixConstPtr_Type getXpetraMatrix() const;
    
    void apply(const MultiVector_Type& X,
               MultiVector_Type& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               SC alpha = Teuchos::ScalarTraits<SC>::one(),
               SC beta = Teuchos::ScalarTraits<SC>::zero() ) const;

    void scale(const SC& alpha);

    void writeMM(std::string fileName="matrix.mm") const;

    void addMatrix(SC alpha, const MatrixPtr_Type &B, SC beta);
    
    void toMV( MultiVectorPtr_Type& mv );
    
    LO getGlobalMaxNumRowEntries() const;

    void insertLocalValues (LO localRow, const Teuchos::ArrayView< const LO > &cols, const Teuchos::ArrayView< const SC > &vals);

    void importFromVector( MatrixPtr_Type mvIn, bool reuseImport = false, std::string combineMode = "Insert", std::string type="Forward" );
    void exportFromVector( MatrixPtr_Type mvIn, bool reuseExport = false, std::string combineMode = "Insert", std::string type="Forward" );
	

private:

    XpetraMatrixPtr_Type matrix_;
    XpetraImportPtr_Type importer_;    
	XpetraExportPtr_Type exporter_;

};
}

#endif
