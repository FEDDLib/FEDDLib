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

	/*! 
		\brief Intertion of values in global row 'globalRow'. Matrix is distributed nodewise. If values are added to same row, they automatically get summed up.

	*/
    void insertGlobalValues (GO globalRow, const Teuchos::ArrayView< const GO > &cols, const Teuchos::ArrayView< const SC > &vals);

	/*!
		\brief Returns the local number of rows.
		
	*/
    LO getNodeNumRows() const;

	/*!
		\brief Returns map of type " ". i.e. row or column map
	*/
    MapConstPtr_Type getMap(std::string map_string="");

	/*!
		\brief Returns map of type " ". i.e. row or column map
	*/
    MapConstPtr_Type getMap(std::string map_string="") const;

	/*!
		\brief Return map in Xpetra Format of type " ".
	*/
    XpetraMapConstPtr_Type getMapXpetra(std::string map_string="");

	/*!
		\brief i.e. for NOX
	*/
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > getThyraLinOp() const;

	/*!
		\brief i.e. for NOX
	*/
    Teuchos::RCP<Thyra::LinearOpBase<SC> > getThyraLinOpNonConst();

	/*!
		\brief printing matrix
	*/
    void print(Teuchos::EVerbosityLevel verbLevel=Teuchos::VERB_EXTREME);

	/*!
		\brief Resuming filling process. But only limited options i.e. scaling remain.
	*/
    void resumeFill();

	/*!
		\brief after inserting global values into matrix. After this step the column map is fixed. Row map is used for filling.
	*/
    void fillComplete();

	/*!
		\brief Filling of Matrix based on specific domainmap (column map) and rangeMap (rowmap).
	*/
    void fillComplete(MapConstPtr_Type domainMap, MapConstPtr_Type rangeMap);

	/*!
		\brief
	*/
    bool isLocallyIndexed();

	/*!
		\brief Check if matrix is already filled complete.
	*/
    bool isFillComplete();
//    void fillComplete(with maps);

//    ThyraLinOpPtr_Type getThyraLinOp();

	/*!
		\brief Extracting single rows of Matrix with global row ID. Indices returns global indices of entries stored in values.
	*/
    void getGlobalRowView(GO globalRow, Teuchos::ArrayView< const GO > &indices, Teuchos::ArrayView< const SC > &values) const;

	/*!
		\brief Extracting single rows of Matrix with local row ID. Indices returns local indices of entries stored in values.
	*/
    void getLocalRowView(LO localRow, Teuchos::ArrayView< const LO > &indices, Teuchos::ArrayView< const SC > &values) const;

	/*!
		\brief Replacing single rows of Matrix with global row ID. Indices returns global indices of entries stored in values.
	*/
    void replaceGlobalValues(GO globalRow, const Teuchos::ArrayView< const GO > &indices, const Teuchos::ArrayView< const SC > &values);

	/*!
		\brief Replacing single rows of Matrix with local row ID. Indices returns local indices of entries stored in values.
	*/
    void replaceLocalValues(LO localRow, const Teuchos::ArrayView< const LO > &indices, const Teuchos::ArrayView< const SC > &values);

	/*!
		\brief Return matrix in Xpetra Format of type " ".
	*/
    XpetraMatrixConstPtr_Type getXpetraMatrix() const;
    
	/*!
		\brief Matrix Vector Operation. Applying MultiVector X to this. Y = alpha * (this)^mode * X + beta * Y. Mode being transposed or not. 
	*/
    void apply(const MultiVector_Type& X,
               MultiVector_Type& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               SC alpha = Teuchos::ScalarTraits<SC>::one(),
               SC beta = Teuchos::ScalarTraits<SC>::zero() ) const;

	/*!
		\brief Scaling this with constant alpha.
	*/
    void scale(const SC& alpha);

	/*!
		\brief Writing Matrix in file.
	*/
    void writeMM(std::string fileName="matrix.mm") const;

	/*!
		\brief B = alpha*this + beta*B.
	*/
    void addMatrix(SC alpha, const MatrixPtr_Type &B, SC beta);
    
	/*!
		\brief Turning Matrix into MultiVector Format
	*/
    void toMV( MultiVectorPtr_Type& mv );
    
	/*!
		\brief Maximum number of entries in any row of the matrix, over all processes.
	*/
    LO getGlobalMaxNumRowEntries() const;

    void insertLocalValues (LO localRow, const Teuchos::ArrayView< const LO > &cols, const Teuchos::ArrayView< const SC > &vals);
	/* !
		Matrix Analogue to MultiVector Import. Based on Row Map of Matrix mvIn. 
	*/

    void importFromVector( MatrixPtr_Type mvIn, bool reuseImport = false, std::string combineMode = "Insert", std::string type="Forward" );

	/* !
		Matrix Analogue to MultiVector Export. Based on Row Map of Matrix mvIn. 
	*/
    void exportFromVector( MatrixPtr_Type mvIn, bool reuseExport = false, std::string combineMode = "Insert", std::string type="Forward" );
	

private:

    XpetraMatrixPtr_Type matrix_;
    XpetraImportPtr_Type importer_;    
	XpetraExportPtr_Type exporter_;
};
}

#endif
