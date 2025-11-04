#ifndef TIMESTEPPINGTOOLS_hpp
#define TIMESTEPPINGTOOLS_hpp

#include "feddlib/problems/problems_config.h"
#include "feddlib/core/FEDDCore.hpp"

#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/core/General/ExporterTxt.hpp"

/*!
 Declaration of TimeSteppingTools
 
 @brief  TimeSteppingTools
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {


class TimeSteppingTools {

public:
    typedef default_sc SC;
    typedef default_lo LO;
    typedef default_go GO;
    typedef default_no NO;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
    typedef typename Matrix_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef Teuchos::Array<BlockMultiVectorPtr_Type> BlockMultiVectorPtrArray_Type;

    typedef ExporterTxt ExporterTxt_Type;
    typedef Teuchos::RCP<ExporterTxt_Type> ExporterTxtPtr_Type;

    typedef typename Matrix_Type::CommConstPtr_Type CommConstPtr_Type;

    enum timeSteppingType   {NON_ADAPTIVE, ADAPTIVE} ;
    enum errorType          {EUCLIDIAN, L2} ;
    /* --------------------------------------------------------------------------------- */
    CommConstPtr_Type		comm_;
    ParameterListPtr_Type parameterList_;
    int                 butcherTableNmb_;
    int                 BDFNmb_;
    timeSteppingType    tsType_;
    double 	tEnd_;
    double 	dt_;
    double 	t_;
    /* adaptive variables*/
    double 	dt_prev_;
    double  dt_adaptive_;
    double 	rho_;
    double 	tolAdaptive_;
    double 	dtmin_;
    double 	dtmax_;
    int		convOrder_;
    errorType     adaptiveError_;
    int 	adaptiveCalculation_;
    double  error_;
    double  error_prev_;
    /* Butcher table information*/
    int 	stages_;
    bool    stifflyAcc_;
    bool    stifflyAccEmbedded_;
    vec2D_dbl_ptr_Type  butcherTable_;
    vec_dbl_ptr_Type    BDFInformation_;
    vec_dbl_ptr_Type 	b_embedded_;
    vec_dbl_ptr_Type 	gamma_vec_;
    bool 	verbose_;
    ExporterTxtPtr_Type exporterTxtTime_;
    ExporterTxtPtr_Type exporterTxtDt_;
    ExporterTxtPtr_Type exporterTxtError_;

    int RKType_;
    // Newmark-Variablen
    double beta_;
    double gamma_;

    /* --------------------------------------------------------------------------------- */

    TimeSteppingTools();

    TimeSteppingTools(ParameterListPtr_Type parameterList , CommConstPtr_Type comm);

    void setParameter();

    void setTableInformationRK();

    void setInformationBDF();

    double getInformationBDF(int i);

    void correctPressure(MultiVectorPtr_Type &newP, MultiVectorConstPtr_Type lastP);
    
    int getBDFNumber();

    double currentTime();

    bool continueTimeStepping();

    void calculateSolution(BlockMultiVectorPtr_Type &sol, BlockMultiVectorPtrArray_Type &rkSolVec, BlockMatrixPtr_Type massSystem, BlockMultiVectorPtr_Type solShort = Teuchos::null);

    void adaptiveTimestep(BlockMultiVectorPtr_Type &sol, BlockMultiVectorPtrArray_Type &rkSolVec, BlockMatrixPtr_Type massSystem, BlockMultiVectorPtr_Type solShort);

    void calculateNewDt(BlockMultiVectorPtr_Type &solDiff, BlockMatrixPtr_Type massSystem);

    void advanceTime(bool printInfo=false);

    void printInfo();

    int getNmbStages();

    double getButcherTableCoefficient(int row , int col);

    double getButcherTableC(int row);

    double get_dt();

    double get_dt_prev();

    void setupTxtExporter();

    double get_beta();

    double get_gamma();
};
}
#endif
