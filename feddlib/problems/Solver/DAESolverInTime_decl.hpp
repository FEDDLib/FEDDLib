#ifndef DAESOLVERINTIME_DECL_hpp
#define DAESOLVERINTIME_DECL_hpp

#include "feddlib/problems/problems_config.h"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/General/ExporterTxt.hpp"

#include "NonLinearSolver.hpp"
#include "TimeSteppingTools.hpp"

/*!
 Declaration of DAESolverInTime

 @brief  DAESolverInTime
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

template <class SC, class LO, class GO, class NO>
class FSI;
template <class SC, class LO, class GO, class NO>
class MeshUnstructured;
template <class SC, class LO, class GO, class NO>
class Domain;

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class DAESolverInTime {

public:
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;
    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;
    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;
    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;
    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;
    typedef Teuchos::Array<MultiVectorPtr_Type> MultiVectorPtrArray_Type;
    typedef Teuchos::Array<MultiVectorConstPtr_Type> MultiVectorConstPtrArray_Type;
    typedef Teuchos::Array<BlockMultiVectorPtr_Type> BlockMultiVectorPtrArray_Type;
    typedef typename Problem_Type::FEFacPtr_Type FEFacPtr_Type;
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    
    typedef NonLinearProblem<SC,LO,GO,NO> NonLinProb_Type;
    typedef Teuchos::RCP<NonLinProb_Type> NonLinProbPtr_Type;
    typedef TimeProblem<SC,LO,GO,NO> TimeProblem_Type;
    typedef Teuchos::RCP<TimeProblem_Type> TimeProblemPtr_Type;

    typedef ExporterParaView<SC,LO,GO,NO> Exporter_Type;
    typedef Teuchos::RCP<Exporter_Type> ExporterPtr_Type;
    typedef Teuchos::RCP<ExporterTxt> ExporterTxtPtr_Type;

    typedef FSI<SC,LO,GO,NO> FSIProblem_Type;
    typedef Teuchos::RCP<FSIProblem_Type> FSIProblemPtr_Type;

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef Teuchos::RCP<Domain_Type > DomainPtr_Type;
    typedef typename Domain_Type::Mesh_Type Mesh_Type;
    typedef typename Domain_Type::MeshPtr_Type MeshPtr_Type;

    
    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef Teuchos::RCP<MeshUnstr_Type> MeshUnstrPtr_Type;

    DAESolverInTime(CommConstPtr_Type comm);

    DAESolverInTime(ParameterListPtr_Type &parameterList, CommConstPtr_Type comm);

    DAESolverInTime(SmallMatrix<int> &timeStepDef, ParameterListPtr_Type &parameterList, CommConstPtr_Type comm);

    ~DAESolverInTime();

//    void setProblem(NonLinearProblem_Type &problem);
    
    void setProblem(Problem_Type &problem);
    
    void defineTimeStepping(SmallMatrix<int> &timeStepDef);

    void setupTimeStepping();

    void advanceInTime();

    void advanceWithLoadStepping();

    void advanceInTimeLinear();

    void advanceInTimeNonLinear();

    // Newmark-Verfahren (fuer lineare Probleme)
    void advanceInTimeLinearNewmark();

    // Newmkar-Verfahren (fuer nichtlineare Probleme)
    void advanceInTimeNonLinearNewmark();

    // Wendet die Zeitdiskretisierung separat auf Fluid und Struktur an und
    // schreibe dann vor dem nlSolve() alles in das FSI-System hinein
    void advanceInTimeFSI();

    void advanceInTimeLinearMultistep();

    void advanceInTimeNonLinearMultistep();
    
    void advanceInTimeLinearExternal();

    void advanceInTimeNonLinearExternal();
    
    void setTimeStep(double dt);

    void setFinalTime(double T);

    void checkTimeSteppingDef();

    void exportTimestep();
    
    void exportTimestep(BlockMultiVectorPtr_Type& solShort);

    void setupExporter();
    
    void setupExporter(BlockMultiVectorPtr_Type& solShort);

    void closeExporter();

    void addRhsDAE(SmallMatrix<double> coeff, BlockMatrixPtr_Type bMat, BlockMultiVectorPtr_Type vec);

    void addRhsDAE(SmallMatrix<double> coeff, BlockMultiVectorPtr_Type vec);
    
//    void addRhsDAE(SmallMatrix<double> coeff, MultiVector_Type &vec);
//
//    void addRhsDAE(SmallMatrix<double> coeff, BE_ptr_vec_Type &bENl, MultiVector_Type &vec, MultiVector_Type &sourceTerm);
//
//    void addRhsDAE(SmallMatrix<double> coeff, MultiVector_Type &vec, MultiVector_Type &sourceTerm);

    void addSourceTermToRHS(double coeff);

    void getMassCoefficients(SmallMatrix<double> &massCoeff);
    
    void getMultiStageCoefficients(SmallMatrix<double> &problemCoeff, int stage, int stagePrior, bool forRhs = false);
    
    void buildMultiStageRhs( int stage, Teuchos::Array<BlockMatrixPtr_Type>& matrixPrevStages, BlockMultiVectorPtrArray_Type& solutionPrevStages );
    
    CommConstPtr_Type comm_;
    bool verbose_;
    ParameterListPtr_Type parameterList_;
    bool isTimeSteppingDefined_;
    
    ProblemPtr_Type problem_;
    TimeProblemPtr_Type problemTime_;

    SmallMatrix<int> timeStepDef_;
    Teuchos::RCP<TimeSteppingTools>	timeSteppingTool_;
    FEFacPtr_Type feFactory_;
    
    std::vector<ExporterPtr_Type> exporter_vector_;
    MultiVectorConstPtrArray_Type export_solution_vector_;
    bool boolExporterSetup_;

private:

#ifdef FEDD_TIMER
    TimePtr_Type solveProblemTimer_;
#endif
#ifdef FEDD_TIMER
    TimePtr_Type reassmbleAddInterfaceRHSTimer_;
    TimePtr_Type reassmbleUpdateMeshDisplacementTimer_;
    TimePtr_Type reassmbleSolveGeometryTimer_;
    TimePtr_Type reassmbleMoveMeshTimer_;
    TimePtr_Type reassmbleSolidMassAndRHSTimer_;
    TimePtr_Type reassmbleForTimeTimer_;
    TimePtr_Type reassmbleUpdateFluidInTimeTimer_;
#endif
};
}
#endif
