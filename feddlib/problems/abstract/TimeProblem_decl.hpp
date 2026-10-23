#ifndef TIMEPROBLEM_DECL_hpp
#define TIMEPROBLEM_DECL_hpp

#include "feddlib/problems/problems_config.h"
#include "feddlib/core/FEDDCore.hpp"
#include "Problem.hpp"

#include <Thyra_StateFuncModelEvaluatorBase.hpp>
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/General/CheckpointFiles.hpp"
#include <set>
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/FE/Domain.hpp"

/*!
 Declaration of TimeProblem

 @brief  TimeProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


/*
 If nonlinear iterative scheme is used, assemble() always assembles the Oseen system (FixedPoint).
 reAssemble("Newton") then assembles the Newton system. We need to adapt the functionality to general systems/problems
 */
namespace FEDD {
template<class SC_, class LO_, class GO_, class NO_>
class Preconditioner;
template<class SC_, class LO_, class GO_, class NO_>
class NonLinearProblem;

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class TimeProblem: public Thyra::StateFuncModelEvaluatorBase<SC>  {

public:

    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef typename Problem_Type::Map_Type Map_Type;
    typedef typename Problem_Type::MapPtr_Type MapPtr_Type;
    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;
    typedef Teuchos::RCP<const Matrix_Type> MatrixConstPtr_Type;

    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;
    typedef typename Problem_Type::BlockMatrixConstPtr_Type BlockMatrixConstPtr_Type;
    typedef Teuchos::Array<BlockMatrixPtr_Type> BlockMatrixPtrArray_Type;

    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;

    typedef typename Problem_Type::Domain_Type::Mesh_Type Mesh_Type;
    typedef typename Problem_Type::Domain_Type::MeshPtr_Type MeshPtr_Type;

    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;
    typedef typename Problem_Type::BlockMultiVectorConstPtr_Type BlockMultiVectorConstPtr_Type;
    typedef Teuchos::Array<BlockMultiVectorPtr_Type> BlockMultiVectorPtrArray_Type;

    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;

    typedef typename Problem_Type::FEFac_Type FEFac_Type;
    typedef Teuchos::RCP<FEFac_Type> FEFacPtr_Type;
    typedef Teuchos::RCP<const FEFac_Type> FEFacConstPtr_Type;
    typedef typename Problem_Type::BCConstPtr_Type BCConstPtr_Type;

    typedef NonLinearProblem<SC,LO,GO,NO> NonLinProb_Type;
    typedef Teuchos::RCP<NonLinProb_Type> NonLinProbPtr_Type;

    typedef typename Problem_Type::Preconditioner_Type Preconditioner_Type;
    typedef typename Problem_Type::PreconditionerPtr_Type PreconditionerPtr_Type;

    typedef typename Problem_Type::LinSolverBuilderPtr_Type LinSolverBuilderPtr_Type;

    using TpetraTypes = TpetraTypedefs<SC,LO,GO,NO>;
    using TpetraMatrix_Type = typename TpetraTypes::TpetraMatrix_Type;
    using TpetraOp_Type = typename TpetraTypes::TpetraOp_Type;

    using ThyraTypes = ThyraTypedefs<SC>;
    using ThyraVecSpace_Type = typename ThyraTypes::ThyraVecSpace_Type;
    using ThyraVec_Type = typename ThyraTypes::ThyraVec_Type;
    typedef Thyra::BlockedLinearOpBase<SC> ThyraBlockOp_Type;

    using ThyraOp_Type = typename ThyraTypes::ThyraOp_Type;    
    typedef ExporterParaView<SC,LO,GO,NO> Exporter_Type;
    typedef Teuchos::RCP<Exporter_Type> ExporterPtr_Type;
    typedef std::vector<ExporterPtr_Type> ExporterPtrVec_Type;
    
        
    TimeProblem(Problem_Type& problem, CommConstPtr_Type comm);
    
    void setTimeDef( SmallMatrix<int>& def );
    
    void assemble(std::string type="") const;
    
    virtual void reAssembleAndFill( BlockMatrixPtr_Type bMat, std::string type="FixedPoint" );

    void updateRhs();

    void updateMultistepRhs(vec_dbl_Type& coeff, int nmbToUse);

    // Mit Massematrix aus vorherigen Zeitschritten
    void updateMultistepRhsFSI(vec_dbl_Type& coeff, int nmbToUse);

    // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
    // Bei Newmark lautet dies:
    // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
    // wobei u' = v (velocity) und u'' = w (acceleration).
    void updateNewmarkRhs(double dt, double beta, double gamma, vec_dbl_Type coeff);

    void combineSystems() const;

    void initializeCombinedSystems() const;

    void assembleMassSystem( ) const;

    void setTimeParameters(SmallMatrix<double> &daeParameters, SmallMatrix<double> &timeParameters);

    int solveAndUpdate( const std::string& criterion, double& criterionValue );

    int solveAndUpdate( const std::string& criterion, double& criterionValue, vec_dbl_Type& criterionValueVec );

    int solveUpdate();

    int update();

    int solve( BlockMultiVectorPtr_Type rhs = Teuchos::null );
    
    double calculateResidualNorm();

    vec_dbl_Type calculateResidualNormVec() const;

    void calculateNonLinResidualVec(std::string type="standard", double time=0.) const;

    bool getVerbose();

    void addBoundaries(BCConstPtr_Type bcFactory);

    void setBoundaries(double time=.0);

    void setBoundariesRHS(double time=.0);

    void setBoundariesSystem() const;

//    void assembleExternal( std::string type ) const;    

    BlockMultiVectorPtr_Type getRhs();

    BlockMultiVectorPtr_Type getRhs() const;
    
    BlockMultiVectorPtr_Type getSolution();

    BlockMultiVectorConstPtr_Type getSolutionConst() const;
    
    BlockMultiVectorPtr_Type getResidual();

    BlockMultiVectorConstPtr_Type getResidualConst() const;
    
    DomainConstPtr_Type getDomain(int i) const;

    std::string getFEType(int i);

    std::string getVariableName(int i);

    int getDofsPerNode(int i);

    BlockMatrixPtr_Type getSystem();

    BlockMatrixPtr_Type getSystemCombined();

    BlockMatrixPtr_Type getSystemCombined() const;

    ProblemPtr_Type getUnderlyingProblem();

    void updateSolutionPreviousStep();

    /// @brief After a restart: read the problem's own restart data (e.g. element history) of the restart time, once.
    void importRestartValues();

    void updateSolutionMultiPreviousStep(int nmbSteps);

    /// Adaptive time stepping: keep the state at the start of a time step (solution, previous
    /// solutions, Newmark velocity and acceleration, time) and go back to it for a repeated step
    void saveState();
    void restoreState();

    void updateSystemMassMultiPreviousStep(int nmbSteps);

    // Verschiebung (u), Geschwindigkeit (u' = v) und Beschleunigung (u'' = w)
    // mit Hilfe der Newmark-Vorschrift aktualisieren.
    // BEACHTE: Wir haben u als primaere Variable (die Variable nach der das GLS geloest wird),
    // anstatt klassischerweise u''. Fuer die Umformungen siehe MA.
    void updateSolutionNewmarkPreviousStep(double dt, double beta, double gamma);

    BlockMultiVectorPtr_Type getSolutionPreviousTimestep();

    BlockMultiVectorPtrArray_Type getSolutionAllPreviousTimestep();

    BlockMatrixPtr_Type getMassSystem();

    ParameterListPtr_Type getParameterList();

    void assembleSourceTerm( double time=0. );

    BlockMultiVectorPtr_Type getSourceTerm( );

    bool hasSourceTerm() const;

    CommConstPtr_Type getComm() const{return  comm_;}

    LinSolverBuilderPtr_Type getLinearSolverBuilder() const;

    void getValuesOfInterest( vec_dbl_Type& values );

    void computeValuesOfInterestAndExport();

    void updateTime( double time ){ time_ = time;}

    void addToRhs(BlockMultiVectorPtr_Type x);
    
    /// @brief Plotting residual Vector of problem
    /// @param problem Timeproblem
    /// @param time current time
    void plotLinResVec(double time =0.) const;

    void initExporterResidual() const;

    mutable ExporterPtrVec_Type exporterResidual_;
    mutable double currentTimeExport_=0.;
    mutable double timeStep_ =0;
    mutable double newtonStep_=0;
    

    ProblemPtr_Type problem_;
    CommConstPtr_Type comm_;
    
    // Exporter for various parts of the solution or vectors that are needed for restarts
    Teuchos::RCP <HDF5Export<SC,LO,GO,NO>> HDF5exporterDsVelocity_; // Verlocity for Newmark
    Teuchos::RCP <HDF5Export<SC,LO,GO,NO>> HDF5exporterDsAcceleration_; // Acceleration for Newmark
    Teuchos::RCP <HDF5Export<SC,LO,GO,NO>> HDF5exporterSolutionNewmark_; // Acceleration for Newmark
    
    std::vector<Teuchos::RCP <HDF5Export<SC,LO,GO,NO>>> HDF5exporterRhs_; // Acceleration for Newmark
    std::vector<Teuchos::RCP <HDF5Export<SC,LO,GO,NO>>> HDF5exporterHistory_; // Solution displacement
    std::vector<Teuchos::RCP <HDF5Export<SC,LO,GO,NO>>> HDF5exporterSolution_; // Solution displacement


    mutable BlockMatrixPtr_Type systemCombined_;
    mutable BlockMatrixPtr_Type systemMass_;
    mutable SmallMatrix<double> timeParameters_;
    SmallMatrix<int>        timeStepDef_;
    SmallMatrix<double>     massParameters_;

    FEFacPtr_Type feFactory_;
//    bool					boolLinearProblem_;
    int                     dimension_;
    bool                    verbose_;
    ParameterListPtr_Type	parameterList_;
    mutable BCConstPtr_Type bcFactory_;
    bool                    massBCSet_;
    BlockMultiVectorPtrArray_Type solutionPreviousTimesteps_;
//    std::vector<MultiVector_ptr_Type>    solutionPreviousTimesteps_;

    // ###########################
    // Fuer das Newmark-Verfahren
    // ###########################
    BlockMultiVectorPtrArray_Type velocityPreviousTimesteps_; // entspricht u' bzw. v
    BlockMultiVectorPtrArray_Type accelerationPreviousTimesteps_; // entspricht u'' bzw. w

    // ###########################
    // Fuer FSI
    // ###########################
    BlockMatrixPtrArray_Type systemMassPreviousTimeSteps_;
    //#################
    void checkForExportAndExport(BlockMultiVectorPtrArray_Type solutionVec, std::string fileName);


    double time_;
protected:
#define TIMEPROBLEM_TIMER
#ifdef TIMEPROBLEM_TIMER
    TimePtr_Type TimeSolveTimer_;
#endif
    
    // Functions used for NOX.
public:
    
    virtual Thyra::ModelEvaluatorBase::InArgs<SC> getNominalValues() const;
    
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > get_x_space() const;
    
    virtual Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > get_f_space() const;
    
    virtual ::Thyra::ModelEvaluatorBase::InArgs<SC> createInArgs() const;
    
    virtual ::Thyra::ModelEvaluatorBase::OutArgs<SC> createOutArgsImpl() const;
    
//    void initNOXParameters();
//    
//    void initVectorSpaces( );
//    
//    void initVectorSpacesMonolithic( );
//    
//    void initVectorSpacesBlock( );
    
    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() ;

    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Monolithic() ;

    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Block() ;

    Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() ;
    
    std::string description() const; //reimplement description function and use description of underlying nonLinearProblem
private:
    
    virtual void evalModelImpl(
                               const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                               const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                               ) const;
    
    void evalModelImplMonolithic(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                 const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;
    
    void evalModelImplBlock(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                            const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;
    
    mutable bool precInitOnly_; //Help variable to signal that we constructed the initial preconditioner for NOX with the Stokes system and we do not need to compute it if fill_W_prec is called for the first time. However, the preconditioner is only correct if a Stokes system is solved in the first nonlinear iteration. This only affects the block preconditioners of Teko

    void initExporter(std::string fileName  );
    void initCheckPoints();
    double getPreviousTimeIncrement(double timeStep =-1.0);
    Teuchos::RCP<HDF5Export<SC, LO, GO, NO>> getExporter(std::string fileName, int i);

    std::vector<std::tuple<double,bool>> checkPointTupel_;
    std::set<std::string> exportsWritten_; // files and times (or checkpoints) written by checkForExportAndExport()
    bool restartValuesImported_ = false; // importRestartValues() ran
    bool restartNewmarkUpdate_ = false; // after a restart: the imported Newmark state is that of t_{r-1} and is updated (see updateSolutionNewmarkPreviousStep)
    BlockMultiVectorPtrArray_Type restartMassSolutions_; // after a restart: the products M_i u_i of the checkpoint (see updateMultistepRhsFSI)
    int timeStepsSinceRestart_ = 0; // time steps of updateMultistepRhsFSI after the one that read restartMassSolutions_
    // State kept by saveState()
    BlockMultiVectorPtr_Type savedSolution_;
    BlockMultiVectorPtrArray_Type savedSolutionPreviousTimesteps_;
    BlockMultiVectorPtrArray_Type savedVelocityPreviousTimesteps_;
    BlockMultiVectorPtrArray_Type savedAccelerationPreviousTimesteps_;
    double savedTime_ = 0.;

};
}

#endif
