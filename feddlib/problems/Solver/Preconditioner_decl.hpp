#ifndef Preconditioner_DECL_hpp
#define Preconditioner_DECL_hpp

#include "feddlib/problems/problems_config.h"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/problems/abstract/MinPrecProblem.hpp"
#include "feddlib/problems/abstract/Problem.hpp"
#include "feddlib/problems/abstract/TimeProblem.hpp"
#include "feddlib/problems/Solver/PrecOpFaCSI.hpp"
#include "feddlib/problems/Solver/PrecBlock2x2.hpp"
#include "feddlib/problems/specific/LaplaceBlocks.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "Xpetra_ThyraUtils.hpp"
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_DefaultPreconditioner_decl.hpp>
#include <Stratimikos_FROSch_def.hpp>
//#include <Stratimikos_FROSch_def.hpp> // hier werden schon alle FROSch VK eingebunden
#ifdef FEDD_HAVE_TEKO
#include "Teko_StratimikosFactory.hpp"
#include <Teko_StaticRequestCallback.hpp>
#endif

/*!
 Declaration of Preconditioner

 @brief  Preconditioner
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC , class LO , class GO , class NO >
class Problem;
template <class SC , class LO , class GO , class NO >
class TimeProblem;
template <class SC , class LO , class GO , class NO >
class MinPrecProblem;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Preconditioner {

public:

    typedef Teuchos::RCP<Teuchos::Comm<int> > CommPtr_Type;
    typedef Teuchos::RCP<const Teuchos::Comm<int> > CommConstPtr_Type;

    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef TimeProblem<SC,LO,GO,NO> TimeProblem_Type;
    typedef MinPrecProblem<SC,LO,GO,NO> MinPrecProblem_Type;
    typedef Teuchos::RCP<Problem_Type> ProblemPtr_Type;
    typedef Teuchos::RCP<TimeProblem_Type> TimeProblemPtr_Type;
    typedef Teuchos::RCP<MinPrecProblem_Type> MinPrecProblemPtr_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map<LO,GO,NO> > MapPtr_Type;
    typedef Teuchos::RCP<const Map<LO,GO,NO> > MapConstPtr_Type;
    typedef Teuchos::RCP<Thyra::PreconditionerBase<SC> > ThyraPrecPtr_Type;
    typedef Teuchos::RCP<const Thyra::PreconditionerBase<SC> > ThyraPrecConstPtr_Type;
    typedef Teuchos::RCP<Thyra::LinearOpBase<SC> > ThyraLinOpPtr_Type;
    typedef Teuchos::RCP<const Thyra::LinearOpBase<SC> > ThyraLinOpConstPtr_Type;
    typedef Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder >  LinSolverBuilderPtr_Type;

    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<Mesh_Type> MeshPtr_Type;
    typedef Teuchos::RCP<const Mesh_Type> MeshConstPtr_Type;

    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;

    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;

    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;

    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;

    typedef typename Problem_Type::BC_Type BC_Type;
    typedef typename Problem_Type::BCConstPtr_Type BCConstPtr_Type;

    Preconditioner(Problem_Type* problem);

    Preconditioner(TimeProblem_Type* problem);

    ~Preconditioner();

    ThyraPrecPtr_Type getThyraPrec();

    ThyraPrecConstPtr_Type getThyraPrecConst() const;

    void setPreconditionerThyraFromLinOp( ThyraLinOpPtr_Type precLinOp );

    void buildPreconditioner( std::string type="Monolithic" );

    void initializePreconditioner( std::string type="Monolithic" );

    void buildPreconditionerMonolithic( );

    void buildPreconditionerMonolithicFSI( );

    void initPreconditionerMonolithic( );

    void initPreconditionerBlock( );

#ifdef FEDD_HAVE_TEKO
    void buildPreconditionerTeko( );

    void setVelocityParameters( ParameterListPtr_Type parameterList, int coarseRanks );

    void setPressureParameters( ParameterListPtr_Type parameterList, int coarseRanks );

    ThyraLinOpConstPtr_Type getTekoOp();

    void setVelocityMassMatrix(MatrixPtr_Type massMatrix) const;
    MatrixPtr_Type getVelocityMassMatrix(){return velocityMassMatrixMatrixPtr_;};

    void setPressureLaplaceMatrix(MatrixPtr_Type matrix) const;
    MatrixPtr_Type getPressureLaplaceMatrix(){return pressureLaplaceMatrixPtr_;};

    void setPCDOperator(MatrixPtr_Type matrix) const;
    MatrixPtr_Type getPCDOperatorMatrix(){return pcdOperatorMatrixPtr_;};
#endif

    void setPressureMassMatrix(MatrixPtr_Type massMatrix) const;
    MatrixPtr_Type getPressureMassMatrix(){return pressureMassMatrixPtr_;};

    void buildPreconditionerFaCSI( std::string type );

    void buildPreconditionerBlock2x2();

    void setFaCSIBCFactory( BCConstPtr_Type bcFactory ){ faCSIBCFactory_ = bcFactory; };

    bool hasFaCSIBCFactory(){ return !faCSIBCFactory_.is_null(); };

    BCConstPtr_Type getFaCSIBCFactory( ){ return faCSIBCFactory_; };

    void exportCoarseBasis( );

    void exportCoarseBasisFSI( );

    bool isPreconditionerComputed() const{return precondtionerIsBuilt_;};


private:
    ThyraPrecPtr_Type thyraPrec_;
    bool precondtionerIsBuilt_;
    ProblemPtr_Type problem_;
    TimeProblemPtr_Type timeProblem_;
    Teuchos::RCP<Thyra::PreconditionerFactoryBase<SC> > precFactory_;
#ifdef FEDD_HAVE_TEKO
    ThyraLinOpConstPtr_Type tekoLinOp_;
    Teuchos::RCP<Teko::RequestHandler> rh_;
    Teuchos::RCP< Teko::StaticRequestCallback<Teko::LinearOp> > callbackPCD_;
#endif
    // For FaCSI precondtioner
    ThyraLinOpConstPtr_Type fsiLinOp_;
    ThyraLinOpPtr_Type precFluid_;
    ThyraLinOpPtr_Type precStruct_;
    ThyraLinOpPtr_Type precGeo_;
    MinPrecProblemPtr_Type probFluid_;
    MinPrecProblemPtr_Type probSolid_;
    MinPrecProblemPtr_Type probGeo_;
    BCConstPtr_Type faCSIBCFactory_;

    // For Stokes-type block diagonal and triangular precondtioner
    ThyraLinOpPtr_Type precVelocity_;
    ThyraLinOpPtr_Type precSchur_;
    MinPrecProblemPtr_Type probVelocity_;
    MinPrecProblemPtr_Type probSchur_;
    // mutable MatrixPtr_Type pressureMassMatrix_;

    // For LSC and PCD preconditioner
    mutable ThyraLinOpConstPtr_Type velocityMassMatrix_; // LSC
    mutable MatrixPtr_Type velocityMassMatrixMatrixPtr_; // LSC

    mutable MatrixPtr_Type pressureLaplaceMatrixPtr_; // PCD
    mutable ThyraLinOpConstPtr_Type pressureLaplace_; // PCD

    mutable MatrixPtr_Type pressureMassMatrixPtr_; // PCD
    mutable ThyraLinOpConstPtr_Type pressureMass_; // PCD
    
    mutable MatrixPtr_Type pcdOperatorMatrixPtr_; // PCD
    mutable ThyraLinOpConstPtr_Type pcdOperator_; // PCD

    // For construction in the FEDDLib
    MinPrecProblemPtr_Type probLaplace_;
    MinPrecProblemPtr_Type probMass_;
    MinPrecProblemPtr_Type probVMass_;
    ThyraLinOpPtr_Type laplaceInverse_;
    ThyraLinOpPtr_Type massMatrixInverse_;
    ThyraLinOpPtr_Type massMatrixVInverse_;  

    ParameterListPtr_Type pListPhiExport_;
#define PRECONDITIONER_TIMER
#ifdef PRECONDITIONER_TIMER
    TimePtr_Type preconditionerTimer_;
#endif
};
}
#endif
