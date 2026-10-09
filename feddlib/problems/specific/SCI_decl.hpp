#ifndef SCI_decl_hpp
#define SCI_decl_hpp
#include "feddlib/problems/abstract/TimeProblem.hpp"
#include "feddlib/problems/specific/DiffusionReaction.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/core/General/ExporterTxt.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/core/General/TimeSteppingTools.hpp"
#include "Xpetra_ThyraUtils.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_ModelEvaluatorBase_decl.hpp>
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
    
    
namespace FEDD{

template <class SC , class LO , class GO , class NO >
class TimeProblem;
template <class SC , class LO , class GO , class NO >
class DiffusionReaction;
template <class SC , class LO , class GO , class NO >
class LinElas;
template <class SC , class LO , class GO , class NO >
class NonLinElasticity;
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class SCI : public NonLinearProblem<SC,LO,GO,NO>  {

public:
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef typename Problem_Type::Matrix_Type Matrix_Type;
    typedef typename Problem_Type::MatrixPtr_Type MatrixPtr_Type;

    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;

    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;
    typedef typename Problem_Type::BlockMultiVector_Type BlockMultiVector_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

    typedef typename Problem_Type::Domain_Type Domain_Type;
    typedef Teuchos::RCP<Domain_Type > DomainPtr_Type;
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    
    typedef typename Problem_Type::Domain_Type::Mesh_Type Mesh_Type;
    typedef typename Problem_Type::Domain_Type::MeshPtr_Type MeshPtr_Type;
    
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef NonLinearProblem<SC,LO,GO,NO> NonLinearProblem_Type;
    
    typedef typename NonLinearProblem_Type::BlockMultiVectorPtrArray_Type BlockMultiVectorPtrArray_Type;

    typedef TimeProblem<SC,LO,GO,NO> TimeProblem_Type;
    typedef Teuchos::RCP<TimeProblem_Type> TimeProblemPtr_Type;

    typedef DiffusionReaction<SC,LO,GO,NO> ChemProblem_Type;
    typedef LinElas<SC,LO,GO,NO> StructureProblem_Type;
    typedef NonLinElasticity<SC,LO,GO,NO> StructureNonLinProblem_Type;

    typedef Teuchos::RCP<ChemProblem_Type> ChemProblemPtr_Type;
    typedef Teuchos::RCP<StructureProblem_Type> StructureProblemPtr_Type;
    typedef Teuchos::RCP<StructureNonLinProblem_Type> StructureNonLinProblemPtr_Type;

    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef typename Problem_Type::BC_Type BC_Type;
    typedef typename Teuchos::RCP<BC_Type> BCPtr_Type;
    
    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef Teuchos::RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    
    typedef ExporterParaView<SC,LO,GO,NO> Exporter_Type;
    typedef Teuchos::RCP<Exporter_Type> ExporterPtr_Type;
    typedef Teuchos::RCP<ExporterTxt> ExporterTxtPtr_Type;
    
    typedef std::vector<GO> vec_GO_Type;
    typedef std::vector<vec_GO_Type> vec2D_GO_Type;
    typedef std::vector<vec2D_GO_Type> vec3D_GO_Type;
    typedef Teuchos::RCP<vec3D_GO_Type> vec3D_GO_ptr_Type;

    // FETypeVelocity muss gleich FETypeStructure sein, wegen Interface.
    // Zudem wird FETypeVelocity auch fuer das Geometrieproblem genutzt.
    SCI(const DomainConstPtr_Type &domainStructure, std::string FETypeStructure,
					const DomainConstPtr_Type &domainChem, std::string FETypeChem,vec2D_dbl_Type diffusionTensor, RhsFunc_Type reactionFunc,
                    ParameterListPtr_Type parameterListStructure, ParameterListPtr_Type parameterListChem,
                    ParameterListPtr_Type parameterListSCI, Teuchos::RCP<SmallMatrix<int> > &defTS);

    ~SCI();

    virtual void info();

    virtual void assemble( std::string type = "" ) const;
    
    // init FSI vectors from partial problems
    void setFromPartialVectorsInit() const;
    
    // Setze die aktuelle Loesung als vergangene Loesung
    void updateMeshDisplacement() const;

    // Berechnet die Massematrix und die daraus resultierende rechte Seite nach BDF2
    // void getFluidMassmatrixAndRHSInTime(BlockMatrixPtr_Type massmatrix, BlockMultiVectorPtr_Type rhs) const;

    // Berechne die Massematrix fuer das FluidProblem.
    void setChemMassmatrix(MatrixPtr_Type& massmatrix) const;

    // Berechnet die Massematrix und die daraus resultierende rechte Seite nach Newmark
    // und macht direkt ein Update. Dies koennen wir bei Struktur machen, da Massematrix
    // innerhalb einer Zeitschleife konstant ist.
    void setSolidMassmatrix( MatrixPtr_Type& massmatrix ) const;

    void computeSolidRHSInTime() const;
    
    void computeChemRHSInTime() const;

    void updateTime() const;

    // Hier wird im Prinzip updateSolution() fuer problemTimeChem_ aufgerufen
    void updateChemInTime() const;

    // Solving chemistry component of problem in case of explicit chemistry
    void solveChemistryProblem() const;

    // Verschiebt die notwendigen Gitter
    void moveMesh() const;

    void initializeCE();
    // Macht setupTimeStepping() auf problemTimeFluid_ und problemTimeStructure_
    void setupSubTimeProblems(ParameterListPtr_Type parameterListFluid, ParameterListPtr_Type parameterListStructure) const;

    void setBoundariesSubProblems() const;
    ChemProblemPtr_Type getChemProblem(){
        return problemChem_;
    }
    
    StructureProblemPtr_Type getStructureProblem(){
        return problemStructure_;
    }
    
    // Berechnet von einer dofID, d.h. dim*nodeID+(0,1,2), die entsprechende nodeID.
    // IN localDofNumber steht dann, ob es die x- (=0), y- (=1) oder z-Komponente (=2) ist.
    void toNodeID(UN dim, GO dofID, GO& nodeID, LO& localDofNumber ) const
    {
        nodeID = (GO) (dofID/dim);
        localDofNumber = (LO) (dofID%dim);
    }

    // Diese Funktion berechnet genau das umgekehrte. Also von einer nodeID die entsprechende dofID
    void toDofID(UN dim, GO nodeID, LO localDofNumber, GO& dofID)  const
    {
        dofID = (GO) ( dim * nodeID + localDofNumber);
    }
    
    void getValuesOfInterest2DBenchmark( vec_dbl_Type& values );

    void getValuesOfInterest3DBenchmark( vec_dbl_Type& values );
    
    virtual void getValuesOfInterest( vec_dbl_Type& values ) {}  ;

    virtual void getValuesOfInterest( BlockMultiVectorPtr_Type& values );

    virtual void exportValuesOfInterest(double time);
    virtual void importValuesOfInterest(double time);

    virtual void computeValuesOfInterestAndExport() {} ;

    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const{};
    
    virtual void reAssemble(std::string type) const;

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions) {};

    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const; //standard or reverse    
    
    BlockMultiVectorPtr_Type getPostProcessingData();

    vec_string_Type getPostprocessingNames();

    MultiVectorPtr_Type getHistoryData();
    /*####################*/

    // Alternativ wie in reAssembleExtrapolation() in NS?

    MultiVectorPtr_Type meshDisplacementOld_rep_;
    MultiVectorPtr_Type meshDisplacementNew_rep_;
    MultiVectorPtr_Type c_rep_;
    MultiVectorPtr_Type d_rep_;

    // stationaere Systeme
    ChemProblemPtr_Type problemChem_;
    StructureProblemPtr_Type problemStructure_;
    StructureNonLinProblemPtr_Type problemStructureNonLin_; // CH: we want to combine both structure models to one general model later

    // zeitabhaengige Systeme
    mutable TimeProblemPtr_Type problemTimeChem_;
    mutable TimeProblemPtr_Type problemTimeStructure_;

    Teuchos::RCP<SmallMatrix<int>> defTS_;
    mutable Teuchos::RCP<TimeSteppingTools>	timeSteppingTool_;

    // Sets the time and time increment of the elements to those of timeSteppingTool_
    void synchronizeElementTime() const { this->feFactory_->synchronizeTime(timeSteppingTool_); }

    // Adaptive time stepping (DAESolverInTime::advanceInTimeSCI): the state at the start of a time
    // step, to repeat a failed one from, and the failures of the elements
    void saveStepState() const;
    void restoreStepState() const;
    /// With adaptive true an element that cannot compute its state records it, and the residual
    /// becomes NaN on every process (the Newton iteration fails) instead of the run stopping. With
    /// acceptElementFailures true as well the residual is left as it is: the time step stands or falls
    /// with the Newton iteration alone.
    void setAdaptiveStep(bool adaptive, bool acceptElementFailures = false) const;
    /// Whether an element of any process failed since the last restoreStepState() or setAdaptiveStep()
    bool elementFailed() const { return elementFailed_; }
    /// The number of elements (of all processes) that failed, as of the last residual
    int numberOfFailedElements() const { return failedElements_; }

private:
    std::string materialModel_;
    vec_dbl_Type valuesForExport_;
    vec_string_Type postProcessingnames_;
    std::vector<std::string> requestedFields_; // The post processing fields that are requested as per the parameters file
    bool geometryExplicit_;
    mutable BlockMatrixPtr_Type systemC_;
    ExporterTxtPtr_Type exporterIterationsChem_;
    mutable ExporterPtr_Type exporterEMod_;
    mutable ExporterPtr_Type exporterChem_;
    mutable bool exportedEMod_ ;
    mutable bool setUpTimeStep_;
    mutable MultiVectorPtr_Type eModVec_;
    bool loadStepping_;
    mutable bool solidMassBuilt_ = false; // the structure mass matrix (setSolidMassmatrix) is built once
    mutable bool adaptiveStep_ = false; // see setAdaptiveStep()
    mutable bool acceptElementFailures_ = false;
    mutable bool elementFailed_ = false;
    mutable int failedElements_ = 0;
    mutable Teuchos::RCP<TimeSteppingTools> savedTimeSteppingTool_; // state kept by saveStepState()
    bool chemistryExplicit_;
    bool externalForce_;
    bool nonlinearExternalForce_;

    /*####################*/

public:
        // NOX and FSI only implement in combination with TimeProblem

private:
    
    virtual void evalModelImpl(
                               const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                               const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                               ) const;
    
//    void evalModelImplMonolithic(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
//                                 const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;
//    
//#ifdef FEDD_HAVE_TEKO
//    void evalModelImplBlock(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
//                            const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;
//#endif
};
}
#endif
