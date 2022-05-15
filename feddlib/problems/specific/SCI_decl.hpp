#ifndef SCI_decl_hpp
#define SCI_decl_hpp
#include "feddlib/problems/abstract/TimeProblem.hpp"
#include "feddlib/problems/specific/DiffusionReaction.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NonLinElasAssFE.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/problems/Solver/TimeSteppingTools.hpp"
#include "Xpetra_ThyraUtils.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_ModelEvaluatorBase_decl.hpp>
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
    typedef NonLinElasAssFE<SC,LO,GO,NO> StructureNonLinProblem_Type;

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

    // Hier wird timeSteppingTool_->t_ inkrementiert
    void updateTime() const;

    // Hier wird im Prinzip updateSolution() fuer problemTimeChem_ aufgerufen
    void updateChemInTime() const;

    // Verschiebt die notwendigen Gitter
    void moveMesh() const;

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

    virtual void computeValuesOfInterestAndExport() {} ;

    virtual void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const{};
    
    virtual void reAssemble(std::string type) const;

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions) {};

    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const; //standard or reverse    
    
    /*####################*/

    // Alternativ wie in reAssembleExtrapolation() in NS?

    MultiVectorPtr_Type meshDisplacementOld_rep_;
    MultiVectorPtr_Type meshDisplacementNew_rep_;
    MultiVectorPtr_Type c_rep_;

    mutable MatrixPtr_Type C2_;

    mutable MatrixPtr_Type 	P_;
    mutable int counterP;
    // stationaere Systeme
    ChemProblemPtr_Type problemChem_;
    StructureProblemPtr_Type problemStructure_;
    StructureNonLinProblemPtr_Type problemStructureNonLin_; // CH: we want to combine both structure models to one general model later

    // zeitabhaengige Systeme
    mutable TimeProblemPtr_Type problemTimeChem_;
    mutable TimeProblemPtr_Type problemTimeStructure_;

    Teuchos::RCP<SmallMatrix<int>> defTS_;
    mutable Teuchos::RCP<TimeSteppingTools>	timeSteppingTool_;

private:
    std::string materialModel_;
    vec_dbl_Type valuesForExport_;
    bool geometryExplicit_;
    ExporterTxtPtr_Type exporterTxtDrag_;
    ExporterTxtPtr_Type exporterTxtLift_;
    mutable ExporterPtr_Type exporterEMod_;
    mutable bool exportedEMod_ ;
    mutable bool setUpTimeStep_;
    mutable MultiVectorPtr_Type eModVec_;
    /*####################*/

public:
        // NOX and FSI only implement in combination with TimeProblem

//    typedef Thyra::VectorSpaceBase<SC> thyra_vec_space;
//    typedef Thyra::VectorBase<SC> thyra_vec;
//    typedef Tpetra::Map<LO, GO, NO> tpetra_map;
//    typedef Tpetra::CrsMatrix<SC, LO, GO, NO> tpetra_matrix;
//    typedef Thyra::LinearOpBase<SC> thyra_op;
//    typedef Tpetra::Operator<SC,LO,GO,NO> tpetra_op;
//

//    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op() const;
//    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Monolithic() const;
//#ifdef FEDD_HAVE_TEKO
//    Teuchos::RCP< Thyra::LinearOpBase<SC> > create_W_op_Block() const;
//#endif
//    Teuchos::RCP<Thyra::PreconditionerBase<SC> > create_W_prec() const;
    
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
