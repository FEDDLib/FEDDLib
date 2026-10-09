#ifndef AssembleFE_SCI_SMC_CMM_Active_Growth_Reorientation_DECL_hpp
#define AssembleFE_SCI_SMC_CMM_Active_Growth_Reorientation_DECL_hpp

#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/AssembleFEBlock.hpp"
#include "feddlib/core/FE/Helper.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/TimeSteppingTools.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#ifdef FEDD_HAVE_ACEGENINTERFACE
#include "aceinterface.hpp"
#endif

/*!
\class AssembleFE_SCI_SMC_CMM_Active_Growth_Reorientation
        Coupled deformation diffusion problem with the constrained-mixture (CMM) smooth-muscle
        model with active response, growth and reorientation (AceGen element
        DeformationDiffusionConstrainedMixtureModelSmoothMuscleActiveGrowthReorientationTetrahedra3D10).
        Mirrors AssembleFE_SCI_SMC_Active_Growth_Reorientation; differs in the history layout
        (39 values per Gauss point) and in the domain-data names.
        Derived from AssembleFE base class
*/

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFE_SCI_SMC_CMM_Active_Growth_Reorientation : public AssembleFE<SC, LO, GO, NO> {
   public:
    typedef Matrix<SC, LO, GO, NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

    typedef MultiVector<SC, LO, GO, NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

    typedef AssembleFE<SC, LO, GO, NO> AssembleFE_Type;

    /*!
    \brief Assemble the element Jacobian matrix.
    */
    void assembleJacobian() override;

    /*!
    \brief Assemble the element right hand side vector.
    */
    void assembleRHS() override;

    /*!
    \brief Assemble block parts of the element Jacobian matrix.
    \return the element Jacobian matrix of block i
    */
    virtual void assembleJacobianBlock(LO i) {};

    void advanceInTime(Teuchos::RCP<TimeSteppingTools> timeSteppingTool) override;

    void synchronizeTime(Teuchos::RCP<TimeSteppingTools> timeSteppingTool) override;

    /// Adaptive time stepping: the base state, the concentrations of the time step and the
    /// switches (and one-time initializations) of the active, growth and reorientation segments
    void saveState() override;
    void restoreState() override;

    void postProcessing() override;

    void getMassMatrix(SmallMatrixPtr_Type& massMatrix) { massMatrix = massMatrix_; };

    void initializeGrowth();

    void initializeActiveResponse();

    void updateDomainData(const std::string& dataName, double dataValue);

    std::vector<std::string> getPostDataNames() { return postDataNames_; }
    std::map<std::string, int> getFieldNameToPosition() { return fieldNameToPosition_; }

    int getNumberOfIntegrationPoints() { return numberOfIntegrationPoints_; }

   protected:
    AssembleFE_SCI_SMC_CMM_Active_Growth_Reorientation(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple);

   private:
    void checkingReorientationActiveGrowth();

    void assemble_SCI_SMC_CMM_Active_Growth_Reorientation(bool computeTangent);

    /*!
    \brief Domain data with the active-acceleration factors applied (if time < "Accelerated Active Until").
    */
    std::vector<double> modifiedDomainData(double time) const;

    friend class AssembleFEFactory<SC, LO, GO, NO>;  // Must have for specfic classes

    int findPosition(const std::string& subString, const std::vector<std::string>& stringArray);

    std::string FEType_;  // FEType of Disk

    SmallMatrixPtr_Type massMatrix_;

    int dofsSolid_;  // Degrees of freedom per node
    int dofsChem_;
    int numNodesSolid_;  // Number of nodes of element
    int numNodesChem_;   // Number of nodes of element

    int dofsElement_;  // "Dimension of return matrix"

    int iCode_;                      // Integration Code
    int numberOfIntegrationPoints_;  // Number of integration points
    int postDataLength_;             // Number of post processing variables
    int domainDataLength_;           // Number of domain data parameters
    double subiterationTolerance_;   // Tolerance for element level NR

    vec_dbl_Type solutionC_n_;
    vec_dbl_Type solutionC_n1_;

    vec_dbl_Type domainData_;
    vec_dbl_Type positions_;
    vec_dbl_Type displacements_;
    vec_dbl_Type accelerations_;
    vec_dbl_Type concentrations_;
    vec_dbl_Type rates_;

    std::vector<std::string> domainDataNames_;
    std::vector<std::string> postDataNames_;
    std::map<std::string, int> fieldNameToPosition_;

    vec_dbl_Type residuumRint_;
    vec_dbl_Type residuumRdyn_;
    vec_dbl_Type residuumRc_;

    vec2D_dbl_Type stiffnessMatrixKuu_;
    vec2D_dbl_Type stiffnessMatrixKuc_;
    vec2D_dbl_Type stiffnessMatrixKcu_;
    vec2D_dbl_Type stiffnessMatrixKcc_;
    vec2D_dbl_Type massMatrixMc_;

    // Timeintervals for Active Response and Reorientation
    vec2D_dbl_Type segmentsActive_;
    vec2D_dbl_Type segmentsGrowth_;
    vec2D_dbl_Type segmentsReorientation_;

    int activeBool_ = 0;
    int growthBool_ = 0;
    int reorientationBool_ = 0;

    // Growth-type switches (e.g. GroundGrowthBool, SMCGrowthBool) driven by the
    // "Timestepping Intervalls Growth" intervals ("Flags" parameter).
    std::vector<std::string> growthFlags_;
    double activeAcceleratedEndTime_ = 0.;
    double activeAcceleratedMultiplier_ = 1.;

    bool activeInitialized_ = false;
    bool growthInitialized_ = false;
    // State kept by saveState(). The element data of the last computation are part of it:
    // initializeActiveResponse() and initializeGrowth() compute with them.
    vec_dbl_Type savedDisplacements_;
    vec_dbl_Type savedConcentrations_;
    vec_dbl_Type savedRates_;
    vec_dbl_Type savedAccelerations_;
    vec_dbl_Type savedSolutionC_n_;
    vec_dbl_Type savedSolutionC_n1_;
    vec_dbl_Type savedDomainData_;
    int savedActiveBool_ = 0;
    int savedGrowthBool_ = 0;
    int savedReorientationBool_ = 0;
    bool savedActiveInitialized_ = false;
    bool savedGrowthInitialized_ = false;

    // Pre-computed indices for domain data modification
    std::vector<int> acceleratedParamIndices_;  // Parameters to multiply
    std::vector<int> deceleratedParamIndices_;  // Parameters to divide

#ifdef FEDD_HAVE_ACEGENINTERFACE
    typedef AceGenInterface::DeformationDiffusionConstrainedMixtureModelSmoothMuscleActiveGrowthReorientationTetrahedra3D10 AceGenElement_Type;

    /*!
    \brief The AceGen element of this element, given the current data (created on the first call).

    Constructing an AceGen element allocates its result arrays and sets up the element
    specification, so the element keeps one and only passes it the data of each use.
    */
    AceGenElement_Type& aceGenElement(double* displacements, double* concentrations, double* accelerations, double* rates, double* domainData, double deltaT, double time);

    Teuchos::RCP<AceGenElement_Type> aceGenElement_;
#endif
};

}  // namespace FEDD
#endif  // AssembleFE_SCI_SMC_CMM_Active_Growth_Reorientation_DECL_hpp
