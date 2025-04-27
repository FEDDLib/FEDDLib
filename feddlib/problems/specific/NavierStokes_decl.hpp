#ifndef NAVIERSTOKES_decl_hpp
#define NAVIERSTOKES_decl_hpp
#include "feddlib/problems/abstract/NonLinearProblem.hpp"
#include <Xpetra_ThyraUtils.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Thyra_ProductVectorBase.hpp>
#include <Thyra_PreconditionerBase.hpp>
#include <Thyra_ModelEvaluatorBase_decl.hpp>
/*!
 Declaration of Navier-Stokes

 @brief Navier-Stokes
 @authors Christian Hochmuth, Lea Saßmannshausen
 @version 1.0
 @copyright CH
 */

namespace FEDD{


template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class NavierStokes : public NonLinearProblem<SC,LO,GO,NO>  {

public:
    //! @name Public Types
    //@{
    typedef Problem<SC,LO,GO,NO> Problem_Type;
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
    typedef Teuchos::RCP<const MatrixPtr_Type> MatrixConstPtr_Type;

    typedef typename Problem_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef typename Problem_Type::BlockMatrix_Type BlockMatrix_Type;
    typedef typename Problem_Type::BlockMatrixPtr_Type BlockMatrixPtr_Type;

    typedef typename Problem_Type::MultiVector_Type MultiVector_Type;
    typedef typename Problem_Type::MultiVectorPtr_Type MultiVectorPtr_Type;
    typedef typename Problem_Type::MultiVectorConstPtr_Type MultiVectorConstPtr_Type;
    typedef typename Problem_Type::BlockMultiVectorPtr_Type BlockMultiVectorPtr_Type;

    typedef typename Problem_Type::Domain_Type Domain_Type;
    typedef typename Problem_Type::DomainPtr_Type DomainPtr_Type;
    typedef typename Problem_Type::DomainConstPtr_Type DomainConstPtr_Type;
    typedef typename Problem_Type::CommConstPtr_Type CommConstPtr_Type;

    typedef NonLinearProblem<SC,LO,GO,NO> NonLinearProblem_Type;
    typedef typename NonLinearProblem_Type::BlockMultiVectorPtrArray_Type BlockMultiVectorPtrArray_Type;

    typedef typename NonLinearProblem_Type::TpetraMatrix_Type TpetraMatrix_Type;

    typedef typename NonLinearProblem_Type::ThyraVecSpace_Type ThyraVecSpace_Type;
    typedef typename NonLinearProblem_Type::ThyraVec_Type ThyraVec_Type;
    typedef typename NonLinearProblem_Type::ThyraOp_Type ThyraOp_Type;
    typedef Thyra::BlockedLinearOpBase<SC> ThyraBlockOp_Type;

    typedef BCBuilder<SC,LO,GO,NO> BC_Type;
    typedef Teuchos::RCP<BC_Type> BCPtr_Type;
    typedef Teuchos::RCP<const BC_Type> BCConstPtr_Type;

    typedef typename NonLinearProblem_Type::TpetraOp_Type TpetraOp_Type;
    //@}

    //! @name Constructor/Destructor
    //@{
    NavierStokes( const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity, const DomainConstPtr_Type &domainPressure, std::string FETypePressure, ParameterListPtr_Type parameterList );
    //@}

    virtual void info();

    virtual void assemble( std::string type = "" ) const;
    
    void assembleConstantMatrices() const;
    
    void assembleDivAndStab() const;

    void updateConvectionDiffusionOperator() const;
    
    void reAssemble( std::string type ) const;
    
    void reAssemble( BlockMultiVectorPtr_Type previousSolution ) const override{}
    
    void reAssembleFSI(std::string type, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const;
    
    virtual void reAssemble(MatrixPtr_Type& massmatrix, std::string type ) const;

    virtual void reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions);

    virtual void calculateNonLinResidualVec(std::string type="standard", double time=0.) const; //standard or reverse
    
    void calculateNonLinResidualVecWithMeshVelo(std::string type, double time, MultiVectorPtr_Type u_minus_w, MatrixPtr_Type P) const;
    // virtual int ComputeDragLift(vec_dbl_ptr_Type &values);

    void getValuesOfInterest( vec_dbl_Type& values ) override {}
    
    void computeValuesOfInterestAndExport() override {}

//    virtual void assembleExternal( std::string type ){}
    /*####################*/

    mutable MatrixPtr_Type 	A_;
    vec_int_ptr_Type pressureIDsLoc;
    MultiVectorPtr_Type u_rep_;

    BCPtr_Type bcFactoryPCD_;
    mutable MatrixPtr_Type 	Mp_;
    mutable MatrixPtr_Type 	BT_Mp_;
    mutable MatrixPtr_Type BT_Mp_B_;
    mutable MatrixPtr_Type 	Ap_;

    bool augmentedLagrange_=false;

private:


};
}
#endif
