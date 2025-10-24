#ifndef PrecBlock2x2_DECL_hpp
#define PrecBlock2x2_DECL_hpp
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/problems/problems_config.h"
#include "PreconditionerOperator.hpp"
#include <Thyra_DefaultProductMultiVector_decl.hpp>
#include <Thyra_DefaultMultiVectorProductVectorSpace_decl.hpp>
#include <Thyra_OperatorVectorTypes.hpp>
#include <Thyra_MultiVectorStdOps_decl.hpp>
#include <Thyra_MultiVectorBase_decl.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_TpetraVector_decl.hpp>
#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_SpmdVectorSpaceBase_decl.hpp>
/*!
 Declaration of PrecBlock2x2
 
 @brief  PrecBlock2x2
 @author Christian Hochmuth/Lea Saßmannshausen
 @version 1.0
 @copyright CH
 */


/*! 
	Applying a 2x2 block preconditioner for the (Navier)-Stokes problem
    | F  B^T |
    | B -C   |
    The matrix C corresponds to a stabilzation, if needed. For stable finite element discretizations,
    the matrix C corresponds to a zero matrix. In PCD there is no special consideration of the 
    stabiization needed. In LSC it would be considered, but that is not implemented here.
    
    Diagonal Preconditioner:
    | \hat{F}^-1  0        |
    | 0         \hat{S}^-1 |

    Triangular Preconditioner:
    | \hat{F}^-1  B^T      |
    | 0         \hat{S}^-1 |

    The inverse of the fluid system is always approximated the same by a Schwarz method

    The Schur complement is replaced by different approximations depending on the strategy:
    -> 'Diagonal' Prec:     the Schur complement is replaced by  -1/nu M_p
    -> 'Triangular' Prec:   the Schur complement is replaced by  -1/nu M_p
    -> 'PCD' Prec:          the Schur complement is replaced by  -M_p F_p^-1 A_p
    -> 'LSC' Prec:          the Schur complement is replaced by  -A_p^-1 (B (M_v^-1) F (M_v^-1) B^T ) A_p^-1

    the arising inverses are again approximated by a Schwarz method

*/

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class PrecBlock2x2 : public PreconditionerOperator<SC,LO,GO,NO> {

public:
    typedef Teuchos::RCP<Thyra::LinearOpBase<SC> > ThyraLinOpPtr_Type;
    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    PrecBlock2x2();
    
    PrecBlock2x2(CommConstPtr_Type comm);
    /*! In the following, we assume that we build a precondtioner for a Stokes-type fluid problem with velocity and pressure variables */

    /// Diagonal preconditioner with \hat{S}= -1/nu M_p schur complement approximation
    void setDiagonal(ThyraLinOpPtr_Type velocityInv,
                     ThyraLinOpPtr_Type pressureInv);
    
    /// Triangular preconditioner with \hat{S}= -1/nu M_p schur complement approximation
    void setTriangular(ThyraLinOpPtr_Type velocityInv,
                       ThyraLinOpPtr_Type pressureInv,
                       ThyraLinOpPtr_Type BT);

    /// Pressure-Convection-Diffusion (PCD) block triangular preconditioner
    void setTriangular(ThyraLinOpPtr_Type velocityInv,
                        ThyraLinOpPtr_Type laplaceInverse,
                        ThyraLinOpPtr_Type convectionDiffusionOperator,
                        ThyraLinOpPtr_Type massMatrixInverse,
                        ThyraLinOpPtr_Type massMatrixVInverse,
                       ThyraLinOpPtr_Type BT);

    /// Least-Squares-Commutator (LSC) block triangular preconditioner
    void setTriangular(ThyraLinOpPtr_Type velocityInv,
                    ThyraLinOpPtr_Type laplaceInverse,
                    ThyraLinOpPtr_Type massMatrixVInverse,
                    ThyraLinOpPtr_Type BT);
    
    /// Setting velocity inverse approximation
    void setVeloctiyInv(ThyraLinOpPtr_Type veloctiyInv);

    /// Setting inverse approximation of Schur complement
    void setPressureInv(ThyraLinOpPtr_Type pressureInv);

    /// Setting inverse approximation of Schur complement by multiple operators. Corresponds to PCD
    void setPressureInvs(ThyraLinOpPtr_Type laplaceInverse, ThyraLinOpPtr_Type convectionDiffusionOperator, ThyraLinOpPtr_Type massMatrixInverse, ThyraLinOpPtr_Type massMatrixVInverse);

    /// Setting inverse approximation of Schur complement by multiple operators. Corresponds to LSC
    void setPressureInvs(ThyraLinOpPtr_Type laplaceInverse, ThyraLinOpPtr_Type massMatrixVInverse);

    /// Setting fluid system matrix B
    void setB(ThyraLinOpPtr_Type B) {B_ = B;}

    /// Setting fluid system matrix F
    void setF(ThyraLinOpPtr_Type F) {F_ = F;}

    /// Setting the preconditioning ype
    void setType(std::string type);
    
    void initialize();
    
    virtual void applyIt(
                           const Thyra::EOpTransp M_trans,
                           const Thyra::MultiVectorBase<SC> &X,
                           const Teuchos::Ptr<Thyra::MultiVectorBase<SC> > &Y,
                           const SC alpha,
                           const SC beta
                           ) const;

protected:
    /// Apply of preconditioner in e.g. GMRES
    virtual void applyImpl(
                           const Thyra::EOpTransp M_trans,
                           const Thyra::MultiVectorBase<SC> &X,
                           const Teuchos::Ptr<Thyra::MultiVectorBase<SC> > &Y,
                           const SC alpha,
                           const SC beta
                           ) const;

    


    
private:    
    
    ThyraLinOpPtr_Type velocityInv_;
    ThyraLinOpPtr_Type pressureInv_;
    ThyraLinOpPtr_Type BT_;
    ThyraLinOpPtr_Type B_;
    ThyraLinOpPtr_Type F_;
    ThyraLinOpPtr_Type laplaceInverse_;
    ThyraLinOpPtr_Type convectionDiffusionOperator_;
    ThyraLinOpPtr_Type massMatrixInverse_;
    ThyraLinOpPtr_Type massMatrixVInverse_;
    
    mutable Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > productRangeFluid_;
    
    CommConstPtr_Type comm_;
    std::string type_;
    
};
}
#endif
