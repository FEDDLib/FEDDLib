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
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
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
    void setDiagonal(ThyraLinOpPtr_Type velocityInv,
                     ThyraLinOpPtr_Type pressureInv);
    
    void setTriangular(ThyraLinOpPtr_Type velocityInv,
                       ThyraLinOpPtr_Type pressureInv,
                       ThyraLinOpPtr_Type BT);
    
    void setVeloctiyInv(ThyraLinOpPtr_Type veloctiyInv);
    
    void setPressureInv(ThyraLinOpPtr_Type pressureInv);

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
    
    mutable Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > productRangeFluid_;
    
    CommConstPtr_Type comm_;
    std::string type_;
    
};
}
#endif
