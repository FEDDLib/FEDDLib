#ifndef PrecOpFaCSI_DECL_hpp
#define PrecOpFaCSI_DECL_hpp
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
 Declaration of PrecOpFaCSI
 
 @brief  PrecOpFaCSI
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class PrecOpFaCSI : public PreconditionerOperator<SC,LO,GO,NO> {

public:
    typedef Teuchos::RCP<Thyra::LinearOpBase<SC> > ThyraLinOpPtr_Type;
    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    PrecOpFaCSI();
    
    PrecOpFaCSI(CommConstPtr_Type comm, bool fluidPrecMonolithic, bool useFluidPreconditioner = true, bool useSolidPreconditioner = true, bool onlyDiagonal=false);
    
    void setGE(ThyraLinOpPtr_Type C1,
          ThyraLinOpPtr_Type C1T,
          ThyraLinOpPtr_Type C2,
          ThyraLinOpPtr_Type sInv,
          ThyraLinOpPtr_Type fInv,
          ThyraLinOpPtr_Type fF,
          ThyraLinOpPtr_Type fBT);

    void setGI(ThyraLinOpPtr_Type C1,
          ThyraLinOpPtr_Type C1T,
          ThyraLinOpPtr_Type C2,
          ThyraLinOpPtr_Type C4,
          ThyraLinOpPtr_Type sInv,
          ThyraLinOpPtr_Type fInv,
          ThyraLinOpPtr_Type fF,
          ThyraLinOpPtr_Type fBT,
          ThyraLinOpPtr_Type gInv);

    void setGIShape(ThyraLinOpPtr_Type C1,
               ThyraLinOpPtr_Type C1T,
               ThyraLinOpPtr_Type C2,
               ThyraLinOpPtr_Type C4,
               ThyraLinOpPtr_Type sInv,
               ThyraLinOpPtr_Type fInv,
               ThyraLinOpPtr_Type fF,
               ThyraLinOpPtr_Type fBT,
               ThyraLinOpPtr_Type gInv,
               ThyraLinOpPtr_Type shape_v,
               ThyraLinOpPtr_Type shape_p);

    void setC1(ThyraLinOpPtr_Type C1);
    void setC1T(ThyraLinOpPtr_Type C1T);
    void setC2(ThyraLinOpPtr_Type C2);
    void setC4(ThyraLinOpPtr_Type C4);
    void setShapeDeriv(ThyraLinOpPtr_Type shape_v, ThyraLinOpPtr_Type shape_p);
    void setStructInv(ThyraLinOpPtr_Type sInv);
    void setFluidInv(ThyraLinOpPtr_Type fInv);
    void setGeoInv(ThyraLinOpPtr_Type gInv);
    void setFluidF(ThyraLinOpPtr_Type fF);
    void setFluidBT(ThyraLinOpPtr_Type fBT);

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

    void copyToMono(Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid) const;
    
    void copyFromMono(Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid) const;
    
    ThyraLinOpPtr_Type C1_;
    ThyraLinOpPtr_Type C1T_;
    ThyraLinOpPtr_Type C2_;
    ThyraLinOpPtr_Type C4_;
    ThyraLinOpPtr_Type sInv_;
    ThyraLinOpPtr_Type fInv_;
    ThyraLinOpPtr_Type fF_;
    ThyraLinOpPtr_Type fBT_;
    ThyraLinOpPtr_Type gInv_;
    ThyraLinOpPtr_Type shape_v_;
    ThyraLinOpPtr_Type shape_p_;
    mutable Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > productRangeFluid_;
    bool fluidPrecMonolithic_;
    mutable Teuchos::RCP< Thyra::MultiVectorBase<SC> > X_fmono_;
    mutable Teuchos::RCP< Thyra::MultiVectorBase<SC> > Y_fmono_;
    CommConstPtr_Type comm_;
    bool useFluidPreconditioner_;
    bool useSolidPreconditioner_;
    bool onlyDiagonal_;
    mutable Teuchos::RCP< Thyra::MultiVectorBase<SC> > tmp_l_;
    mutable Teuchos::RCP< Thyra::MultiVectorBase<SC> > Z_fv_;
};
}
#endif
