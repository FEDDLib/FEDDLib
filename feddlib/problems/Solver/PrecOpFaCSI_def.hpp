#ifndef PrecOpFaCSI_DEF_hpp
#define PrecOpFaCSI_DEF_hpp

#include <Thyra_TpetraMultiVector_decl.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_Range1D.hpp>

/*!
 Definition of PrecOpFaCSI

 @brief  PrecOpFaCSI
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
        
// Constructors

template<class SC, class LO, class GO, class NO>
PrecOpFaCSI<SC,LO,GO,NO>::PrecOpFaCSI()
:PreconditionerOperator<SC,LO,GO,NO>(),
fluidPrecMonolithic_(false),
useSolidPreconditioner_(true)
{
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSI, cant be used yet.");
}

template<class SC, class LO, class GO, class NO>
PrecOpFaCSI<SC,LO,GO,NO>::PrecOpFaCSI(CommConstPtr_Type comm, bool fluidPrecMonolithic, bool useFluidPreconditioner, bool useSolidPreconditioner, bool onlyDiagonal)
:PreconditionerOperator<SC,LO,GO,NO>(),
fluidPrecMonolithic_(fluidPrecMonolithic),
useFluidPreconditioner_(useFluidPreconditioner),
useSolidPreconditioner_(useSolidPreconditioner),
onlyDiagonal_(onlyDiagonal)
{
    comm_=comm;
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSI, cant be used yet.");
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGE(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type sInv,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT){

    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setStructInv(sInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);

    initialize();
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGI(ThyraLinOpPtr_Type C1,
                                     ThyraLinOpPtr_Type C1T,
                                     ThyraLinOpPtr_Type C2,
                                     ThyraLinOpPtr_Type C4,
                                     ThyraLinOpPtr_Type sInv,
                                     ThyraLinOpPtr_Type fInv,
                                     ThyraLinOpPtr_Type fF,
                                     ThyraLinOpPtr_Type fBT,
                                     ThyraLinOpPtr_Type gInv){
    
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setC4(C4);
    setStructInv(sInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);
    setGeoInv(gInv);

    initialize();
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGIShape(ThyraLinOpPtr_Type C1,
                                          ThyraLinOpPtr_Type C1T,
                                          ThyraLinOpPtr_Type C2,
                                          ThyraLinOpPtr_Type C4,
                                          ThyraLinOpPtr_Type sInv,
                                          ThyraLinOpPtr_Type fInv,
                                          ThyraLinOpPtr_Type fF,
                                          ThyraLinOpPtr_Type fBT,
                                          ThyraLinOpPtr_Type gInv,
                                          ThyraLinOpPtr_Type shape_v,
                                          ThyraLinOpPtr_Type shape_p){
    
    setC1(C1);
    setC1T(C1T);
    setC2(C2);
    setC4(C4);
    setStructInv(sInv);
    setFluidInv(fInv);
    setFluidF(fF);
    setFluidBT(fBT);
    setGeoInv(gInv);
    setShapeDeriv(shape_v, shape_p);
 
    initialize();
}
    
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC1(ThyraLinOpPtr_Type C1){
    C1_ = C1;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC1T(ThyraLinOpPtr_Type C1T){
    C1T_ = C1T;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC2(ThyraLinOpPtr_Type C2){
    C2_ = C2;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setC4(ThyraLinOpPtr_Type C4){
    C4_ = C4;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setShapeDeriv(ThyraLinOpPtr_Type shape_v, ThyraLinOpPtr_Type shape_p){
    shape_v_ = shape_v;
    shape_p_ = shape_p;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setStructInv(ThyraLinOpPtr_Type sInv){
    sInv_ = sInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setFluidInv(ThyraLinOpPtr_Type fInv){
    fInv_ = fInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setGeoInv(ThyraLinOpPtr_Type gInv){
    gInv_ = gInv;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setFluidF(ThyraLinOpPtr_Type fF){
    fF_ = fF;
}
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::setFluidBT(ThyraLinOpPtr_Type fBT){
    fBT_ = fBT;
}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::initialize(){
    TEUCHOS_TEST_FOR_EXCEPTION(fInv_.is_null(), std::runtime_error,"Can not initialize FaCSI preconditioner: Fluid preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(sInv_.is_null(), std::runtime_error,"Can not initialize FaCSI preconditioner: Structure preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(C1_.is_null(), std::runtime_error,"Can not initialize FaCSI preconditioner: C1 not set.");
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRange( 4 );
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesDomain( 4 );
    vectorSpacesRange[0] = fF_->range();
    vectorSpacesRange[1] = fBT_->domain();
    vectorSpacesRange[2] = sInv_->range();
    vectorSpacesRange[3] = C1_->range();
    
    vectorSpacesDomain[0] = fF_->domain();
    vectorSpacesDomain[1] = fBT_->domain();
    vectorSpacesDomain[2] = sInv_->domain();
    vectorSpacesDomain[3] = C1T_->domain();
    if ( !gInv_.is_null() ) {
        vectorSpacesRange.push_back( gInv_->range() );
        vectorSpacesDomain.push_back( gInv_->domain() );
    }
    //     defaultProductRange_;
    //    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > defaultProductDomain_;
    
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pR = Thyra::productVectorSpace<SC>( vectorSpacesRange );
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pD = Thyra::productVectorSpace<SC>( vectorSpacesDomain );
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSR = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pR);
//    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > pVSD = Teuchos::rcp_dynamic_cast<const Thyra::VectorSpaceBase<SC> >(pD);
//    
//    this->defaultProductRange_ = Thyra::multiVectorProductVectorSpace<SC>( pVSR , 1);
//    this->defaultProductDomain_ = Thyra::multiVectorProductVectorSpace<SC>( pVSD, 1);
    
    this->defaultProductRange_ = pR;
    this->defaultProductDomain_ = pD;
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::applyIt(
                                         const Thyra::EOpTransp M_trans,
                                         const Thyra::MultiVectorBase<SC> &X_in,
                                         const Thyra::Ptr<Thyra::MultiVectorBase<SC> > &Y_inout,
                                         const SC alpha,
                                         const SC beta
                                         ) const
{
    applyImpl(M_trans, X_in, Y_inout, alpha, beta);

}
    
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::applyImpl(
                                   const Thyra::EOpTransp M_trans,
                                   const Thyra::MultiVectorBase<SC> &X_in,
                                   const Thyra::Ptr<Thyra::MultiVectorBase<SC> > &Y_inout,
                                   const SC alpha,
                                   const SC beta
                                   ) const
{
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    typedef Teuchos::ScalarTraits<SC> ST;
    typedef Teuchos::RCP<Thyra::MultiVectorBase<SC> > MultiVectorPtr;
    typedef Teuchos::RCP<const Thyra::MultiVectorBase<SC> > ConstMultiVectorPtr;
    typedef Teuchos::RCP<const Thyra::LinearOpBase<SC> > ConstLinearOpPtr;

    int rank = comm_->getRank();
    if (onlyDiagonal_) {
        
        Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > X
        = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC> > ( Teuchos::rcpFromRef(X_in) );
        
        Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > Y
        = Teuchos::rcp_dynamic_cast< Thyra::ProductMultiVectorBase<SC> > ( rcpFromPtr(Y_inout) );
        
        Y_inout->assign(0.);

        Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_s = X->getMultiVectorBlock(2);
        Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_s = Y->getNonconstMultiVectorBlock(2);

        
        Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_fv = X->getMultiVectorBlock(0);
        Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_fv = Y->getNonconstMultiVectorBlock(0);
        Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_fp = X->getMultiVectorBlock(1);
        Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_fp = Y->getNonconstMultiVectorBlock(1);
        
        
        Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_l = Y->getNonconstMultiVectorBlock(3);
        Teuchos::RCP<const Thyra::MultiVectorBase< SC > > X_l = X->getMultiVectorBlock(3);

        assign(Y_fv.ptr(), *X_fv);
        assign(Y_fp.ptr(), *X_fp);
        assign(Y_s.ptr(), *X_s);
        assign(Y_l.ptr(), *X_l);
        
//        Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > YsTpetra =
//        Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_s );
//
//        
////        std::cout << "Ys pri:" << std::endl;
////        YsTpetra->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
////        comm_->barrier();    comm_->barrier();    comm_->barrier();
////
//        Teuchos::RCP< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > XsTpetra =
//            Teuchos::rcp_dynamic_cast< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( X_s );
//        Teuchos::RCP< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > XfvT =
//            Teuchos::rcp_dynamic_cast< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( X_fv );
//        Teuchos::RCP< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > XfpT =
//            Teuchos::rcp_dynamic_cast< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( X_fp );
//        Teuchos::RCP< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > XlT =
//            Teuchos::rcp_dynamic_cast< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( X_l );
//        std::cout << "Xs:" << std::endl;
//        XsTpetra->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//        std::cout << "Xfv:" << std::endl;
//        XfvT->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//        std::cout << "Xfp:" << std::endl;
//        XfpT->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//        std::cout << "Xl:" << std::endl;
//        XlT->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//
//
//        
//        if (useSolidPreconditioner_)
//            sInv_->apply(Thyra::NOTRANS, *X_s, Y_s.ptr(), 1., 0.);
//        else
//            assign(Y_s.ptr(), *X_s);
//        
//        
////        std::cout << "Ys after:" << std::endl;
////        YsTpetra->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
////        comm_->barrier();    comm_->barrier();    comm_->barrier();
//
//        
//        if (productRangeFluid_.is_null()) {
//            Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRangeFluid( 2 );
//            
//            vectorSpacesRangeFluid[0] = X_fv->range();
//            vectorSpacesRangeFluid[1] = X_fp->range();
//            
//            productRangeFluid_ = Thyra::productVectorSpace<SC>( vectorSpacesRangeFluid );
//        }
//        
//        Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid( 2 );
//        Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid( 2 );
//        
//        X_fluid[0] = Y_fv;
//        X_fluid[1] = Y_fp;
//        
//        Y_fluid[0] = Y_fv;
//        Y_fluid[1] = Y_fp;
//        
//
//        
//        //    std::cout << "Y_fp before mono" << std::endl;
//        //    Y_fp->describe(*out,Teuchos::VERB_EXTREME);
//        //    comm_->barrier();    comm_->barrier();    comm_->barrier();
//        
//        if (useFluidPreconditioner_){
//            
//            if (fluidPrecMonolithic_) {
//                Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > fMonoVS = fInv_->domain();
//                
//                if ( X_fmono_.is_null() ){
//                    X_fmono_ = createMembers( fMonoVS, X_fluid[0]->domain()->dim() );
//                    Y_fmono_ = createMembers( fMonoVS, Y_fluid[0]->domain()->dim() );
//                }
//                //We can/should speedup this process
//                copyToMono(X_fluid);
//                fInv_->apply(Thyra::NOTRANS, *X_fmono_, Y_fmono_.ptr(), 1., 0.);
//                copyFromMono(Y_fluid);
//            }
//            else{
//                Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodX_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, X_fluid );
//                Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodY_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, Y_fluid );
//                
//                Teuchos::RCP< Thyra::MultiVectorBase<SC> > X_f = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase< SC > >(prodX_f);
//                Teuchos::RCP< Thyra::MultiVectorBase<SC> > Y_f = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase< SC > >(prodY_f);
//                
//                fInv_->apply(Thyra::NOTRANS, *X_f, Y_f.ptr(), 1., 0.);
//            }
//        }
//        else{
//            assign(Y_fv.ptr(), *X_fv);
//            assign(Y_fp.ptr(), *X_fp);
//        }
//        
//        Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > YfvTpetra =
//        Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_fv );
//        Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > YfpTpetra =
//        Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_fp );
//
////        std::cout << "Yfv after:" << std::endl;
////        YfvTpetra->getTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
////        comm_->barrier();    comm_->barrier();    comm_->barrier();
////
////        std::cout << "Yfp after:" << std::endl;
////        YfpTpetra->getTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
////        comm_->barrier();    comm_->barrier();    comm_->barrier();
//
//        
//        assign(Y_l.ptr(), *X_l);
//        
////        Y_fv->describe(*out,Teuchos::VERB_EXTREME);
////        Y_fp->describe(*out,Teuchos::VERB_EXTREME);
//
//        
////        Y_l->describe(*out,Teuchos::VERB_EXTREME);
        
    }
    else{
    Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > X
        = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC> > ( Teuchos::rcpFromRef(X_in) );

    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > Y
        = Teuchos::rcp_dynamic_cast< Thyra::ProductMultiVectorBase<SC> > ( rcpFromPtr(Y_inout) );

    
//    Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> >
//    X = Thyra::castOrCreateSingleBlockProductMultiVector<SC>( this->defaultProductDomain_, Teuchos::rcpFromRef(X_in) );
//    Teuchos::RCP<Thyra::ProductMultiVectorBase<SC> >
//    Y = Thyra::nonconstCastOrCreateSingleBlockProductMultiVector<SC>( this->defaultProductRange_, rcpFromPtr(Y_inout) );
    Y_inout->assign(0.);
//    std::cout << "Xin" << std::endl;
//    X_in.describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
//    
//    std::cout << "Xrcp" << std::endl;
//    X->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
//    std::cout << "Y_inoutin" << std::endl;
//    Y_inout->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_s = X->getMultiVectorBlock(2);
    Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_s = Y->getNonconstMultiVectorBlock(2);
        
        
    Teuchos::RCP< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > XsTpetra =
    Teuchos::rcp_dynamic_cast< const Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( X_s );
//    XsTpetra->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);

//    std::cout << "Xs" << std::endl;
//    X_s->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
//    typedef Thyra::TpetraVector< SC, LO, GO, NO > TpetraVector_Type;
//    typedef Teuchos::RCP<TpetraVector_Type> TpetraVectorPtr_Type;
//    typedef Teuchos::RCP<const TpetraVector_Type> TpetraVectorConstPtr_Type;
//    
//    TpetraVectorConstPtr_Type tpetraThyraXs = Teuchos::rcp_dynamic_cast<const TpetraVector_Type>( X_s );
//    
//    std::cout << "Xs tpetra:" << std::endl;
//    tpetraThyraXs->getConstTpetraVector()->describe(*out,Teuchos::VERB_EXTREME);
//
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
    
    // apply solid preconditioner
    if (useSolidPreconditioner_)
        sInv_->apply(Thyra::NOTRANS, *X_s, Y_s.ptr(), 1., 0.);
    else
        assign(Y_s.ptr(), *X_s);
    
//    std::cout << "Ys" << std::endl;
//    Y_s->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_g;
    Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_g;
    // apply geometry preconditioner
    if ( !gInv_.is_null() ) {
        X_g = X->getMultiVectorBlock(4);
        Y_g = Y->getNonconstMultiVectorBlock(4);

        assign(Y_g.ptr(), *X_g);

        C4_->apply(Thyra::NOTRANS, *Y_s, Y_g.ptr(), -1., 1.);
        
        gInv_->apply(Thyra::NOTRANS, *Y_g, Y_g.ptr(), 1., 0.);
    }
    
    Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_l = X->getMultiVectorBlock(3);
    Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_l = Y->getNonconstMultiVectorBlock(3);

    assign(Y_l.ptr(), *X_l);
//    std::cout << "Yl pre c2 apply" << std::endl;
//    Y_l->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
    C2_->apply( Thyra::NOTRANS, *Y_s, Y_l.ptr(), -1., 1. );
//    C2_->apply( Thyra::NOTRANS, *Y_s, Y_l.ptr(), 1., 0. );
//    std::cout << "Yl after c2 apply" << std::endl;
//    Y_l->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
//    scale(-1., Y_l.ptr());
//    update(1., *X_l, Y_l.ptr());
//    std::cout << "Yl2" << std::endl;
//    Y_l->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
    Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_fv = X->getMultiVectorBlock(0);
    Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_fv = Y->getNonconstMultiVectorBlock(0);
    assign(Y_fv.ptr(), *X_fv);
    //    if (Z_fv_.is_null())
//        Z_fv_ = X_fv->clone_mv();
//    else
    
    

    
    Teuchos::RCP< const Thyra::MultiVectorBase< SC > > X_fp = X->getMultiVectorBlock(1);
    Teuchos::RCP< Thyra::MultiVectorBase< SC > > Y_fp = Y->getNonconstMultiVectorBlock(1);
    assign(Y_fp.ptr(), *X_fp);
//    if (Z_fp_.is_null())
//        Z_fp_ = X_fv->clone_mv();
//    else
//        assign(Z_fp_.ptr(), *X_fp);
    
//    std::cout << "zp" << std::endl;
//    Z_fp->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//    std::cout << "Y-fp set:" << std::endl;
//    Y_fp->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    if (!shape_v_.is_null() && !shape_p_.is_null()) {
        shape_v_->apply(Thyra::NOTRANS, *Y_g, Y_fv.ptr(), -1., 1.);
        shape_p_->apply(Thyra::NOTRANS, *Y_g, Y_fp.ptr(), -1., 1.);
    }
    
    
    if (Z_fv_.is_null())
        Z_fv_ = Y_fv->clone_mv();
    else
        assign(Z_fv_.ptr(), *Y_fv);
    
//    std::cout << "Z_fv_ set:" << std::endl;
//    Z_fv_->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
//
//    std::cout << "Z_fp_ set:" << std::endl;
//    Y_fp->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();

    
    // fluid condensation
    if (tmp_l_.is_null())
        tmp_l_ = Y_l->clone_mv();

    C1_->apply(Thyra::NOTRANS, *Y_fv, tmp_l_.ptr(), 1., 0);
    C1T_->apply(Thyra::NOTRANS, *tmp_l_, Y_fv.ptr(), -1., 1.);
    
    C1T_->apply(Thyra::NOTRANS, *Y_l, Y_fv.ptr(), 1., 1.);
    
//    std::cout << "W_fv" << std::endl;
//    Y_fv->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
    if (productRangeFluid_.is_null()) {
        Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRangeFluid( 2 );
        
        vectorSpacesRangeFluid[0] = X_fv->range();
        vectorSpacesRangeFluid[1] = X_fp->range();
        
        productRangeFluid_ = Thyra::productVectorSpace<SC>( vectorSpacesRangeFluid );
    }
    
    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid( 2 );
    Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid( 2 );

    X_fluid[0] = Y_fv;
    X_fluid[1] = Y_fp;

    Y_fluid[0] = Y_fv;
    Y_fluid[1] = Y_fp;

    
//    std::cout << "Y_fp before mono" << std::endl;
//    Y_fp->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
    if (useFluidPreconditioner_){
    
        if (fluidPrecMonolithic_) {
            Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > fMonoVS = fInv_->domain();

            if ( X_fmono_.is_null() ){
                X_fmono_ = createMembers( fMonoVS, X_fluid[0]->domain()->dim() );
                Y_fmono_ = createMembers( fMonoVS, Y_fluid[0]->domain()->dim() );
            }
            //We can/should speedup this process
            copyToMono(X_fluid);
            fInv_->apply(Thyra::NOTRANS, *X_fmono_, Y_fmono_.ptr(), 1., 0.);
            copyFromMono(Y_fluid);
        }
        else{
            Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodX_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, X_fluid );
            Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > prodY_f = Thyra::defaultProductMultiVector<SC>( productRangeFluid_, Y_fluid );
            
            Teuchos::RCP< Thyra::MultiVectorBase<SC> > X_f = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase< SC > >(prodX_f);
            Teuchos::RCP< Thyra::MultiVectorBase<SC> > Y_f = Teuchos::rcp_dynamic_cast<Thyra::MultiVectorBase< SC > >(prodY_f);
            
            fInv_->apply(Thyra::NOTRANS, *X_f, Y_f.ptr(), 1., 0.);
        }
    }
    else{
        assign(Y_fv.ptr(), *X_fv);
        assign(Y_fp.ptr(), *X_fp);
    }

//    std::cout << "Y_fv after mono" << std::endl;
//    Y_fv->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
//    std::cout << "Y_fp after mono" << std::endl;
//    Y_fp->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();


    fBT_->apply(Thyra::NOTRANS, *Y_fp, Z_fv_.ptr(), -1., 1.);

//    std::cout << "W_fv after BT" << std::endl;
//    W_fv->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
    fF_->apply(Thyra::NOTRANS, *Y_fv, Z_fv_.ptr(), -1., 1.);

//    std::cout << "W_fv after F" << std::endl;
//    W_fv->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
//    std::cout << "W_fv after update" << std::endl;
//    Z_fv_->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
    
    C1_->apply(Thyra::NOTRANS, *Z_fv_, Y_l.ptr(), 1., 0.);
    
//    std::cout << "Y_l after C1" << std::endl;
//    Y_l->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
//    std::cout << "Y_inoutout" << std::endl;
//    Y_inout->describe(*out,Teuchos::VERB_EXTREME);
//    comm_->barrier();    comm_->barrier();    comm_->barrier();
//    std::cout << "FaCSI done!" << std::endl;
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,"End of FaCSI apply.");
        
        
        
//        Teuchos::RCP< MultiVectorBase< SC > > Y_fv = Y->getNonconstMultiVectorBlock(0);
//        Teuchos::RCP< MultiVectorBase< SC > > Y_fp = Y->getNonconstMultiVectorBlock(1);
//        Teuchos::RCP< MultiVectorBase< SC > > Y_s = Y->getNonconstMultiVectorBlock(2);
//        Teuchos::RCP< MultiVectorBase< SC > > Y_l = Y->getNonconstMultiVectorBlock(3);

        Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_fvT =
        Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_fv );
        Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_fpT =
        Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_fp );
        Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_sT =
            Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_s );
        Teuchos::RCP< Thyra::TpetraMultiVector< SC, LO, GO, NO > > Y_lT =
        Teuchos::rcp_dynamic_cast< Thyra::TpetraMultiVector< SC, LO, GO, NO > > ( Y_l );

//        std::cout << "Yfv:" << std::endl;
//        Y_fvT->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//        std::cout << "Yfp:" << std::endl;
//        Y_fpT->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//        std::cout << "Ys:" << std::endl;
//        Y_sT->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();
//        std::cout << "Yl:" << std::endl;
//        Y_lT->getConstTpetraMultiVector()->describe(*out,Teuchos::VERB_EXTREME);
//        comm_->barrier();    comm_->barrier();    comm_->barrier();

        
        
        
    }
}
// private
template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::copyToMono( Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > X_fluid ) const{
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fluid[1]->range());
    
    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[0] ,Teuchos::Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fluid[1] ,Teuchos::Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(X_fmono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( X_fmono_ ,Teuchos::Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );

    for (int j=0; j<X_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_fluid)(i,j) = (*thyData_v)(i,j);
        
        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_fluid)( i + localSubDim_v, j ) = (*thyData_p)(i,j);

    }
}

template<class SC, class LO, class GO, class NO>
void PrecOpFaCSI<SC,LO,GO,NO>::copyFromMono(Teuchos::Array< Teuchos::RCP< Thyra::MultiVectorBase< SC > > > Y_fluid) const{
    

    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_v = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[0]->range());
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_p = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fluid[1]->range());
    
    const LO localOffset_v = mpiVS_v->localOffset();
    const LO localSubDim_v = mpiVS_v->localSubDim();
    
    const LO localOffset_p = mpiVS_p->localOffset();
    const LO localSubDim_p = mpiVS_p->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_v =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[0] ,Teuchos::Range1D(localOffset_v,localOffset_v+localSubDim_v-1) ) );
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_p =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fluid[1] ,Teuchos::Range1D(localOffset_p,localOffset_p+localSubDim_p-1) ) );
    
    Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<SC> > mpiVS_mono = Teuchos::rcp_dynamic_cast<const Thyra::SpmdVectorSpaceBase<SC> >(Y_fmono_->range());
    
    const LO localOffset_mono = mpiVS_mono->localOffset();
    const LO localSubDim_mono = mpiVS_mono->localSubDim();
    
    Teuchos::RCP<Thyra::DetachedMultiVectorView<SC> > thyData_fluid =
    Teuchos::rcp(new Thyra::DetachedMultiVectorView<SC>( Y_fmono_ ,Teuchos::Range1D(localOffset_mono,localOffset_mono+localSubDim_mono-1) ) );
    
    for (int j=0; j<Y_fluid[0]->domain()->dim(); j++) {
        
        for (LO i=0; i < localSubDim_v; i++)
            (*thyData_v)(i,j) = (*thyData_fluid)(i,j);
        
        for (LO i=0; i<localSubDim_p; i++)
            (*thyData_p)(i,j) = (*thyData_fluid)( i + localSubDim_v, j );
        
    }    
}

}

#endif
