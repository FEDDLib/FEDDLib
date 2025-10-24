#ifndef PrecBlock2x2_DEF_hpp
#define PrecBlock2x2_DEF_hpp
#include "PrecBlock2x2_decl.hpp"
#include <Thyra_TpetraMultiVector_decl.hpp>
#include <Teuchos_VerboseObject.hpp>
/*!
 Definition of PrecBlock2x2

 @brief  PrecBlock2x2
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */



namespace FEDD {
using namespace Thyra;
        
// Constructors

template<class SC, class LO, class GO, class NO>
PrecBlock2x2<SC,LO,GO,NO>::PrecBlock2x2()
:PreconditionerOperator<SC,LO,GO,NO>()
{
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSI, cant be used yet.");
}

template<class SC, class LO, class GO, class NO>
PrecBlock2x2<SC,LO,GO,NO>::PrecBlock2x2( CommConstPtr_Type comm )
:PreconditionerOperator<SC,LO,GO,NO>()
{
    comm_=comm;
//    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Still a problem for FaCSI, cant be used yet.");
}
    
template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setDiagonal(
                                            ThyraLinOpPtr_Type velocityInv,
                                            ThyraLinOpPtr_Type pressureInv
                                            ){
    setVeloctiyInv(velocityInv);
    
    setPressureInv(pressureInv);

    setType("Diagonal");
    
    initialize();
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setTriangular(
                                              ThyraLinOpPtr_Type velocityInv,
                                              ThyraLinOpPtr_Type pressureInv,
                                              ThyraLinOpPtr_Type BT
                                              ){
    setVeloctiyInv(velocityInv);
    
    setPressureInv(pressureInv);

    setType("Triangular");

    BT_ = BT;
    
    initialize();
}


template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setTriangular(ThyraLinOpPtr_Type velocityInv,
                        ThyraLinOpPtr_Type laplaceInverse,
                        ThyraLinOpPtr_Type convectionDiffusionOperator,
                        ThyraLinOpPtr_Type massMatrixInverse,
                        ThyraLinOpPtr_Type massMatrixVInverse,
                       ThyraLinOpPtr_Type BT){

    setVeloctiyInv(velocityInv);
    
    setPressureInvs(laplaceInverse,convectionDiffusionOperator,massMatrixInverse,massMatrixVInverse);

    setType("PCD");

    BT_ = BT;
    
    initialize();
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setTriangular(ThyraLinOpPtr_Type velocityInv,
                        ThyraLinOpPtr_Type laplaceInverse,
                        ThyraLinOpPtr_Type massMatrixVInverse,
                       ThyraLinOpPtr_Type BT){

    setVeloctiyInv(velocityInv);
    
    setPressureInvs(laplaceInverse,massMatrixVInverse);

    setType("LSC");

    BT_ = BT;
    
    initialize();
}


template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setVeloctiyInv(ThyraLinOpPtr_Type velocityInv){
    velocityInv_ = velocityInv;
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setPressureInv(ThyraLinOpPtr_Type pressureInv){
    pressureInv_ = pressureInv;
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setPressureInvs(ThyraLinOpPtr_Type laplaceInverse,
                        ThyraLinOpPtr_Type convectionDiffusionOperator,
                        ThyraLinOpPtr_Type massMatrixInverse,
                        ThyraLinOpPtr_Type massMatrixVInverse){

    laplaceInverse_ = laplaceInverse;
    convectionDiffusionOperator_=convectionDiffusionOperator;
    massMatrixInverse_=massMatrixInverse;
    massMatrixVInverse_=massMatrixVInverse;
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setPressureInvs(ThyraLinOpPtr_Type laplaceInverse,
                        ThyraLinOpPtr_Type massMatrixVInverse){

    laplaceInverse_ = laplaceInverse;
    massMatrixVInverse_=massMatrixVInverse;
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::setType(std::string type){
    type_ = type;
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::initialize(){
   TEUCHOS_TEST_FOR_EXCEPTION(velocityInv_.is_null(), std::runtime_error,"Can not initialize Block2x2 preconditioner: 1 preconditioner not set.");
    TEUCHOS_TEST_FOR_EXCEPTION(pressureInv_.is_null() && laplaceInverse_.is_null(), std::runtime_error,"Can not initialize Block2x2 preconditioner: 2 preconditioner not set.");
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRange( 2 );
    Teuchos::Array< Teuchos::RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesDomain( 2 );
    vectorSpacesRange[0] = velocityInv_->range();
    
    vectorSpacesDomain[0] = velocityInv_->domain();

    if(!pressureInv_.is_null()){
        vectorSpacesRange[1] = pressureInv_->domain();
        vectorSpacesDomain[1] = pressureInv_->domain();
    }
    else{
        vectorSpacesRange[1] = laplaceInverse_->domain();
        vectorSpacesDomain[1] = laplaceInverse_->domain();
    }
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pR = Thyra::productVectorSpace<SC>( vectorSpacesRange );
    Teuchos::RCP<const Thyra::DefaultProductVectorSpace<SC> > pD = Thyra::productVectorSpace<SC>( vectorSpacesDomain );

    this->defaultProductRange_ = pR;
    this->defaultProductDomain_ = pD;
}

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::applyIt(
                                         const EOpTransp M_trans,
                                         const MultiVectorBase<SC> &X_in,
                                         const Ptr<MultiVectorBase<SC> > &Y_inout,
                                         const SC alpha,
                                         const SC beta
                                         ) const
{
    applyImpl(M_trans, X_in, Y_inout, alpha, beta);

}
    
/// Function called within e.g. GMRES to apply the input vector X_in (corresponding to the RHS = (f_u,f_p)) to the preconditioner
/// Output is the contents of vector Y_inout
// Block Preconditioners:
// - 'Diagonal' Prec: where the Schur complement is replaced by - 1/nu M_p
// - 'Triangular' Prec: where the Schur complement is replaced by - 1/nu M_p
// - 'PCD' Prec: where the Schur complement is replaced by -M_p F_p^-1 A_p
// - 'LSC' Prec: where the Schur complement is replaced by  -A_p^-1 (B (M_v^-1) F (M_v^-1) B^T ) A_p^-1

template<class SC, class LO, class GO, class NO>
void PrecBlock2x2<SC,LO,GO,NO>::applyImpl(
                                   const EOpTransp M_trans,
                                   const MultiVectorBase<SC> &X_in,
                                   const Ptr<MultiVectorBase<SC> > &Y_inout,
                                   const SC alpha,
                                   const SC beta
                                   ) const
{
    // alpha and beta are ignored!
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    using Teuchos::rcpFromRef;
    typedef Teuchos::ScalarTraits<SC> ST;
    typedef RCP<MultiVectorBase<SC> > MultiVectorPtr;
    typedef RCP<const MultiVectorBase<SC> > ConstMultiVectorPtr;
    typedef RCP<const LinearOpBase<SC> > ConstLinearOpPtr;

    int rank = comm_->getRank();
        
    Teuchos::RCP<const Thyra::ProductMultiVectorBase<SC> > X
        = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC> > ( rcpFromRef(X_in) );
    
    Teuchos::RCP< Thyra::ProductMultiVectorBase<SC> > Y
        = Teuchos::rcp_dynamic_cast< Thyra::ProductMultiVectorBase<SC> > ( rcpFromPtr(Y_inout) );

    // First component correspondig to velocity block
    Teuchos::RCP< const MultiVectorBase< SC > > X_0 = X->getMultiVectorBlock(0); 
    Teuchos::RCP< MultiVectorBase< SC > > Y_0 = Y->getNonconstMultiVectorBlock(0);

    // Second component corresponding to pressure block
    Teuchos::RCP< const MultiVectorBase< SC > > X_1 = X->getMultiVectorBlock(1);
    Teuchos::RCP< MultiVectorBase< SC > > Y_1 = Y->getNonconstMultiVectorBlock(1);
    
    // We distinguish between the different preconditioning types here
    if (type_ == "Diagonal"){
        velocityInv_->apply(NOTRANS, *X_0, Y_0.ptr(), 1., 0.); // Apply input vector to fluid inverse approximation

        pressureInv_->apply(NOTRANS, *X_1, Y_1.ptr(), 1., 0.); // Apply input vector to Schur complement inverse approximation
                                                               // Here: 1/nu Mp
    }
    else if (type_ == "Triangular"){

        pressureInv_->apply(NOTRANS, *X_1, Y_1.ptr(), 1., 0.); // Apply input vector pressure component to Schur complement inverse approximation
        
        Teuchos::RCP< MultiVectorBase< SC > > Z_0 = X_0->clone_mv();
        
        BT_->apply(NOTRANS, *Y_1, Z_0.ptr(), -1., 1.); //Z0= BT*Y1 + X0
        
        velocityInv_->apply(NOTRANS, *Z_0, Y_0.ptr(), 1., 0.);
                        
    }
    else if (type_ == "PCD"){
        // PCD block triangular preconditioner
        TEUCHOS_TEST_FOR_EXCEPTION(laplaceInverse_.is_null(), std::runtime_error,"laplaceInverse_ not set.");
        TEUCHOS_TEST_FOR_EXCEPTION(convectionDiffusionOperator_.is_null(), std::runtime_error,"convectionDiffusionOperator_ not set.");
        TEUCHOS_TEST_FOR_EXCEPTION(massMatrixInverse_.is_null(), std::runtime_error,"massMatrixInverse_ not set.");
        // For PCD we need apply the 'pressure inverse' differently, as it is made up of three components.
        // X_1->describe(*out,Teuchos::VERB_EXTREME); // Examplary print 

        // Apply input vector to Schur complement inverse approximation
        massMatrixInverse_->apply(NOTRANS, *X_1, Y_1.ptr(), 1., 0.);  // 1. pressure mass matrix inverse approximation

        convectionDiffusionOperator_->apply(NOTRANS, *Y_1, Y_1.ptr(), 1., 0.); // 2. pressure convection-diffusion operator

        laplaceInverse_->apply(NOTRANS, *Y_1, Y_1.ptr(), 1., 0.);  // 3. pressure laplace inverse approximation

        // Alternative option to use the Laplace constructed via B M_v B^T
        // else{ // We operate in different dimensions here
        //     Teuchos::RCP< MultiVectorBase< SC > > X_res_0 = X_0->clone_mv();
        //     BT_->apply(NOTRANS, *Y_1, X_res_0.ptr(), 1., 0.); //BT*y
        //     massMatrixVInverse_->apply(NOTRANS, *X_res_0, X_res_0.ptr(), 1., 0.);
        //     BT_->apply(TRANS, *X_res_0, Y_1.ptr(), 1., 0.);  
        // }

        Teuchos::RCP< MultiVectorBase< SC > > Z_0 = X_0->clone_mv();
        
        BT_->apply(NOTRANS, *Y_1, Z_0.ptr(), -1., 1.); //Z0= - BT*Y1 + X0
        
        velocityInv_->apply(NOTRANS, *Z_0, Y_0.ptr(), 1., 0.); // Apply input vector to fluid inverse approximation
            
    }
    else if (type_ == "LSC"){
        // LSC block triangular preconditioner
        TEUCHOS_TEST_FOR_EXCEPTION(laplaceInverse_.is_null(), std::runtime_error,"laplaceInverse_ not set.");
        TEUCHOS_TEST_FOR_EXCEPTION(massMatrixVInverse_.is_null(), std::runtime_error,"massMatrixVInverse_ not set.");
        
        laplaceInverse_->apply(NOTRANS, *X_1, Y_1.ptr(), 1., 0.);   // 1. pressure laplace inverse approximation

        // 2. Component: B (M_v^-1) F (M_v^-1) B^T 
        Teuchos::RCP< MultiVectorBase< SC > > X_res_0 = X_0->clone_mv();
        BT_->apply(NOTRANS, *Y_1, X_res_0.ptr(), 1., 0.);
        massMatrixVInverse_->apply(NOTRANS, *X_res_0, X_res_0.ptr(), 1., 0.);
        F_->apply(NOTRANS, *X_res_0, X_res_0.ptr(), 1., 0.);
        massMatrixVInverse_->apply(NOTRANS, *X_res_0, X_res_0.ptr(), 1., 0.);
        BT_->apply(TRANS, *X_res_0, Y_1.ptr(), 1., 0.);

        laplaceInverse_->apply(NOTRANS, *Y_1, Y_1.ptr(), -1., 0.);  // 3. pressure laplace inverse approximation
 
        Teuchos::RCP< MultiVectorBase< SC > > Z_0 = X_0->clone_mv();
        BT_->apply(NOTRANS, *Y_1, Z_0.ptr(), -1., 1.); //Z0= BT*Y1 + X0
        
        velocityInv_->apply(NOTRANS, *Z_0, Y_0.ptr(), 1., 0.);
            
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,"Unknow 2x2 block preconditioner type. Select Diagonal or Triangular.");
    }
}
}

#endif
