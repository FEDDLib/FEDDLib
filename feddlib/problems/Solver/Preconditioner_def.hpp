#ifndef PRECONDITIONER_START
#define PRECONDITIONER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Preconditioner: ") + std::string(S))));
#endif

#ifndef PRECONDITIONER_STOP
#define PRECONDITIONER_STOP(A) A.reset();
#endif

#ifndef Preconditioner_DEF_hpp
#define Preconditioner_DEF_hpp
#include "Preconditioner_decl.hpp"
#include <Thyra_DefaultZeroLinearOp_decl.hpp>
#ifdef FEDD_HAVE_IFPACK2
#include <Thyra_Ifpack2PreconditionerFactory_def.hpp>
#endif

/*!
 Definition of Preconditioner

 @brief  Preconditioner
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC,class LO,class GO,class NO>
Preconditioner<SC,LO,GO,NO>::Preconditioner(Problem_Type* problem):
thyraPrec_(),
precondtionerIsBuilt_(false),
problem_(),
timeProblem_(),
precFactory_()
#ifdef FEDD_HAVE_TEKO
,tekoLinOp_()
,velocityMassMatrix_()
,rh_()
#endif
,fsiLinOp_()
,precFluid_()
,precStruct_()
,precGeo_()
,probFluid_()
,probSolid_()
,probGeo_()
#ifdef PRECONDITIONER_TIMER
,preconditionerTimer_ (Teuchos::TimeMonitor::getNewCounter("Preconditioner: Setup"))
#endif
{
    problem_.reset( problem, false );

    Stratimikos::enableFROSch<LO,GO,NO>( *problem_->getLinearSolverBuilder() );

#ifdef FEDD_HAVE_IFPACK2
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> CRSMat;
    typedef Tpetra::Operator<SC,LO,GO,NO> TP_op;
    typedef Thyra::PreconditionerFactoryBase<SC>        Base;
    typedef Thyra::Ifpack2PreconditionerFactory<CRSMat> Impl;

    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix_Type;
    typedef Teuchos::RCP<TpetraCrsMatrix_Type> TpetraCrsMatrixPtr_Type;

    problem_->getLinearSolverBuilder()->setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
#ifdef FEDD_HAVE_TEKO
    Teko::addTekoToStratimikosBuilder( *problem_->getLinearSolverBuilder() );
#endif



}

template <class SC,class LO,class GO,class NO>
Preconditioner<SC,LO,GO,NO>::Preconditioner(TimeProblem_Type* problem):
thyraPrec_(),
precondtionerIsBuilt_(false),
problem_(),
timeProblem_(),
precFactory_()
#ifdef FEDD_HAVE_TEKO
,tekoLinOp_()
,velocityMassMatrix_()
,rh_()
#endif
,fsiLinOp_()
,precFluid_()
,precStruct_()
,precGeo_()
,probFluid_()
,probSolid_()
,probGeo_()
#ifdef PRECONDITIONER_TIMER
,preconditionerTimer_ (Teuchos::TimeMonitor::getNewCounter("Preconditioner: Setup"))
#endif
{
    timeProblem_.reset( problem, false );

    // We need to ensure that already set information is upheld
    if(!problem->getUnderlyingProblem()->preconditioner_->getPressureProjection().is_null()){
        setPressureProjection(problem->getUnderlyingProblem()->preconditioner_->getPressureProjection());
    }

    if(!problem->getUnderlyingProblem()->preconditioner_->getVelocityMassMatrix().is_null()){
        setVelocityMassMatrix(problem->getUnderlyingProblem()->preconditioner_->getVelocityMassMatrix());
    }

    if(!problem->getUnderlyingProblem()->preconditioner_->getPressureLaplaceMatrix().is_null()){
        setPressureLaplaceMatrix(problem->getUnderlyingProblem()->preconditioner_->getPressureLaplaceMatrix());
    }
    if(!problem->getUnderlyingProblem()->preconditioner_->getPressureMassMatrix().is_null()){
        setPressureMassMatrix(problem->getUnderlyingProblem()->preconditioner_->getPressureMassMatrix());
    
    }
    if(!problem->getUnderlyingProblem()->preconditioner_->getPCDOperatorMatrix().is_null()){
        setPCDOperator(problem->getUnderlyingProblem()->preconditioner_->getPCDOperatorMatrix());
    }

}

template <class SC,class LO,class GO,class NO>
Preconditioner<SC,LO,GO,NO>::~Preconditioner()
{
}

template<class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setPreconditionerThyraFromLinOp( ThyraLinOpPtr_Type precLinOp ){
    TEUCHOS_TEST_FOR_EXCEPTION( thyraPrec_.is_null(), std::runtime_error, "thyraPrec_ is null." );
    Teuchos::RCP<Thyra::DefaultPreconditioner<SC> > defaultPrec = Teuchos::rcp_dynamic_cast<Thyra::DefaultPreconditioner<SC> > (thyraPrec_);
    TEUCHOS_TEST_FOR_EXCEPTION( defaultPrec.is_null(), std::runtime_error, "Cast to DefaultPrecondtioner failed." );
    TEUCHOS_TEST_FOR_EXCEPTION( precLinOp.is_null(), std::runtime_error, "precLinOp is null." );
    defaultPrec->initializeUnspecified( precLinOp );
}
    
template <class SC,class LO,class GO,class NO>
typename Preconditioner<SC,LO,GO,NO>::ThyraPrecPtr_Type Preconditioner<SC,LO,GO,NO>::getThyraPrec(){
    return thyraPrec_;
}

template <class SC,class LO,class GO,class NO>
typename Preconditioner<SC,LO,GO,NO>::ThyraPrecConstPtr_Type Preconditioner<SC,LO,GO,NO>::getThyraPrecConst() const{
    return thyraPrec_;
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::initializePreconditioner( std::string type )
{
    if ( type == "Monolithic" || type == "FaCSI" || type == "Diagonal" || type == "Triangular"){
        if (type == "Monolithic")
            initPreconditionerMonolithic( );
        else if (type == "FaCSI" || type == "Diagonal" || type == "Triangular" || type == "PCD" || type == "LSC")
            initPreconditionerBlock( );
        
    }
    else if (type == "Teko" || type == "FaCSI-Teko" || type == "FaCSI-Block"){
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Please construct the Teko precondtioner completely.");
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner type for initialization.");
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setPressureProjection(BlockMultiVectorPtr_Type pressureProjection) const{
    pressureProjection_ = pressureProjection;
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::initPreconditionerMonolithic( )
{

    LinSolverBuilderPtr_Type solverBuilder;
    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > thyraRangeSpace;
    Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > thyraDomainSpace;

    if (!problem_.is_null()){
        solverBuilder = problem_->getLinearSolverBuilder();
        thyraRangeSpace = Thyra::tpetraVectorSpace<SC,LO,GO,NO>( problem_->getSystem()->getMap()->getMergedMap()->getTpetraMap()); //Tpetra::ThyraUtils<SC,LO,GO,NO>::toThyra( problem_->getSystem()->getMap()->getMergedMap()->getTpetraMap() );
        thyraDomainSpace = Thyra::tpetraVectorSpace<SC,LO,GO,NO>( problem_->getSystem()->getMap()->getMergedMap()->getTpetraMap()); //Tpetra::ThyraUtils<SC,LO,GO,NO>::toThyra( problem_->getSystem()->getMap()->getMergedMap()->getTpetraMap() );
    }
    else if(!timeProblem_.is_null()){
        solverBuilder = timeProblem_->getUnderlyingProblem()->getLinearSolverBuilder();
        thyraRangeSpace = Thyra::tpetraVectorSpace<SC,LO,GO,NO>( timeProblem_->getSystem()->getMap()->getMergedMap()->getTpetraMap() );
        thyraDomainSpace = Thyra::tpetraVectorSpace<SC,LO,GO,NO>( timeProblem_->getSystem()->getMap()->getMergedMap()->getTpetraMap() );

    }
    
    Teuchos::RCP<Thyra::LinearOpBase<SC> > thyraLinOp = Teuchos::rcp( new Thyra::DefaultZeroLinearOp<SC>( thyraRangeSpace, thyraDomainSpace ) );

    Teuchos::RCP<Thyra::PreconditionerFactoryBase<SC> > precFactory = solverBuilder->createPreconditioningStrategy("");
    thyraPrec_ = precFactory->createPrec();

    Teuchos::RCP<Thyra::DefaultPreconditioner<SC> > defaultPrec = Teuchos::rcp_dynamic_cast<Thyra::DefaultPreconditioner<SC> >(thyraPrec_);

    defaultPrec->initializeUnspecified(thyraLinOp);
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::initPreconditionerBlock( )
{
    using Teuchos::Array;
    using Teuchos::RCP;
    using Teuchos::rcp;
    BlockMultiVectorPtr_Type sol;
    BlockMultiVectorPtr_Type rhs;
    LinSolverBuilderPtr_Type solverBuilder;
    int size = 0;
    if (!problem_.is_null()){
        solverBuilder = problem_->getLinearSolverBuilder();
        size = problem_->getSystem()->size();
        sol = problem_->getSolution();
        rhs = problem_->getRhs();
    }
    else if(!timeProblem_.is_null()){
        solverBuilder = timeProblem_->getUnderlyingProblem()->getLinearSolverBuilder();
        size = timeProblem_->getSystem()->size();
        sol = timeProblem_->getSolution();
        rhs = timeProblem_->getRhs();
    }
    
    Array< RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesRange( size );
    Array< RCP< const Thyra::VectorSpaceBase< SC > > > vectorSpacesDomain( size );

    for (int i=0; i<size; i++) {
        vectorSpacesRange[i] = rhs->getBlock(i)->getMap()->getThyraVectorSpaceBase();
        vectorSpacesDomain[i] = sol->getBlock(i)->getMap()->getThyraVectorSpaceBase();
    }

    RCP<const Thyra::DefaultProductVectorSpace<SC> > pR = Thyra::productVectorSpace<SC>( vectorSpacesRange );
    RCP<const Thyra::DefaultProductVectorSpace<SC> > pD = Thyra::productVectorSpace<SC>( vectorSpacesDomain );
    
    RCP<Thyra::LinearOpBase<SC> > thyraLinOp = rcp( new Thyra::DefaultZeroLinearOp<SC>( pR, pD ) );
    
    Teuchos::RCP<Thyra::PreconditionerFactoryBase<SC> > precFactory = solverBuilder->createPreconditioningStrategy("");
    thyraPrec_ = precFactory->createPrec();
    
    Teuchos::RCP<Thyra::DefaultPreconditioner<SC> > defaultPrec = Teuchos::rcp_dynamic_cast<Thyra::DefaultPreconditioner<SC> >(thyraPrec_);
    
    defaultPrec->initializeUnspecified(thyraLinOp);
}
    
template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::buildPreconditioner( std::string type )
{

#ifdef PRECONDITIONER_TIMER
    CommConstPtr_Type comm;
    if (!problem_.is_null())
        comm = problem_->getComm();
    else if(!timeProblem_.is_null())
        comm = timeProblem_->getComm();
    comm->barrier();
    Teuchos::TimeMonitor preconditionerTimeMonitor(*preconditionerTimer_);
#endif

    if (!type.compare("Monolithic")){
        buildPreconditionerMonolithic( );
    }
    else if (!type.compare("Teko")){
#ifdef FEDD_HAVE_TEKO
        buildPreconditionerTeko( );
#else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Teko not found! Build Trilinos with Teko.");
#endif
    }
    else if( type == "FaCSI" || type == "FaCSI-Teko" || type == "FaCSI-Block" ){
        buildPreconditionerFaCSI( type );
    }
    else if(type == "Triangular" || type == "Diagonal" || type == "PCD" || type == "LSC"){
        buildPreconditionerBlock2x2( );
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner type.");

#ifdef PRECONDITIONER_TIMER
    comm->barrier();
#endif
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::buildPreconditionerMonolithic( )
{

    CommConstPtr_Type comm;
    if (!problem_.is_null())
        comm = problem_->getComm();
    else if(!timeProblem_.is_null())
        comm = timeProblem_->getComm();
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Preconditioner can not be used without a problem.");

    bool verbose ( comm->getRank() == 0 );
    
    ParameterListPtr_Type parameterList;

    if (!problem_.is_null())
        parameterList = problem_->getParameterList();
    else if(!timeProblem_.is_null())
        parameterList = timeProblem_->getParameterList();

    std::string precType = parameterList->sublist("ThyraPreconditioner").get("Preconditioner Type", "FROSch");

    bool useRepeatedMaps = parameterList->get( "Use repeated maps", true );
    bool useNodeLists = parameterList->get( "Use node lists", true );
    ParameterListPtr_Type pListThyraPrec = sublist( parameterList, "ThyraPreconditioner" );
    ParameterListPtr_Type plFrosch = sublist( sublist( pListThyraPrec, "Preconditioner Types" ), "FROSch");
    
    ThyraLinOpConstPtr_Type thyraMatrix;
    if (!problem_.is_null())
        thyraMatrix = problem_->getSystem()->getThyraLinOp();
    else if(!timeProblem_.is_null())
        thyraMatrix = timeProblem_->getSystemCombined()->getThyraLinOp();

    UN numberOfBlocks = parameterList->get("Number of blocks",1);

    typedef Tpetra::MultiVector<SC,LO,GO,NO> TMultiVector;
    typedef Teuchos::RCP<TMultiVector> TMultiVectorPtr;
    typedef Teuchos::ArrayRCP<TMultiVectorPtr> TMultiVectorPtrVecPtr;
    
    // ------------------------
    // Defs to cast back from tpetra to xpetra
    typedef Xpetra::MultiVector<SC,LO,GO,NO> XMultiVector;
    typedef Teuchos::RCP<XMultiVector> XMultiVectorPtr;
    typedef Teuchos::ArrayRCP<XMultiVectorPtr> XMultiVectorPtrVecPtr;

    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const XpetraMapConstPtr_Type XpetraMapConstPtrConst_Type;
    // ---------
    // XMapVecPtrVecPtr
    //Teuchos::ArrayRCP<Teuchos::RCP<Tpetra::Map<LO,GO,NO> > > repeatedMaps(numberOfBlocks);
    Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > repeatedMaps(numberOfBlocks);
    // --------

    XMultiVectorPtrVecPtr nodeListVec( numberOfBlocks );
    if (!useNodeLists)
        nodeListVec = Teuchos::null;

   // timeProblem_->getSystemCombined()->print();
    //Set Precondtioner lists
    if (!precondtionerIsBuilt_) {
        if ( !precType.compare("FROSch") ){
            Teuchos::ArrayRCP<FROSch::DofOrdering> dofOrderings(numberOfBlocks);
            Teuchos::ArrayRCP<UN> dofsPerNodeVector(numberOfBlocks);
            for (UN i = 0; i < numberOfBlocks; i++) {
                dofsPerNodeVector[i] = (UN) pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get( "DofsPerNode" + std::to_string(i+1), 1);

                if (useRepeatedMaps) {
                    if (dofsPerNodeVector[i] > 1){

                        if (!problem_.is_null()) {
                            if (problem_->getDomain(i)->getFEType() == "P0") {
                                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Vector field map not implemented for P0 elements.");
                            }

                            //Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = problem_->getDomain(i)->getMapVecFieldRepeated()->getTpetraMap();
                            //Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);

                            MapConstPtr_Type mapConstTmp = problem_->getDomain(i)->getMapVecFieldRepeated();//->getTpetraMap();
                            XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);
                            
                            repeatedMaps[i] = mapX;


                        }
                        else if(!timeProblem_.is_null()){
                            if (timeProblem_->getDomain(i)->getFEType() == "P0") {
                                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Vector field map not implemented for P0 elements.");
                            }
                            // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = timeProblem_->getDomain(i)->getMapVecFieldRepeated()->getTpetraMap();
                            // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                            
                            MapConstPtr_Type mapConstTmp = timeProblem_->getDomain(i)->getMapVecFieldRepeated();//->getTpetraMap();
                            XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);
                            
                            repeatedMaps[i] = mapX;
                        }
                    }
                    else{
                        if (!problem_.is_null()){
                            if (problem_->getDomain(i)->getFEType() == "P0") {
                                 // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = problem_->getDomain(i)->getElementMap()->getTpetraMap();
                                // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                MapConstPtr_Type mapConstTmp = problem_->getDomain(i)->getElementMap();//->getTpetraMap();
                                XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

                                repeatedMaps[i] = mapX;
                            }
                            else{
                                // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = problem_->getDomain(i)->getMapRepeated()->getTpetraMap();
                                // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                MapConstPtr_Type mapConstTmp = problem_->getDomain(i)->getMapRepeated();//->getTpetraMap();
                                XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);
                                
                                repeatedMaps[i] = mapX;
                            }
                        }
                        else if (!timeProblem_.is_null()){
                            if (timeProblem_->getDomain(i)->getFEType() == "P0") {
                                 // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = timeProblem_->getDomain(i)->getElementMap()->getTpetraMap();
                                // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                MapConstPtr_Type mapConstTmp = timeProblem_->getDomain(i)->getElementMap();//->getTpetraMap();
                                XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

                                repeatedMaps[i] = mapX;
                            }
                            else{
                                // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = timeProblem_->getDomain(i)->getMapRepeated()->getTpetraMap();
                                // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                MapConstPtr_Type mapConstTmp = timeProblem_->getDomain(i)->getMapRepeated();//->getTpetraMap();
                                XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

                                repeatedMaps[i] = mapX;
                            }
                        }
                    }
                }
                
                if (useNodeLists) {
                    if (!problem_.is_null()){
                        TEUCHOS_TEST_FOR_EXCEPTION( problem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Node lists cannot be used for P0 elements." );
                        Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > nodeListTpetra =  problem_->getDomain(i)->getNodeListMV()->getTpetraMultiVectorNonConst();
                        Teuchos::RCP< Xpetra::TpetraMultiVector<SC,LO,GO,NO> > nodeListXpetraTpetra = Teuchos::rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(nodeListTpetra));
                        Teuchos::RCP< Xpetra::MultiVector<SC,LO,GO,NO> > nodeListXpetra = Teuchos::rcp_dynamic_cast<Xpetra::MultiVector<SC,LO,GO,NO>>(nodeListXpetraTpetra);
         
                        nodeListVec[i] = nodeListXpetra; //problem_->getDomain(i)->getNodeListMV()->getTpetraMultiVectorNonConst();
                    }
                    else if (!timeProblem_.is_null()){
                        TEUCHOS_TEST_FOR_EXCEPTION( timeProblem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Node lists cannot be used for P0 elements." );
                        Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > nodeListTpetra =  timeProblem_->getDomain(i)->getNodeListMV()->getTpetraMultiVectorNonConst();
                        Teuchos::RCP< Xpetra::TpetraMultiVector<SC,LO,GO,NO> > nodeListXpetraTpetra = Teuchos::rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(nodeListTpetra));
                        Teuchos::RCP< Xpetra::MultiVector<SC,LO,GO,NO> > nodeListXpetra = Teuchos::rcp_dynamic_cast<Xpetra::MultiVector<SC,LO,GO,NO>>(nodeListXpetraTpetra);
         
                        nodeListVec[i] =  nodeListXpetra;//timeProblem_->getDomain(i)->getNodeListMV()->getTpetraMultiVectorNonConst();
                    }
                    
                }
                
                
                if (!pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get( "DofOrdering" + std::to_string(i+1), "NodeWise" ).compare("DimensionWise"))
                    dofOrderings[i] = FROSch::DimensionWise;
                else if (!pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get( "DofOrdering" + std::to_string(i+1), "NodeWise" ).compare("NodeWise"))
                    dofOrderings[i] = FROSch::NodeWise;
                else
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Please chose a valid  DOF ordering for ThyraPreconditioner(FROSch).");
            }
            int dim;
            if (!problem_.is_null())
                dim = problem_->getDomain(0)->getDimension();
            else if(!timeProblem_.is_null())
                dim = timeProblem_->getDomain(0)->getDimension();

            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Dimension", dim);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Repeated Map Vector",repeatedMaps);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Coordinates List Vector",nodeListVec);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("DofOrdering Vector",dofOrderings);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("DofsPerNode Vector",dofsPerNodeVector);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Mpi Ranks Coarse",parameterList->sublist("General").get("Mpi Ranks Coarse",0) );

            // This a pressure projection is only used for saddle point problems. We check here if we have a pressure projection set and if we have more than one block or one block with dim dof per node (i.e. fluid problem)
            // This is unfortunately called pressure correction in FROSch, but is a projection!
            if(!pressureProjection_.is_null() && ( dofsPerNodeVector.size() > 1 || dofsPerNodeVector[0] == 1) ){
                pressureProjection_->merge(); // We merge the projection vector, as FROSch does not distinguish between blocks

                Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > vectorTpetra =  pressureProjection_->getMergedVectorNonConst()->getTpetraMultiVectorNonConst();
                Teuchos::RCP< Xpetra::TpetraMultiVector<SC,LO,GO,NO> > vectorXpetraTpetra = Teuchos::rcp(new Xpetra::TpetraMultiVector<SC,LO,GO,NO>(vectorTpetra));
                Teuchos::RCP< Xpetra::MultiVector<SC,LO,GO,NO> > vectorXpetra = Teuchos::rcp_dynamic_cast<Xpetra::MultiVector<SC,LO,GO,NO>>(vectorXpetraTpetra);

                pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").sublist("AlgebraicOverlappingOperator").set("Projection",vectorXpetra);
                // In case of pressure correction we set the parameter in the parameterlist to true
                pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").sublist("AlgebraicOverlappingOperator").set("Use Pressure Correction", true);
                pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").sublist("AlgebraicOverlappingOperator").set("Use Local Pressure Correction", true);
            }

            /*  We need to set the ranges of local problems and the coarse problem here.
                When using an unstructured decomposition of, e.g., FSI, with 2 domains, which might be on a different set of ranks, we need to set the following parameters for FROSch here. Similarly we need to set a coarse rank problem range. For now, we use extra coarse ranks only for structured decompositions
             */
            
            int lowerBound = 100000000;
            int upperBound = -1;
            for (UN i = 0; i < numberOfBlocks; i++) {
                tuple_intint_Type rankRange;
                if (!problem_.is_null())
                    rankRange = problem_->getDomain(i)->getMesh()->getRankRange();
                else if(!timeProblem_.is_null())
                    rankRange = timeProblem_->getDomain(i)->getMesh()->getRankRange();
                
                if (std::get<0>(rankRange) < lowerBound)
                    lowerBound = std::get<0>(rankRange);
                if (std::get<1>(rankRange) > upperBound)
                    upperBound = std::get<1>(rankRange);
            }
            
            int lowerBoundCoarse = lowerBound;
            int upperBoundCoarse = upperBound;
            // For now only for structured decompositions
            if ( parameterList->sublist("General").get("Mpi Ranks Coarse",0) > 0 ){
                lowerBoundCoarse = upperBound + 1;
                upperBoundCoarse = comm->getSize() - 1;
            }
            if (verbose) {
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
                std::cout << "\t --- Range for local problems of preconditioner from " << lowerBound << " to " << upperBound << std::endl;
                std::cout << "\t --- Range for coarse problem of preconditioner from " << lowerBoundCoarse << " to " << upperBoundCoarse << std::endl;
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
            }

            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Local problem ranks lower bound", lowerBound );
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Local problem ranks upper bound", upperBound );
            
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Coarse problem ranks lower bound", lowerBoundCoarse );
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Coarse problem ranks upper bound", upperBoundCoarse );
            
    //        if ( !pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive").compare("Multiplicative") ){
    //            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);
    //        }

            // Here we set the parameterlist which is a member variable of class Problem
            ParameterListPtr_Type pListThyraSolver = sublist( parameterList, "ThyraSolver" ); //ch 11.02.19 To we need full Solver details
            pListThyraSolver->setParameters( *pListThyraPrec );

            LinSolverBuilderPtr_Type solverBuilder;
            if (!problem_.is_null())
                solverBuilder = problem_->getLinearSolverBuilder();
            else if(!timeProblem_.is_null())
                solverBuilder = timeProblem_->getLinearSolverBuilder();

            // We save the pointer to the coarse matrix in this parameter list inside FROSch
            pListPhiExport_ = pListThyraSolver;
            
            solverBuilder->setParameterList(pListThyraSolver);
            precFactory_ = solverBuilder->createPreconditioningStrategy("");

            if ( thyraPrec_.is_null() ){
                thyraPrec_ = precFactory_->createPrec();
           	}     
           
            Thyra::initializePrec<SC>(*precFactory_, thyraMatrix, thyraPrec_.ptr()); // (precfactory, fwdOp, prec) Problem: PreconditionerBase<SC>* thyraPrec_
            precondtionerIsBuilt_ = true;
            
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unknown preconditioner type; You can only compute a FROSch preconditioner here." );
    }
    else {
        TEUCHOS_TEST_FOR_EXCEPTION( precFactory_.is_null(), std::runtime_error, "precFactory_ is null.");
        Thyra::initializePrec<SC>(*precFactory_, thyraMatrix, thyraPrec_.ptr());
    }

    if (parameterList->sublist( "Exporter" ).get("Export coarse functions",false))
            exportCoarseBasis();
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::buildPreconditionerMonolithicFSI( )
{

    CommConstPtr_Type comm;
    if (!problem_.is_null())
        comm = problem_->getComm();
    else if(!timeProblem_.is_null())
        comm = timeProblem_->getComm();
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Preconditioner can not be used without a problem.");

    bool verbose ( comm->getRank() == 0 );
    
    ParameterListPtr_Type parameterList;

    if (!problem_.is_null())
        parameterList = problem_->getParameterList();
    else if(!timeProblem_.is_null())
        parameterList = timeProblem_->getParameterList();

    std::string precType = parameterList->sublist("ThyraPreconditioner").get("Preconditioner Type", "FROSch");

    bool useRepeatedMaps = parameterList->get( "Use repeated maps", true );
    bool useNodeLists = parameterList->get( "Use node lists", true );
    ParameterListPtr_Type pListThyraPrec = sublist( parameterList, "ThyraPreconditioner" );
    ParameterListPtr_Type plFrosch = sublist( sublist( pListThyraPrec, "Preconditioner Types" ), "FROSch");
//    setRecyclingParameter( plFrosch );
    
    
    ThyraLinOpConstPtr_Type thyraMatrix;
    if (!problem_.is_null())
        thyraMatrix = problem_->getSystem()->getThyraLinOp();
    else if(!timeProblem_.is_null())
        thyraMatrix = timeProblem_->getSystemCombined()->getThyraLinOp();

    UN numberOfBlocks = parameterList->get("Number of blocks",1);
    TEUCHOS_TEST_FOR_EXCEPTION( numberOfBlocks<4 || numberOfBlocks>5, std::logic_error, "Unknown FSI size." );
    

    typedef Tpetra::MultiVector<SC,LO,GO,NO> TMultiVector;
    typedef Teuchos::RCP<TMultiVector> TMultiVectorPtr;
    typedef Teuchos::ArrayRCP<TMultiVectorPtr> TMultiVectorPtrVecPtr;
    
    // ------------------------
    // Defs to cast back from tpetra to xpetra
    typedef Xpetra::MultiVector<SC,LO,GO,NO> XMultiVector;
    typedef Teuchos::RCP<XMultiVector> XMultiVectorPtr;
    typedef Teuchos::ArrayRCP<XMultiVectorPtr> XMultiVectorPtrVecPtr;

    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const XpetraMapConstPtr_Type XpetraMapConstPtrConst_Type;

    // XMapPtrVecPtr
    Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > repeatedMaps(numberOfBlocks);
    // ------------------------

    TMultiVectorPtrVecPtr nodeListVec( numberOfBlocks );
    nodeListVec = Teuchos::null;

    //Set Precondtioner lists
    if (!precondtionerIsBuilt_) {
        if ( !precType.compare("FROSch") ){
            Teuchos::ArrayRCP<FROSch::DofOrdering> dofOrderings(numberOfBlocks);
            Teuchos::ArrayRCP<UN> dofsPerNodeVector(numberOfBlocks);
            for (UN i = 0; i < numberOfBlocks; i++) {
                dofsPerNodeVector[i] = (UN) pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get( "DofsPerNode" + std::to_string(i+1), 1);
                if (i==3) { //interface coupling
                    TEUCHOS_TEST_FOR_EXCEPTION( timeProblem_.is_null(), std::logic_error, "FSI time problem is null!" );
                    TEUCHOS_TEST_FOR_EXCEPTION( timeProblem_->getDomain(i)->getFEType() == "P0", std::logic_error, "We should not be able to use P0 for interface coupling." );
                    //Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp timeProblem_->getDomain(i)->getInterfaceMapUnique()->getTpetraMap();
                    //Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                    MapConstPtr_Type mapConstTmp =timeProblem_->getDomain(i)->getInterfaceMapUnique();//->getTpetraMap();
                    XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

                    repeatedMaps[i] = mapX;
                }
                else {
                    if (useRepeatedMaps) {
                        if (dofsPerNodeVector[i] > 1){

                            if (!problem_.is_null()) {
                                if (problem_->getDomain(i)->getFEType() == "P0") {
                                    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Vector field map not implemented for P0 elements.");
                                }

                                //Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = problem_->getDomain(i)->getMapVecFieldRepeated()->getTpetraMap();
                                //Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);

                                MapConstPtr_Type mapConstTmp = problem_->getDomain(i)->getMapVecFieldRepeated();//->getTpetraMap();
                                XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);
                               
                                repeatedMaps[i] = mapX;

                            }
                            else if(!timeProblem_.is_null()){
                                if (timeProblem_->getDomain(i)->getFEType() == "P0") {
                                    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Vector field map not implemented for P0 elements.");
                                }
                                // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = timeProblem_->getDomain(i)->getMapVecFieldRepeated()->getTpetraMap();
                                // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                
                                MapConstPtr_Type mapConstTmp = timeProblem_->getDomain(i)->getMapVecFieldRepeated();//->getTpetraMap();
                                XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);
                                
                                repeatedMaps[i] = mapX;
                            }
                        }
                        else{
                            if (!problem_.is_null()){
                                if (problem_->getDomain(i)->getFEType() == "P0") {
                                    // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = problem_->getDomain(i)->getElementMap()->getTpetraMap();
                                    // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                    MapConstPtr_Type mapConstTmp = problem_->getDomain(i)->getElementMap();//->getTpetraMap();
                                    XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

                                   
                                    
                                    repeatedMaps[i] = mapX;
                                }
                                else{
                                    // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = problem_->getDomain(i)->getMapRepeated()->getTpetraMap();
                                    // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                    MapConstPtr_Type mapConstTmp = problem_->getDomain(i)->getMapRepeated();//->getTpetraMap();
                                    XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);
                                    
                                    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
                                    mapX->describe(*out,Teuchos::VERB_EXTREME);
                                    repeatedMaps[i] = mapX;
                                }
                            }
                            else if (!timeProblem_.is_null()){
                                if (timeProblem_->getDomain(i)->getFEType() == "P0") {
                                    // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = timeProblem_->getDomain(i)->getElementMap()->getTpetraMap();
                                    // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                    MapConstPtr_Type mapConstTmp = timeProblem_->getDomain(i)->getElementMap();//->getTpetraMap();
                                    XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

                                    repeatedMaps[i] = mapX;
                                }
                                else{
                                    // Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp = timeProblem_->getDomain(i)->getMapRepeated()->getTpetraMap();
                                    // Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
                                    MapConstPtr_Type mapConstTmp = timeProblem_->getDomain(i)->getMapRepeated();//->getTpetraMap();
                                    XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
                                    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

                                    repeatedMaps[i] = mapX;
                                }
                            }
                        }
                    }
                    
                    if (!pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get( "DofOrdering" + std::to_string(i+1), "NodeWise" ).compare("DimensionWise"))
                        dofOrderings[i] = FROSch::DimensionWise;
                    else if (!pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get( "DofOrdering" + std::to_string(i+1), "NodeWise" ).compare("NodeWise"))
                        dofOrderings[i] = FROSch::NodeWise;
                    else
                        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Please chose a valid  DOF ordering for ThyraPreconditioner(FROSch).");
                        
                }
            }
            int dim;
            if (!problem_.is_null())
                dim = problem_->getDomain(0)->getDimension();
            else if(!timeProblem_.is_null())
                dim = timeProblem_->getDomain(0)->getDimension();

            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Dimension", dim);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Repeated Map Vector",repeatedMaps);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Coordinates List Vector",nodeListVec);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("DofOrdering Vector",dofOrderings);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("DofsPerNode Vector",dofsPerNodeVector);
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Mpi Ranks Coarse",parameterList->sublist("General").get("Mpi Ranks Coarse",0) );

            /*  We need to set the ranges of local problems and the coarse problem here.
                When using an unstructured decomposition of, e.g., FSI, with 2 domains, which might be on a different set of ranks, we need to set the following parameters for FROSch here. Similarly we need to set a coarse rank problem range. For now, we use extra coarse ranks only for structured decompositions
             */
            
            int lowerBound = 100000000;
            int upperBound = -1;
            for (UN i = 0; i < numberOfBlocks; i++) {
                tuple_intint_Type rankRange;
                if (!problem_.is_null())
                    rankRange = problem_->getDomain(i)->getMesh()->getRankRange();
                else if(!timeProblem_.is_null())
                    rankRange = timeProblem_->getDomain(i)->getMesh()->getRankRange();
                
                if (std::get<0>(rankRange) < lowerBound)
                    lowerBound = std::get<0>(rankRange);
                if (std::get<1>(rankRange) > upperBound)
                    upperBound = std::get<1>(rankRange);
            }
            
            int lowerBoundCoarse = lowerBound;
            int upperBoundCoarse = upperBound;
            // For now only for structured decompositions
            if ( parameterList->sublist("General").get("Mpi Ranks Coarse",0) > 0 ){
                lowerBoundCoarse = upperBound + 1;
                upperBoundCoarse = comm->getSize() - 1;
            }
            if (verbose) {
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
                std::cout << "\t --- Range for local problems of preconditioner from " << lowerBound << " to " << upperBound << std::endl;
                std::cout << "\t --- Range for coarse problem of preconditioner from " << lowerBoundCoarse << " to " << upperBoundCoarse << std::endl;
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
            }

            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Local problem ranks lower bound", lowerBound );
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Local problem ranks upper bound", upperBound );
            
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Coarse problem ranks lower bound", lowerBoundCoarse );
            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set( "Coarse problem ranks upper bound", upperBoundCoarse );
            
    //        if ( !pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").get("Level Combination","Additive").compare("Multiplicative") ){
    //            pListThyraPrec->sublist("Preconditioner Types").sublist("FROSch").set("Only apply coarse",true);
    //        }

            // Here we set the parameterlist which is a member variable of class Problem
            ParameterListPtr_Type pListThyraSolver = sublist( parameterList, "ThyraSolver" ); //ch 11.02.19 To we need full Solver details
            pListThyraSolver->setParameters( *pListThyraPrec );

            LinSolverBuilderPtr_Type solverBuilder;
            if (!problem_.is_null())
                solverBuilder = problem_->getLinearSolverBuilder();
            else if(!timeProblem_.is_null())
                solverBuilder = timeProblem_->getUnderlyingProblem()->getLinearSolverBuilder();

            // We save the pointer to the coarse matrix in this parameter list inside FROSch
            pListPhiExport_ = pListThyraSolver;
            
            solverBuilder->setParameterList(pListThyraSolver);
            precFactory_ = solverBuilder->createPreconditioningStrategy("");

            if ( thyraPrec_.is_null() )
                thyraPrec_ = precFactory_->createPrec();


            Thyra::initializePrec<SC>(*precFactory_, thyraMatrix, thyraPrec_.ptr());
            precondtionerIsBuilt_ = true;
            
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unknown preconditioner type; You can only compute a FROSch preconditioner here." );
    }
    else {
        TEUCHOS_TEST_FOR_EXCEPTION( precFactory_.is_null(), std::runtime_error, "precFactory_ is null.");
        Thyra::initializePrec<SC>(*precFactory_, thyraMatrix, thyraPrec_.ptr());
    }

    if (parameterList->sublist( "Exporter" ).get("Export coarse functions",false))
        exportCoarseBasisFSI();

}


#ifdef FEDD_HAVE_TEKO
template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::buildPreconditionerTeko( )
{

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    ParameterListPtr_Type parameterList;
    if (!problem_.is_null())
        parameterList = problem_->getParameterList();
    else if(!timeProblem_.is_null())
        parameterList = timeProblem_->getParameterList();
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Preconditioner can not be used without a problem.");

    
    LinSolverBuilderPtr_Type solverBuilder;
    if (!problem_.is_null())
        solverBuilder = problem_->getLinearSolverBuilder();
    else if(!timeProblem_.is_null())
        solverBuilder = timeProblem_->getUnderlyingProblem()->getLinearSolverBuilder();

    ParameterListPtr_Type tekoPList= sublist( parameterList, "Teko Parameters" );

    if (precFactory_.is_null()) {
        ParameterListPtr_Type tmpSubList = sublist( sublist( sublist( sublist( parameterList, "Teko Parameters" ) , "Preconditioner Types" ) , "Teko" ) , "Inverse Factory Library" );

        //only sets repeated maps in parameterlist if FROSch is used for both block
        if ( !tmpSubList->sublist("FROSch-Pressure").get("Type","FROSch").compare("FROSch") &&
             !tmpSubList->sublist("FROSch-Velocity").get("Type","FROSch").compare("FROSch") ) {
            setVelocityParameters( tekoPList, parameterList->sublist("General").get("Mpi Ranks Coarse",0) );
            setPressureParameters( tekoPList, parameterList->sublist("General").get("Mpi Ranks Coarse",0) );
        }
    }

    BlockMatrixPtr_Type system;
    if (!problem_.is_null())
        system = problem_->getSystem();
    else if(!timeProblem_.is_null())
        system = timeProblem_->getSystemCombined();

    TEUCHOS_TEST_FOR_EXCEPTION( system->size()!=2, std::logic_error, "Wrong size of system for Teko-Block-Preconditioners.");

    Teko::LinearOp thyraF = system->getBlock(0,0)->getThyraLinOp();
    Teko::LinearOp thyraBT = system->getBlock(0,1)->getThyraLinOp();
    Teko::LinearOp thyraB = system->getBlock(1,0)->getThyraLinOp();

    if (!system->blockExists(1,1)){
        MatrixPtr_Type dummy;
        dummy.reset( new Matrix_Type( system->getBlock(1,0)->getMap(), 1 ) );
        dummy->fillComplete();
        system->addBlock( dummy, 1, 1 );
    }
    Teko::LinearOp thyraC = system->getBlock(1,1)->getThyraLinOp();

    tekoLinOp_ = Thyra::block2x2(thyraF,thyraBT,thyraB,thyraC);

    if (!precondtionerIsBuilt_) {
        if ( precFactory_.is_null() ){
            ParameterListPtr_Type pListThyraSolver = sublist( parameterList, "ThyraSolver" );

            pListThyraSolver->setParameters( *tekoPList );

            solverBuilder->setParameterList( pListThyraSolver );
            precFactory_ = solverBuilder->createPreconditioningStrategy("");//createPreconditioningStrategy(*solverBuilder);
            
            rh_.reset(new Teko::RequestHandler());
            
            if(!tekoPList->sublist("Preconditioner Types").sublist("Teko").get("Inverse Type", "SIMPLE").compare("LSC") || 
                !tekoPList->sublist("Preconditioner Types").sublist("Teko").get("Inverse Type", "SIMPLE").compare("LSC-Pressure-Laplace")  || 
                !tekoPList->sublist("Preconditioner Types").sublist("Teko").get("Inverse Type", "SIMPLE").compare("SIMPLE")){
                Teko::LinearOp thyraMass = velocityMassMatrix_;
                Teuchos::RCP< Teko::StaticRequestCallback<Teko::LinearOp> > callbackMass = Teuchos::rcp(new Teko::StaticRequestCallback<Teko::LinearOp> ( "Velocity Mass Matrix", thyraMass ) );
                rh_->addRequestCallback( callbackMass );

                Teko::LinearOp thyraLaplace = pressureLaplace_;

                Teuchos::RCP< Teko::StaticRequestCallback<Teko::LinearOp> > callbackLaplace = Teuchos::rcp(new Teko::StaticRequestCallback<Teko::LinearOp> ( "Pressure Laplace Operator", thyraLaplace ) );
                rh_->addRequestCallback( callbackLaplace );
            }
            else if(!tekoPList->sublist("Preconditioner Types").sublist("Teko").get("Inverse Type", "SIMPLE").compare("PCD")){

                // Velocity Mass Matrix
                Teko::LinearOp thyraMass = velocityMassMatrix_;
                Teuchos::RCP< Teko::StaticRequestCallback<Teko::LinearOp> > callbackMass = Teuchos::rcp(new Teko::StaticRequestCallback<Teko::LinearOp> ( "Velocity Mass Matrix", thyraMass ) );
                rh_->addRequestCallback( callbackMass );

                // Pressure Laplace
                Teko::LinearOp thyraLaplace = pressureLaplace_;
                Teuchos::RCP< Teko::StaticRequestCallback<Teko::LinearOp> > callbackLaplace = Teuchos::rcp(new Teko::StaticRequestCallback<Teko::LinearOp> ( "Pressure Laplace Operator", thyraLaplace ) );
                rh_->addRequestCallback( callbackLaplace );

                // Pressure Mass
                Teko::LinearOp thyraPressureMass = pressureMass_;
                Teuchos::RCP< Teko::StaticRequestCallback<Teko::LinearOp> > callbackPressureMass = Teuchos::rcp(new Teko::StaticRequestCallback<Teko::LinearOp> ( "Pressure Mass Matrix", thyraPressureMass ) );
                rh_->addRequestCallback( callbackPressureMass );

                // PCD
                if (!timeProblem_.is_null()){
                    timeProblem_->assemble("UpdateConvectionDiffusionOperator");
                }
                else{
                    problem_->assemble("UpdateConvectionDiffusionOperator");
                }
                Teko::LinearOp thyraPCD = pcdOperator_;
                callbackPCD_ = Teuchos::rcp(new Teko::StaticRequestCallback<Teko::LinearOp> ( "PCD Operator", thyraPCD ) );
                rh_->addRequestCallback( callbackPCD_ );

            }
            Teuchos::RCP< Teko::StratimikosFactory > tekoFactory = Teuchos::rcp_dynamic_cast<Teko::StratimikosFactory>(precFactory_);
            tekoFactory->setRequestHandler( rh_ );
        }
        

        if ( thyraPrec_.is_null() ){
            thyraPrec_ = precFactory_->createPrec();
        }   

        Teuchos::RCP< const Thyra::DefaultLinearOpSource< SC > > thyraMatrixSourceOp =  defaultLinearOpSource (tekoLinOp_);
        //    Thyra::initializePrec<SC>(*precFactory, thyraMatrixSourceOp, thyraPrec_.ptr());

        precFactory_->initializePrec(thyraMatrixSourceOp, thyraPrec_.get());
        precondtionerIsBuilt_ = true;
        
    }
    else{
        if(!tekoPList->sublist("Preconditioner Types").sublist("Teko").get("Inverse Type", "SIMPLE").compare("PCD")){
            // PCD: As part of the pcd preconditioner depends on the current velocity, we need to update it in each iteration.
            //pcdOperatorMatrixPtr_->print();
                // PCD
            if (!timeProblem_.is_null()){
                timeProblem_->assemble("UpdateConvectionDiffusionOperator");
            }
            else{
                problem_->assemble("UpdateConvectionDiffusionOperator");
            }      
            Teuchos::RCP<Teko::RequestHandler> rh = Teuchos::rcp(new Teko::RequestHandler());

            // PCD
            Teko::LinearOp thyraPCD;
            
            // When we deal with a time problem only the underlying problem holds the updated pcd matrix. Why: When we assemble the PCD operator it within i.e. the Navier Stokes class. 
            // The Navier stokes class is derived from a nonlinear problem originally derived from a problem which has a precondidioner object to which the PCD operator is added in each reassembly.
            // Since the Navier Stokes problem does not know the timeProblem, we need to extract the information from the problem which is added to the time problem / which the timeProblem is based on.
            if(!timeProblem_.is_null())
            {
              // timeProblem_->getUnderlyingProblem()->preconditioner_->getPCDOperatorMatrix()->print(); 
              pcdOperator_ = timeProblem_->getUnderlyingProblem()->preconditioner_->getPCDOperatorMatrix()->getThyraLinOp();
              thyraPCD = timeProblem_->getUnderlyingProblem()->preconditioner_->getPCDOperatorMatrix()->getThyraLinOp();     
            }
            else
                thyraPCD= pcdOperator_;

            // Updating matrix in the pointer
            Teuchos::RCP< Teko::StaticRequestCallback<Teko::LinearOp> > callbackTmp = Teuchos::rcp(new Teko::StaticRequestCallback<Teko::LinearOp> ( "PCD Operator", thyraPCD ) );
            *callbackPCD_ = *callbackTmp;
           
         }

        Teuchos::RCP< const Thyra::DefaultLinearOpSource< SC > > thyraMatrixSourceOp =  defaultLinearOpSource (tekoLinOp_);
        //    Thyra::initializePrec<SC>(*precFactory, thyraMatrixSourceOp, thyraPrec_.ptr());

        precFactory_->initializePrec(thyraMatrixSourceOp, thyraPrec_.get());
    }

}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::buildPreconditionerFaCSI( std::string type )
{
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef Teuchos::RCP<const Domain_Type> DomainConstPtr_Type;
    typedef std::vector<DomainConstPtr_Type> DomainConstPtr_vec_Type;

    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    // here we assume that FSI is always a time problem, this can be
    ParameterListPtr_Type parameterList;
    if (!timeProblem_.is_null())
        parameterList = timeProblem_->getParameterList();
    else
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Preconditioner can not be used without a time problem.");

    // Get FSI problem
    ProblemPtr_Type steadyProblem = timeProblem_->getUnderlyingProblem();
    Teuchos::RCP< FSI<SC,LO,GO,NO> > steadyFSI = Teuchos::rcp_dynamic_cast<FSI<SC,LO,GO,NO> >(steadyProblem);
    BlockMatrixPtr_Type fsiSystem = timeProblem_->getSystemCombined();

    ParameterListPtr_Type pLFluid = steadyFSI->getFluidProblem()->getParameterList();
    
 
    std::string precTypeFluid;
    if (type == "FaCSI")
        precTypeFluid = "Monolithic";
    else if (type == "FaCSI-Teko"){
        precTypeFluid = "Teko";
    }
    else if (type == "FaCSI-Block"){
        precTypeFluid = parameterList->sublist("Parameter Fluid").get("Preconditioner Type", "PCD");
    }

    CommConstPtr_Type comm = timeProblem_->getComm();
    bool useFluidPreconditioner = parameterList->sublist("General").get("Use Fluid Preconditioner", true);
    bool useSolidPreconditioner = parameterList->sublist("General").get("Use Solid Preconditioner", true);
    bool onlyDiagonal = parameterList->sublist("General").get("Only Diagonal", false);
    Teuchos::RCP< PrecOpFaCSI<SC,LO,GO,NO> > facsi
        = Teuchos::rcp(new PrecOpFaCSI<SC,LO,GO,NO> ( comm, precTypeFluid == "Monolithic", useFluidPreconditioner, useSolidPreconditioner, onlyDiagonal) );
    
    if (comm->getRank() == 0) {
        if (onlyDiagonal)
            std::cout << "\t### No preconditioner will be used! ###" << std::endl;
        else
            std::cout << "\t### FaCSI standard ###" << std::endl;
    }

    BlockMatrixPtr_Type fluidSystem = Teuchos::rcp( new BlockMatrix_Type(2) );
    
    // We build copies of the fluid system with homogenous Dirichlet boundary conditions on the interface
    // MatrixPtr_Type f = Teuchos::rcp(new Matrix_Type( fsiSystem->getBlock(0,0) ) );
   
    // MatrixPtr_Type bt = Teuchos::rcp(new Matrix_Type( fsiSystem->getBlock(0,1) ) );
    // MatrixPtr_Type b = Teuchos::rcp(new Matrix_Type( fsiSystem->getBlock(1,0) ) );
    // MatrixPtr_Type c;
    // if ( fsiSystem->blockExists(1,1) )
    //     c = Teuchos::rcp(new Matrix_Type( fsiSystem->getBlock(1,1) ) );
    // fluidSystem->addBlock( f, 0, 0 );
    // fluidSystem->addBlock( bt, 0, 1 );
    // fluidSystem->addBlock( b, 1, 0 );
    // if ( fsiSystem->blockExists(1,1) )
    //     fluidSystem->addBlock( c, 1, 1 );


   // We want to use the underlying Navier-Stokes Fluid Problem to build the preconditioner
    // We start with the fluid time problem
    Teuchos::RCP< TimeProblem<SC,LO,GO,NO> > fluidProblem = steadyFSI->problemTimeFluid_;
    fluidProblem->combineSystems(); // Build combined system || check if even is neccesary
    fluidProblem->setBoundariesSystem(); // Set boundaries || might also need fsi bc
    // The we cast the timeproblem to original Navier-Stokes problem and use it to build preconditioner
    Teuchos::RCP< NavierStokes<SC,LO,GO,NO> > fluidProblemSteady = Teuchos::rcp_dynamic_cast<NavierStokes<SC,LO,GO,NO> >(fluidProblem->getUnderlyingProblem());

    faCSIBCFactory_->setSystem( fluidProblem->getSystemCombined() );

    fluidProblemSteady->setupPreconditioner( precTypeFluid );
    precFluid_ = fluidProblemSteady->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();

    //Setup structure problem
    bool nonlinearStructure = false;
    if (steadyFSI->getStructureProblem().is_null())
        nonlinearStructure = true;

    ParameterListPtr_Type pLStructure;
    if (nonlinearStructure)
        pLStructure = steadyFSI->getNonLinStructureProblem()->getParameterList();
    else
        pLStructure = steadyFSI->getStructureProblem()->getParameterList();

    if (probSolid_.is_null()){
        probSolid_ = Teuchos::rcp( new MinPrecProblem_Type(  pLStructure, timeProblem_->getComm() ) );
        DomainConstPtr_vec_Type structDomain;
        if (nonlinearStructure)
            structDomain = steadyFSI->getNonLinStructureProblem()->getDomainVector();
        else
            structDomain = steadyFSI->getStructureProblem()->getDomainVector();
        
        probSolid_->initializeDomains( structDomain );
        probSolid_->initializeLinSolverBuilder( timeProblem_->getLinearSolverBuilder() );
    }
    
    BlockMatrixPtr_Type structSystem = Teuchos::rcp( new BlockMatrix_Type(1) );

    structSystem->addBlock( fsiSystem->getBlock(2,2), 0, 0 );

    probSolid_->initializeSystem( structSystem );
    
    probSolid_->setupPreconditioner( );

    precStruct_ = probSolid_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();


    //Setup geometry problem

    if (timeProblem_->getSystem()->size()>4) {
        ParameterListPtr_Type pLGeometry = steadyFSI->getGeometryProblem()->getParameterList();
        if (probGeo_.is_null()) {
            probGeo_ = Teuchos::rcp( new MinPrecProblem_Type( pLGeometry, timeProblem_->getComm() ) );
            DomainConstPtr_vec_Type geoDomain = steadyFSI->getGeometryProblem()->getDomainVector();
            probGeo_->initializeDomains( geoDomain );
            probGeo_->initializeLinSolverBuilder( timeProblem_->getLinearSolverBuilder() );
        }
        
        BlockMatrixPtr_Type geoSystem = Teuchos::rcp( new BlockMatrix_Type(1) );

        geoSystem->addBlock( fsiSystem->getBlock(4,4), 0, 0 );

        probGeo_->initializeSystem( geoSystem );

        probGeo_->setupPreconditioner( );

        precGeo_ = probGeo_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();
    }
    bool shape = false;
    if (fsiSystem->size()>4) {
        if (shape){
            facsi->setGIShape(   fsiSystem->getBlock(3,0)->getThyraLinOpNonConst()/*C1*/,
                                 fsiSystem->getBlock(0,3)->getThyraLinOpNonConst()/*C1T*/,
                                 fsiSystem->getBlock(3,2)->getThyraLinOpNonConst()/*C2*/,
                                 fsiSystem->getBlock(4,2)->getThyraLinOpNonConst()/*C4*/,
                                 precStruct_,
                                 precFluid_,
                                 fsiSystem->getBlock(0,0)->getThyraLinOpNonConst()/*fF*/,
                                 fsiSystem->getBlock(0,1)->getThyraLinOpNonConst()/*fBT*/,
                                 precGeo_,
                                 fsiSystem->getBlock(0,4)->getThyraLinOpNonConst()/*shape v*/,
                                 fsiSystem->getBlock(1,4)->getThyraLinOpNonConst()/*shape p*/ );
        }
        else {
            facsi->setGI(   fsiSystem->getBlock(3,0)->getThyraLinOpNonConst()/*C1*/,
                            fsiSystem->getBlock(0,3)->getThyraLinOpNonConst()/*C1T*/,
                            fsiSystem->getBlock(3,2)->getThyraLinOpNonConst()/*C2*/,
                            fsiSystem->getBlock(4,2)->getThyraLinOpNonConst()/*C4*/,
                            precStruct_,
                            precFluid_,
                            fsiSystem->getBlock(0,0)->getThyraLinOpNonConst()/*fF*/,
                            fsiSystem->getBlock(0,1)->getThyraLinOpNonConst()/*fBT*/,
                            precGeo_ );
        }
    }
    else {
        facsi->setGE(   fsiSystem->getBlock(3,0)->getThyraLinOpNonConst()/*C1*/,
                        fsiSystem->getBlock(0,3)->getThyraLinOpNonConst()/*C1T*/,
                        fsiSystem->getBlock(3,2)->getThyraLinOpNonConst()/*C2*/,
                        precStruct_,
                        precFluid_,
                        fsiSystem->getBlock(0,0)->getThyraLinOpNonConst()/*fF*/,
                        fsiSystem->getBlock(0,1)->getThyraLinOpNonConst()/*fBT*/ );
    }

    LinSolverBuilderPtr_Type solverBuilder = timeProblem_->getUnderlyingProblem()->getLinearSolverBuilder();

    if ( precFactory_.is_null() )
        precFactory_ = solverBuilder->createPreconditioningStrategy("");

    if ( thyraPrec_.is_null() )
        thyraPrec_ = precFactory_->createPrec();

    Teuchos::RCP< Thyra::DefaultPreconditioner<SC> > defaultPrec =
        Teuchos::rcp_dynamic_cast< Thyra::DefaultPreconditioner<SC> > (thyraPrec_);
    ThyraLinOpPtr_Type linOp =
        Teuchos::rcp_dynamic_cast< Thyra::LinearOpBase<SC> > (facsi);

    defaultPrec->initializeUnspecified( linOp );

    precondtionerIsBuilt_ = true;

}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setPressureMassMatrix(MatrixPtr_Type massMatrix) const{
    pressureMassMatrixPtr_ = massMatrix;
    pressureMass_= massMatrix->getThyraLinOp();
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setPressureLaplaceMatrix(MatrixPtr_Type matrix) const{
    pressureLaplace_ =matrix->getThyraLinOp();
    pressureLaplaceMatrixPtr_ = matrix; 
}

// template <class SC,class LO,class GO,class NO>
// void Preconditioner<SC,LO,GO,NO>::setPressureMass(MatrixPtr_Type matrix) const{
//     pressureMass_ = matrix->getThyraLinOp();
//     pressureMassMatrixPtr_ = matrix;
// }

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setPCDOperator(MatrixPtr_Type matrix) const{
    pcdOperator_ = matrix->getThyraLinOp();
    pcdOperatorMatrixPtr_ = matrix;
}

// Function to build a general 2 x 2 Block preconditioner
// Currently only used for (Navier-)Stokes type problems
// This includes the
// - Diagonal Prec, where the Schur complement is replaced by - 1/nu M_p
// - Triangular Prec, where the Schur complement is replaced by - 1/nu M_p
// - PCD Prec, where the Schur complement is replaced by -M_p F_p^-1 A_p
// - LSC Prec, where the Schur complement is replaced by  -A_p^-1 (B (M_v^-1) F (M_v^-1) B^T ) A_p^-1

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::buildPreconditionerBlock2x2( )
{
    PRECONDITIONER_START(buildPreconditionerBlock2x2, " buildPreconditionerBlock2x2");
   
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef Teuchos::RCP<const Domain_Type> DomainConstPtr_Type;
    typedef std::vector<DomainConstPtr_Type> DomainConstPtr_vec_Type;
    
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    
    ParameterListPtr_Type parameterList;
    BlockMatrixPtr_Type system;
    CommConstPtr_Type comm;
    ProblemPtr_Type steadyProblem;
    if (!timeProblem_.is_null()){
        std::cout << "buildPreconditionerBlock2x2 timeProblem_ " << std::endl;
        parameterList = timeProblem_->getParameterList();
        system = timeProblem_->getSystemCombined();
        comm = timeProblem_->getComm();
        steadyProblem = timeProblem_->getUnderlyingProblem();
    }
    else{
        parameterList = problem_->getParameterList();
        system = problem_->getSystem();
        comm = problem_->getComm();
        steadyProblem = problem_;
    }
    
    bool verbose( comm->getRank() == 0 );

    if(verbose){
        std::cout << " ######################## " << std::endl;
        std::cout << " Build 2x2 Preconditioner " << std::endl;
        std::cout << " ######################## " << std::endl;
    }
    ParameterListPtr_Type plVelocity( new Teuchos::ParameterList( parameterList->sublist("Velocity preconditioner") ) );
    ParameterListPtr_Type plSchur( new Teuchos::ParameterList( parameterList->sublist("Schur complement preconditioner") ) );
    
    //Addig General information to both lists
    plVelocity->sublist("General").setParameters(parameterList->sublist("General"));
    plSchur->sublist("General").setParameters(parameterList->sublist("General"));
    
    Teuchos::RCP< PrecBlock2x2<SC,LO,GO,NO> > blockPrec2x2
        = Teuchos::rcp(new PrecBlock2x2<SC,LO,GO,NO> ( comm ) );
    
    // The velocity problem is always treated the same
    if (probVelocity_.is_null()){
        probVelocity_ = Teuchos::rcp( new MinPrecProblem_Type( plVelocity, comm ) );
        DomainConstPtr_vec_Type domain1(0);
        domain1.push_back( steadyProblem->getDomain(0) );
        probVelocity_->initializeDomains( domain1 );
        probVelocity_->initializeLinSolverBuilder( steadyProblem->getLinearSolverBuilder() );
    }
    
    BlockMatrixPtr_Type system1 = Teuchos::rcp( new BlockMatrix_Type(1) );
    
    MatrixPtr_Type F = system->getBlock(0,0);//Teuchos::rcp(new Matrix_Type( system->getBlock(0,0) ) );

    system1->addBlock( F, 0, 0 );
    
    probVelocity_->initializeSystem( system1 );
    
    PRECONDITIONER_START(setupFInv, " Setup Preconditioner for F");
    probVelocity_->setupPreconditioner( "Monolithic" ); // single matrix
    PRECONDITIONER_STOP(setupFInv);

    precVelocity_ = probVelocity_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();
    
    if (probSchur_.is_null()) {
        if(!timeProblem_.is_null()){
            probSchur_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
            DomainConstPtr_vec_Type domain2(0);
            domain2.push_back( timeProblem_->getDomain(1) );
            probSchur_->initializeDomains( domain2 );
            probSchur_->initializeLinSolverBuilder( timeProblem_->getLinearSolverBuilder() );
        }
        else
        {
            probSchur_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
            DomainConstPtr_vec_Type domain2(0);
            domain2.push_back( problem_->getDomain(1) );
            probSchur_->initializeDomains( domain2 );
            probSchur_->initializeLinSolverBuilder( problem_->getLinearSolverBuilder() );
            
        }
    }
   
    std::string type = parameterList->sublist("General").get("Preconditioner Method","Diagonal");

    // We distinguish for the Schur complement component
    // Setup additional things
    PRECONDITIONER_START(setupSInv, " Setup Preconditioner for S");

    if(type == "Diagonal" || type == "Triangular"){
        BlockMatrixPtr_Type Mp = Teuchos::rcp( new BlockMatrix_Type(1) );
            
        Mp->addBlock( pressureMassMatrixPtr_, 0, 0 );
        
        probSchur_->initializeSystem( Mp );
        
        probSchur_->setupPreconditioner( "Monolithic" ); // single matrix
        
        precSchur_ = probSchur_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();
    }
    else if (type == "PCD") {
        // For PCD we additionally need to setup the monolithic preconditioner for the 
        // Laplace operator and the pressure mass matrix
        // We include a diagonal inverse approximation of the mass matrix
        if (verbose) {
            std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
            std::cout << "\t --- Building PCD Operator Components " << std::endl;
            std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
        }

        if (probLaplace_.is_null()) {
            if (!timeProblem_.is_null()){
                probLaplace_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
                DomainConstPtr_vec_Type domain2(0);
                domain2.push_back( timeProblem_->getDomain(1) );
                probLaplace_->initializeDomains( domain2 );
                probLaplace_->initializeLinSolverBuilder( timeProblem_->getLinearSolverBuilder() );
            }
            else{
                probLaplace_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
                DomainConstPtr_vec_Type domain2(0);
                domain2.push_back( problem_->getDomain(1) );
                probLaplace_->initializeDomains( domain2 );
                probLaplace_->initializeLinSolverBuilder( problem_->getLinearSolverBuilder() );
            }
        }

        if(laplaceInverse_.is_null()){// The Schwarz approximation of Laplace operator only needs to be build once
            if (verbose) {
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
                std::cout << "\t --- PCD: Setup A_p Schwarz Approximation " << std::endl;
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
            }
            BlockMatrixPtr_Type Ap = Teuchos::rcp( new BlockMatrix_Type(1) );
            Ap->addBlock(pressureLaplaceMatrixPtr_,0,0);

            probLaplace_->initializeSystem( Ap );
         
            probLaplace_->setupPreconditioner( "Monolithic" ); // single matrix
            laplaceInverse_ = probLaplace_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();
        }
        
        if (probMass_.is_null()) {
            if (!timeProblem_.is_null()){
                probMass_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
                DomainConstPtr_vec_Type domain2(0);
                domain2.push_back( timeProblem_->getDomain(1) );
                probMass_->initializeDomains( domain2 );
                probMass_->initializeLinSolverBuilder( timeProblem_->getLinearSolverBuilder() );
            }
            else{
                probMass_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
                DomainConstPtr_vec_Type domain2(0);
                domain2.push_back( problem_->getDomain(1) );
                probMass_->initializeDomains( domain2 );
                probMass_->initializeLinSolverBuilder( problem_->getLinearSolverBuilder() );
            }
        }

        if(massMatrixInverse_.is_null()){// The Schwarz approximation of pressure mass matrix only needs to be build once
            if (verbose) {
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
                std::cout << "\t --- PCD: Setup M_p Schwarz Approximation " << std::endl;
                std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
            }
            BlockMatrixPtr_Type Qp = Teuchos::rcp( new BlockMatrix_Type(1) );       
            Qp->addBlock(pressureMassMatrixPtr_,0,0);

            // Approximation of Mp is either done by Monolithic preconditioner or by a diagonal
            // Approximation with 'Diagonal' or TODO: AbsRowSum
            bool explicitInverse = parameterList->sublist("General").get("Mu Explicit Inverse",true);
            std::string typeDiag = parameterList->sublist("General").get("Diagonal Approximation","Diagonal");

            if(explicitInverse)
            {
                probMass_->initializeSystem( Qp );
                probMass_->setupPreconditioner( "Monolithic" ); // single matrix
                massMatrixInverse_ = probMass_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();
            }
            else
            {
                massMatrixInverse_ = pressureMassMatrixPtr_->buildDiagonalInverse(typeDiag)->getThyraLinOpNonConst() ;
            }
        }

    }
    else if (type == "LSC") {

        if (verbose) {
            std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
            std::cout << "\t --- Building LSC Operator Components " << std::endl;
            std::cout << "\t --- -------------------------------------------------------- ---"<< std::endl;
        }

        if (probLaplace_.is_null()) {
            if (!timeProblem_.is_null()){
                probLaplace_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
                DomainConstPtr_vec_Type domain2(0);
                domain2.push_back( timeProblem_->getDomain(1) );
                probLaplace_->initializeDomains( domain2 );
                probLaplace_->initializeLinSolverBuilder( timeProblem_->getLinearSolverBuilder() );
            }
            else{
                probLaplace_ = Teuchos::rcp( new MinPrecProblem_Type( plSchur, comm ) );
                DomainConstPtr_vec_Type domain2(0);
                domain2.push_back( problem_->getDomain(1) );
                probLaplace_->initializeDomains( domain2 );
                probLaplace_->initializeLinSolverBuilder( problem_->getLinearSolverBuilder() );
            }
        }

        BlockMatrixPtr_Type Ap = Teuchos::rcp( new BlockMatrix_Type(1) );
        Ap->addBlock(pressureLaplaceMatrixPtr_,0,0);

        probLaplace_->initializeSystem( Ap );

        probLaplace_->setupPreconditioner( "Monolithic" ); // single matrix
        laplaceInverse_ = probLaplace_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();

        if (probVMass_.is_null()) {
            if (!timeProblem_.is_null()){
                probVMass_ = Teuchos::rcp( new MinPrecProblem_Type( plVelocity, comm ) );
                DomainConstPtr_vec_Type domain1(0);
                domain1.push_back( timeProblem_->getDomain(0) );
                probVMass_->initializeDomains( domain1 );
                probVMass_->initializeLinSolverBuilder( timeProblem_->getLinearSolverBuilder() );
            }
            else{
                probVMass_ = Teuchos::rcp( new MinPrecProblem_Type( plVelocity, comm ) );
                DomainConstPtr_vec_Type domain1(0);
                domain1.push_back( problem_->getDomain(0) );
                probVMass_->initializeDomains( domain1 );
                probVMass_->initializeLinSolverBuilder( problem_->getLinearSolverBuilder() );
            }
        }
        BlockMatrixPtr_Type Qv = Teuchos::rcp( new BlockMatrix_Type(1) );       
        Qv->addBlock(velocityMassMatrixMatrixPtr_,0,0);

        // Approximation of Mp is either done by Monolithic preconditioner or by a diagonal
        // Approximation with 'Diagonal' or TODO: AbsRowSum
        bool explicitInverse = parameterList->sublist("General").get("Mu Explicit Inverse",true);
        std::string typeDiag = parameterList->sublist("General").get("Diagonal Approximation","Diagonal");

        if(explicitInverse)
        {
            probVMass_->initializeSystem( Qv );
            probVMass_->setupPreconditioner( "Monolithic" ); // single matrix
            massMatrixVInverse_ = probVMass_->getPreconditioner()->getThyraPrec()->getNonconstUnspecifiedPrecOp();
        }
        else
        {
           massMatrixVInverse_ = velocityMassMatrixMatrixPtr_->buildDiagonalInverse(typeDiag)->getThyraLinOpNonConst() ;
        }
            
    }
    PRECONDITIONER_STOP(setupSInv);

    // Building block Prec and passing along the different operators
    // that are required to build the respective preconditioners
    if (type == "Diagonal") {
        blockPrec2x2->setDiagonal(precVelocity_,
                                  precSchur_);
    }
    else if (type == "Triangular") {
        ThyraLinOpPtr_Type BT = system->getBlock(0,1)->getThyraLinOpNonConst();
        blockPrec2x2->setTriangular(precVelocity_,
                                    precSchur_,
                                    BT);
    }
    else if (type == "PCD") {
        if (!timeProblem_.is_null()){
            timeProblem_->assemble("UpdateConvectionDiffusionOperator");
        }
        else{
            problem_->assemble("UpdateConvectionDiffusionOperator");
        }
        MatrixPtr_Type pcdOperatorScaled = Teuchos::rcp( new Matrix_Type( pcdOperatorMatrixPtr_ ) );
        pcdOperatorScaled->resumeFill();
        pcdOperatorScaled->scale(-1.0);
        pcdOperatorScaled->fillComplete();
        ThyraLinOpPtr_Type BT = system->getBlock(0,1)->getThyraLinOpNonConst();
        blockPrec2x2->setTriangular(precVelocity_,
                                    laplaceInverse_,
                                    pcdOperatorScaled->getThyraLinOpNonConst(),
                                    massMatrixInverse_,
                                    massMatrixVInverse_,
                                    BT);
        ThyraLinOpPtr_Type B = system->getBlock(1,0)->getThyraLinOpNonConst();
        blockPrec2x2->setB(B);        
    }
    else if (type == "LSC") {
        ThyraLinOpPtr_Type BT = system->getBlock(0,1)->getThyraLinOpNonConst();
        blockPrec2x2->setTriangular(precVelocity_,
                                    laplaceInverse_,
                                    massMatrixVInverse_,
                                    BT);

        
        ThyraLinOpPtr_Type B = system->getBlock(1,0)->getThyraLinOpNonConst();
        blockPrec2x2->setB(B); 
        
        ThyraLinOpPtr_Type F = system->getBlock(0,0)->getThyraLinOpNonConst();
        blockPrec2x2->setF(F);    
            
       
    }
    
    LinSolverBuilderPtr_Type solverBuilder;
    if (!timeProblem_.is_null())
        solverBuilder = timeProblem_->getLinearSolverBuilder();
    else
        solverBuilder = problem_->getLinearSolverBuilder();

    if ( precFactory_.is_null() )
        precFactory_ = solverBuilder->createPreconditioningStrategy("");
    
    if ( thyraPrec_.is_null() )
        thyraPrec_ = precFactory_->createPrec();
    
    Teuchos::RCP< Thyra::DefaultPreconditioner<SC> > defaultPrec =
    Teuchos::rcp_dynamic_cast< Thyra::DefaultPreconditioner<SC> > (thyraPrec_);
    ThyraLinOpPtr_Type linOp =
        Teuchos::rcp_dynamic_cast< Thyra::LinearOpBase<SC> > (blockPrec2x2);
    
    defaultPrec->initializeUnspecified( linOp );
    
    precondtionerIsBuilt_ = true;

    PRECONDITIONER_STOP(buildPreconditionerBlock2x2);

    
}

    

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setVelocityParameters( ParameterListPtr_Type parameterList, int coarseRanks )
{
    CommConstPtr_Type comm;
    if (!problem_.is_null())
        comm = problem_->getComm();
    else if(!timeProblem_.is_null())
        comm = timeProblem_->getComm();

    bool verbose( comm->getRank() == 0 );
    
    // Xpetra for now
    Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > repeatedMaps(1);
    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const XpetraMapConstPtr_Type XpetraMapConstPtrConst_Type;

    Teuchos::ArrayRCP<FROSch::DofOrdering> dofOrderings(1);
    Teuchos::ArrayRCP<UN> dofsPerNodeVector(1);
    ParameterListPtr_Type velocitySubList = sublist( sublist( sublist( sublist( parameterList, "Preconditioner Types" ) , "Teko" ) , "Inverse Factory Library" ) , "FROSch-Velocity" );
    dofsPerNodeVector[0] = (UN) velocitySubList->get( "DofsPerNode", 2);
    TEUCHOS_TEST_FOR_EXCEPTION(dofsPerNodeVector[0]<2, std::logic_error, "DofsPerNode for velocity must be atleast 2.");

    /*Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp;
    if (!problem_.is_null())
        mapConstTmp = problem_->getDomain(0)->getMapVecFieldRepeated()->getTpetraMap();
    else if(!timeProblem_.is_null())
        mapConstTmp = timeProblem_->getDomain(0)->getMapVecFieldRepeated()->getTpetraMap();

    Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);*/

    MapConstPtr_Type mapConstTmp;
    if (!problem_.is_null())
        mapConstTmp = problem_->getDomain(0)->getMapVecFieldRepeated();
    else if(!timeProblem_.is_null())
        mapConstTmp = timeProblem_->getDomain(0)->getMapVecFieldRepeated();
    
    XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);

    repeatedMaps[0] = mapX;

    if (!velocitySubList->get( "DofOrdering", "NodeWise" ).compare("DimensionWise"))
        dofOrderings[0] = FROSch::DimensionWise;
    else if (!velocitySubList->get( "DofOrdering", "NodeWise" ).compare("NodeWise"))
        dofOrderings[0] = FROSch::NodeWise;
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Please chose a valid  DOF ordering for ThyraPreconditioner(FROSch).");

    velocitySubList->set("Repeated Map Vector",repeatedMaps);
    velocitySubList->set("DofOrdering Vector",dofOrderings);
    velocitySubList->set("DofsPerNode Vector",dofsPerNodeVector);
    velocitySubList->set( "Mpi Ranks Coarse", coarseRanks );
    
    int lowerBound = 100000000;
    int upperBound = -1;

    tuple_intint_Type rankRange;
    if (!problem_.is_null())
        rankRange = problem_->getDomain(0)->getMesh()->getRankRange();
    else if(!timeProblem_.is_null())
        rankRange = timeProblem_->getDomain(0)->getMesh()->getRankRange();
    
    if (std::get<0>(rankRange) < lowerBound)
        lowerBound = std::get<0>(rankRange);
    if (std::get<1>(rankRange) > upperBound)
        upperBound = std::get<1>(rankRange);
    
    int lowerBoundCoarse = lowerBound;
    int upperBoundCoarse = upperBound;
    // For now only for structured decompositions
    if ( coarseRanks > 0 ){
        lowerBoundCoarse = upperBound + 1;
        upperBoundCoarse = comm->getSize() - 1;
    }
    if (verbose) {
        std::cout << "\t --- ------------------------------------------------------------------- ---"<< std::endl;
        std::cout << "\t --- Range for local problems of preconditioner Teko velocity from " << lowerBound << " to " << upperBound << std::endl;
        std::cout << "\t --- Range for coarse problem of preconditioner Teko velocity from " << lowerBoundCoarse << " to " << upperBoundCoarse << std::endl;
        std::cout << "\t --- ------------------------------------------------------------------- ---"<< std::endl;
    }
    
    velocitySubList->set( "Local problem ranks lower bound", lowerBound );
    velocitySubList->set( "Local problem ranks upper bound", upperBound );
    
    velocitySubList->set( "Coarse problem ranks lower bound", lowerBoundCoarse );
    velocitySubList->set( "Coarse problem ranks upper bound", upperBoundCoarse );

    

}


template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setPressureParameters( ParameterListPtr_Type parameterList, int coarseRanks )
{
    CommConstPtr_Type comm;
    if (!problem_.is_null())
        comm = problem_->getComm();
    else if(!timeProblem_.is_null())
        comm = timeProblem_->getComm();

    bool verbose( comm->getRank() == 0 );
    
    Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > repeatedMaps(1);
    typedef Xpetra::Map<LO,GO,NO> XpetraMap_Type;
    typedef Teuchos::RCP<XpetraMap_Type> XpetraMapPtr_Type;
    typedef Teuchos::RCP<const XpetraMap_Type> XpetraMapConstPtr_Type;
    typedef const XpetraMapConstPtr_Type XpetraMapConstPtrConst_Type;

    Teuchos::ArrayRCP<FROSch::DofOrdering> dofOrderings(1);
    Teuchos::ArrayRCP<UN> dofsPerNodeVector(1);
    ParameterListPtr_Type pressureSubList = sublist( sublist( sublist( sublist( parameterList, "Preconditioner Types" ) , "Teko" ) , "Inverse Factory Library" ) , "FROSch-Pressure" );
    dofsPerNodeVector[0] = (UN) pressureSubList->get( "DofsPerNode", 1);
    TEUCHOS_TEST_FOR_EXCEPTION(dofsPerNodeVector[0]!=1, std::logic_error, "DofsPerNode for pressure must be  1.");

    /*Teuchos::RCP<const Tpetra::Map<LO,GO,NO> > mapConstTmp;

    if (!problem_.is_null()){
        if ( problem_->getDomain(1)->getFEType()=="P0" )
            mapConstTmp = problem_->getDomain(1)->getElementMap()->getTpetraMap();
        else
            mapConstTmp = problem_->getDomain(1)->getMapRepeated()->getTpetraMap();
    }
    else if(!timeProblem_.is_null()){
        if ( timeProblem_->getDomain(1)->getFEType()=="P0" )
            mapConstTmp = timeProblem_->getDomain(1)->getElementMap()->getTpetraMap();
        else
            mapConstTmp = timeProblem_->getDomain(1)->getMapRepeated()->getTpetraMap();
    }*/
    MapConstPtr_Type mapConstTmp;
    if (!problem_.is_null()){
        if ( problem_->getDomain(1)->getFEType()=="P0" )
            mapConstTmp = problem_->getDomain(1)->getElementMap();
        else
            mapConstTmp = problem_->getDomain(1)->getMapRepeated();
    }
    else if(!timeProblem_.is_null()){
        if ( timeProblem_->getDomain(1)->getFEType()=="P0" )
            mapConstTmp = timeProblem_->getDomain(1)->getElementMap();
        else
            mapConstTmp = timeProblem_->getDomain(1)->getMapRepeated();
    }

    //Teuchos::RCP<Tpetra::Map<LO,GO,NO> > mapTmp = Teuchos::rcp_const_cast<Tpetra::Map<LO,GO,NO> > (mapConstTmp);
    
    XpetraMapConstPtr_Type mapConstX = Xpetra::MapFactory<LO,GO,NO>::Build( Xpetra::UseTpetra, mapConstTmp->getGlobalNumElements(), mapConstTmp->getNodeElementList(), mapConstTmp->getIndexBase(), mapConstTmp->getComm() );
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapX= Teuchos::rcp_const_cast<Xpetra::Map<LO,GO,NO> > (mapConstX);
                           
    repeatedMaps[0] = mapX;

    if (!pressureSubList->get( "DofOrdering", "NodeWise" ).compare("DimensionWise"))
        dofOrderings[0] = FROSch::DimensionWise;
    else if (!pressureSubList->get( "DofOrdering", "NodeWise" ).compare("NodeWise"))
        dofOrderings[0] = FROSch::NodeWise;
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Please chose a valid  DOF ordering for ThyraPreconditioner(FROSch).");

    pressureSubList->set("Repeated Map Vector",repeatedMaps);
    pressureSubList->set("DofOrdering Vector",dofOrderings);
    pressureSubList->set("DofsPerNode Vector",dofsPerNodeVector);
    pressureSubList->set( "Mpi Ranks Coarse", coarseRanks );
    
    int lowerBound = 100000000;
    int upperBound = -1;
    

    tuple_intint_Type rankRange;
    if (!problem_.is_null())
        rankRange = problem_->getDomain(1)->getMesh()->getRankRange();
    else if(!timeProblem_.is_null())
        rankRange = timeProblem_->getDomain(1)->getMesh()->getRankRange();
    
    if (std::get<0>(rankRange) < lowerBound)
        lowerBound = std::get<0>(rankRange);
    if (std::get<1>(rankRange) > upperBound)
        upperBound = std::get<1>(rankRange);
    
    int lowerBoundCoarse = lowerBound;
    int upperBoundCoarse = upperBound;
    // For now only for structured decompositions
    if ( coarseRanks > 0 ){
        lowerBoundCoarse = upperBound + 1;
        upperBoundCoarse = comm->getSize() - 1;
    }
    if (verbose) {
        std::cout << "\t --- ------------------------------------------------------------------- ---"<< std::endl;
        std::cout << "\t --- Range for local problems of preconditioner Teko pressure from " << lowerBound << " to " << upperBound << std::endl;
        std::cout << "\t --- Range for coarse problem of preconditioner Teko pressure from " << lowerBoundCoarse << " to " << upperBoundCoarse << std::endl;
        std::cout << "\t --- ------------------------------------------------------------------- ---"<< std::endl;
    }
    
    pressureSubList->set( "Local problem ranks lower bound", lowerBound );
    pressureSubList->set( "Local problem ranks upper bound", upperBound );
    
    pressureSubList->set( "Coarse problem ranks lower bound", lowerBoundCoarse );
    pressureSubList->set( "Coarse problem ranks upper bound", upperBoundCoarse );

    
}

template <class SC,class LO,class GO,class NO>
typename Preconditioner<SC,LO,GO,NO>::ThyraLinOpConstPtr_Type Preconditioner<SC,LO,GO,NO>::getTekoOp(){
    return tekoLinOp_;
}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::setVelocityMassMatrix(MatrixPtr_Type massMatrix) const{
    velocityMassMatrix_ = massMatrix->getThyraLinOp();
    velocityMassMatrixMatrixPtr_ = massMatrix;
}
#endif

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::exportCoarseBasis( ){
    
    CommConstPtr_Type comm;
    if (!problem_.is_null())
        comm = problem_->getComm();
    else if(!timeProblem_.is_null())
        comm = timeProblem_->getComm();
    
    TEUCHOS_TEST_FOR_EXCEPTION( pListPhiExport_.is_null(), std::runtime_error, "No parameterlist to extract Phi pointer.");
    ParameterListPtr_Type pLPrec = sublist( sublist( pListPhiExport_, "Preconditioner Types" ) ,"FROSch" );
    std::string coarseType = pLPrec->get( "CoarseOperator Type", "RGDSWCoarseOperator" );
    ParameterListPtr_Type pLCoarse = sublist( pLPrec, coarseType );
//    TEUCHOS_TEST_FOR_EXCEPTION( coarseType!="RGDSWCoarseOperator" && coarseType!="GDSWCoarseOperator", std::runtime_error, "Export Phi only for GDSWCoarseOperator and RGDSWCoarseOperator.");

    TEUCHOS_TEST_FOR_EXCEPTION( !pLCoarse->isParameter("RCP(Phi)"), std::runtime_error, "No parameter to extract Phi pointer.");
    
    Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > phiTpetra;
    
    TEUCHOS_TEST_FOR_EXCEPTION( !pLCoarse->isType<decltype(phiTpetra)>("RCP(Phi)"), std::runtime_error, "Wrong type of pointer to extract Phi.");
    
    phiTpetra = pLCoarse->get<decltype(phiTpetra)>("RCP(Phi)");
    
    MatrixPtr_Type phiMatrix = Teuchos::rcp( new Matrix_Type( phiTpetra ) );
    int numberOfBlocks;
    {
        ParameterListPtr_Type parameterList;
        if (!problem_.is_null())
            parameterList = problem_->getParameterList();
        else if(!timeProblem_.is_null())
            parameterList = timeProblem_->getParameterList();
        numberOfBlocks = parameterList->get("Number of blocks",1);        
    }

    std::vector<MapConstPtr_Type> blockMaps(numberOfBlocks);
    MultiVectorPtr_Type phiMV;
    phiMatrix->toMV( phiMV );
    
    for (UN i = 0; i < numberOfBlocks; i++) {
        int dofsPerNode =  pLPrec->get( "DofsPerNode" + std::to_string(i+1), 1);
        MapConstPtr_Type map;
        if (dofsPerNode>1) {
            if (!problem_.is_null()) {
                TEUCHOS_TEST_FOR_EXCEPTION( problem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                map = problem_->getDomain(i)->getMapVecFieldUnique();
            }
            else if(!timeProblem_.is_null()){
                TEUCHOS_TEST_FOR_EXCEPTION( timeProblem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                map = timeProblem_->getDomain(i)->getMapVecFieldUnique();
            }
        }
        else{
            if (!problem_.is_null()) {
                TEUCHOS_TEST_FOR_EXCEPTION( problem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                map = problem_->getDomain(i)->getMapUnique();
            }
            else if(!timeProblem_.is_null()){
                TEUCHOS_TEST_FOR_EXCEPTION( timeProblem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                map = timeProblem_->getDomain(i)->getMapUnique();
            }
        }
        blockMaps[i] = map;
    }
    
    
    BlockMultiVectorPtr_Type phiBlockMV = Teuchos::rcp( new BlockMultiVector_Type ( blockMaps, phiMV->getNumVectors() ) );
    phiBlockMV->setMergedVector( phiMV );
    phiBlockMV->split();
    
    ParameterListPtr_Type pLProblem;
    if (!problem_.is_null())
        pLProblem = problem_->getParameterList();
    else if(!timeProblem_.is_null())
        pLProblem = timeProblem_->getParameterList();
    
    BlockMultiVectorPtr_Type phiSumBlockMV = phiBlockMV->sumColumns();
    
    for (int i=0; i<phiBlockMV->size(); i++) {
        bool exportThisBlock = !pLProblem->sublist("Exporter").get( "Exclude coarse functions block" + std::to_string(i+1), false);
        
        if (exportThisBlock){
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exporter(new ExporterParaView<SC,LO,GO,NO>());

            DomainConstPtr_Type domain;
            if (!problem_.is_null())
                domain = problem_->getDomain(i);
            else if(!timeProblem_.is_null())
                domain = timeProblem_->getDomain(i);

            int dofsPerNode = domain->getDimension();
            std::string fileName = pLProblem->sublist("Exporter").get( "Name coarse functions block" + std::to_string(i+1), "phi");

            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( domain->getMesh() );
            exporter->setup(fileName, meshNonConst, domain->getFEType());
            
            for (int j=0; j<phiBlockMV->getNumVectors(); j++) {
                
                MultiVectorConstPtr_Type exportPhi = phiBlockMV->getBlock( i )->getVector( j );
                
                std::string varName = fileName + std::to_string(j);
                if ( pLPrec->get( "DofsPerNode" + std::to_string(i+1), 1) > 1 )
                    exporter->addVariable( exportPhi, varName, "Vector", dofsPerNode, domain->getMapUnique() );
                else
                    exporter->addVariable( exportPhi, varName, "Scalar", 1, domain->getMapUnique() );
                
            }
            MultiVectorConstPtr_Type exportSumPhi = phiSumBlockMV->getBlock(i);
            std::string varName = fileName + "Sum";
            if ( pLPrec->get( "DofsPerNode" + std::to_string(i+1), 1) > 1 )
                exporter->addVariable( exportSumPhi, varName, "Vector", dofsPerNode, domain->getMapUnique() );
            else
                exporter->addVariable( exportSumPhi, varName, "Scalar", 1, domain->getMapUnique() );
            
            exporter->save(0.0);
        }
    }

}

template <class SC,class LO,class GO,class NO>
void Preconditioner<SC,LO,GO,NO>::exportCoarseBasisFSI( ){
    
    CommConstPtr_Type comm;
    if (!problem_.is_null())
        comm = problem_->getComm();
    else if(!timeProblem_.is_null())
        comm = timeProblem_->getComm();
    
    TEUCHOS_TEST_FOR_EXCEPTION( pListPhiExport_.is_null(), std::runtime_error, "No parameterlist to extract Phi pointer.");
    ParameterListPtr_Type pLPrec = sublist( sublist( pListPhiExport_, "Preconditioner Types" ) ,"FROSch" );
    std::string coarseType = pLPrec->get( "CoarseOperator Type", "RGDSWCoarseOperator" );
    ParameterListPtr_Type pLCoarse = sublist( pLPrec, coarseType );
//    TEUCHOS_TEST_FOR_EXCEPTION( coarseType!="RGDSWCoarseOperator" && coarseType!="GDSWCoarseOperator", std::runtime_error, "Export Phi only for GDSWCoarseOperator and RGDSWCoarseOperator.");

    TEUCHOS_TEST_FOR_EXCEPTION( !pLCoarse->isParameter("Phi Pointer"), std::runtime_error, "No parameter to extract Phi pointer.");
    
    Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > phiTpetra;
    
    TEUCHOS_TEST_FOR_EXCEPTION( !pLCoarse->isType<decltype(phiTpetra)>("Phi Pointer"), std::runtime_error, "Wrong type of pointer to extract Phi.");
    
    phiTpetra = pLCoarse->get<decltype(phiTpetra)>("Phi Pointer");
    
    MatrixPtr_Type phiMatrix = Teuchos::rcp( new Matrix_Type( phiTpetra ) );
    int numberOfBlocks;
    {
        ParameterListPtr_Type parameterList;
        if (!problem_.is_null())
            parameterList = problem_->getParameterList();
        else if(!timeProblem_.is_null())
            parameterList = timeProblem_->getParameterList();
        numberOfBlocks = parameterList->get("Number of blocks",1);
    }

    std::vector<MapConstPtr_Type> blockMaps(numberOfBlocks);
    MultiVectorPtr_Type phiMV;
    phiMatrix->toMV( phiMV );
    
    for (UN i = 0; i < numberOfBlocks; i++) {
        int dofsPerNode =  pLPrec->get( "DofsPerNode" + std::to_string(i+1), 1);
        MapConstPtr_Type map;
        if (i==3) {
            map = timeProblem_->getDomain(i)->getInterfaceMapVecFieldUnique();
        }
        else {
            if (dofsPerNode>1) {
                if (!problem_.is_null()) {
                    TEUCHOS_TEST_FOR_EXCEPTION( problem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                    map = problem_->getDomain(i)->getMapVecFieldUnique();
                }
                else if(!timeProblem_.is_null()){
                    TEUCHOS_TEST_FOR_EXCEPTION( timeProblem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                    map = timeProblem_->getDomain(i)->getMapVecFieldUnique();
                }
            }
            else{
                if (!problem_.is_null()) {
                    TEUCHOS_TEST_FOR_EXCEPTION( problem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                    map = problem_->getDomain(i)->getMapUnique();
                }
                else if(!timeProblem_.is_null()){
                    TEUCHOS_TEST_FOR_EXCEPTION( timeProblem_->getDomain(i)->getFEType() == "P0", std::logic_error, "Vector field map not implemented for P0 elements.");
                    map = timeProblem_->getDomain(i)->getMapUnique();
                }
            }
        }
        blockMaps[i] = map;
    }
    
    
    BlockMultiVectorPtr_Type phiBlockMV = Teuchos::rcp( new BlockMultiVector_Type ( blockMaps, phiMV->getNumVectors() ) );
    phiBlockMV->setMergedVector( phiMV );
    phiBlockMV->split();
    
    ParameterListPtr_Type pLProblem;
    if (!problem_.is_null())
        pLProblem = problem_->getParameterList();
    else if(!timeProblem_.is_null())
        pLProblem = timeProblem_->getParameterList();
    
    BlockMultiVectorPtr_Type phiSumBlockMV = phiBlockMV->sumColumns();
    
    for (int i=0; i<phiBlockMV->size(); i++) {
        bool exportThisBlock = !pLProblem->sublist("Exporter").get( "Exclude coarse functions block" + std::to_string(i+1), false);
        
        if (exportThisBlock){
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exporter(new ExporterParaView<SC,LO,GO,NO>());

            DomainConstPtr_Type domain;
            if (!problem_.is_null())
                domain = problem_->getDomain(i);
            else if(!timeProblem_.is_null())
                domain = timeProblem_->getDomain(i);

            int dofsPerNode = domain->getDimension();
            std::string fileName = pLProblem->sublist("Exporter").get( "Name coarse functions block" + std::to_string(i+1), "phi");

            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( domain->getMesh() );
            exporter->setup(fileName, meshNonConst, domain->getFEType());
            
            for (int j=0; j<phiBlockMV->getNumVectors(); j++) {
                
                MultiVectorConstPtr_Type exportPhi = phiBlockMV->getBlock( i )->getVector( j );
                
                std::string varName = fileName + std::to_string(j);
                if ( pLPrec->get( "DofsPerNode" + std::to_string(i+1), 1) > 1 )
                    exporter->addVariable( exportPhi, varName, "Vector", dofsPerNode, domain->getMapUnique() );
                else
                    exporter->addVariable( exportPhi, varName, "Scalar", 1, domain->getMapUnique() );
                
            }
            MultiVectorConstPtr_Type exportSumPhi = phiSumBlockMV->getBlock(i);
            std::string varName = fileName + "Sum";
            if ( pLPrec->get( "DofsPerNode" + std::to_string(i+1), 1) > 1 )
                exporter->addVariable( exportSumPhi, varName, "Vector", dofsPerNode, domain->getMapUnique() );
            else
                exporter->addVariable( exportSumPhi, varName, "Scalar", 1, domain->getMapUnique() );
            
            exporter->save(0.0);
        }
    }

}


}

#endif
