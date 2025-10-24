#ifndef TIMEPROBLEM_DEF_hpp
#define TIMEPROBLEM_DEF_hpp
#include "TimeProblem_decl.hpp"
/*!
 Definition of TimeProblem

 @brief  TimeProblem
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {

template<class SC,class LO,class GO,class NO>
TimeProblem<SC,LO,GO,NO>::TimeProblem(Problem_Type& problem, CommConstPtr_Type comm):
comm_(comm),
systemCombined_(),
systemMass_(),
timeParameters_(),
timeStepDef_(),
massParameters_(),
feFactory_(),
dimension_(problem.getDomain(0)->getDimension()),
verbose_(comm->getRank()==0),
parameterList_(problem.getParameterList()),
bcFactory_(problem.getBCFactory()),
massBCSet_(false),
solutionPreviousTimesteps_(0),
velocityPreviousTimesteps_(0),
accelerationPreviousTimesteps_(0),
systemMassPreviousTimeSteps_(0),
time_(0.)
#ifdef TIMEPROBLEM_TIMER
,TimeSolveTimer_ (Teuchos::TimeMonitor::getNewCounter("TT Problem: Solve"))
#endif
,precInitOnly_(false)
{
    problem_ = Teuchos::rcpFromRef(problem);
    BCConstPtr_Type bcFaCSI;
    if( problem_->preconditioner_->hasFaCSIBCFactory() )
        bcFaCSI = problem_->preconditioner_->getFaCSIBCFactory();
    
    problem_->preconditioner_ = Teuchos::rcp( new Preconditioner_Type( this ) ); //we reset the preconditioner of underlying problem!
    
    if( !bcFaCSI.is_null() )
        problem_->preconditioner_->setFaCSIBCFactory( bcFaCSI );
    FEFacConstPtr_Type  facConst = problem.getFEFactory();
    feFactory_ = Teuchos::rcp_const_cast<FEFac_Type>(facConst);
    
    systemMass_.reset(new BlockMatrix_Type(1));
    systemCombined_.reset( new BlockMatrix_Type( 1 ) );

}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::setTimeDef( SmallMatrix<int>& def ){
    timeStepDef_ = def;
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::assemble( std::string type ) const{
    // If timestepping class is external, it is assumed that the full timedependent problem matrix and rhs are assembled during the assemble call(s)
    std::string timestepping = parameterList_->sublist("Timestepping Parameter").get("Class","Singlestep");

    if (type == "MassSystem"){
        // is not used in FSI
        // systemMass_ wird gebaut (Massematrix), welche schon mit der Dichte \rho skaliert wurde
        assembleMassSystem();
        initializeCombinedSystems();
        NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
        if (!nonLinProb.is_null()){// we combine the nonlinear system with the mass matrix in the NonLinearSolver after the reassembly of each linear system
        }
        else{
            if (timestepping=="External" )
                this->systemCombined_ = problem_->getSystem();
            else
                this->combineSystems();
        }
    }
    else{
        //we need to tell the problem about the last solution if we use extrapolation!
        problem_->assemble(type);

        if (timestepping=="External") 
            this->systemCombined_ = problem_->getSystem();
        else
            this->combineSystems();
    }
}


template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::reAssembleAndFill( BlockMatrixPtr_Type bMat, std::string type ){

    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    nonLinProb->reAssembleAndFill( bMat, type );
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::combineSystems() const{
    BlockMatrixPtr_Type tmpSystem = problem_->getSystem();
    int size = tmpSystem->size();
    systemCombined_.reset( new BlockMatrix_Type ( size ) );

    for (int i=0; i<size; i++) {
        DomainConstPtr_Type dom = this->getDomain(i);
        for (int j=0; j<size; j++) {
            if ( tmpSystem->blockExists(i,j) ) {
                LO maxNumEntriesPerRow = tmpSystem->getBlock(i,j)->getGlobalMaxNumRowEntries();
                MatrixPtr_Type matrix = Teuchos::rcp( new Matrix_Type( tmpSystem->getBlock(i,j)->getMap(), maxNumEntriesPerRow ) );
                
                systemCombined_->addBlock( matrix, i, j );

            }
            else if (systemMass_->blockExists(i,j)) {
                LO maxNumEntriesPerRow = systemMass_->getBlock(i,j)->getGlobalMaxNumRowEntries();
                MatrixPtr_Type matrix = Teuchos::rcp( new Matrix_Type( systemMass_->getBlock(i,j)->getMap(), maxNumEntriesPerRow) );
                systemCombined_->addBlock( matrix, i, j );
            }
        }
    }
    SmallMatrix<SC> ones( size , Teuchos::ScalarTraits<SC>::one());
    SmallMatrix<SC> zeros( size , Teuchos::ScalarTraits<SC>::zero());
    systemMass_->addMatrix( massParameters_, systemCombined_, zeros );
    tmpSystem->addMatrix( timeParameters_, systemCombined_, ones );
    
    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            if ( systemCombined_->blockExists(i,j) ) {
                if ( tmpSystem->blockExists(i,j) ) {
                    MatrixPtr_Type matrix = tmpSystem->getBlock(i,j);
                    MatrixConstPtr_Type matrixConstComb = systemCombined_->getBlock(i,j);
                    MatrixPtr_Type matrixComb = Teuchos::rcp_const_cast<Matrix_Type>(matrixConstComb);
                    matrixComb->fillComplete( matrix->getMap("domain"), matrix->getMap("range") );
                }
                else if ( systemMass_->blockExists(i,j) ) {
                    MatrixPtr_Type matrix = systemMass_->getBlock(i,j);
                    MatrixConstPtr_Type matrixConstComb = systemCombined_->getBlock(i,j);
                    MatrixPtr_Type matrixComb = Teuchos::rcp_const_cast<Matrix_Type>(matrixConstComb);
                    matrixComb->fillComplete( matrix->getMap("domain"), matrix->getMap("range") );
                }
                else
                    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error,"TimeProblem: Combined systems possess a block which does not exist in partial systems.");
            }
        }
    }
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateRhs(){

    systemMass_->apply( *solutionPreviousTimesteps_[0], *problem_->getRhs(), massParameters_ );
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateMultistepRhs(vec_dbl_Type& coeff, int nmbToUse){

    problem_->getRhs()->putScalar(0.);

    int size = massParameters_.size();

    for (int i=0; i<nmbToUse; i++) {
        SmallMatrix<SC> tmpMassParameter(size);
        for (int r=0; r<size; r++) {
            for (int s=0; s<size; s++) {
                if (massParameters_[r][s]!=0.)
                    tmpMassParameter[r][s] = coeff[i];
            }
        }
        BlockMultiVectorPtr_Type tmpVector
            = Teuchos::rcp( new BlockMultiVector_Type( problem_->getRhs() ) );
        BlockMultiVectorPtr_Type tmpBlockVector = solutionPreviousTimesteps_[i];
        systemMass_->apply( *tmpBlockVector, *tmpVector, tmpMassParameter );
        problem_->getRhs()->update( 1., tmpVector, 1. );
    }

}


template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateMultistepRhsFSI(vec_dbl_Type& coeff, int nmbToUse){

    problem_->getRhs()->putScalar(0.);

    int size = massParameters_.size();

    for (int i=0; i<nmbToUse; i++) {
        SmallMatrix<SC> tmpMassParameter(size);
        for (int r=0; r<size; r++) {
            for (int s=0; s<size; s++) {
                if (massParameters_[r][s]!=0.) {
                    tmpMassParameter[r][s] = coeff[i];
                }

            }
        }
        BlockMultiVectorPtr_Type tmpVector
            = Teuchos::rcp( new BlockMultiVector_Type( problem_->getRhs() ) );
        BlockMultiVectorPtr_Type tmpBlockVector = solutionPreviousTimesteps_[i];
        systemMassPreviousTimeSteps_[i]->apply( *tmpBlockVector, *tmpVector, tmpMassParameter );
        problem_->getRhs()->update( 1., tmpVector, 1. );
    }

}


// Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
// Bei Newmark lautet dies:
// M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
// wobei u' = v (velocity) und u'' = w (acceleration).
template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateNewmarkRhs(double dt, double beta, double gamma, vec_dbl_Type coeff)
{
    // ######################
    // Wir reseten die rechte Seite (beim ersten Schritt ist dies die rechte Seite der DGL [= SourceTerm] und
    // spaeter ist es dann die rechte Seite des zeitintegrierten System aus dem vorherigen Zeitschritt)
    // und fuegen alle noetigen Terme bis auf den SourceTerm hinzu.
    // Der SourceTerm wir in AdvanceInTime() hinzuaddiert.
    // BEACHTE: Bei Struktur gibt es nur einen Block, z.B. tempVector1[0].
    // ######################

    // Temperaerer Vektor fuer die Vektoradditionen in der rechte Seite.
    // ACHTUNG: BlockMultiVector_Type tempVector1(*(solutionPreviousTimesteps_.at(0)))
    // veraendert solutionPreviousTimesteps_.at(0)!!! NICHT NUTZEN!!!
    BlockMultiVectorPtrArray_Type tempVector1(1);
    // tempVector1 = u_n (in der neuen Zeiteration)
    tempVector1.at(0) = Teuchos::rcp(new BlockMultiVector_Type(solutionPreviousTimesteps_.at(0)));

    // tempVector1 = \frac{1}{dt^2*beta}*u_n
    tempVector1.at(0)->scale(1.0/(dt*dt*beta));

    // tempVector1 = tempVector1 + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n
    tempVector1.at(0)->update(1.0/(dt*beta), *(velocityPreviousTimesteps_[0]), (0.5 - beta)/beta, *(accelerationPreviousTimesteps_[0]), 1.0);


    // In tempVector2 steht dann das Endergebniss am Ende.
    // Siehe ACHTUNG von oben.
    BlockMultiVectorPtrArray_Type tempVector2(1);
    tempVector2.at(0) = Teuchos::rcp(new BlockMultiVector_Type(solutionPreviousTimesteps_.at(0)));

    // Koeffizient, um die mit \frac{rho}{dt^2*beta} skalierte Massematrix wieder zur Normalen zur machen
    // Skalierung mit \rho ist allerdings notwendig!
    int size = massParameters_.size();
    SmallMatrix<double> tmpMassParameter(size);
    for(int r = 0; r < size; r++)
    {
        for(int s = 0; s < size; s++)
        {
            if(massParameters_[r][s] != 0.0)
            {
                tmpMassParameter[r][s] = coeff.at(0);
            }

        }
    }

    // tempVector2 = tmpMassParameter.*M*tempVector1
    systemMass_->apply(*(tempVector1.at(0)), *(tempVector2.at(0)), tmpMassParameter);

    // RHS = 0*RHS + 1.0*tempVector2
    problem_->getRhs()->update(1.0, *(tempVector2.at(0)), .0);

}


template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorPtr_Type TimeProblem<SC,LO,GO,NO>::getSolution(){

    return problem_->getSolution();
}
    
template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorConstPtr_Type TimeProblem<SC,LO,GO,NO>::getSolutionConst() const {
    return problem_->getSolution();
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorPtr_Type TimeProblem<SC,LO,GO,NO>::getResidual(){
    
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    return nonLinProb->getResidualVector ();
}
    
template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorConstPtr_Type TimeProblem<SC,LO,GO,NO>::getResidualConst() const{
    
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    return nonLinProb->getResidualVector ();
}
    
template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorPtr_Type TimeProblem<SC,LO,GO,NO>::getSolutionPreviousTimestep(){

    return solutionPreviousTimesteps_[0];
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorPtrArray_Type TimeProblem<SC,LO,GO,NO>::getSolutionAllPreviousTimestep(){

    return solutionPreviousTimesteps_;
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::initializeCombinedSystems() const{

    BlockMatrixConstPtr_Type blockMatrixProblem = problem_->getSystem();

    int size = blockMatrixProblem->size();

    for (int i=0; i<size; i++) {
        for (int j=0; j<size; j++) {
            MatrixPtr_Type matrix;
            bool foundMatrix = false;
            if ( blockMatrixProblem->blockExists(i,j) ) {
                MatrixConstPtr_Type matrixProblem = blockMatrixProblem->getBlockConst(i,j);
                MapConstPtr_Type map = matrixProblem->getMap();
                MapPtr_Type mapNonConst = Teuchos::rcp_const_cast<Map_Type> (map);
                matrix = Teuchos::rcp( new Matrix_Type( mapNonConst, matrixProblem->getGlobalMaxNumRowEntries() ) ); //We should build matrices with graphs!
                foundMatrix = true;
            }
            else if( systemMass_->blockExists(i,j) && !foundMatrix ){
                MatrixConstPtr_Type matrixMass = systemMass_->getBlockConst(i,j);
                MapConstPtr_Type map = matrixMass->getMap();
                MapPtr_Type mapNonConst = Teuchos::rcp_const_cast<Map_Type> (map);
                matrix = Teuchos::rcp( new Matrix_Type( mapNonConst, matrixMass->getGlobalMaxNumRowEntries() ) ); //We should build matrices with graphs!
                foundMatrix = true;
            }
            if (foundMatrix)
                systemCombined_->addBlock( matrix, i, j );

        }
    }
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::assembleMassSystem( ) const {


    ProblemPtr_Type tmpProblem;
    SC eps100 = 100.*Teuchos::ScalarTraits<SC>::eps();

    // Die Massematrix wird noch mit der Dichte \rho skaliert
    double density = parameterList_->sublist("Parameter").get("Density",1.);

    int size = problem_->getSystem()->size();
    systemMass_->resize( size );
    int dofsPerNode;
    for (int i=0; i<size; i++ ) {
        dofsPerNode = problem_->getDofsPerNode(i);
        MatrixPtr_Type M;
        if ( timeStepDef_[i][i]>0 ) {
            if (dofsPerNode>1) {
                M = Teuchos::rcp(new Matrix_Type( this->getDomain(i)->getMapVecFieldUnique(), dimension_*this->getDomain(i)->getApproxEntriesPerRow() ) );
                feFactory_->assemblyMass(dimension_, problem_->getFEType(i), "Vector", M, true);

                M->resumeFill();
                M->scale(density);
                M->fillComplete( this->getDomain(i)->getMapVecFieldUnique(), this->getDomain(i)->getMapVecFieldUnique());

                systemMass_->addBlock(M, i, i);
            }
            else{
                M = Teuchos::rcp(new Matrix_Type( this->getDomain(i)->getMapUnique(), this->getDomain(i)->getApproxEntriesPerRow() ) );
                feFactory_->assemblyMass(dimension_, problem_->getFEType(i), "Scalar", M, true);

                M->resumeFill();
                M->scale(density);
                M->fillComplete( this->getDomain(i)->getMapUnique(), this->getDomain(i)->getMapUnique());

                systemMass_->addBlock(M, i, i);
            }
        }
        else{
            for (int j=0; j<size; j++) {
                if (timeStepDef_[i][j]==2) {
                    if (dofsPerNode>1) {
                        M = Teuchos::rcp(new Matrix_Type( this->getDomain(i)->getMapVecFieldUnique(), this->getDomain(i)->getApproxEntriesPerRow() ) );
                        feFactory_->assemblyMass(dimension_, problem_->getFEType(i), "Vector", M, true);
                        
                        M->resumeFill();
                        M->scale(density);
                        M->fillComplete( this->getDomain(i)->getMapVecFieldUnique(), this->getDomain(i)->getMapVecFieldUnique());
                        
                        systemMass_->addBlock(M, i, j);
                    }
                    else{
                        M = Teuchos::rcp(new Matrix_Type( this->getDomain(i)->getMapUnique(), this->getDomain(i)->getApproxEntriesPerRow() ) );
                        feFactory_->assemblyMass(dimension_, problem_->getFEType(i), "Scalar", M, true);
                        
                        M->resumeFill();
                        M->scale(density);
                        M->fillComplete( this->getDomain(i)->getMapUnique(), this->getDomain(i)->getMapUnique());
                        
                        systemMass_->addBlock(M, i, j);
                    }
                }
            }
        }
    }
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::setTimeParameters(SmallMatrix<double> &massParameters, SmallMatrix<double> &timeParameters){

    timeParameters_ = timeParameters;

    massParameters_ = massParameters;

}

template<class SC,class LO,class GO,class NO>  
bool TimeProblem<SC,LO,GO,NO>::getVerbose(){

    return verbose_;
}

template<class SC,class LO,class GO,class NO>
double TimeProblem<SC,LO,GO,NO>::calculateResidualNorm(){

    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    Teuchos::Array<SC> res(1);
    nonLinProb->getResidualVector()->norm2(res);
    return res[0];
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::calculateNonLinResidualVec( std::string type, double time ) const{

    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    
    if (type == "external") { // we get the complete residual from AceGEN
        nonLinProb->calculateNonLinResidualVec( "external" );
    }
    else{
        // rhs and sourceterm is accounted for in calculateNonLinResidualVec.
        // rhs must bet assembled (computed) correctly in DAESolver (e.g. M/dt*u_t). sourceterm aswell
        nonLinProb->calculateNonLinResidualVec( timeParameters_ ,type, time );

        // for FSI we need to reassemble the massmatrix system if the mesh was moved for geometry implicit computations
        if (this->parameterList_->sublist("Parameter").get("FSI",false) ){
            bool geometryExplicit = this->parameterList_->sublist("Parameter").get("Geometry Explicit",true);
            if( !geometryExplicit ) {
                typedef FSI<SC,LO,GO,NO> FSI_Type;
                typedef Teuchos::RCP<FSI_Type> FSIPtr_Type;
                
                MatrixPtr_Type massmatrix;                
                FSIPtr_Type fsi = Teuchos::rcp_dynamic_cast<FSI_Type>( this->problem_ );
                fsi->setFluidMassmatrix( massmatrix );
                this->systemMass_->addBlock( massmatrix, 0, 0 );
            }
        }
        // we need to add M/dt*u_(t+1)^k (the last results of the nonlinear method) to the residualVec_
        //Copy
        BlockMultiVectorPtr_Type tmpMV = Teuchos::rcp(new BlockMultiVector_Type( nonLinProb->getSolution() ) );
        tmpMV->putScalar(0.);
        
        systemMass_->apply( *nonLinProb->getSolution(), *tmpMV, massParameters_ );

        if (type=="reverse")// reverse: b-Ax
            nonLinProb->getResidualVector()->update(-1,*tmpMV,1.);//this=1.*this + -1.*tmpMV
        else if(type=="standard")// standard: Ax-b
            nonLinProb->getResidualVector()->update(1,*tmpMV,1.);
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Unkown type to compute the residual for a time problem.");
        
        if (type=="reverse")// we set the Dirichlet BC to the residual
            this->bcFactory_->setBCMinusVector( nonLinProb->getResidualVector(), nonLinProb->getSolution(), time );
        else if(type=="standard")// we set the negative Dirichlet BC to the residual
            this->bcFactory_->setVectorMinusBC( nonLinProb->getResidualVector(), nonLinProb->getSolution(), time );
            
    }

    
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::setBoundaries(double time){

    BlockMultiVectorPtr_Type tmpRhs = problem_->getRhs();

    bcFactory_->set( systemCombined_, tmpRhs, time );
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::setBoundariesRHS(double time){

    BlockMultiVectorPtr_Type tmpRhs = problem_->getRhs();

    bcFactory_->setRHS( tmpRhs, time );
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::setBoundariesSystem() const{

    bcFactory_->setSystem(systemCombined_);

}

template<class SC,class LO,class GO,class NO>
int TimeProblem<SC,LO,GO,NO>::solveUpdate(  ){

    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    
    *nonLinProb->previousSolution_ = *nonLinProb->getSolution();
    int its = this->solve( nonLinProb->residualVec_ );

    return its;
}

template<class SC,class LO,class GO,class NO>
int TimeProblem<SC,LO,GO,NO>::solveAndUpdate( const std::string& criterion, double& criterionValue ){

    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    

    int its = solveUpdate(  );

    if (criterion=="Update") {
        Teuchos::Array<SC> updateNorm(1);
        nonLinProb->getSolution()->norm2(updateNorm());
        criterionValue = updateNorm[0];
    }


    if (criterion=="ResidualAceGen") {
        nonLinProb->getSolution()->update( 1., *nonLinProb->previousSolution_, -1. );
        nonLinProb->assemble( "SetSolutionNewton" );
        //nonLinProb->assembleExternal( "OnlyUpdate" ); // CH 04.06.2020: Ist das richtig?
    }
    else
        nonLinProb->getSolution()->update( 1., *nonLinProb->previousSolution_, 1. );

    return its;
}

template<class SC,class LO,class GO,class NO>
int TimeProblem<SC,LO,GO,NO>::solve( BlockMultiVectorPtr_Type rhs ){

    int its;
    if (verbose_)
        std::cout << "-- Solve System ..." << std::endl;
    {
        std::string type = parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
        LinearSolver<SC,LO,GO,NO> linSolver;
        its = linSolver.solve( this, rhs, type ); // if rhs is null. Then the rhs_ of this is used in the linear solver
    }
    if (verbose_)
        std::cout << " done. -- " << std::endl;

    return its;
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateSolutionPreviousStep(){

    if (solutionPreviousTimesteps_.size()==0)
        solutionPreviousTimesteps_.resize(1);

    solutionPreviousTimesteps_[0] = Teuchos::rcp( new BlockMultiVector_Type( problem_->getSolution() ) );
    
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateSolutionMultiPreviousStep(int nmbSteps){

    int size = solutionPreviousTimesteps_.size();
    if (size<nmbSteps &&  size > 0) {
        BlockMultiVectorPtr_Type toAddMVreset = Teuchos::rcp( new BlockMultiVector_Type( solutionPreviousTimesteps_[size-1] ) );
        solutionPreviousTimesteps_.push_back( toAddMVreset );
    }
    else if(size == 0)
        solutionPreviousTimesteps_.resize(1);
    else{
        for (int i=size-1; i>0; i--)
            solutionPreviousTimesteps_[i] = Teuchos::rcp( new BlockMultiVector_Type( solutionPreviousTimesteps_[i-1] ) );
    }
    
    solutionPreviousTimesteps_[0] = Teuchos::rcp( new BlockMultiVector_Type( problem_->getSolution() ) );
    
}


template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateSystemMassMultiPreviousStep(int nmbSteps){

    int size = systemMassPreviousTimeSteps_.size();
    if (size<nmbSteps &&  size > 0) {// only works for BDF2, for BDF3 we would have to adjust below else case
        BlockMatrixPtr_Type toAddMatrixreset = Teuchos::rcp( new BlockMatrix_Type( systemMassPreviousTimeSteps_[size-1] ) );

        systemMassPreviousTimeSteps_.push_back( toAddMatrixreset );
    }
    else if(size == 0)
        systemMassPreviousTimeSteps_.resize(1);
    else{
        for (int i=size-1; i>0; i--)
            systemMassPreviousTimeSteps_[i] = Teuchos::rcp( new BlockMatrix_Type( systemMassPreviousTimeSteps_[i-1] ) );
    }

    systemMassPreviousTimeSteps_[0] = Teuchos::rcp( new BlockMatrix_Type( systemMass_ ) );
    
}


// Beachte: Diese Funktion wird zu Beginn jeder time-loop aufgerufen!!!
template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::updateSolutionNewmarkPreviousStep(double dt, double beta, double gamma)
{
    // ########################
    // Fuer die Newmark-Variable u (displacement)
    // ########################
    // Zur Berechnung von u'_{n+1} und u''_{n+1} wird die letzte Loesung u_{n+1} und die Vorherige u_n gebraucht.
    // Da wir aber bereits im neuen Zeitschritt sind (vgl. Zeitpunkt des Aufrufs der Funktion), ist u_n = u_{n+1} und u_{n-1} = u_n.
    // Wir benoetigen solutionPreviousTimesteps_ mit zwei Eintraegen.
    int size = solutionPreviousTimesteps_.size();

    if(size < 2 &&  size > 0) // Sofern der Vektor solutionPreviousTimesteps_ noch nicht komplett belegt ist (bei Newmark benoetigen wir als vergangene Loesung noch u_n)
    {
        BlockMultiVectorPtr_Type toAddMVreset = Teuchos::rcp(new BlockMultiVector_Type(solutionPreviousTimesteps_.at(size-1)));
        // Fuege die z.B. zweite Loesung u_2 hinten drauf, sofern eine vergangene Loesung bei der Zeitintegration gebraucht wird
        solutionPreviousTimesteps_.push_back(toAddMVreset);
    }
    else if(size == 0) // Falls noch keine alte Loesung vorhanden ist
    {
        solutionPreviousTimesteps_.resize(1);
    }
    else // Verschiebe alle vergangenen Loesungen
    {
        for(int i = size-1; i > 0; i--)
        {
            solutionPreviousTimesteps_.at(i) = Teuchos::rcp(new BlockMultiVector_Type(solutionPreviousTimesteps_.at(i-1)));
        }

    }

    solutionPreviousTimesteps_.at(0) = Teuchos::rcp(new BlockMultiVector_Type(problem_->getSolution()));

    // Bzgl. der neuen Zeititeration steht am Ende nun an der Stelle 0 die vergangene Loesung u_n und
    // an der Stelle die 1 davor, also u_{n-1}. Dies entspricht u_{n+1} und u_n, wenn man u'_{n+1} und
    // u''_{n+1} bereits am Ende der letzten Zeititeration berechnet.

    // ########################
    // Fuer die Newmark-Variablen u' = v (velocity) und u'' = w (acceleration)
    // ########################
    if(velocityPreviousTimesteps_.size() == 0)
    {
        // Hier steht noch kein MultiVector drin
        velocityPreviousTimesteps_.resize(1);
        accelerationPreviousTimesteps_.resize(1);

        // Am Anfang haben wir als Startloesung die Nulloesung.
        // Wir wollen als Startwerte ebenso u' = u'' = 0 setzen.
        // Schreibe die Startloesung als MultiVector hinein (damit dort irgendetwas steht) und setze diese dann auf Null
        velocityPreviousTimesteps_.at(0) = Teuchos::rcp(new BlockMultiVector_Type(problem_->getSolution()));
        accelerationPreviousTimesteps_.at(0) = Teuchos::rcp(new BlockMultiVector_Type(problem_->getSolution()));
        velocityPreviousTimesteps_.at(0)->putScalar(0.0);
        accelerationPreviousTimesteps_.at(0)->putScalar(0.0);
    }
    else
    {
        // ########################
        // u'_{n+1} = \frac{gamma}{dt*beta}*(u_{n+1} - u_n) + (1 - \frac{gamma}{beta})*u'_n + dt*\frac{beta - 0.5*gamma}{beta}*u''_n, vgl. MA
        // ########################

        // Speichere die alte (nicht geupdatete) velocity u'_n fuer die Berechnung von u''_{n+1} ab.
        // Siehe hier drunter fuer ACHTUNG.
        BlockMultiVectorPtrArray_Type velocityOld;
        velocityOld.resize(1);
        velocityOld.at(0) = Teuchos::rcp(new BlockMultiVector_Type(velocityPreviousTimesteps_.at(0)));


        // Berechne \frac{gamma}{dt*beta}*(u_{n+1} - u_n):
        // ACHTUNG: BlockMultiVector_Type tmpVector1(*(solutionPreviousTimesteps_.at(0)))
        // veraendert solutionPreviousTimesteps_.at(0)!!! NICHT NUTZEN!!!
        // BlockMultiVector_Type tmpVector1(*(solutionPreviousTimesteps_.at(0)));
        BlockMultiVectorPtrArray_Type tmpVector1;
        tmpVector1.resize(1);
        // tmpVector1 = u_{n+1}
        tmpVector1.at(0) = Teuchos::rcp(new BlockMultiVector_Type(solutionPreviousTimesteps_.at(0)));

        // Funktionsaufruf: Update(ScalarA, A, ScalarThis) => this = ScalarThis*this + ScalarA*A
        tmpVector1.at(0)->update(-gamma/(dt*beta), *(solutionPreviousTimesteps_.at(1)), gamma/(dt*beta));

        // Berechne (1 - \frac{gamma}{beta})*u'_n
        velocityPreviousTimesteps_.at(0)->scale(1.0 - (gamma/beta));

        // Addiere noch \tmpVector1 + dt*\frac{beta - 0.5*gamma}{beta}*u''_n HINZU
        // Funktionsaufruf: Update(ScalarA, A, ScalarB, B, ScalarThis) => this = ScalarThis*this + ScalarA*A + ScalarB*B
        velocityPreviousTimesteps_.at(0)->update(1.0, *(tmpVector1.at(0)), dt*((beta - 0.5*gamma)/(beta)), *(accelerationPreviousTimesteps_.at(0)), 1.0);


        // ########################
        // u''_{n+1} = \frac{1}{dt*dt*beta}*(u_{n+1} - u_n) - \frac{1}{dt*beta})*u'_n - \frac{0.5 - beta}{beta}*u''_n, vgl. MA
        // ########################

        // Berechne \frac{1}{dt*dt*beta}*(u_{n+1} - u_n):
        // Beachte ACHTUNG von oben.
        BlockMultiVectorPtrArray_Type tmpVector2;
        tmpVector2.resize(1);
         // tmpVector2 = u_{n+1}
        tmpVector2.at(0) = Teuchos::rcp(new BlockMultiVector_Type(solutionPreviousTimesteps_.at(0)));

        // Funktionsaufruf: Update(ScalarA, A, ScalarThis) => this = ScalarThis*this + ScalarA*A
        tmpVector2.at(0)->update(-1.0/(dt*dt*beta), *(solutionPreviousTimesteps_.at(1)), 1.0/(dt*dt*beta));

        // Berechne -\frac{0.5 - beta}{beta}*u''_n. ACHTUNG VORZEICHEN!!!
        accelerationPreviousTimesteps_.at(0)->scale(-(0.5 - beta)/beta);

        // Addiere noch \tmpVector2 - \frac{1}{dt*beta})*u'_n HINZU
        // Funktionsaufruf: Update(ScalarA, A, ScalarB, B, ScalarThis) => this = ScalarThis*this + ScalarA*A + ScalarB*B
        accelerationPreviousTimesteps_.at(0)->update(1.0, *(tmpVector2.at(0)), -1.0/(dt*beta), *(velocityOld.at(0)), 1.0);
    }
}
template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::assembleSourceTerm( double time ){
    
    problem_->assembleSourceTerm(time); 

}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorPtr_Type TimeProblem<SC,LO,GO,NO>::getSourceTerm( ) {
    
    return problem_->getSourceTerm();
    
}

template<class SC,class LO,class GO,class NO>
bool TimeProblem<SC,LO,GO,NO>::hasSourceTerm( ) const{
    
    return problem_->hasSourceTerm();
    
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorPtr_Type TimeProblem<SC,LO,GO,NO>::getRhs(){
    
    return problem_->getRhs();
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMultiVectorPtr_Type TimeProblem<SC,LO,GO,NO>::getRhs() const{

    return problem_->getRhs();
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::DomainConstPtr_Type TimeProblem<SC,LO,GO,NO>::getDomain(int i) const {
    
    return problem_->getDomain(i);
}

template<class SC,class LO,class GO,class NO>
std::string TimeProblem<SC,LO,GO,NO>::getFEType(int i){

    return problem_->getFEType(i);
}

template<class SC,class LO,class GO,class NO>
int TimeProblem<SC,LO,GO,NO>::getDofsPerNode(int i){

    return problem_->getDofsPerNode(i);
}

template<class SC,class LO,class GO,class NO>
std::string TimeProblem<SC,LO,GO,NO>::getVariableName(int i){

    return problem_->getVariableName(i);
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMatrixPtr_Type TimeProblem<SC,LO,GO,NO>::getSystem(){

    return problem_->getSystem();
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMatrixPtr_Type TimeProblem<SC,LO,GO,NO>::getSystemCombined(){

    TEUCHOS_TEST_FOR_EXCEPTION(systemCombined_.is_null(), std::logic_error,"system combined of TimeProblem is null.");

    return systemCombined_;
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMatrixPtr_Type TimeProblem<SC,LO,GO,NO>::getSystemCombined() const{

    TEUCHOS_TEST_FOR_EXCEPTION(systemCombined_.is_null(), std::logic_error,"system combined of TimeProblem is null.");

    return systemCombined_;
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::ProblemPtr_Type TimeProblem<SC,LO,GO,NO>::getUnderlyingProblem(){

    return problem_;
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::addToRhs(BlockMultiVectorPtr_Type x){
    
    problem_->addToRhs(x);
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::BlockMatrixPtr_Type TimeProblem<SC,LO,GO,NO>::getMassSystem(){

    return systemMass_;
}

template<class SC,class LO,class GO,class NO>
ParameterListPtr_Type TimeProblem<SC,LO,GO,NO>::getParameterList(){

   return parameterList_;
}

template<class SC,class LO,class GO,class NO>
typename TimeProblem<SC,LO,GO,NO>::LinSolverBuilderPtr_Type TimeProblem<SC,LO,GO,NO>::getLinearSolverBuilder() const{

    return problem_->getLinearSolverBuilder();
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::getValuesOfInterest( vec_dbl_Type& values ){
    
    problem_->getValuesOfInterest( values );
}
 
template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::computeValuesOfInterestAndExport( ){
    
    problem_->computeValuesOfInterestAndExport( );
}
    
    
template<class SC, class LO, class GO, class NO>
Thyra::ModelEvaluatorBase::InArgs<SC> TimeProblem<SC,LO,GO,NO>::getNominalValues() const
{
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    
    return nonLinProb->getNominalValues();
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > TimeProblem<SC,LO,GO,NO>::get_x_space() const{
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    
    return nonLinProb->get_x_space();
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC> > TimeProblem<SC,LO,GO,NO>::get_f_space() const{
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    return nonLinProb->get_f_space();
}

template<class SC,class LO,class GO,class NO>
Thyra::ModelEvaluatorBase::InArgs<SC> TimeProblem<SC,LO,GO,NO>::createInArgs() const
{
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    return nonLinProb->createInArgs();
}

// NOX
template<class SC,class LO,class GO,class NO>
Thyra::ModelEvaluatorBase::OutArgs<SC> TimeProblem<SC,LO,GO,NO>::createOutArgsImpl() const
{
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    return nonLinProb->createOutArgsImpl();
}
    
template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > TimeProblem<SC,LO,GO,NO>::create_W_op()
{

    this->calculateNonLinResidualVec( "standard", time_ );
    this->assemble("Newton");
    
//    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "we need to fix the operators and preconditioners for TimeProblem NOX setup.");
    
    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( type  == "Monolithic" )
        return create_W_op_Monolithic( );
    else
        return create_W_op_Block( );
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > TimeProblem<SC,LO,GO,NO>::create_W_op_Monolithic()
{
    Teuchos::RCP<const Thyra::LinearOpBase<SC> > W_opConst = this->getSystemCombined()->getThyraLinOp();
    Teuchos::RCP<Thyra::LinearOpBase<SC> > W_op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC> >(W_opConst);
    return W_op;
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::LinearOpBase<SC> > TimeProblem<SC,LO,GO,NO>::create_W_op_Block()
{
    
    BlockMatrixPtr_Type system = this->getSystemCombined();
    
    Teuchos::RCP<const ThyraBlockOp_Type> W_opBlocksConst = system->getThyraLinBlockOp();
    Teuchos::RCP<ThyraBlockOp_Type> W_opBlocks = Teuchos::rcp_const_cast<ThyraBlockOp_Type >(W_opBlocksConst);
    Teuchos::RCP<ThyraOp_Type> W_op = Teuchos::rcp_dynamic_cast<ThyraOp_Type >(W_opBlocks);
    
    return W_op;
}

template<class SC,class LO,class GO,class NO>
Teuchos::RCP<Thyra::PreconditionerBase<SC> > TimeProblem<SC,LO,GO,NO>::create_W_prec()
{
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    
    if (!nonLinProb->getPreconditionerConst()->isPreconditionerComputed()) {
        nonLinProb->initializeSolverBuilder();
        
        std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
        this->setBoundariesSystem();
        
        // if ( type == "Teko" || type == "FaCSI-Teko" || || type == "FaCSI-Block" || type =="Diagonal" ) { //we need to construct the whole preconditioner if Teko is used
        //     nonLinProb->setupPreconditioner( type );
        //     precInitOnly_ = false;
        // }
        // else{
            nonLinProb->setupPreconditioner( type ); //nonLinProb->initializePreconditioner( type );
            precInitOnly_ = false;
        // }
    }
    
    Teuchos::RCP<const Thyra::PreconditionerBase<SC> > thyraPrec =  nonLinProb->getPreconditionerConst()->getThyraPrecConst();
    Teuchos::RCP<Thyra::PreconditionerBase<SC> > thyraPrecNonConst = Teuchos::rcp_const_cast<Thyra::PreconditionerBase<SC> >(thyraPrec);
    
    return thyraPrecNonConst;
    
}
    
template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::evalModelImpl( const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                              const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                              ) const
{
    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
    if ( !type.compare("Monolithic"))
        evalModelImplMonolithic( inArgs, outArgs );
    else
        evalModelImplBlock( inArgs, outArgs );
}

template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::evalModelImplMonolithic( const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                        const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs ) const
{
    
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;
    using Teuchos::ArrayView;
    using Teuchos::Array;

    TEUCHOS_TEST_FOR_EXCEPTION( inArgs.get_x().is_null(), std::logic_error, "inArgs.get_x() is null.");
    
    RCP< const Thyra::VectorBase< SC > > vecThyra = inArgs.get_x();
    
    RCP< Thyra::VectorBase< SC > > vecThyraNonConst = rcp_const_cast<Thyra::VectorBase< SC > >(vecThyra);
    
    BlockMultiVectorConstPtr_Type solConst = this->getSolutionConst();
    BlockMultiVectorPtr_Type sol = rcp_const_cast<BlockMultiVector_Type >(solConst);
    sol->fromThyraMultiVector(vecThyraNonConst);
    
    const RCP<Thyra::MultiVectorBase<SC> > f_out = outArgs.get_f();
    const RCP<Thyra::LinearOpBase<SC> > W_out = outArgs.get_W_op();
    const RCP<Thyra::PreconditionerBase<SC> > W_prec_out = outArgs.get_W_prec();
    
    typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
    typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;
    typedef RCP<const TpetraMatrix_Type> TpetraMatrixConstPtr_Type;
    
    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);
    
    if ( fill_f || fill_W || fill_W_prec ) {
        
        // ****************
        // Get the underlying xpetra objects
        // ****************
        if (fill_f) {
            
            this->calculateNonLinResidualVec( "standard", time_ );
                        
            BlockMultiVectorConstPtr_Type resConst = this->getResidualConst();
            BlockMultiVectorPtr_Type res = rcp_const_cast<BlockMultiVector_Type >(resConst);
            
            Teuchos::RCP<Thyra::MultiVectorBase<SC> > f_thyra = res->getThyraMultiVector();
            f_out->assign(*f_thyra);
        }
        
        TpetraMatrixPtr_Type W;
        if (fill_W) {
            
            this->assemble("Newton");

            this->setBoundariesSystem();
            
            Teuchos::RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
            Teuchos::RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);
            
            TpetraMatrixConstPtr_Type W_systemTpetra = this->getSystemCombined()->getMergedMatrix()->getTpetraMatrix();           
            TpetraMatrixPtr_Type W_systemTpetraNonConst = rcp_const_cast<TpetraMatrix_Type>(W_systemTpetra);
            
            //Tpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_systemXpetraNonConst);
            //Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
            
            Teuchos::RCP<TpetraMatrix_Type> tpetraMatTpetra = W_systemTpetraNonConst; //xTpetraMat.getTpetra_CrsMatrixNonConst();
            
            W_tpetraMat->resumeFill();
            
            for (auto i=0; i<tpetraMatTpetra->getMap()->getLocalNumElements(); i++) {
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
                typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;  //ArrayView< const LO > indices
                tpetraMatTpetra->getLocalRowView( i, indices, values);
                W_tpetraMat->replaceLocalValues( i,  indices, values);
            }
            W_tpetraMat->fillComplete();
            
        }
        
        if (fill_W_prec) {
            int newtonLimit = this->parameterList_->sublist("Parameter").get("newtonLimit",2);
            NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);

            if (precInitOnly_){

                if(nonLinProb->newtonStep_ < newtonLimit || this->parameterList_->sublist("Parameter").get("Rebuild Preconditioner every Newton Iteration",true) )
                {
                    this->problem_->setupPreconditioner( "Monolithic" );
                }
                else{
                    if (this->verbose_)
                        std::cout << " \t\t TimeProblem<SC,LO,GO,NO>::evalModelImplMonolithic: Skipping preconditioner reconstruction" << std::endl;
                }
            }
            else
                precInitOnly_ = true; // If a Monolithic preconditioner was constructed for the first time this variable is false. Because the preconditioner was not only initialized but already constructed. We can now set this variable to true to always setup all following preconditioners in the above if case
           
            nonLinProb->newtonStep_++;
            
        }
    }
}


template<class SC,class LO,class GO,class NO>
void TimeProblem<SC,LO,GO,NO>::evalModelImplBlock( const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                                   const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs ) const
{
    
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::rcp_const_cast;
    using Teuchos::ArrayView;
    using Teuchos::Array;
    
    TEUCHOS_TEST_FOR_EXCEPTION( inArgs.get_x().is_null(), std::logic_error, "inArgs.get_x() is null.");
    
    RCP< const Thyra::VectorBase< SC > > vecThyra = inArgs.get_x();
    
    RCP< Thyra::VectorBase< SC > > vecThyraNonConst = rcp_const_cast<Thyra::VectorBase< SC > >(vecThyra);
    
    RCP< Thyra::ProductVectorBase< SC > > vecThyraBlock = rcp_dynamic_cast<Thyra::ProductVectorBase< SC > > (vecThyraNonConst);
    BlockMultiVectorConstPtr_Type solConst = this->getSolutionConst();
    BlockMultiVectorPtr_Type sol = rcp_const_cast<BlockMultiVector_Type >(solConst);
    for (int i=0; i<sol->size(); i++)
        sol->getBlockNonConst(i)->fromThyraMultiVector( vecThyraBlock->getNonconstVectorBlock(i) );
    
    const RCP<Thyra::MultiVectorBase<SC> > f_out = outArgs.get_f();
    const RCP<Thyra::LinearOpBase<SC> > W_out = outArgs.get_W_op();
    const RCP<Thyra::PreconditionerBase<SC> > W_prec_out = outArgs.get_W_prec();
  
    typedef Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraMatrix_Type;
    typedef RCP<TpetraMatrix_Type> TpetraMatrixPtr_Type;
    typedef RCP<const TpetraMatrix_Type> TpetraMatrixConstPtr_Type;
     
    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);
    
    if ( fill_f || fill_W || fill_W_prec ) {
        
        // ****************
        // Get the underlying xpetra objects
        // ****************
        if (fill_f) {
            
            this->calculateNonLinResidualVec("standard" , time_);
            
            RCP< Thyra::ProductMultiVectorBase<SC> > f_outBlocks
                = rcp_dynamic_cast< Thyra::ProductMultiVectorBase<SC> > ( f_out );

            BlockMultiVectorConstPtr_Type resConst = this->getResidualConst();
            BlockMultiVectorPtr_Type res = rcp_const_cast<BlockMultiVector_Type >(resConst);

            for (int i=0; i<res->size(); i++) {
                RCP< Thyra::MultiVectorBase< SC > > res_i =
                    res->getBlockNonConst(i)->getThyraMultiVector();

                RCP< Thyra::MultiVectorBase< SC > > f_i = f_outBlocks->getNonconstMultiVectorBlock(i);
                f_i->assign(*res_i);
            }
        }
        
        TpetraMatrixPtr_Type W;
        if (fill_W) {
            
            typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;

            this->assemble("Newton");
            
            this->setBoundariesSystem();
            
            RCP<ThyraBlockOp_Type> W_blocks = rcp_dynamic_cast<ThyraBlockOp_Type>(W_out);

            BlockMatrixPtr_Type system = this->getSystemCombined();
            
            for (int i=0; i<system->size(); i++) {
                for (int j=0; j<system->size(); j++) {
                    if ( system->blockExists(i,j) ) {
                        RCP<const ThyraOp_Type> W = W_blocks->getBlock(i,j);
                        RCP<ThyraOp_Type> W_NonConst = rcp_const_cast<ThyraOp_Type>( W );
                        RCP<TpetraOp_Type> W_tpetra = tpetra_extract::getTpetraOperator( W_NonConst );
                        RCP<TpetraMatrix_Type> W_tpetraMat = Teuchos::rcp_dynamic_cast<TpetraMatrix_Type>(W_tpetra);
                        
                        TpetraMatrixConstPtr_Type W_matrixTpetra = system->getBlock(i,j)->getTpetraMatrix();   
                        TpetraMatrixPtr_Type W_matrixTpetraNonConst = rcp_const_cast<TpetraMatrix_Type>(W_matrixTpetra);
                        
                        //Xpetra::CrsMatrixWrap<SC,LO,GO,NO>& crsOp = dynamic_cast<Xpetra::CrsMatrixWrap<SC,LO,GO,NO>&>(*W_matrixXpetraNonConst);
                        //Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>& xTpetraMat = dynamic_cast<Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>&>(*crsOp.getCrsMatrix());
                        
                        RCP<TpetraMatrix_Type> tpetraMatTpetra = W_matrixTpetraNonConst;
                        
                        W_tpetraMat->resumeFill();
                        
                        for (auto i=0; i<tpetraMatTpetra->getMap()->getLocalNumElements(); i++) {
						    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_inds_host_view_type indices;  //ArrayView< const LO > indices
						    typename Tpetra::CrsMatrix<SC,LO,GO,NO>::values_host_view_type values;  //ArrayView< const LO > indices
						    tpetraMatTpetra->getLocalRowView( i, indices, values);
						    W_tpetraMat->replaceLocalValues( i,  indices, values);
						}
                        W_tpetraMat->fillComplete( W_tpetraMat->getDomainMap(), W_tpetraMat->getRangeMap() );
                    }
                }
            }
        }
        
        if (fill_W_prec) {
            std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
            NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);

            if (precInitOnly_){
                int newtonLimit = this->parameterList_->sublist("Parameter").get("newtonLimit",2);

                if(nonLinProb->newtonStep_ < newtonLimit || this->parameterList_->sublist("Parameter").get("Rebuild Preconditioner every Newton Iteration",true) )
                {
                    this->problem_->setupPreconditioner( type );
                }
                else{
                    if (this->verbose_)
                        std::cout << " \t\t TimeProblem<SC,LO,GO,NO>::evalModelImplBlock Skipping preconditioner reconstruction" << std::endl;
                }
            }
            else
                precInitOnly_ = true; // If a Teko preconditioner was constructed for the first time this variable is false. Because the preconditioner was not only initialized but already constructed. We can now set this variable to true to always setup all following preconditioners in the above if case
            
            nonLinProb->newtonStep_++;
        }
    }
}
template<class SC,class LO,class GO,class NO>
std::string TimeProblem<SC,LO,GO,NO>::description() const{ //reimplement description function and use description of underlying nonLinearProblem
    NonLinProbPtr_Type nonLinProb = Teuchos::rcp_dynamic_cast<NonLinProb_Type>(problem_);
    TEUCHOS_TEST_FOR_EXCEPTION(nonLinProb.is_null(), std::runtime_error, "Nonlinear problem is null.");
    return nonLinProb->description();
}
}
#endif
