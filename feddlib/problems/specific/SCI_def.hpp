#ifndef SCI_def_hpp
#define SCI_def_hpp
#include "SCI_decl.hpp"
#include <limits>
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/General/BCBuilder.hpp"


namespace FEDD {
// Funktionen fuer die rechte Seite der Struktur/ Chem/ Geometrie sind im jeweiligen Problem

template<class SC,class LO,class GO,class NO>
SCI<SC,LO,GO,NO>::SCI(const DomainConstPtr_Type &domainStructure, std::string FETypeStructure,
					const DomainConstPtr_Type &domainChem, std::string FETypeChem, vec2D_dbl_Type diffusionTensor, RhsFunc_Type reactionFunc,
                    ParameterListPtr_Type parameterListStructure, ParameterListPtr_Type parameterListChem,
                    ParameterListPtr_Type parameterListSCI, Teuchos::RCP<SmallMatrix<int> > &defTS):
NonLinearProblem<SC,LO,GO,NO>( parameterListSCI, domainChem->getComm() ),
// hasSourceTerm = drittes Arguement. assembleSourceTerm() fuer NS nicht programmiert.
// Deswegen hier erstmal false (default Parameter).
// Fuer Struktur hingegen ist default Parameter true, da programmiert.
problemStructure_(),
problemChem_(),
meshDisplacementOld_rep_(),
meshDisplacementNew_rep_(),
c_rep_(),
defTS_(defTS),
timeSteppingTool_(),
exporterEMod_(),
materialModel_( parameterListSCI->sublist("Parameter").get("Structure Model","SCI_NH") )
{

    this->initNOXParameters();
    
    this->addVariable( domainStructure, FETypeStructure, "d_s", domainStructure->getDimension() ); // Structure
    this->addVariable( domainChem, FETypeChem, "c", 1); // Chemistry scalar valued problem

    this->dim_ = this->getDomain(0)->getDimension();

    // The structure's own time problem keeps the Newmark state (u, u', u''): it
    // writes and reads it with the checkpoints of the SCI problem, so it gets the
    // same timestepping settings (Checkpointing, Restart, directories, times).
    parameterListStructure->sublist("Timestepping Parameter").setParameters(parameterListSCI->sublist("Timestepping Parameter"));

    if (materialModel_=="SCI_Linear"){
        problemStructure_ = Teuchos::rcp( new StructureProblem_Type( domainStructure, FETypeStructure, parameterListStructure ) );
        problemStructure_->initializeProblem();
    }
    else {
        problemStructureNonLin_ = Teuchos::rcp( new StructureNonLinProblem_Type( domainStructure, FETypeStructure, parameterListStructure) );
        problemStructureNonLin_->initializeProblem();
    }
    problemChem_ = Teuchos::rcp( new ChemProblem_Type( domainChem, FETypeChem, parameterListChem, diffusionTensor, reactionFunc ) );
    problemChem_->initializeProblem();

    //We initialize the subproblems. In the main routine, we need to call initializeFSI(). There, we first initialize the vectors of the FSI problem and then we set the pointers of the subproblems to the vectors of the full monolithic FSI system. This way all values are only saved once in the subproblems and can be used by the monolithic FSI system.
    meshDisplacementNew_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    meshDisplacementOld_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    
    d_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    c_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );
    exportedEMod_ = false;
    setUpTimeStep_=false;
    eModVec_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getElementMap() ) );
    chemistryExplicit_ =    parameterListSCI->sublist("Parameter").get("Chemistry Explicit",false);
    loadStepping_ =    parameterListSCI->sublist("Parameter").get("Load Stepping",false);
    externalForce_ =   parameterListSCI->sublist("Parameter").get("External Force",false);
    nonlinearExternalForce_ = parameterListSCI->sublist("Parameter").get("Nonlinear External Force",false);
    this->info();

    if ( chemistryExplicit_){
        exporterIterationsChem_ = Teuchos::rcp(new ExporterTxt());
        exporterIterationsChem_->setup( "linearIterations_chem", this->comm_ );
    }

    timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));

}


template<class SC,class LO,class GO,class NO>
SCI<SC,LO,GO,NO>::~SCI()
{
    if (!exporterEMod_.is_null()) {
       exporterEMod_->closeExporter();
    }
    if (!exporterChem_.is_null()) {
       exporterChem_->closeExporter();
    }
    if ( chemistryExplicit_){
        exporterIterationsChem_->closeExporter();
    }
}   

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::info()
{
    this->infoProblem();
    if (this->verbose_)
    {
        std::cout << "\t SCI:: info() ... " << std::endl;
        std::cout << "\t SCI:: loadStepping_ " << loadStepping_  << " ... " << std::endl;
        std::cout << "\t SCI:: chemistryExplicit_ " << chemistryExplicit_  << " ... " << std::endl;
        std::cout << "\t SCI:: externalForce_ " << externalForce_  << " ... " << std::endl;
        std::cout << "\t SCI:: nonlinearExternalForce_ " << nonlinearExternalForce_  << " ... " << std::endl;
    }
    
}

/*! 
    Generally in assemble() all linear equations are assembled. Suppose we have a linear elasticity and reaction-diffusion equation,
    we get a system with linear Diagonal entries.
    The Offdiagonal parts of the system Matrix can then either be treated as nonliarities or eliminated by treating them explicitly 
    in the time stepping scheme at hand.

    We start with the latter approach.

    type:
    " " -> assembly of constant matrices
    " EMod " -> assembly E module with concentration of previous timestep
    " MoveMesh " -> the move mesh operation moves the mesh, then the reaction diffusion includes the diplacement

*/
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::assemble( std::string type ) const
{

    if (type == "") {
        if (this->verbose_)
        {
            std::cout << "SCI_def: assemble() ... " << std::endl;
        }
    
        // Maybe nothing should happen here as there are no constant matrices
        this->system_.reset(new BlockMatrix_Type(2));
        if ( chemistryExplicit_){
            this->system_.reset(new BlockMatrix_Type(1));
            this->systemC_.reset(new BlockMatrix_Type(1));
        }


        MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        
        if (this->residualVec_.is_null())
            this->residualVec_.reset(new BlockMultiVector_Type(2));
        
        MatrixPtr_Type B(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type BT(new Matrix_Type(this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type C(new Matrix_Type( this->getDomain(1)->getMapUnique(),this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ));

        // For implicit the system is ordered differently with solid block in 0,0 and diffusion in 1,1
        
        // As the implementation is based on a 2x2 System we use a tmp system for the actual assembly routine.
        BlockMatrixPtr_Type systemTmp(new BlockMatrix_Type(2)); 
        systemTmp->addBlock(A,0,0);
        systemTmp->addBlock(BT,0,1);
        systemTmp->addBlock(B,1,0);
        systemTmp->addBlock(C,1,1);

        this->setupSubTimeProblems(this->problemChem_->getParameterList(), this->problemStructureNonLin_->getParameterList());

        this->problemTimeStructure_->assembleSourceTerm( 0. );
        this->problemTimeStructure_->addToRhs( this->problemTimeStructure_->getSourceTerm() );       
        this->problemTimeStructure_->setBoundariesRHS();

        this->setFromPartialVectorsInit();

        MultiVectorConstPtr_Type c; 
        if(chemistryExplicit_)
            c= this->problemTimeChem_->getSolution()->getBlock(0);
        else
            c = this->solution_->getBlock(1);

        MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
        d_rep_->importFromVector(d, true); 
        c_rep_->importFromVector(c, true);
    
        BlockMultiVectorPtr_Type blockSol = Teuchos::rcp( new BlockMultiVector_Type(2) );
        blockSol->addBlock(d_rep_,0);
        blockSol->addBlock(c_rep_,1);

        this->feFactory_->initAssembleFEAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 1,this->dim_,this->parameterList_);
        this->feFactory_->synchronizeTime(this->timeSteppingTool_);
        this->feFactory_->assemblyAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,c_rep_,d_rep_,systemTmp,this->residualVec_, this->parameterList_, "Jacobian", true/*call fillComplete*/);

        if(chemistryExplicit_){
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->systemC_->addBlock(systemTmp->getBlock(1,1),0,0);
        }
        else {
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->system_->addBlock(systemTmp->getBlock(0,1),0,1);
            this->system_->addBlock(systemTmp->getBlock(1,0),1,0);
            this->system_->addBlock(systemTmp->getBlock(1,1),1,1);
        }
  
        if ( chemistryExplicit_ && this->parameterList_->sublist("General").get("ParaViewExport",false)){
            exporterChem_ = Teuchos::rcp(new Exporter_Type());
            
            DomainConstPtr_Type dom = this->getDomain(1);

            int exportEveryXTimesteps = this->parameterList_->sublist("Exporter").get( "Export every X timesteps", 1 );
            std::string varName = "c";
            
            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( dom->getMesh() );
            exporterChem_->setup(varName, meshNonConst, dom->getFEType(), exportEveryXTimesteps, this->parameterList_);
            
            MultiVectorConstPtr_Type exportVector = this->problemTimeChem_->getSolution()->getBlock(0);
            
            exporterChem_->addVariable( exportVector, varName, "Scalar", 1, dom->getMapUnique() );
        }
        

    

        if (this->verbose_)
        {
            std::cout << "done -- " << std::endl;
        }
    }
    else
        reAssemble(type);

}
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::solveChemistryProblem() const
{
    // Move Mesh operation was already called for chemistry domain
    SmallMatrix<double> massCoeffChem(1);
    SmallMatrix<double> problemCoeffChem(1);
    SmallMatrix<int> defChem(1);
    double dt = timeSteppingTool_->get_dt();

    // Always bdf 1
    defChem[0][0] = 1;
    massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt;
    problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1);
    double  coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1);

    this->problemTimeChem_->setTimeDef(defChem);
    this->problemTimeChem_->setTimeParameters(massCoeffChem,problemCoeffChem);
    // 1. Assemble Mass System
    MatrixPtr_Type massmatrix;
    this->setChemMassmatrix(massmatrix);
    // 1. Assemble Chemisty Problem
    this->problemTimeChem_->assemble();
    // 2. Rhs
    this->computeChemRHSInTime();

    this->problemTimeChem_->combineSystems();

    this->problemTimeChem_->setBoundaries(timeSteppingTool_->currentTime()); 

    double iter= 0.;

    iter = this->problemTimeChem_->solve(this->problemTimeChem_->getRhs());

    exporterIterationsChem_->exportData(  timeSteppingTool_->currentTime() , iter);

    if (!exporterChem_.is_null())
            this->exporterChem_->save( this->timeSteppingTool_->currentTime() );

    // If we want to safe/export the solution, we need to call the function here, as the chemistry problem is decoupled from the rest
    // and does not encounter the usual export point within the time stepping loop
    // BlockMultiVectorPtrArray_Type solution; solution.resize(1);
    // solution[0] = Teuchos::rcp( new BlockMultiVector_Type( this->problemTimeChem_->getSolution()->getMap() ) );
    // solution[0]->addBlock(this->problemTimeChem_->getSolution()->getBlock(0),0);


}
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::reAssemble(std::string type) const
{

    double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);

    if (type == "MoveMesh"){
        moveMesh();
        return;
    }
    else if (type == "UpdateMeshDisplacement"){
        updateMeshDisplacement();
        return;
    }

    else if (type == "ComputeSolidRHSInTime"){
        if (this->verbose_)
            std::cout << "-- Assembly (ComputeSolidRHSInTime)" << '\n';
          
        computeSolidRHSInTime();
        return;
    }
    else if (type == "ComputeChemRHSInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (ComputeChemRHSInTime)" << '\n';
          
        computeChemRHSInTime();
        return;
    }
    else if (type == "UpdateChemInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (UpdateChemInTime)" << '\n';
        updateChemInTime();
        return;
    }

    else if(type == "UpdateTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateTime)" << '\n';

        updateTime();
        return;
    }
     else if(type == "SolveChemistryProblem")
    {
        if(this->verbose_)
            std::cout << "-- SolveChemistryProblem" << '\n';

        solveChemistryProblem();
        return;
    }
    else if(type == "SetBoundaries")
    {
        if(this->verbose_)
            std::cout << "-- set Boundaries" << '\n';

        setBoundariesSubProblems();
        return;
    }
    else if(type == "UpdateCoupling")
    {
        if(this->verbose_)
            std::cout << "-- Update Coupling" << '\n';

       
        MultiVectorPtr_Type solChemRep = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );
        solChemRep->importFromVector(this->problemTimeChem_->getSolution()->getBlock(0));

        if (materialModel_=="SCI_Linear"){
            int dim = this->getDomain(0)->getDimension();
            this->feFactory_->determineEMod(this->getDomain(0)->getFEType(),solChemRep,eModVec_,this->getDomain(1),this->parameterList_);
            double nu = this->parameterList_->sublist("Parameter").get("PoissonRatio",0.4);
            MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) ); // Structure-Matrix

            this->feFactory_->assemblyLinElasXDimE(this->dim_, this->getDomain(0)->getFEType(), A, eModVec_, nu, true);
            this->problemStructure_->system_->addBlock(A,0,0);// assemble(); //
            this->system_->addBlock( this->problemStructure_->system_->getBlock(0,0), 0, 0);
        }
        else{
            this->problemStructureNonLin_->updateConcentration(solChemRep);          

        }
        this->moveMesh();

        this->problemChem_->assemble();     


        return;
    }
    

       
    if(type == "FixedPoint")
    {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "should always be called by the residual and never here.");
    }
    else if(type == "Newton")
    {

        if(this->verbose_)
            std::cout << "-- reassemble Newton " << '\n';

        // Maybe nothing should happen here as there are no constant matrices
        this->system_.reset(new BlockMatrix_Type(2));
        if (chemistryExplicit_)
            this->system_.reset(new BlockMatrix_Type(1));


        MatrixPtr_Type A(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );       
        MatrixPtr_Type B(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type BT(new Matrix_Type(this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) );
        MatrixPtr_Type C(new Matrix_Type( this->getDomain(1)->getMapUnique(),this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ));

        // For implicit the system is ordered differently with solid block in 0,0 and diffusion in 1,1
        
        // As the implementation is based on a 2x2 System we use a tmp system for the actual assembly routine.
        BlockMatrixPtr_Type systemTmp(new BlockMatrix_Type(2)); 
        systemTmp->addBlock(A,0,0);
        systemTmp->addBlock(BT,0,1);
        systemTmp->addBlock(B,1,0);
        systemTmp->addBlock(C,1,1);
        
        MultiVectorConstPtr_Type c; 
        if(chemistryExplicit_)
            c= this->problemTimeChem_->getSolution()->getBlock(0);
        else
            c = this->solution_->getBlock(1);

        c_rep_->importFromVector(c, true);

        MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
        d_rep_->importFromVector(d, true); 
    
        this->feFactory_->assemblyAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,c_rep_,d_rep_,systemTmp,this->residualVec_, this->parameterList_, "Jacobian", true/*call fillComplete*/);

        if(chemistryExplicit_){
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->systemC_->addBlock(systemTmp->getBlock(1,1),0,0);
        }
        else{
            this->system_->addBlock(systemTmp->getBlock(0,0),0,0);
            this->system_->addBlock(systemTmp->getBlock(0,1),0,1);
            this->system_->addBlock(systemTmp->getBlock(1,0),1,0);
            this->system_->addBlock(systemTmp->getBlock(1,1),1,1);
        }

        if(nonlinearExternalForce_)
            computeSolidRHSInTime();

        if(this->verbose_)
            std::cout << " done -- " << '\n';   
                    
        
    }
      // ########################
        // Prec def -- experimental -- if explicit system we use 'diffusion matrix' as pressure matrix
       /* string precType = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
        if ( precType == "Diagonal" || precType == "Triangular" ) {
            
            this->problemChem_->assemble();
            this->getPreconditionerConst()->setPressureMassMatrix( this->problemChem_->system_->getBlock(0,0) );
        }*/

}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const
{
    
    MultiVectorConstPtr_Type c; 
    if(chemistryExplicit_)
        c= this->problemTimeChem_->getSolution()->getBlock(0);
    else
        c = this->solution_->getBlock(1);

    c_rep_->importFromVector(c, true);

    MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
    d_rep_->importFromVector(d, true); 

    BlockMultiVectorPtr_Type resTmp(new BlockMultiVector_Type(2));
    MultiVectorPtr_Type res_d = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    MultiVectorPtr_Type res_c = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );
    resTmp->addBlock(res_d,0);
    resTmp->addBlock(res_c,1);

    this->feFactory_->assemblyAceDeformDiffu(this->dim_, this->getDomain(1)->getFEType(), this->getDomain(0)->getFEType(), 2, 1,this->dim_,c_rep_,d_rep_,this->system_,resTmp, this->parameterList_, "Rhs", true/*call fillComplete*/);
    
    this->residualVec_->addBlock(Teuchos::rcp_const_cast<MultiVector_Type>(resTmp->getBlock(0)),0);

    if(!chemistryExplicit_){
        this->residualVec_->addBlock(Teuchos::rcp_const_cast<MultiVector_Type>(resTmp->getBlock(1)),1);
    }

    if(nonlinearExternalForce_)
        computeSolidRHSInTime();

    if (!type.compare("standard")){
        if(this->verbose_)
            std::cout << " Residual Type : " << type  << std::endl;
        this->residualVec_->getBlockNonConst(0)->update(-1.,*this->rhs_->getBlockNonConst(0),1.);
    
    }
    else if(!type.compare("reverse")){
        if(!chemistryExplicit_)
            this->residualVec_->getBlockNonConst(1)->scale(-1.0);

        if(this->verbose_)
            std::cout << " Residual Type : " << type  << std::endl;

        this->residualVec_->getBlockNonConst(0)->update(1.,*this->rhs_->getBlockNonConst(0),-1.);
        
    }

    /*this->problemChem_->assemble();

    MultiVectorPtr_Type residualChemSCI =  Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(1) );
    this->problemChem_->getSystem()->getBlock(0,0)->apply( *this->problemChem_->getSolution()->getBlock(0), *residualChemSCI, Teuchos::NO_TRANS, -1. ); // y= -Ax + 0*y
    MultiVectorPtr_Type resChemNonConst = Teuchos::rcp_const_cast<MultiVector_Type> ( this->residualVec_->getBlock(1) );
    resChemNonConst->update(1., *this->problemChem_->getRhs()->getBlock(0), 1.);        
    */
    if(!chemistryExplicit_){
        MultiVectorPtr_Type resChemNonConst = Teuchos::rcp_const_cast<MultiVector_Type> ( this->residualVec_->getBlock(1) );
        resChemNonConst->update(1., *this->problemChem_->getRhs()->getBlock(0), 1.);
    }
        

   
   // might also be called in the sub calculateNonLinResidualVec() methods which were used above
    if (type == "standard")
        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    else if (type == "reverse"){
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
    }

    this->setBoundariesRHS(this->timeSteppingTool_->currentTime());

    // Adaptive time stepping: an element of any process that cannot compute its state fails the
    // time step. A NaN residual on every process ends the Newton iteration (its finite value test).
    // With acceptElementFailures_ the failures are only counted.
    if (adaptiveStep_){
        std::string message;
        int failedHere = this->feFactory_->assemblyFEElementsFailed(message);
        int failedAnywhere = 0;
        Teuchos::reduceAll<int, int>( *this->comm_, Teuchos::REDUCE_SUM, failedHere, Teuchos::outArg(failedAnywhere) );
        if (failedAnywhere > 0){
            if (failedHere > 0 && !elementFailed_)
                std::cout << " [adaptive] process " << this->comm_->getRank() << ": " << message << std::endl;
            elementFailed_ = true;
            failedElements_ = failedAnywhere;
            if (!acceptElementFailures_)
                this->residualVec_->putScalar( std::numeric_limits<SC>::quiet_NaN() );
        }
    }


    Teuchos::Array<SC> norm_d(1); 
    this->residualVec_->getBlock(0)->norm2(norm_d);
    
    if(this->verbose_)
        std::cout << "2-Norm of residual of displacement: " << norm_d[0] << std::endl;

    if(!chemistryExplicit_){
        Teuchos::Array<SC> norm_c(1); 
        this->residualVec_->getBlock(1)->norm2(norm_c);
        
        if(this->verbose_)
            std::cout << "2-Norm of residual of concentration: " << norm_c[0] << std::endl;

    }


}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateMeshDisplacement() const
{

     *meshDisplacementOld_rep_ = *meshDisplacementNew_rep_;
}
// Muss derzeit nur am Anfang jeder Zeititeration aufgerufen werden, damit
// problemTimeFluid_ und problemTimeStructure_ die aktuelle Loesung haben.
// ACHTUNG: Wenn wir irgendwann einmal anfangen reAssemble() auf problemFluid_ und
// problemStructure_ aufzurufen, dann muessen wir in jeder nichtlinearen Iteration
// diese setPartialSolutions() aufrufen, damit problemFluid_ und problemStructure_
// den korrekten nichtlinearen Term ausrechnen koennen.
// CH: Ist das noch relevant?
// We need to build FSI so this method is not needed anymore
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setFromPartialVectorsInit() const
{
    
    //Chem 
    if(!chemistryExplicit_){
        this->solution_->addBlock( this->problemChem_->getSolution()->getBlockNonConst(0), 1);
        this->rhs_->addBlock( this->problemChem_->getRhs()->getBlockNonConst(0), 1 );
        this->sourceTerm_->addBlock( this->problemChem_->getSourceTerm()->getBlockNonConst(0), 1 );
    }

    if (materialModel_=="SCI_Linear"){
        this->solution_->addBlock( this->problemStructure_->getSolution()->getBlockNonConst(0), 0);
        // we dont have a residual vector for linear problems
        this->rhs_->addBlock( this->problemStructure_->getRhs()->getBlockNonConst(0), 0 );
        this->sourceTerm_->addBlock( this->problemStructure_->getSourceTerm()->getBlockNonConst(0), 0 );
    }
    else{
        this->solution_->addBlock( this->problemStructureNonLin_->getSolution()->getBlockNonConst(0), 0);
        this->residualVec_->addBlock( this->problemStructureNonLin_->getResidualVector()->getBlockNonConst(0), 0 );
        this->rhs_->addBlock( this->problemStructureNonLin_->getRhs()->getBlockNonConst(0), 0 );
        this->previousSolution_->addBlock( this->problemStructureNonLin_->getPreviousSolution()->getBlockNonConst(0), 0 );
        this->sourceTerm_->addBlock( this->problemStructureNonLin_->getSourceTerm()->getBlockNonConst(0), 0 );
    }
      
}

// We need to add 'External' as an option

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setupSubTimeProblems(ParameterListPtr_Type parameterListChem, ParameterListPtr_Type parameterListStructure) const
{
    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblems " << std::endl;

    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();

    int sizeChem = this->problemChem_->getSystem()->size();
    int sizeStructure;
    if (materialModel_=="SCI_Linear")
        sizeStructure = this->problemStructure_->getSystem()->size();
    else
        sizeStructure = this->problemStructureNonLin_->getSystem()->size();
    
    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblem for Chem ";

    problemTimeChem_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemChem_, this->comm_));
        
    if(this->verbose_)
        std::cout << "-- done " << std::endl;

    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblem for Elasticity " ;

    if (materialModel_=="SCI_Linear")
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructure_, this->comm_));
    else
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructureNonLin_, this->comm_));

    if(this->verbose_)
        std::cout << "-- done " << std::endl;

    // ######################
    // Chem: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffChem(sizeChem);
    SmallMatrix<double> problemCoeffChem(sizeChem);
    SmallMatrix<int> defChem(sizeChem);

    double coeffSourceTermChem = 0.0;
    
    if ( this->getParameterList()->sublist("Timestepping Parameter").get("Class","Multistep") == "Multistep" ) {
       /* for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++) {
                if ((*defTS_)[i+sizeStructure][j+sizeStructure]==1 && i==j) {
                    defChem[i][j] = 1;
                    massCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffChem[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++){
                if ((*defTS_)[i+sizeStructure][j+sizeStructure]==1){
                    problemCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(1);
                    coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffChem[i][j] = 1.;
                }
            }
        }*/
        defChem[0][0] = 1;
        massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt;
        problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1);
        coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1);

        this->problemTimeChem_->setTimeDef(defChem);
        this->problemTimeChem_->setTimeParameters(massCoeffChem,problemCoeffChem);
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Implement other SCI Chem time stepping than BDF.");
    }
    // ######################
    // Struktur: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems

    SmallMatrix<double> massCoeffStructure(sizeStructure);
    SmallMatrix<double> problemCoeffStructure(sizeStructure);
    SmallMatrix<int> defStructure(sizeStructure);
    double coeffSourceTermStructure = 0.0; // Koeffizient fuer den Source-Term (= rechte Seite der DGL); mit Null initialisieren

    // Koeffizient vor der Massematrix
    for(int i = 0; i < sizeStructure; i++)
    {
        for(int j = 0; j < sizeStructure; j++)
        {
            // Falls in dem Block von timeStepDef_ zeitintegriert werden soll.
            // i == j, da vektorwertige Massematrix blockdiagonal ist
            if((*defTS_)[i][j] == 1 && i == j) // Weil: (d_s,c) und timeStepDef_ von SCI
            {
                defStructure[i][j] = 1;
            // Vorfaktor der Massematrix in der LHS
                massCoeffStructure[i][j] = 1.0/(dt*dt*beta);
            }
            else
            {
                massCoeffStructure[i][j] = 0.;
            }
        }
    }

    // Die anderen beiden Koeffizienten
    for(int i = 0; i < sizeStructure; i++)
    {
        for(int j = 0; j < sizeStructure; j++)
        {
            if((*defTS_)[i ][j] == 1)
            {
                problemCoeffStructure[i][j] =  1.0;
                // Der Source Term ist schon nach der Assemblierung mit der Dichte \rho skaliert worden
                coeffSourceTermStructure = 1.0; // ACHTUNG FUER SOURCE TERM, DER NICHT IN DER ZEIT DISKRETISIERT WIRD!
            }
            else // Die steady-Systemmatrix ist nicht zwingend blockdiagonal
            {
                problemCoeffStructure[i][j] = 1.0;
            }
        }
    }
    this->problemTimeStructure_->setTimeDef(defStructure);
    this->problemTimeStructure_->setTimeParameters(massCoeffStructure,problemCoeffStructure);

    this->problemTimeChem_->assemble( "MassSystem" );
    this->problemTimeStructure_->assemble( "MassSystem" );
  
    if(this->verbose_)
        std::cout << "done -- \n" << std::endl;

    // TIMER ERSTELLEN AN DER

}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setChemMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix fuer SCI combineSystems(), ggf nichtlinear. Mass Matrix same as for Chem.
    //######################
    int size = this->problemTimeChem_->getSystem()->size();

    this->problemTimeChem_->systemMass_.reset(new BlockMatrix_Type(size));
    {
        BlockMatrixPtr_Type massSystem = Teuchos::rcp(new BlockMatrix_Type(1));

        massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeChem_->getDomain(0)->getMapUnique(), this->getDomain(1)->getApproxEntriesPerRow() ) );
        // 1 = Chem
        this->feFactory_->assemblyMass( this->dim_, this->problemTimeChem_->getFEType(0), "Scalar",  massmatrix, 1, true );

        /*if(chemistryExplicit_){
            massmatrix->resumeFill();
            massmatrix->scale(0.); 
            massmatrix->fillComplete( this->problemTimeChem_->getDomain(0)->getMapUnique(), this->problemTimeChem_->getDomain(0)->getMapUnique() );
        }*/
        this->problemTimeChem_->systemMass_->addBlock(massmatrix, 0, 0);
    }
}

// For now: Leave it like that. ProblemChem. We use BDF2 for the chemistry as well
// TODO: updateMultistepRhsFSI() einbauen!
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::computeChemRHSInTime( ) const
{
    //######################
    // RHS nach BDF2
    //######################
    int sizeChem = this->problemChem_->getSystem()->size();
    int sizeStructure =1; 
    
    double dt = timeSteppingTool_->get_dt();
    double dt_prev = timeSteppingTool_->get_dt_prev();
    int nmbBDF = timeSteppingTool_->getBDFNumber();

    vec_dbl_Type coeffPrevSteps(nmbBDF);
    for(int i = 0; i < coeffPrevSteps.size(); i++)
    {
        coeffPrevSteps.at(i) = timeSteppingTool_->getInformationBDF(i+2) / dt_prev;
    }

    if (timeSteppingTool_->currentTime()==0.) {
        SmallMatrix<double> tmpmassCoeff(sizeChem);
        SmallMatrix<double> tmpproblemCoeff(sizeChem);
        tmpmassCoeff[0][0] = 1. / dt;
        tmpproblemCoeff[0][0] =  1.; // ist das richtig? Vermutlich schon, da BDF so geschrieben ist, dass zu berechnende Lsg den Koeffizienten 1 hat

        /*for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++) {
                if ((*defTS_)[i+sizeStructure][j+sizeStructure]==1 && i==j) {
                    tmpmassCoeff[i][j] = 1. / dt;
                }
                else{
                    tmpmassCoeff[i][j] = 0.;
                }
            }
        }
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++){
                if ((*defTS_)[i+sizeStructure][j+sizeStructure]==1){
                    tmpproblemCoeff[i][j] =  1.; // ist das richtig? Vermutlich schon, da BDF so geschrieben ist, dass zu berechnende Lsg den Koeffizienten 1 hat
                }
                else{
                    tmpproblemCoeff[i][j] = 1.;
                }
            }
        }*/
        this->problemTimeChem_->setTimeParameters(tmpmassCoeff, tmpproblemCoeff);
    }
    if (timeSteppingTool_->currentTime()==0.) {
        vec_dbl_Type tmpcoeffPrevSteps(1, 1. / dt);
        this->problemTimeChem_->updateMultistepRhsFSI(coeffPrevSteps,1);/*apply (mass matrix_t / dt) to u_t and more*/
    }
    else{
        this->problemTimeChem_->updateMultistepRhsFSI(coeffPrevSteps,nmbBDF);/*apply (mass matrix_t / dt) to u_t and more*/
    }

    // TODO
    if (this->problemTimeChem_->hasSourceTerm()) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Check sourceterm.");
    }

    // Wieder zu den eigentlichen Parametern zuruecksetzen, nachdem die temporaeren
    // genommen wurden.
    if (timeSteppingTool_->currentTime()==0.) {
        SmallMatrix<double> massCoeffChem(sizeChem);
        SmallMatrix<double> problemCoeffChem(sizeChem);
        massCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(0) / dt;
        problemCoeffChem[0][0] = timeSteppingTool_->getInformationBDF(1);

        /*for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++) {
                if ((*defTS_)[i+sizeStructure][j+sizeStructure]==1 && i==j) {
                    massCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffChem[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++){
                if ((*defTS_)[i+sizeStructure][j+sizeStructure]==1){
                    problemCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffChem[i][j] = 1.;
                }
            }
        }*/

        this->problemTimeChem_->setTimeParameters(massCoeffChem, problemCoeffChem);
    }
}


// This is equivalent to the FSI Structure part.
// We add one additionale feature where we can distinguish between load or timestep
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::computeSolidRHSInTime() const { // TODO: Rename because it computes only the forcing term and not the internal forces
    //######################
    // RHS nach Newmark
    //######################
    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();

    // Temporaerer Koeffizienten fuer die Skalierung der Massematrix in der rechten Seite des Systems in UpdateNewmarkRhs()
    vec_dbl_Type coeffTemp(1);
    coeffTemp.at(0) = 1.0;
    
    if( loadStepping_ == false){
        // Update u und berechne u' und u'' mit Hilfe der Newmark-Vorschrift
        this->problemTimeStructure_->updateSolutionNewmarkPreviousStep(dt, beta, gamma);
        
        // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
        // Bei Newmark lautet dies:
        // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
        // wobei u' = v (velocity) und u'' = w (acceleration).
        this->problemTimeStructure_->updateNewmarkRhs(dt, beta, gamma, coeffTemp);
    }
    else{
        // The if condition resets the rhs. If we skip it when we attemp loadstepping, the rhs will be updated continously and wrongly increase with each timestep
        this->problemTimeStructure_->getRhs()->scale(0.0);
    }

    //can we get rid of this?
   
    // TODO: SourceTerm wird in jedem Zeitschritt neu berechnet; auch wenn konstant!!!
    // if(time == 0){nur dann konstanten SourceTerm berechnen}

    if (this->problemTimeStructure_->hasSourceTerm())
    {
        if(externalForce_){

            MultiVectorPtr_Type FERhs = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ));
            /*vec_dbl_Type funcParameter(6,0.);
            funcParameter[0] = timeSteppingTool_->t_;            
            // how can we use different parameters for different blocks here?
            funcParameter[1] =this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Volume force",0.00211);
            funcParameter[2]= this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Load Step Size",1.);
            funcParameter[3] = this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Load Ramp End",1.);
    
            funcParameter[4] =this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Heart Beat Start",70.);
            funcParameter[5] = 0.;*/
            vec_dbl_Type funcParameter(1,0.);
            funcParameter[0] = timeSteppingTool_->t_;            
            // how can we use different parameters for different blocks here?
            for (int j = 0; j < this->problemTimeStructure_->getUnderlyingProblem()->getParameterCount(); j++)
                funcParameter.push_back(this->problemTimeStructure_->getUnderlyingProblem()->getParameterRhs(j));
            funcParameter.push_back(0.);
            
            
            if(nonlinearExternalForce_){

                MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
                d_rep_->importFromVector(d, true); 
                MatrixPtr_Type A( new Matrix_Type (this->system_->getBlock(0,0)));
                MatrixPtr_Type AKext(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
                MatrixPtr_Type Kext(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow()*2 ) );
                MultiVectorPtr_Type Kext_vec;
                this->feFactory_->assemblyNonlinearSurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, d_rep_,Kext, funcParameter, this->problemTimeStructure_->getUnderlyingProblem()->rhsFuncVec_[0],this->parameterList_);
                A->addMatrix(1.,AKext,0.);
                Kext->addMatrix(1.,AKext,1.);

                AKext->fillComplete(this->getDomain(0)->getMapVecFieldUnique(),this->getDomain(0)->getMapVecFieldUnique());
                this->system_->addBlock(AKext,0,0);

            }
            else            
                this->feFactory_->assemblySurfaceIntegralExternal(this->dim_, this->getDomain(0)->getFEType(),FERhs, d_rep_, funcParameter, this->problemTimeStructure_->getUnderlyingProblem()->rhsFuncVec_[0],this->parameterList_);


            this->sourceTerm_->getBlockNonConst(0)->exportFromVector( FERhs, false, "Add" );

           

        }
        else
            this->problemTimeStructure_->assembleSourceTerm( timeSteppingTool_->t_ );
        //this->problemTimeStructure_->getSourceTerm()->scale(density);
        // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
        // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!

        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
        BlockMultiVectorPtr_Type tmpSourceterm = Teuchos::rcp(new BlockMultiVector_Type(1)) ;
        tmpSourceterm->addBlock(this->sourceTerm_->getBlockNonConst(0),0);

        BlockMultiVectorPtr_Type tmpPtr ; 
        if(externalForce_)
            tmpPtr =  tmpSourceterm;
        else
            tmpPtr = this->problemTimeStructure_->getSourceTerm();
        


        this->problemTimeStructure_->getRhs()->update(coeffSourceTermStructure, *tmpPtr, 1.);
        this->rhs_->addBlock( this->problemTimeStructure_->getRhs()->getBlockNonConst(0), 0 );


    }

}

// Can stay the same
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setSolidMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix
    //######################
    double density = this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Density",1.0);

    int size = this->problemTimeStructure_->getSystem()->size();

    // Built once, in the first time step of the run (after a restart, the
    // first time step after the restart time).
    if(!solidMassBuilt_)
    {
        solidMassBuilt_ = true;
        this->problemTimeStructure_->systemMass_.reset(new BlockMatrix_Type(size));
        {

            massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // 0 = Struktur
            this->feFactory_->assemblyMass(this->dim_, this->problemTimeStructure_->getFEType(0), "Vector", massmatrix, 0, true);
            massmatrix->resumeFill();
            massmatrix->scale(density);
            if(loadStepping_ == true)
                massmatrix->scale(0.0);

            massmatrix->fillComplete( this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique());

            this->problemTimeStructure_->systemMass_->addBlock( massmatrix, 0, 0 );
        }
    }
}

// Can stay the same
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setBoundariesSubProblems( ) const
{
    //######################
    // Boundaries 
    //######################
  
    this->problemTimeStructure_->setBoundaries();
    this->problemTimeChem_->setBoundaries();

        
}


// Damit die richtige timeSteppingTool_->currentTime() genommen wird.
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::saveStepState() const
{
    TEUCHOS_TEST_FOR_EXCEPTION( chemistryExplicit_, std::logic_error, "Adaptive time stepping needs the implicit coupling (Chemistry Explicit false).");
    savedTimeSteppingTool_ = Teuchos::rcp( new TimeSteppingTools( *timeSteppingTool_ ) );
    if (!problemTimeStructure_.is_null())
        problemTimeStructure_->saveState();
    this->feFactory_->saveStateAssemblyFEElements();
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::restoreStepState() const
{
    TEUCHOS_TEST_FOR_EXCEPTION( savedTimeSteppingTool_.is_null(), std::logic_error, "SCI::restoreStepState() without saveStepState().");
    *timeSteppingTool_ = *savedTimeSteppingTool_;
    if (!problemTimeStructure_.is_null())
        problemTimeStructure_->restoreState();
    this->feFactory_->restoreStateAssemblyFEElements();
    elementFailed_ = false;
    failedElements_ = 0;
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setAdaptiveStep(bool adaptive, bool acceptElementFailures) const
{
    adaptiveStep_ = adaptive;
    acceptElementFailures_ = adaptive && acceptElementFailures;
    elementFailed_ = false;
    failedElements_ = 0;
    this->feFactory_->setRecordFailureAssemblyFEElements(adaptive);
    this->feFactory_->resetFailureAssemblyFEElements();
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateTime() const
{
    timeSteppingTool_->advanceTime(); // Now SCI time stepper has t_{n+1} and corresponding dt(already applied)

    MultiVectorConstPtr_Type c; 
    if(chemistryExplicit_)
        c= this->problemTimeChem_->getSolution()->getBlock(0);
    else
        c = this->solution_->getBlock(1);
    
    MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
    c_rep_->importFromVector(c, true);
    d_rep_->importFromVector(d, true); 
    this->feFactory_->advanceInTimeAssemblyFEElements(timeSteppingTool_, d_rep_, c_rep_ ); // Syncs time t_{n+1} and dt to all AceGen elements

    this->problemTimeChem_->updateTime(timeSteppingTool_->t_); // ProblemChem has t_{n+1}
    this->problemTimeStructure_->updateTime(timeSteppingTool_->t_); // ProblemStructure has t_{n+1}

}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                     const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                    ) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "implement NOX for steady SCI.");
    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
}   
    
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::moveMesh() const
{

    MultiVectorConstPtr_Type displacementUniqueConst;

    displacementUniqueConst = this->solution_->getBlock(0);
    MultiVectorPtr_Type displacementRepeated = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );

    displacementRepeated->importFromVector( displacementUniqueConst );
    MultiVectorPtr_Type displacementUnique = Teuchos::rcp_const_cast<MultiVector_Type>(displacementUniqueConst);


    // Verschiebe die Gitter fuer Chemistry
    // ACHTUNG: Klappt nur, weil die P2-Knoten hinter den P1-Knoten kommen.
    // Sonst muessen fuer den Druck die P1-Knoten extrahiert werden.
    // TODO: Wahrscheinlich reicht nur FSI-Domain, da selbes Objekt bei problemChem_ und problemTimeChem_.
   // ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
   // ( Teuchos::rcp_const_cast<Domain_Type>(this->problemChem_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeChem_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);

}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::initializeCE(){
/*All vectors (solution, rhs, previousSolution,...) should be initialized at this point (initializeProblem() was called)
 Therefore, all these BlockMVectors should have a length of 5 since the geometry problem is included in the general setup. We need to resize these BlockMVs here if the Geometry Explicit (GE) system is used.
*/
    if (chemistryExplicit_) {
        this->solution_->resize( 1 );
        this->rhs_->resize( 1 );
        this->sourceTerm_->resize( 1 );
        this->rhsFuncVec_.resize( 1 );
        this->previousSolution_->resize( 1 );
        this->residualVec_->resize( 1 );
        this->initVectorSpaces();  //reinitialize NOX vector spaces
    }
}


// Am Anfang der Zeititeration erst updateSolutionMultiPreviousStep() aufrufen und dann erst updateMultistepRhs(),
// damit die previousSolution_ initialisiert sind. Genauso fuer SystemMass
// TODO: updateSystemMassMultiPreviousStep() fertig programmieren
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateChemInTime() const
{
    int nmbBDF = 1; //timeSteppingTool_->getBDFNumber();

    if(nmbBDF<2 && !this->parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation")) {
        if (timeSteppingTool_->currentTime()!=0.){
            this->problemTimeChem_->updateSolutionMultiPreviousStep(2);
            this->problemTimeChem_->updateSystemMassMultiPreviousStep(2);
        }
        else{
            this->problemTimeChem_->updateSolutionMultiPreviousStep(1);
            this->problemTimeChem_->updateSystemMassMultiPreviousStep(1);
        }
    }
    else{
        this->problemTimeChem_->updateSolutionMultiPreviousStep(nmbBDF);
        this->problemTimeChem_->updateSystemMassMultiPreviousStep(nmbBDF);
    }
}
template<class SC,class LO,class GO,class NO>
typename SCI<SC,LO,GO,NO>::BlockMultiVectorPtr_Type SCI<SC,LO,GO,NO>::getPostProcessingData()
{
        
    // /*
    // 0 -- "Volume","
    // 1 -- Sxx",
    // 2 -- "Sxy"
    // 3 -- "Sxz" 
    // 4 -- "Syx"
    // 5 -- "Syy"
    // 6 -- "Syz" 
    // 7 -- "Szx" 
    // 8 -- "Szy" 
    // 9 -- "Szz" 
    // 10 -- "MisesStress" 
    // 11 -- "SCirc"
    // 12 -- "SAxial",
    // 13 -- "SRadial"
    // 14 -- "Exx"
    // 15 -- "Exy"
    // 16 -- "Exz" 
    // 17 -- "Eyx"
    // 18 -- "Eyy"
    // 19 -- "Eyz" 
    // 20 -- "Ezx"
    // 21 -- "Ezy"
    // 22 -- "Ezz"
    // 23 -- "W"
    // 24 -- "Growth1"
    // 25 -- "Growth2" 
    // 26 -- "Growth3"
    // 27 -- "Stretch1"
    // 28 -- "Stretch2"
    // 29 -- "DetF",
    // 30 - 38      "Ag1n1","Ag1n2","Ag1n3","Ag2n1","Ag2n2","Ag2n3",
    //                    "Ag3n1","Ag3n2","Ag3n3"
    // 39 -- "a11"
    // 40 -- "a12"
    // 41 -- "a13"
    // 42 -- "a21"
    // 43 -- "a22"
    // 44 -- "a23"
    // 45 -- "nC1" <----- !!
    // 46 -- "nC2" <----- !!
    // 47 -- "nD1" <----- !! 
    // 48 -- "nD2" <----- !!
    // 49 -- "ScDir1"
    // 50 -- "ScDir2" 
    // 51 -- "ScDir3" 
    // 52 -- "SaDir1"
    // 53 -- "SaDir2"
    // 54 -- "SaDir3"
    // 55 -- "SrDir1"
    // 56 -- "SrDir2"
    // 57 -- "SrDir3"*/

    // Initialize the post-processing names if not already done

    if (postProcessingnames_.empty()) {
        if(this->verbose_)
            std::cout << "Initializing post-processing names..." << std::endl;
        postProcessingnames_ = this->feFactory_->getPostDataNames();
    }

    // Read the desired post-processing fields from the parameter list
    Teuchos::Array<std::string> requestedField;
    if (this->parameterList_->sublist("Parameter").isParameter("Post Processing Fields")) {
        requestedField = this->parameterList_->sublist("Parameter").get("Post Processing Fields", 
                                                                        Teuchos::Array<std::string>());
    } else {
        // Fallback to all available fields if not specified
        requestedField.resize(postProcessingnames_.size());
        for (int i = 0; i < postProcessingnames_.size(); i++) {
            requestedField[i] = postProcessingnames_[i];
        }
    }


    // Create BlockMultiVector with the number of requested fields
    BlockMultiVectorPtr_Type postProcess = Teuchos::rcp(new BlockMultiVector_Type(requestedField.size()));

    // Populate the BlockMultiVector with the requested fields
    for(int i = 0; i < requestedField.size(); i++) {
        MultiVectorPtr_Type fieldData = Teuchos::rcp(new MultiVector_Type(this->getDomain(0)->getMapUnique()));
        this->feFactory_->postProcessing(requestedField[i], fieldData);
        postProcess->addBlock(fieldData, i);
    }

    // This runs after the above block because it ensures that the requested field is available in the postprocessing data
    if(requestedFields_.empty()) {
        requestedFields_ = vec_string_Type(requestedField.begin(), requestedField.end());
    }

    return postProcess;
}

template<class SC,class LO,class GO,class NO>
vec_string_Type SCI<SC,LO,GO,NO>::getPostprocessingNames()
{
    return requestedFields_;
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::getValuesOfInterest(BlockMultiVectorPtr_Type& historyMultiVector)
{
    historyMultiVector = this->feFactory_->getHistoryValues();

}

/// Names of the element history values of one Gauss point in the History checkpoint
/// files: those of the SCI_SMC_Active_Growth_Reorientation element (34 values;
/// checkpoints written before the history layout was taken from the element use
/// them), or H0, H1, ... for an element with a different history length.
inline vec_string_Type sciHistoryVariableNames(int valuesPerGaussPoint)
{
    vec_string_Type names = {"LambdaBarC1", "LambdaBarC2", "nA1", "nA2", "nB1", "nB2", "nC1",
                             "nC2", "nD1", "nD2", "LambdaA1", "LambdaA2", "k251", "k252", "LambdaBarP1",
                             "LambdaBarP2", "Theta1", "Theta2", "Theta3", "Ag11", "Ag12", "Ag13", "Ag21",
                             "Ag22", "Ag23", "Ag31", "Ag32", "Ag33", "a11", "a12", "a13", "a21", "a22", "a23"};
    if (valuesPerGaussPoint == static_cast<int>(names.size()))
        return names;

    names.clear();
    for (int k = 0; k < valuesPerGaussPoint; k++)
        names.push_back("H" + std::to_string(k));
    return names;
}

// Element history of the checkpoint at 'time', in History<time> in the checkpoint directory
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::exportValuesOfInterest(double time)
{
    bool exportHistory = this->parameterList_->sublist("Timestepping Parameter").get("Export history", false);

    if(exportHistory)
    {
        // One block per Gauss point, one vector per history value of a Gauss point
        BlockMultiVectorPtr_Type historyValues;
        this->getValuesOfInterest(historyValues);
        int numGaussPoints = historyValues->size();
        vec_string_Type historyNames = sciHistoryVariableNames(numGaussPoints > 0 ? historyValues->getBlock(0)->getNumVectors() : 0);

        Teuchos::RCP<HDF5Export<SC,LO,GO,NO>> exporter =Teuchos::rcp(new HDF5Export<SC,LO,GO,NO>(this->getDomain(0)->getElementMap(),checkpointFile(this->parameterList_, "History"+std::to_string(time))));

        for(int gp =0; gp<numGaussPoints; gp++){
            for(int k=0; k < historyNames.size(); k++){
                std::string varName = historyNames[k]+"_"+std::to_string(gp);
                exporter->writeVariablesHDF5(varName,historyValues->getBlock(gp)->getVector(k));
            }
        }
    }
}

// Element history and element solution of the restart time 'time'
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::importValuesOfInterest(double time)
{

    bool importHistory = this->parameterList_->sublist("Timestepping Parameter").get("Import history", false);

    if(importHistory)
    {
        // The layout (Gauss points, values per Gauss point) is the elements'
        BlockMultiVectorPtr_Type historyValues;
        this->getValuesOfInterest(historyValues);
        int numGaussPoints = historyValues->size();
        int valuesPerGaussPoint = numGaussPoints > 0 ? historyValues->getBlock(0)->getNumVectors() : 0;
        vec_string_Type historyNames = sciHistoryVariableNames(valuesPerGaussPoint);

        MapConstPtr_Type elementMap = this->getDomain(0)->getElementMap();
        Teuchos::RCP<HDF5Import<SC,LO,GO,NO>> importer =Teuchos::rcp(new HDF5Import<SC,LO,GO,NO>(elementMap,restartFile(this->parameterList_, "History"+std::to_string(time))));

        vec2D_dbl_Type myHistory (elementMap->getNodeNumElements(),vec_dbl_Type(valuesPerGaussPoint*numGaussPoints,-1.));

        for(int gp =0; gp<numGaussPoints; gp++){
            for(int k=0; k < valuesPerGaussPoint; k++){
                std::string varName = historyNames[k]+"_"+std::to_string(gp);
                MultiVectorConstPtr_Type history  = importer->readVariablesHDF5(varName);
                Teuchos::ArrayRCP<const SC>  historyArray = history->getData(0);
                for(int T =0;T<historyArray.size(); T++){
                   myHistory[T][k+gp*valuesPerGaussPoint] = historyArray[T];
                }
            }
        }

        for(int T=0; T< elementMap->getNodeNumElements() ; T++ ){
            this->feFactory_->setHistoryValues(T,myHistory[T]);
        }

    }

    // The elements' solution: the restarted solution
    MultiVectorConstPtr_Type c;
    if(chemistryExplicit_)
        c= this->problemTimeChem_->getSolution()->getBlock(0);
    else
        c = this->solution_->getBlock(1);

    MultiVectorConstPtr_Type d = this->solution_->getBlock(0);
    d_rep_->importFromVector(d, true);
    c_rep_->importFromVector(c, true);
    this->feFactory_->updateSolutionAssemblyFEElements(d_rep_, c_rep_ );
}

}
#endif

