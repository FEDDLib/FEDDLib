#ifndef SCI_def_hpp
#define SCI_def_hpp
#include "SCI_decl.hpp"


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
//problemStructureNonLin_(),
meshDisplacementOld_rep_(),
meshDisplacementNew_rep_(),
c_rep_(),
defTS_(defTS),
timeSteppingTool_(),
exporterEMod_(),
materialModel_( parameterListStructure->sublist("Parameter").get("Material model","linear") )
{
    //this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);

    //this->initNOXParameters();
    counterP = 0;
    
    //std::string linearization = parameterListSCI->sublist("General").get("Linearization","FixedPoint");
    
    //TEUCHOS_TEST_FOR_EXCEPTION( !(linearization == "Newton" || linearization == "NOX")  && materialModel_ != "linear", std::runtime_error, "Nonlinear material models can only be used with Newton's method or FixedPoint (nonlinear material Jacobian will still be used).");
    
    this->addVariable( domainChem, FETypeChem, "c", 1); // Chemistry scalar valued problem
    this->addVariable( domainStructure, FETypeStructure, "d_s", domainStructure->getDimension() ); // Structure

    this->dim_ = this->getDomain(0)->getDimension();
    
  
    
    if (materialModel_=="linear"){
        problemStructure_ = Teuchos::rcp( new StructureProblem_Type( domainStructure, FETypeStructure, parameterListStructure ) );
        problemStructure_->initializeProblem();
    }
    else{
        problemStructureNonLin_ = Teuchos::rcp( new StructureNonLinProblem_Type( domainStructure, FETypeStructure, parameterListStructure) );
        problemStructureNonLin_->initializeProblem();
    }
    
    problemChem_ = Teuchos::rcp( new ChemProblem_Type( domainChem, FETypeChem, parameterListChem, diffusionTensor, reactionFunc ) );
    problemChem_->initializeProblem();

    //We initialize the subproblems. In the main routine, we need to call initializeFSI(). There, we first initialize the vectors of the FSI problem and then we set the pointers of the subproblems to the vectors of the full monolithic FSI system. This way all values are only saved once in the subproblems and can be used by the monolithic FSI system.
    
    meshDisplacementNew_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    meshDisplacementOld_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    
    c_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    //w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    //u_minus_w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    exportedEMod_ = false;
    setUpTimeStep_=false;
    eModVec_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getElementMap() ) );

}


template<class SC,class LO,class GO,class NO>
SCI<SC,LO,GO,NO>::~SCI()
{
    if (!exporterEMod_.is_null()) {
       exporterEMod_->closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::info()
{
    this->infoProblem();
    //this->infoNonlinProblem();
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
            std::cout << "-- Assembly SCI ... " << std::endl;
        }

        // First Assumption: Almost incompressible Material
        this->problemChem_->assemble();
        
        // Elementwise determined E Module
        MultiVectorPtr_Type solChemRep = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapRepeated() ) );
        solChemRep->importFromVector(this->problemChem_->getSolution()->getBlock(0));

        int dim = this->getDomain(0)->getDimension();

        this->feFactory_->determineEMod(this->getDomain(1)->getFEType(),solChemRep,eModVec_,this->getDomain(1),this->parameterList_);
                
        double nu = this->parameterList_->sublist("Parameter").get("PoissonRatio",0.4);

        MatrixPtr_Type A(new Matrix_Type( this->getDomain(1)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) ); // Structure-Matrix

        // steady rhs wird hier assembliert.
        // rhsFunc auf 0 (=x) und 0 (=y) abaendern bei LinElas!
        if (materialModel_=="linear"){
            this->feFactory_->assemblyLinElasXDimE(this->dim_, this->getDomain(0)->getFEType(), A, eModVec_, nu, true);
            this->problemStructure_->system_->addBlock(A,0,0);// assemble(); //

            //this->problemStructure_->assemble();
        }
        else{  
            MultiVectorConstPtr_Type eModVecConst = eModVec_;
            this->problemStructureNonLin_->updateEMod(eModVecConst);        
            this->problemStructureNonLin_->assemble(); //system_->addBlock(A,0,0);// assemble(); //                               
        }

       
        //MatrixPtr_Type C1_T(new Matrix_Type( this->getDomain(1)->getMapVecFieldUnique(), 1 ) ); // Chem-Kopplung
        //MatrixPtr_Type C1(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) ); // Struktur-Kopplung
       
        // ###########################
        // Korrekte Skalierung der entsprechenden Bloecke
        // ###########################
        double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);
        
        // ###########################
        // Bloecke hinzufuegen
        // ###########################
        
        this->system_.reset(new BlockMatrix_Type(2));

         
        // Struktur
        if (materialModel_=="linear"){
            this->system_->addBlock( this->problemStructure_->system_->getBlock(0,0), 1, 1);
        }
        else{       
            this->system_->addBlock( this->problemStructureNonLin_->system_->getBlock(0,0), 1, 1);
        }    

        // Chemistry
        this->system_->addBlock( this->problemChem_->system_->getBlock(0,0), 0, 0 );
       

        // Fuer die Zeitprobleme
        timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_)); 
        ParameterListPtr_Type plStructure;
        if (materialModel_=="linear")
            plStructure = this->problemStructure_->getParameterList();
        else
            plStructure = this->problemStructureNonLin_->getParameterList();

        setupSubTimeProblems(this->problemChem_->getParameterList(), plStructure);
        // We set the vector from the partial problems
        this->setFromPartialVectorsInit();
        
        
        if ( this->parameterList_->sublist("Exporter").get("Export EMod", true)){
            if(exportedEMod_ == false){
                exporterEMod_ = Teuchos::rcp(new Exporter_Type());
                
                DomainConstPtr_Type dom = this->getDomain(1);
                MultiVectorConstPtr_Type exportVector = eModVec_;

                MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( dom->getMesh() );
                exporterEMod_->setup("EModuleValues", meshNonConst, "P0",  this->parameterList_);
                
                
                exporterEMod_->addVariable( exportVector, "EModule", "Scalar", 1, dom->getElementMap() );
                exporterEMod_->save( 0.0);
                exportedEMod_ = true;
            }
            else {
                exporterEMod_->save( timeSteppingTool_->t_);
            }
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
void SCI<SC,LO,GO,NO>::reAssemble(std::string type) const
{

    double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);

    if (type == "MoveMesh")
         moveMesh();
    else if (type == "UpdateMeshDisplacement")
         updateMeshDisplacement();
    else if (type == "ComputeSolidRHSInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (ComputeSolidRHSInTime)" << '\n';
          
         computeSolidRHSInTime();
    }
    else if (type == "ComputeChemRHSInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (ComputeChemRHSInTime)" << '\n';
          
         computeChemRHSInTime();
    }
    else if (type == "UpdateChemInTime"){
        if(this->verbose_)
            std::cout << "-- Assembly (UpdateChemInTime)" << '\n';
          updateChemInTime();}

    else if(type == "UpdateTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateTime)" << '\n';

        updateTime();
        return;
    }
    else if(type == "SetBoundaries")
    {
        if(this->verbose_)
            std::cout << "-- set Boundaries" << '\n';

        setBoundariesSubProblems();
        return;
    }
    else if(type == "UpdateEMod")
    {
        if(this->verbose_)
            std::cout << "-- set Boundaries" << '\n';

        MultiVectorPtr_Type solChemRep = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapRepeated() ) );
        solChemRep->importFromVector(this->problemTimeChem_->getSolution()->getBlock(0));

        int dim = this->getDomain(0)->getDimension();

        this->feFactory_->determineEMod(this->getDomain(1)->getFEType(),solChemRep,eModVec_,this->getDomain(1),this->parameterList_);
                
        double nu = this->parameterList_->sublist("Parameter").get("PoissonRatio",0.4);

        MatrixPtr_Type A(new Matrix_Type( this->getDomain(1)->getMapVecFieldUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) ); // Structure-Matrix

        // steady rhs wird hier assembliert.
        // rhsFunc auf 0 (=x) und 0 (=y) abaendern bei LinElas!
        //if (materialModel_=="linear"){
        if (materialModel_=="linear"){
            this->feFactory_->assemblyLinElasXDimE(this->dim_, this->getDomain(0)->getFEType(), A, eModVec_, nu, true);
            this->problemStructure_->system_->addBlock(A,0,0);// assemble(); //
            this->system_->addBlock( this->problemStructure_->system_->getBlock(0,0), 1, 1);
            //this->problemStructure_->assemble();
        }
        else{
            MultiVectorConstPtr_Type eModVecConst = eModVec_;
            this->problemStructureNonLin_->updateEMod(eModVecConst);                
            this->system_->addBlock( this->problemStructureNonLin_->getSystem()->getBlock(0,0), 1, 1 );                                
        }
      
        exporterEMod_->save( timeSteppingTool_->t_);


        return;
    }
    

       
    if(type == "FixedPoint")
    {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "should always be called by the residual and never here.");
    }
    else if(type == "Newton")
    {

        if (materialModel_ != "linear"){
            this->problemStructureNonLin_->reAssemble("Newton");
        }
        
    }
        
    if (materialModel_ != "linear")
        this->system_->addBlock( this->problemStructureNonLin_->getSystem()->getBlock(0,0), 1, 1 );
}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const
{
      
    //this->problemChem_->calculateNonLinResidualVec( "reverse", time );
    //this->residualVec_->addBlock( this->problemChem_->getResidualVector()->getBlockNonConst(0) , 0);
    // we need to account for the coupling in the residuals
    if (materialModel_!="linear"){
        this->problemStructureNonLin_->calculateNonLinResidualVec( "reverse", time );
        this->residualVec_->addBlock( this->problemStructureNonLin_->getResidualVector()->getBlockNonConst(0) , 1);
        // we need to add a possible source term
    }
    else{
        MultiVectorPtr_Type residualSolidSCI =  Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(1) );
        this->problemStructure_->getSystem()->getBlock(0,0)->apply( *this->problemStructure_->getSolution()->getBlock(0), *residualSolidSCI, Teuchos::NO_TRANS, -1. ); // y= -Ax + 0*y
        MultiVectorPtr_Type resSolidNonConst = Teuchos::rcp_const_cast<MultiVector_Type> ( this->residualVec_->getBlock(1) );
        resSolidNonConst->update(1., *this->problemStructure_->getRhs()->getBlock(0), 1.);
        // we need to add a possible source term
    }
    
    MultiVectorPtr_Type residualChemSCI =  Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(0) );
    this->problemChem_->getSystem()->getBlock(0,0)->apply( *this->problemChem_->getSolution()->getBlock(0), *residualChemSCI, Teuchos::NO_TRANS, -1. ); // y= -Ax + 0*y
    MultiVectorPtr_Type resChemNonConst = Teuchos::rcp_const_cast<MultiVector_Type> ( this->residualVec_->getBlock(0) );
    resChemNonConst->update(1., *this->problemChem_->getRhs()->getBlock(0), 1.);

    // might also be called in the sub calculateNonLinResidualVec() methods which were used above
    if (type == "reverse")
        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    else if (type == "standard"){
        this->residualVec_->scale(-1.);
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
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
    this->solution_->addBlock( this->problemChem_->getSolution()->getBlockNonConst(0), 0);
    //this->residualVec_->addBlock( this->problemChem_->getResidualVector()->getBlockNonConst(0), 0 );
    this->rhs_->addBlock( this->problemChem_->getRhs()->getBlockNonConst(0), 0 );
    this->sourceTerm_->addBlock( this->problemChem_->getSourceTerm()->getBlockNonConst(0), 0 );
   
    if (materialModel_=="linear"){
        this->solution_->addBlock( this->problemStructure_->getSolution()->getBlockNonConst(0), 1 );
        // we dont have a residual vector for linear problems
        this->rhs_->addBlock( this->problemStructure_->getRhs()->getBlockNonConst(0), 1 );
        this->sourceTerm_->addBlock( this->problemStructure_->getSourceTerm()->getBlockNonConst(0), 1 );
    }
    else{
        this->solution_->addBlock( this->problemStructureNonLin_->getSolution()->getBlockNonConst(0), 1 );
        this->residualVec_->addBlock( this->problemStructureNonLin_->getResidualVector()->getBlockNonConst(0), 1 );
        this->rhs_->addBlock( this->problemStructureNonLin_->getRhs()->getBlockNonConst(0), 1 );
        this->previousSolution_->addBlock( this->problemStructureNonLin_->getPreviousSolution()->getBlockNonConst(0), 1 );
        this->sourceTerm_->addBlock( this->problemStructureNonLin_->getSourceTerm()->getBlockNonConst(0), 1 );
    }
      
}



template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setupSubTimeProblems(ParameterListPtr_Type parameterListChem, ParameterListPtr_Type parameterListStructure) const
{
    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblems \n" << std::flush;

    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();

    int sizeChem = this->problemChem_->getSystem()->size();
    int sizeStructure;
    if (materialModel_=="linear")
        sizeStructure = this->problemStructure_->getSystem()->size();
    else
        sizeStructure = this->problemStructureNonLin_->getSystem()->size();
    
    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblem for Chem \n" << std::flush;

    problemTimeChem_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemChem_, this->comm_));
        
    if(this->verbose_)
        std::cout << "-- done \n" << std::flush;

    if(this->verbose_)
        std::cout << "-- Setup SCI Sub-TimeProblem for Elasticity \n" << std::flush;

    if (materialModel_=="linear")
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructure_, this->comm_));
    else
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructureNonLin_, this->comm_));

    if(this->verbose_)
        std::cout << "-- done \n" << std::flush;

    // ######################
    // Chem: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffChem(sizeChem);
    SmallMatrix<double> problemCoeffChem(sizeChem);
    SmallMatrix<int> defChem(sizeChem);

    double coeffSourceTermChem = 0.0;
    if ( this->getParameterList()->sublist("Timestepping Parameter").get("Class","Multistep") == "Multistep" ) {
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++) {
                if ((*defTS_)[i][j]==1 && i==j) {
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
                if ((*defTS_)[i][j]==1){
                    problemCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(1);
                    coeffSourceTermChem = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffChem[i][j] = 1.;
                }
            }
        }
        this->problemTimeChem_->setTimeDef(defChem);
        this->problemTimeChem_->setTimeParameters(massCoeffChem,problemCoeffChem);
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Implement other FSI Chem time stepping than BDF.");
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
            if((*defTS_)[i + sizeChem][j + sizeChem] == 1 && i == j) // Weil: (u_f, p, d_s,...) und timeStepDef_ von FSI
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
            if((*defTS_)[i + sizeChem][j + sizeChem] == 1)
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
}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setChemMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix fuer SCI combineSystems(), ggf nichtlinear. Mass Matrix same as for Chem.
    //######################
    double density = this->problemTimeChem_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
    int size = this->problemTimeChem_->getSystem()->size();

    this->problemTimeChem_->systemMass_.reset(new BlockMatrix_Type(size));
    {
        massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeChem_->getDomain(0)->getMapUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
        // 0 = Chem
        this->feFactory_->assemblyMass( this->dim_, this->problemTimeChem_->getFEType(0), "Scalar",  massmatrix, 0, true );
        massmatrix->resumeFill();
        massmatrix->scale(density);
        massmatrix->fillComplete( this->problemTimeChem_->getDomain(0)->getMapUnique(), this->problemTimeChem_->getDomain(0)->getMapUnique() );

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
    double dt = timeSteppingTool_->get_dt();
    int nmbBDF = timeSteppingTool_->getBDFNumber();

    vec_dbl_Type coeffPrevSteps(nmbBDF);
    for(int i = 0; i < coeffPrevSteps.size(); i++)
    {
        coeffPrevSteps.at(i) = timeSteppingTool_->getInformationBDF(i+2) / dt;
    }

    if (timeSteppingTool_->currentTime()==0.) {
        SmallMatrix<double> tmpmassCoeff(sizeChem);
        SmallMatrix<double> tmpproblemCoeff(sizeChem);
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++) {
                if ((*defTS_)[i][j]==1 && i==j) {
                    tmpmassCoeff[i][j] = 1. / dt;
                }
                else{
                    tmpmassCoeff[i][j] = 0.;
                }
            }
        }
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++){
                if ((*defTS_)[i][j]==1){
                    tmpproblemCoeff[i][j] =  1.; // ist das richtig? Vermutlich schon, da BDF so geschrieben ist, dass zu berechnende Lsg den Koeffizienten 1 hat
                }
                else{
                    tmpproblemCoeff[i][j] = 1.;
                }
            }
        }
        this->problemTimeChem_->setTimeParameters(tmpmassCoeff, tmpproblemCoeff);
    }
    if (timeSteppingTool_->currentTime()==0.) {
        vec_dbl_Type tmpcoeffPrevSteps(1, 1. / dt);
        this->problemTimeChem_->updateMultistepRhs(tmpcoeffPrevSteps,1);/*apply (mass matrix_t / dt) to u_t*/
    }
    else{
        this->problemTimeChem_->updateMultistepRhs(coeffPrevSteps,nmbBDF);/*apply (mass matrix_t / dt) to u_t and more*/
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

        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++) {
                if ((*defTS_)[i][j]==1 && i==j) {
                    massCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffChem[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeChem; i++) {
            for (int j=0; j<sizeChem; j++){
                if ((*defTS_)[i][j]==1){
                    problemCoeffChem[i][j] = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffChem[i][j] = 1.;
                }
            }
        }

        this->problemTimeChem_->setTimeParameters(massCoeffChem, problemCoeffChem);
    }
}


// This is equivalent to the FSI Structure part.
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::computeSolidRHSInTime() const {
    //######################
    // RHS nach Newmark
    //######################
    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();
    double gamma = timeSteppingTool_->get_gamma();
    
    // Temporaerer Koeffizienten fuer die Skalierung der Massematrix in der rechten Seite des Systems in UpdateNewmarkRhs()
    vec_dbl_Type coeffTemp(1);
    coeffTemp.at(0) = 1.0;
    
    // Update u und berechne u' und u'' mit Hilfe der Newmark-Vorschrift
    //this->problemTimeStructure_->updateSolutionNewmarkPreviousStep(dt, beta, gamma);
    
    // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
    // Bei Newmark lautet dies:
    // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
    // wobei u' = v (velocity) und u'' = w (acceleration).
    //this->problemTimeStructure_->updateNewmarkRhs(dt, beta, gamma, coeffTemp);
    
    //can we get rid of this?
    double time = timeSteppingTool_->currentTime() + dt;
    
    // TODO: SourceTerm wird in jedem Zeitschritt neu berechnet; auch wenn konstant!!!
    // if(time == 0){nur dann konstanten SourceTerm berechnen}
    if (this->problemTimeStructure_->hasSourceTerm())
    {
        //this->problemTimeStructure_->getUnderlyingProblem()->addRhsFunction( rhsX3D );

        this->problemTimeStructure_->assembleSourceTerm( time );
        
        // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
        // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!
        
        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
        BlockMultiVectorPtr_Type tmpPtr = this->problemTimeStructure_->getSourceTerm();
        this->problemTimeStructure_->getRhs()->update(coeffSourceTermStructure, *tmpPtr, 1.);

    }

}

// Can stay the same
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setSolidMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix
    //######################
    double density = this->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
    int size = this->problemTimeStructure_->getSystem()->size();

    if(timeSteppingTool_->currentTime() == 0.0)
    {
        this->problemTimeStructure_->systemMass_.reset(new BlockMatrix_Type(size));
        {

            massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // 1 = Struktur
            //this->feFactory_->assemblyMass(this->dim_, this->problemTimeStructure_->getFEType(0), "Vector", massmatrix, 1, true);
            massmatrix->resumeFill();
            massmatrix->scale(density);
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
void SCI<SC,LO,GO,NO>::updateTime() const
{
    timeSteppingTool_->t_ = timeSteppingTool_->t_ + timeSteppingTool_->dt_prev_;
}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                     const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                    ) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "implement NOX for steady FSI.");
    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
//    if ( !type.compare("Monolithic"))
//        evalModelImplMonolithic( inArgs, outArgs );
//    else if ( !type.compare("FaCSI")){
//        evalModelImplBlock( inArgs, outArgs );
//    }
//    else
//        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner/solver type.");
}   
    
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::moveMesh() const
{

    MultiVectorConstPtr_Type displacementUniqueConst;

    displacementUniqueConst = this->solution_->getBlock(1);
    MultiVectorPtr_Type displacementRepeated = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapVecFieldRepeated() ) );

    displacementRepeated->importFromVector( displacementUniqueConst );
    MultiVectorPtr_Type displacementUnique = Teuchos::rcp_const_cast<MultiVector_Type>(displacementUniqueConst);


    // Verschiebe die Gitter fuer Chemistry
    // ACHTUNG: Klappt nur, weil die P2-Knoten hinter den P1-Knoten kommen.
    // Sonst muessen fuer den Druck die P1-Knoten extrahiert werden.
    // TODO: Wahrscheinlich reicht nur FSI-Domain, da selbes Objekt bei problemChem_ und problemTimeChem_.
    ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemChem_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeChem_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
}


// Am Anfang der Zeititeration erst updateSolutionMultiPreviousStep() aufrufen und dann erst updateMultistepRhs(),
// damit die previousSolution_ initialisiert sind. Genauso fuer SystemMass
// TODO: updateSystemMassMultiPreviousStep() fertig programmieren
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateChemInTime() const
{
    int nmbBDF = timeSteppingTool_->getBDFNumber();

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
}
#endif
