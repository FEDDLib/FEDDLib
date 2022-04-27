#ifndef SCI_def_hpp
#define SCI_def_hpp
#include "SCI_decl.hpp"

double OneFunc(double* x, int* parameter)
{
    return 1.0;
}

namespace FEDD {
// Funktionen fuer die rechte Seite der Struktur/ Fluid/ Geometrie sind im jeweiligen Problem

template<class SC,class LO,class GO,class NO>
SCI<SC,LO,GO,NO>::SCI(const DomainConstPtr_Type &domainStructure, std::string FETypeStructure,
					const DomainConstPtr_Type &domainChem, std::string FETypeChem,
                    ParameterListPtr_Type parameterListChem, ParameterListPtr_Type parameterListStructure,
                    ParameterListPtr_Type parameterListSCI, Teuchos::RCP<SmallMatrix<int> > &defTS):
Problem<SC,LO,GO,NO>( parameterListFSI, domainVelocity->getComm() ),
// hasSourceTerm = drittes Arguement. assembleSourceTerm() fuer NS nicht programmiert.
// Deswegen hier erstmal false (default Parameter).
// Fuer Struktur hingegen ist default Parameter true, da programmiert.
problemStructure_(),
problemChem_(),
//problemStructureNonLin_(),
meshDisplacementOld_rep_(),
meshDisplacementNew_rep_(),
u_rep_(),
w_rep_(),
u_minus_w_rep_(),
defTS_(defTS),
timeSteppingTool_(),
materialModel_( parameterListStructure->sublist("Parameter").get("Material model","linear") ),
valuesForExport_(0),
exporterGeo_()
{
    //this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);

    //this->initNOXParameters();
    counterP = 0;
    
    //std::string linearization = parameterListSCI->sublist("General").get("Linearization","FixedPoint");
    
    //TEUCHOS_TEST_FOR_EXCEPTION( !(linearization == "Newton" || linearization == "NOX")  && materialModel_ != "linear", std::runtime_error, "Nonlinear material models can only be used with Newton's method or FixedPoint (nonlinear material Jacobian will still be used).");
    
    this->addVariable( domainStructure, FETypeStructure, "d_s", domainStructure->getDimension() ); // Struktur
    this->addVariable( domainChem, FETypeChem, "c", domainChem->getDimension()); // Chemistry

    this->dim_ = this->getDomain(0)->getDimension();
    
    problemChem_ = Teuchos::rcp( new ChemProblem_Type( domainVelocity, FETypeVelocity, domainPressure, FETypePressure, parameterListFluid ) );
    problemChem_->initializeProblem();
    
    //if (materialModel_=="linear"){
        problemStructure_ = Teuchos::rcp( new StructureProblem_Type( domainStructure, FETypeStructure, parameterListStructure ) );
        problemStructure_->initializeProblem();
    /*}
    else{
        problemStructureNonLin_ = Teuchos::rcp( new StructureNonLinProblem_Type( domainStructure, FETypeStructure, parameterListStructure) );
        problemStructureNonLin_->initializeProblem();
    }*/
    
    //We initialize the subproblems. In the main routine, we need to call initializeFSI(). There, we first initialize the vectors of the FSI problem and then we set the pointers of the subproblems to the vectors of the full monolithic FSI system. This way all values are only saved once in the subproblems and can be used by the monolithic FSI system.
    
    meshDisplacementNew_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(4)->getMapVecFieldRepeated() ) );
    meshDisplacementOld_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(4)->getMapVecFieldRepeated() ) );
    u_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    u_minus_w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
}


template<class SC,class LO,class GO,class NO>
SCI<SC,LO,GO,NO>::~SCI()
{
    if (!exporterGeo_.is_null()) {
       exporterGeo_->closeExporter();
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
        MultiVectorPtr_Type eModVec = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getElementMap() ) );

        MultiVectorPtr_Type solChemRep = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) )
        solChemRep->importFromVector(this->problemChem_->getSolution());

        this->feFactory_->determineEMod(FEType,solChemRep,eModVec,this->getDomain(1));

        /*

        MultiVectorPtr_Type solChemRep = Teuchos::rcp(new MultiVector_Type ( this->getDomain(0)->getMapRepeated()));

        Teuchos::ArrayRCP< SC > c = solChemRep->getDataNonConst(0);
        vec2D_dbl_ptr_Type pointsRep = this->getDomain(0)->getPointsRepeated();
        int dim = this->getDomain(0)->getDimension();
        for(int i=0; i< pointsRep->size(); i++){
            c[i] = 1.; //(1.-pointsRep->at(i).at(dim-1));        
        }

        MultiVectorPtr_Type eModVec = Teuchos::rcp(new MultiVector_Type ( this->getDomain(0)->getElementMap()));
        this->feFactory_->determineEMod(this->getDomain(0)->getFEType(),solChemRep,eModVec,this->getDomain(0));
        
        //this->feFactory_->assemblyLinElasXDimE(this->dim_,this->getDomain(0)->getFEType(), K, eModVec, poissonRatio, true);
        */


        // steady rhs wird hier assembliert.
        // rhsFunc auf 0 (=x) und 0 (=y) abaendern bei LinElas!
        //if (materialModel_=="linear"){
            this->feFactory_->assemblyLinElasXDimE(int dim, FEType,tr_Type &A, eModVec, nu, true);

            //this->problemStructure_->assemble();
        /*}
        else
            this->problemStructureNonLin_->assemble();*/
      
       
        MatrixPtr_Type C1_T(new Matrix_Type( this->getDomain(1)->getMapVecUnique(), 1 ) ); // Chem-Kopplung
        MatrixPtr_Type C1(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) ); // Struktur-Kopplung
       
        // ###########################
        // Korrekte Skalierung der entsprechenden Bloecke
        // ###########################
        double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);
        
        // ###########################
        // Bloecke hinzufuegen
        // ###########################
        
        this->system_.reset(new BlockMatrix_Type(2));

         
        // Struktur
        this->system_->addBlock( this->problemStructure_->system_->getBlock(0,0), 0, 0);
        
        // Chemistry
        this->system_->addBlock( this->probleChem_->system_->getBlock(0,0), 1, 1 );
       
        // We set the vector from the partial problems
        this->setFromPartialVectorsInit();
        
        // Fuer die Zeitprobleme
        timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));
        ParameterListPtr_Type plStructure;
        
        plStructure = this->problemStructure_->getParameterList();
        
        setupSubTimeProblems(this->problemFluid_->getParameterList(), plStructure);
        
        if (this->verbose_)
        {
            std::cout << "done -- " << std::endl;
        }
    }
    else if (type == "MoveMesh")
         moveMesh();
    else if (type == "UpdateMeshDisplacement")
         updateMeshDisplacement();
    else if (type == "ComputeSolidRHSInTime")
         computeSolidRHSInTime();
    else if (type == "UpdateChemInTime")
          updateChemInTime();

}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateMeshDisplacement() const
{

     *meshDisplacementOld_rep_ = *meshDisplacementNew_rep_;

}




template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setupSubTimeProblems(ParameterListPtr_Type parameterListFluid, ParameterListPtr_Type parameterListStructure) const
{
    if(this->verbose_)
        std::cout << "-- Setup FSI Sub-TimeProblems \n" << std::flush;

    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();

    int sizeFluid = this->problemFluid_->getSystem()->size();
    int sizeStructure;
    if (materialModel_=="linear")
        sizeStructure = this->problemStructure_->getSystem()->size();
    else
        sizeStructure = this->problemStructureNonLin_->getSystem()->size();
    
    problemTimeFluid_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemFluid_, this->comm_));
    if (materialModel_=="linear")
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructure_, this->comm_));
    else
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructureNonLin_, this->comm_));

    // ######################
    // Fluid: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    SmallMatrix<double> massCoeffFluid(sizeFluid);
    SmallMatrix<double> problemCoeffFluid(sizeFluid);
    SmallMatrix<int> defFluid(sizeFluid);

    double coeffSourceTermFluid = 0.0;
    if ( this->getParameterList()->sublist("Timestepping Parameter").get("Class","Multistep") == "Multistep" ) {
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++) {
                if ((*defTS_)[i][j]==1 && i==j) {
                    defFluid[i][j] = 1;
                    massCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffFluid[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++){
                if ((*defTS_)[i][j]==1){
                    problemCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(1);
                    coeffSourceTermFluid = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffFluid[i][j] = 1.;
                }
            }
        }
        this->problemTimeFluid_->setTimeDef(defFluid);
        this->problemTimeFluid_->setTimeParameters(massCoeffFluid,problemCoeffFluid);
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Implement other FSI fluid time stepping than BDF.");
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
            if((*defTS_)[i + sizeFluid][j + sizeFluid] == 1 && i == j) // Weil: (u_f, p, d_s,...) und timeStepDef_ von FSI
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
            if((*defTS_)[i + sizeFluid][j + sizeFluid] == 1)
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

    this->problemTimeFluid_->assemble( "MassSystem" );
    this->problemTimeStructure_->assemble( "MassSystem" );
}


template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::setChemMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix fuer FSI combineSystems(), ggf nichtlinear.
    //######################
    double density = this->problemTimeChem_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
    int size = this->problemTimeChem_->getSystem()->size();

    this->problemTimeChem_->systemMass_.reset(new BlockMatrix_Type(size));
    {
        massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeChem_->getDomain(0)->getMapdUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
        // 0 = Fluid
        this->feFactory_->assemblyMass( this->dim_, this->problemTimeChem_->getFEType(0), "Vector",  massmatrix, 0, true );
        massmatrix->resumeFill();
        massmatrix->scale(density);
        massmatrix->fillComplete( this->problemTimeFluid_->getDomain(0)->getMapVecFieldUnique(), this->problemTimeFluid_->getDomain(0)->getMapVecFieldUnique() );

        this->problemTimeChem_->systemMass_->addBlock(massmatrix, 0, 0);
    }
}


// TODO: updateMultistepRhsFSI() einbauen!
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::computeFluidRHSInTime( ) const
{
    //######################
    // RHS nach BDF2
    //######################
    int sizeFluid = this->problemFluid_->getSystem()->size();
    double dt = timeSteppingTool_->get_dt();
    int nmbBDF = timeSteppingTool_->getBDFNumber();

    vec_dbl_Type coeffPrevSteps(nmbBDF);
    for(int i = 0; i < coeffPrevSteps.size(); i++)
    {
        coeffPrevSteps.at(i) = timeSteppingTool_->getInformationBDF(i+2) / dt;
    }

    if (timeSteppingTool_->currentTime()==0.) {
        SmallMatrix<double> tmpmassCoeff(sizeFluid);
        SmallMatrix<double> tmpproblemCoeff(sizeFluid);
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++) {
                if ((*defTS_)[i][j]==1 && i==j) {
                    tmpmassCoeff[i][j] = 1. / dt;
                }
                else{
                    tmpmassCoeff[i][j] = 0.;
                }
            }
        }
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++){
                if ((*defTS_)[i][j]==1){
                    tmpproblemCoeff[i][j] =  1.; // ist das richtig? Vermutlich schon, da BDF so geschrieben ist, dass zu berechnende Lsg den Koeffizienten 1 hat
                }
                else{
                    tmpproblemCoeff[i][j] = 1.;
                }
            }
        }
        this->problemTimeFluid_->setTimeParameters(tmpmassCoeff, tmpproblemCoeff);
    }
    if (timeSteppingTool_->currentTime()==0.) {
        vec_dbl_Type tmpcoeffPrevSteps(1, 1. / dt);
        this->problemTimeFluid_->updateMultistepRhsFSI(tmpcoeffPrevSteps,1);/*apply (mass matrix_t / dt) to u_t*/
    }
    else{
        this->problemTimeFluid_->updateMultistepRhsFSI(coeffPrevSteps,nmbBDF);/*apply (mass matrix_t / dt) to u_t and more*/
    }

    // TODO
    if (this->problemTimeFluid_->hasSourceTerm()) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Check sourceterm.");
    }

    // Wieder zu den eigentlichen Parametern zuruecksetzen, nachdem die temporaeren
    // genommen wurden.
    if (timeSteppingTool_->currentTime()==0.) {
        SmallMatrix<double> massCoeffFluid(sizeFluid);
        SmallMatrix<double> problemCoeffFluid(sizeFluid);

        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++) {
                if ((*defTS_)[i][j]==1 && i==j) {
                    massCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffFluid[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++){
                if ((*defTS_)[i][j]==1){
                    problemCoeffFluid[i][j] = timeSteppingTool_->getInformationBDF(1);
                }
                else{
                    problemCoeffFluid[i][j] = 1.;
                }
            }
        }

        this->problemTimeFluid_->setTimeParameters(massCoeffFluid, problemCoeffFluid);
    }
}



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
    this->problemTimeStructure_->updateSolutionNewmarkPreviousStep(dt, beta, gamma);
    
    // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
    // Bei Newmark lautet dies:
    // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
    // wobei u' = v (velocity) und u'' = w (acceleration).
    this->problemTimeStructure_->updateNewmarkRhs(dt, beta, gamma, coeffTemp);
    
    //can we get rid of this?
    double time = timeSteppingTool_->currentTime() + dt;
    
    // TODO: SourceTerm wird in jedem Zeitschritt neu berechnet; auch wenn konstant!!!
    // if(time == 0){nur dann konstanten SourceTerm berechnen}
    if (this->problemTimeStructure_->hasSourceTerm())
    {
        this->problemTimeStructure_->assembleSourceTerm( time );
        
        // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
        // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!
        
        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
        BlockMultiVectorPtr_Type tmpPtr = this->problemTimeStructure_->getSourceTerm();
        this->problemTimeStructure_->getRhs()->update(coeffSourceTermStructure, *tmpPtr, 1.);
    }

}

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
            // 2 = Struktur
            this->feFactory_->assemblyMass(this->dim_, this->problemTimeStructure_->getFEType(0), "Vector", massmatrix, 2, true);
            massmatrix->resumeFill();
            massmatrix->scale(density);
            massmatrix->fillComplete( this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique());

            this->problemTimeStructure_->systemMass_->addBlock( massmatrix, 0, 0 );
        }
    }
}


// Damit die richtige timeSteppingTool_->currentTime() genommen wird.
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::updateTime() const
{
    timeSteppingTool_->t_ = timeSteppingTool_->t_ + timeSteppingTool_->dt_prev_;
}




template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::addInterfaceBlockRHS() const
{
    MultiVectorPtr_Type vectorToAdd = Teuchos::rcp( new MultiVector_Type( this->rhs_->getBlock(3) ) );

    C2_->apply(*(this->solution_->getBlock(2)), *vectorToAdd);
    this->rhs_->addBlock(vectorToAdd, 3);
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
void SCI<SC,LO,GO,NO>::getValuesOfInterest3DBenchmark( vec_dbl_Type& values ){
    
    if ( valuesForExport_[0] >= 0.   ) { // we set the displacement of the 3D Richter Benchmark in x, y, and z direction
        LO loc = valuesForExport_[0] + 10*Teuchos::ScalarTraits<SC>::eps();
        values[0] = this->getSolution()->getBlock(2)->getData(0)[3*loc];
        values[1] = this->getSolution()->getBlock(2)->getData(0)[3*loc+1];
        values[2] = this->getSolution()->getBlock(2)->getData(0)[3*loc+2];
    }
}
    
template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::computeValuesOfInterestAndExport(){
    if ( this->getParameterList()->sublist("General").get("Export drag and lift",false) ) {
        
        int dim = this->dim_;
        TEUCHOS_TEST_FOR_EXCEPTION( this->parameterList_->sublist("Parameter").get("Criterion","Residual") == "Update",  std::runtime_error, "Wrong nonlinear criterion to calculate the drag coefficient. The last system is the Newton system but we need the fixed point system. Either use Criterion=Residual or implement for Criterion=Update." );
        
        TEUCHOS_TEST_FOR_EXCEPTION( this->problemFluid_->hasSourceTerm(),  std::runtime_error, "We need to substract the additional source term: drag = < F*u + B_T*p + C1_T*lamba - f, v >" );
        
        Teuchos::Array<SC> drag(1);
        Teuchos::Array<SC> lift(1);
        
        BlockMultiVectorPtr_Type uDrag = Teuchos::rcp( new BlockMultiVector_Type( this->problemFluid_->getSolution() ) );
        BlockMultiVectorPtr_Type uLift = Teuchos::rcp( new BlockMultiVector_Type( this->problemFluid_->getSolution() ) );
        // should be the last fixed point system without boundary conditions or the last extrapolation system without boundary values.
        // We need to reassemble B and BT, because we might have set Dirichlet boundary conditions in BT (less likely in B)
        this->problemFluid_->assembleDivAndStab();
        
        this->problemFluid_->getSystem()->apply( *this->problemFluid_->getSolution(), *uDrag );
        this->problemFluid_->getSystem()->apply( *this->problemFluid_->getSolution(), *uLift );
        
        MultiVectorPtr_Type C1T_lambda = Teuchos::rcp( new MultiVector_Type( this->getSolution()->getBlock(0) ) );
        this->system_->getBlock(0,3)->apply( *this->getSolution()->getBlock(3), *C1T_lambda );
        
        uDrag->getBlockNonConst(0)->update( 1., *C1T_lambda, 1. ); // velocity + C1_T * lambda
        uLift->getBlockNonConst(0)->update( 1., *C1T_lambda, 1. ); // velocity + C1_T * lambda
        
        BCPtr_Type bcFactoryDrag = Teuchos::rcp( new BC_Type( ) );
        BCPtr_Type bcFactoryLift = Teuchos::rcp( new BC_Type( ) );
        
        DomainConstPtr_Type domainVelocityConst = this->problemFluid_->getDomain(0);
        DomainPtr_Type domainVelocity = Teuchos::rcp_const_cast<Domain_Type>(domainVelocityConst);
        if( dim == 2 ){
            bcFactoryDrag->addBC(drag2D, 4, 0, domainVelocity, "Dirichlet", dim); // obstacle
            bcFactoryDrag->addBC(drag2D, 5, 0, domainVelocity, "Dirichlet", dim); // interface; check main fsi for matching flags at the obstacle and interface
            bcFactoryLift->addBC(lift2D, 4, 0, domainVelocity, "Dirichlet", dim);
            bcFactoryLift->addBC(lift2D, 5, 0, domainVelocity, "Dirichlet", dim);
        }
        else if( dim == 3 ){
            bcFactoryDrag->addBC(drag3D, 3, 0, domainVelocity, "Dirichlet", dim); // check main fsi for matching
            bcFactoryDrag->addBC(drag3D, 6, 0, domainVelocity, "Dirichlet", dim); // check main fsi for matching flags at the obstacle and interface
            bcFactoryLift->addBC(lift3D, 3, 0, domainVelocity, "Dirichlet", dim);
            bcFactoryLift->addBC(lift3D, 6, 0, domainVelocity, "Dirichlet", dim);
        }
        
        BlockMultiVectorPtr_Type vD = Teuchos::rcp( new BlockMultiVector_Type( this->problemFluid_->getSolution() ) );
        BlockMultiVectorPtr_Type vL = Teuchos::rcp( new BlockMultiVector_Type( this->problemFluid_->getSolution() ) );
        
        vD->putScalar(0.);
        vL->putScalar(0.);
        
        bcFactoryDrag->setRHS( vD );
        bcFactoryLift->setRHS( vL );
        
        uDrag->dot( vD, drag() );
        uLift->dot( vL, lift() );
        
//        double density = this->problemFluid_->getParameterList()->sublist("Parameter").get("Density",1.);
//        double uMean = this->getParameterList()->sublist("Parameter").get("MeanVelocity",2.0);
//        double L = 0.;
//        if ( dim == 2)
//            L = 2.*0.05;
//        else
//            L = 1.;
//        
//        drag[0] *= -(2./(density*uMean*uMean*L));
//        lift[0] *= -(2./(density*uMean*uMean*L));

        drag[0] *= -1.;
        lift[0] *= -1.;
        
        exporterTxtDrag_->exportData( drag[0] );
        exporterTxtLift_->exportData( lift[0] );
    }
}

template<class SC,class LO,class GO,class NO>
void SCI<SC,LO,GO,NO>::moveMesh() const
{

    MultiVectorConstPtr_Type displacementUniqueConst;

    displacementUniqueConst = this->solution_->getBlock(0);
    
    MultiVectorPtr_Type displacementRepeated = Teuchos::rcp( new MultiVector_Type( this->problemGeometry_->getDomain(0)->getMapVecFieldRepeated() ) );

    displacementRepeated->importFromVector( displacementUniqueConst );
    MultiVectorPtr_Type displacementUnique = Teuchos::rcp_const_cast<MultiVector_Type>(displacementUniqueConst);


    // Verschiebe die Gitter fuer Chemistry
    // ACHTUNG: Klappt nur, weil die P2-Knoten hinter den P1-Knoten kommen.
    // Sonst muessen fuer den Druck die P1-Knoten extrahiert werden.
    // TODO: Wahrscheinlich reicht nur FSI-Domain, da selbes Objekt bei problemFluid_ und problemTimeFluid_.
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
