#ifndef FSI_def_hpp
#define FSI_def_hpp
#include "FSI_decl.hpp"

double OneFunc(double* x, int* parameter)
{
    return 1.0;
}

void drag2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 1.;
    res[1] = 0.;
    
    return;
}

void drag3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 1.;
    res[1] = 0.;
    res[2] = 0.;
    
    return;
}

void lift2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 1.;
    
    return;
}

void lift3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 1.;
    res[2] = 0.;
    
    return;
}
void pressureRamp(double* x, double* res, double* parameters){

    double pressureValue = parameters[1];
    double flag = parameters[3];
    double ramp = parameters[2];
  	res[0] =0.;
    
    if(parameters[0]+1.e-12 < ramp)
        pressureValue = parameters[0]*pressureValue/ramp;
    else
        pressureValue = parameters[1];

    res[0] = pressureValue;  // Usually we check here for the correct flag. But as the boundary condition is limited to the outlet anyway, we have no issues
    

    return;
}

namespace FEDD {
// Funktionen fuer die rechte Seite der Struktur/ Fluid/ Geometrie sind im jeweiligen Problem

template<class SC,class LO,class GO,class NO>
FSI<SC,LO,GO,NO>::FSI(const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity,
                    const DomainConstPtr_Type &domainPressure, std::string FETypePressure,
                    const DomainConstPtr_Type &domainStructure, std::string FETypeStructure,
                    const DomainConstPtr_Type &domainInterface, std::string FETypeInterface,
                    const DomainConstPtr_Type &domainGeometry, std::string FETypeGeometry,
                    ParameterListPtr_Type parameterListFluid, ParameterListPtr_Type parameterListStructure,
                    ParameterListPtr_Type parameterListFSI, ParameterListPtr_Type parameterListGeometry,
                    Teuchos::RCP<SmallMatrix<int> > &defTS):
NonLinearProblem<SC,LO,GO,NO>( parameterListFSI, domainVelocity->getComm() ),
// hasSourceTerm = drittes Arguement. assembleSourceTerm() fuer NS nicht programmiert.
// Deswegen hier erstmal false (default Parameter).
// Fuer Struktur hingegen ist default Parameter true, da programmiert.
P_(),
problemFluid_(),
problemStructure_(),
problemStructureNonLin_(),
problemGeometry_(),
meshDisplacementOld_rep_(),
meshDisplacementNew_rep_(),
u_rep_(),
w_rep_(),
u_minus_w_rep_(),
p_rep_(),
defTS_(defTS),
timeSteppingTool_(),
materialModel_( parameterListStructure->sublist("Parameter").get("Material model","linear") ),
valuesForExport_(0),
exporterTxtDrag_(),
exporterGeo_()
{
    this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    geometryExplicit_ = this->parameterList_->sublist("Parameter").get("Geometry Explicit",true);

    this->initNOXParameters();
    counterP = 0;
    
    std::string linearization = parameterListFSI->sublist("General").get("Linearization","FixedPoint");
    
    TEUCHOS_TEST_FOR_EXCEPTION( !(linearization == "Newton" || linearization == "NOX")  && materialModel_ != "linear", std::runtime_error, "Nonlinear material models can only be used with Newton's method or FixedPoint (nonlinear material Jacobian will still be used).");
    
    this->addVariable( domainVelocity, FETypeVelocity, "u_f", domainVelocity->getDimension() ); // Fluid-Geschw.
    this->addVariable( domainPressure, FETypePressure, "p", 1); // Fluid-Druck
    this->addVariable( domainStructure, FETypeStructure, "d_s", domainStructure->getDimension() ); // Struktur
    this->addVariable( domainInterface, FETypeInterface, "lambda", domainInterface->getDimension() ); // Interface
    this->addVariable( domainGeometry, FETypeGeometry, "d_f", domainGeometry->getDimension() ); // Geometrie
    this->dim_ = this->getDomain(0)->getDimension();
    
    problemFluid_ = Teuchos::rcp( new FluidProblem_Type( domainVelocity, FETypeVelocity, domainPressure, FETypePressure, parameterListFluid ) );
    problemFluid_->initializeProblem();
    problemFluid_->infoParameter(false,"Fluid");
    
    if (materialModel_=="linear"){
        problemStructure_ = Teuchos::rcp( new StructureProblem_Type( domainStructure, FETypeStructure, parameterListStructure ) );
        problemStructure_->initializeProblem();
        problemStructure_->infoParameter(false,"Solid");
    }
    else{
        problemStructureNonLin_ = Teuchos::rcp( new StructureNonLinProblem_Type( domainStructure, FETypeStructure, parameterListStructure) );
        problemStructureNonLin_->initializeProblem();
        problemStructureNonLin_->infoParameter(false, "Solid");
    }
    
    problemGeometry_ = Teuchos::rcp( new GeometryProblem_Type( domainGeometry, FETypeGeometry, parameterListGeometry ) );
    problemGeometry_->initializeProblem();
    problemGeometry_->infoParameter(false,"Geometry");
    //We initialize the subproblems. In the main routine, we need to call initializeFSI(). There, we first initialize the vectors of the FSI problem and then we set the pointers of the subproblems to the vectors of the full monolithic FSI system. This way all values are only saved once in the subproblems and can be used by the monolithic FSI system.
    
    meshDisplacementNew_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(4)->getMapVecFieldRepeated() ) );
    meshDisplacementOld_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(4)->getMapVecFieldRepeated() ) );
    u_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    u_minus_w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    
    // Exporting benchmark values
    if ( parameterListFSI->sublist("General").get("Export Extra Data",false) ){
        if (this->dim_==2)
            findDisplacementTurek2DBenchmark();
        else if (this->dim_==3)
            findDisplacementRichter3DBenchmark();
    }
    if ( parameterListFSI->sublist("General").get("Export drag and lift",false) ){
        exporterTxtDrag_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtDrag_->setup( "drag_force", this->comm_ );
        exporterTxtLift_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtLift_->setup( "lift_force", this->comm_ );
    }
    p_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );
    
    // if ( this->parameterList_->sublist("Timestepping Parameter").get("Checkpointing", false)){
    //     exporterBoundaryCondition_ = Teuchos::rcp(new ExporterTxt () );
    //     exporterBoundaryCondition_->setup( "boundaryConditionFluid", this->comm_ );
    // }
}

template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::findDisplacementTurek2DBenchmark(){
    valuesForExport_.resize(1);
    vec_dbl_Type x = {0.6,0.2};
    //we save the local index as a double
    valuesForExport_[0] = this->getDomain(2)->findInPointsUnique( x );
}
    
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::findDisplacementRichter3DBenchmark(){
    valuesForExport_.resize(1);
    vec_dbl_Type x = {0.45,0.15,0.15};
    //we save the local index as a double
    valuesForExport_[0] = this->getDomain(2)->findInPointsUnique( x );
}
    
template<class SC,class LO,class GO,class NO>
FSI<SC,LO,GO,NO>::~FSI()
{
    if (!exporterGeo_.is_null()) {
       exporterGeo_->closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::info()
{
    this->infoProblem();
    this->infoNonlinProblem();
}

// Alle linearen Probleme
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::assemble( std::string type ) const
{
    if (type == "") {
        if (this->verbose_)
        {
            std::cout << "-- Assembly FSI ... " << std::endl;
        }

    //    P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), 10 ) );
        
        this->problemFluid_->assemble();
        
        // steady rhs wird hier assembliert.
        // rhsFunc auf 0 (=x) und 0 (=y) abaendern bei LinElas!
        if (materialModel_=="linear")
            this->problemStructure_->assemble();
        else
            this->problemStructureNonLin_->assemble();
        
        this->problemGeometry_->assemble();
        
        if ( geometryExplicit_ && this->parameterList_->sublist("Exporter").get("Export GE geometry solution",false)){
            exporterGeo_ = Teuchos::rcp(new Exporter_Type());
            
            DomainConstPtr_Type dom = this->getDomain(4);

            int exportEveryXTimesteps = this->parameterList_->sublist("Exporter").get( "Export every X timesteps", 1 );
            std::string suffix = this->parameterList_->sublist("Exporter").get("Geometry Suffix", "" );
            std::string varName = "d_f" + suffix;
            
            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( dom->getMesh() );
            exporterGeo_->setup(varName, meshNonConst, dom->getFEType(), exportEveryXTimesteps, this->parameterList_);
            
            MultiVectorConstPtr_Type exportVector = this->problemGeometry_->getSolution()->getBlock(0);
            
            exporterGeo_->addVariable( exportVector, varName, "Vector", this->dim_, dom->getMapUnique() );
        }
        
        
        if (!geometryExplicit_) {
            // RW fuer die Geometrie-Matrix H setzen, weil wenn wir das spaeter im FSI-System
            // machen, dann wird aufgrund der RW auf dem Interface auch der Kopplungsblock
            // zur Struktur C4 ausgenullt, was wir nicht wollen.
            this->problemGeometry_->setBoundariesSystem(); // RW im System setzen (mit den Flags von oben)
            this->problemGeometry_->getRhs()->putScalar(0.0);
        }

        // ###########################
        // Kopplungsbloecke
        // ###########################
        // ACHTUNG: Anders als im Matlab-Code!!!
        // Matlab: C1_T hat so viele Zeilen wie das komplette Fluid-System (u + p)
        // FEDDLib: C1_T hat so viele Zeilen wie das Geschwindigkeitssystem (nur u)
        // Bemerkung: Verteilung der Zeilen wird angegeben (Range-Map).
        // Interface wird von der Fluid-Seite aus gehalten, deswegen auch
        // getDomain(0) bei C2, obwohl es in der Struktur-Spalte ist.
        MatrixPtr_Type C1_T(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), 1 ) ); // Fluid-Kopplung
        MatrixPtr_Type C1(new Matrix_Type( this->getDomain(0)->getInterfaceMapVecFieldUnique(), 1 ) ); // Fluid-Spalte
        MatrixPtr_Type C2(new Matrix_Type( this->getDomain(0)->getInterfaceMapVecFieldUnique(), 1 ) ); // Struktur-Spalte
        MatrixPtr_Type C3_T(new Matrix_Type( this->getDomain(2)->getMapVecFieldUnique(), 1 ) ); // Struktur-Kopplung
        
        
        // Nur vorhanden, falls geometrisch implizit
        MatrixPtr_Type C4(new Matrix_Type( this->getDomain(4)->getMapVecFieldUnique(), 1 ) ); // Geometrie (=Fluid)

        // Fluid-Bloecke
        this->feFactory_->assemblyFSICoupling(this->dim_, this->domain_FEType_vec_.at(0), C1, C1_T, 0, 0,
            this->getDomain(0)->getInterfaceMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique(), true);
        
        // Struktur-Bloecke
        // ACHTUNG: Die Interface-Variable \lambda wird eindeutig von der Fluid-Seite gehalten.
        // Deswegen auch getDomain(0) fuer die Spalten der Kopplungsbloecke.
        this->feFactory_->assemblyFSICoupling(this->dim_, this->domain_FEType_vec_.at(2), C2, C3_T, 0, 2,
            this->getDomain(0)->getInterfaceMapVecFieldUnique(), this->getDomain(2)->getMapVecFieldUnique(), true);

        
        // Falls geometrisch implizit
        if(!geometryExplicit_)
        {
            // TODO: Wegen IndicesGlobalMatched_ vlt .sicherheitshalber FEloc = 0 nehmen.
            // TODO: Check C4
            this->feFactory_->assemblyGeometryCoupling(this->dim_, this->domain_FEType_vec_.at(4), C4, 4,
                                                       this->getDomain(0)->getGlobalInterfaceMapUnique(),
                                                       this->getDomain(2)->getMapVecFieldUnique(),
                                                       this->getDomain(4)->getMapVecFieldUnique(), true);
        }
        
        MatrixPtr_Type dummyC;
        // we need to set the dummy coupling conditions for monolithic preconditioning with FROSch
        if ( !this->getDomain(0)->getGlobalInterfaceMapVecFieldPartial().is_null() ) {
            dummyC.reset(new Matrix_Type( this->getDomain(0)->getInterfaceMapVecFieldUnique(), 1 ) );
            this->feFactory_->assemblyDummyCoupling(this->dim_, this->domain_FEType_vec_.at(0), dummyC, 0,true);
        }

        // ###########################
        // Korrekte Skalierung der entsprechenden Bloecke
        // ###########################
        double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);

        C2->resumeFill();
        C3_T->resumeFill();

        C2->scale( -(1.0/dt) ); // this will be used in a first order approximation of the solid velocity
        C3_T->scale( -1.0 );
        
        // ACHTUNG: Die Interface-Variable \lambda wird eindeutig von der Fluid-Seite gehalten.
        // Deswegen auch getDomain(0) fuer die Spalten von C3_T.
        C2->fillComplete(this->getDomain(2)->getMapVecFieldUnique(), this->getDomain(0)->getInterfaceMapVecFieldUnique());
        C3_T->fillComplete(this->getDomain(0)->getInterfaceMapVecFieldUnique(), this->getDomain(2)->getMapVecFieldUnique());

        // C2 in Membervariable C2_ speichern, fuer rechte Seite im Interface-Block:
        // C2*d_s^n
        C2_ = C2;

        if(!geometryExplicit_)
        {
            C4->resumeFill();
            C4->scale(-1.0);
            // Domain = Spalten = Struktur; Range = Zeilen = Geometrie
            C4->fillComplete(this->getDomain(2)->getMapVecFieldUnique(), this->getDomain(4)->getMapVecFieldUnique());
        }
        

        
        // ###########################
        // Bloecke hinzufuegen
        // ###########################
        if(geometryExplicit_)
            this->system_.reset(new BlockMatrix_Type(4));
        else
            this->system_.reset(new BlockMatrix_Type(5));
        
        // Fluid
        this->system_->addBlock( this->problemFluid_->system_->getBlock(0,0), 0, 0 );
        this->system_->addBlock( this->problemFluid_->system_->getBlock(0,1), 0, 1 );
        this->system_->addBlock( this->problemFluid_->system_->getBlock(1,0), 1, 0 );
        if (this->getDomain(0)->getFEType()=="P1")
            this->system_->addBlock( this->problemFluid_->system_->getBlock(1,1), 1, 1 );
        
        // Struktur
        if (materialModel_=="linear")
            this->system_->addBlock( this->problemStructure_->system_->getBlock(0,0), 2, 2 );
        else
            this->system_->addBlock( this->problemStructureNonLin_->system_->getBlock(0,0), 2, 2 );
        // Kopplung
        this->system_->addBlock( C1_T, 0, 3 );
        this->system_->addBlock( C3_T, 2, 3 );
        this->system_->addBlock( C1, 3, 0 );
        this->system_->addBlock( C2, 3, 2 );

        if (!dummyC.is_null())
            this->system_->addBlock( dummyC, 3, 3 );
        
        if(!geometryExplicit_)
        {
            // Geometrie
            this->system_->addBlock( this->problemGeometry_->system_->getBlock(0,0), 4, 4 );
            // Kopplung
            this->system_->addBlock( C4, 4, 2 );
        }

        // Sollte (bzw. muss) erst aufgerufen werden, wenn alle Bloecke aus assemble()
        // in das ganze System hineingeschrieben worden ist. Ansonsten findet
        // blockExists() keinen Block
    //    this->initializeVectors();
    //    this->initializeVectorsNonLinear();
        //NOX

        // We set the vector from the partial problems
        this->setFromPartialVectorsInit();
        
        // Fuer die Zeitprobleme
        timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));
        ParameterListPtr_Type plStructure;
        if (materialModel_=="linear")
            plStructure = this->problemStructure_->getParameterList();
        else
            plStructure = this->problemStructureNonLin_->getParameterList();
        
        setupSubTimeProblems(this->problemFluid_->getParameterList(), plStructure);
        
        if (this->verbose_)
        {
            std::cout << "done -- " << std::endl;
        }
    }
    else
        reAssemble(type);
}


template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::reAssemble(std::string type) const
{

    double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);

    // Fluid-Dichte
    double density = this->problemFluid_->getParameterList()->sublist("Parameter").get("Density",1.);
    double viscosity = this->problemFluid_->getParameterList()->sublist("Parameter").get("Viscosity",1.);

    if(type == "UpdateMeshDisplacement")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateMeshDisplacement)" << '\n';

        // Da dieser Abschnitt zu Beginn der neuen Zeititeration aufgerufen wird, muss
        // auch die alte Gitterverschiebung durch die neue Loesung aktualisiert werden.
        updateMeshDisplacement();
        return;
    }

    if(type == "SolveGeometryProblem")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (SolveGeometryProblem)" << '\n';

        solveGeometryProblem();
        return;
    }

    if(type == "UpdateTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateTime)" << '\n';

        updateTime();
        return;
    }

    if(type == "UpdateFluidInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateFluidInTime)" << '\n';

        updateFluidInTime();
        return;
    }

    if(type == "MoveMesh")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (MoveMesh)" << '\n';

        moveMesh();
        return;
    }

    if(type == "AddInterfaceBlockRHS")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (AddInterfaceBlockRHS)" << '\n';

        addInterfaceBlockRHS();
        return;
    }

    if(type == "ComputeFluidRHSInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (ComputeFluidRHSInTime)" << '\n';
        
        computeFluidRHSInTime( );
        return;
    }
    if(type == "ComputeSolidRHSInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (ComputeSolidRHSInTime)" << '\n';
        
        computeSolidRHSInTime( );
        return;
    }
    if(type == "ComputePressureRHSInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (ComputePressureRHSInTime)" << '\n';
        
        computePressureRHSInTime();
        return;
    }

    
    // ###############
    // Berechne die Gittergeschwindigkeit w und FSI-Konvektion-Anteil (u-w)
    // ###############
    MultiVectorConstPtr_Type fluidSolution = this->solution_->getBlock(0);
    u_rep_->importFromVector(fluidSolution, true);
    u_minus_w_rep_->importFromVector(fluidSolution, true); //    u_minus_w_rep_ = *u_rep_;

    MultiVectorConstPtr_Type geometrySolution;
    if(geometryExplicit_)
    {
        geometrySolution = this->problemGeometry_->getSolution()->getBlock(0);
    }
    else
    {
        geometrySolution = this->solution_->getBlock(4);
    }
    meshDisplacementNew_rep_->importFromVector(geometrySolution, true);

    *w_rep_ = *meshDisplacementNew_rep_;
    w_rep_->update(-1.0, *meshDisplacementOld_rep_, 1.0);
    w_rep_->scale( 1.0/dt );

    u_minus_w_rep_->update( -1.0, *w_rep_, 1.0 );

    // Selbiges fuer den Druck
    MultiVectorConstPtr_Type pressureSolution = this->solution_->getBlock(1);
    p_rep_->importFromVector(pressureSolution, true);


    // ###############
    // Neu-Assemblierung zu Beginn der neuen Zeititeration im Falle von geometrisch explizit,
    // da das Geometrieproblem zu Beginn der neuen Zeititeration geloest wird und sich somit
    // das Gitter danach bewegt.
    // ###############
    if(type == "ForTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (ForTime)" << '\n';

        // Do we need ForTime at this point? see DAESolverInTime, is it only needed for extrapolation?
        if(geometryExplicit_)
        {
            // ACHTUNG: Fluid-Loesung wird hier auf Null gesetzt, wegen initializeVectors().
            // Somit dann auch die problemTimeFluid_ wodurch eine falsche BDF2-RHS entsteht.
            // Rufe im DAESolverInTime deswegen erneut setPartialSolutions() auf.
            
            
            this->problemFluid_->assembleConstantMatrices(); // Die Steifikeitsmatrix wird weiter unten erst genutzt
            
            // Es ist P = P_
            P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
//            counterP++;
//            std::string outNameMeshVelo = "meshVelo" + to_string(counterP) + ".mm";
//            w_rep_->writeMM(outNameMeshVelo);
//            Teuchos::ArrayRCP<SC> values = w_rep_->getDataNonConst(0);
//            for (int i=0; i<values.size()/2; i++) {
//                values[2*i] = i;
//            }
            this->feFactory_->assemblyAdditionalConvection( this->dim_, this->domain_FEType_vec_.at(0), P_, w_rep_, true );
            P_->resumeFill();
            P_->scale(density);
            P_->scale(-1.0);
//            P_->scale(.0);
            P_->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
//            std::string outNameP = "P" + to_string(counterP) + ".mm";
//            std::cout << "write P:"<<std::endl;
//            P_->writeMM( "P.mm" );
            // Druck und Divergenz hinzufuegen
            this->system_->addBlock( this->problemFluid_->system_->getBlock(0,1), 0, 1 );
            this->system_->addBlock( this->problemFluid_->system_->getBlock(1,0), 1, 0 );
            if (this->problemFluid_->system_->blockExists(1,1))
                this->system_->addBlock( this->problemFluid_->system_->getBlock(1,1), 1, 1 );            
        }

        return; // Damit wir nicht den ganzen Rest auch noch machen!
    }


    // ###############
    // Neu-Assemblierung
    // ###############
    // Fluid: A+P+N+W
    MatrixPtr_Type APNW = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

    if(geometryExplicit_)
    {
        
        if(type == "FixedPoint")
        {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "should always be called by the residual and never here.");
        }
        else if(type == "Newton")
        {
            if(this->verbose_){
                std::cout << "-- Reassembly GE (Newton)" << '\n';
                std::cout << "-- Reassembly GE (Newton) ... full reassembly" << '\n';
            }
            
            problemFluid_->reAssemble( "Newton" );
            
            if (materialModel_ != "linear")
                this->problemStructureNonLin_->reAssemble("Newton");
            
        }
    }
    else
    {
        if(type == "FixedPoint")
        {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "should always be called by the residual and never here.");
        }
        else if(type == "Newton")
        {
            if(this->verbose_){
                std::cout << "-- Reassembly GI (Newton)" << '\n';
                std::cout << "-- Reassembly GI (Newton) ... only W" << '\n';
            }

            problemFluid_->reAssemble( "Newton" );
            
            // TODO: Shape
            // domain(0) = domain(4)
            MatrixPtr_Type shapeVelocity = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
            MatrixPtr_Type shapeDiv = Teuchos::rcp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) ); // shape fuer div-Nebenbedingung

            this->feFactory_->assemblyShapeDerivativeVelocity(this->dim_, this->domain_FEType_vec_.at(4), this->domain_FEType_vec_.at(1),
                            shapeVelocity, 4, u_rep_, w_rep_, p_rep_, dt, density, viscosity, true);
            this->feFactory_->assemblyShapeDerivativeDivergence(this->dim_, this->domain_FEType_vec_.at(4), this->domain_FEType_vec_.at(1),
                            shapeDiv, 1, 4, this->getDomain(1)->getMapUnique(), this->getDomain(4)->getMapVecFieldUnique(), u_rep_, true);
            shapeDiv->resumeFill();
            shapeDiv->scale(-1.0);
            shapeDiv->fillComplete(this->getDomain(4)->getMapVecFieldUnique(), this->getDomain(1)->getMapUnique());

            // Shape Reinschreiben
            this->system_->addBlock(shapeVelocity, 0, 4);
            this->system_->addBlock(shapeDiv, 1, 4);
            
            if (materialModel_ != "linear")
                this->problemStructureNonLin_->reAssemble("Newton");

        }
    }

    this->system_->addBlock( problemFluid_->getSystem()->getBlock( 0, 0 ), 0, 0 );
    
    if (materialModel_ != "linear")
        this->system_->addBlock( this->problemStructureNonLin_->getSystem()->getBlock(0,0), 2, 2 );
}

template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions)
{
    std::cout << "FSI<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions)"<< std::endl;
    // double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);
    // // Fluid-Dichte
    // double density = this->problemFluid_->getParameterList()->sublist("Parameter").get("Density",1.);

    // // ###############
    // // w bestimmen
    // // ###############
    // MultiVectorConstPtr_Type geometrySolution;
    // if(geometryExplicit_)
    // {
    //     geometrySolution = this->problemGeometry_->getSolution()->getBlock(0);
    // }
    // else
    // {
    //     geometrySolution = this->solution_->getBlock(4);
    // }
    // meshDisplacementNew_rep_->importFromVector(geometrySolution, true);


    // *w_rep_ = *meshDisplacementNew_rep_;
    // // update(): this = alpha*A + beta*this
    // w_rep_->update(-1.0, *meshDisplacementOld_rep_, 1.0);
    // w_rep_->scale(1.0/dt);


    // // ###############
    // // u extrapolieren
    // // ###############
    // // Beachte: getBlock(0) ist hier der Richtige, da dies der u-Variable in FSI entspricht.
    // if (previousSolutions.size() >= 2)
    // {
    //     MultiVectorPtr_Type extrapolatedVector = Teuchos::rcp( new MultiVector_Type( previousSolutions[0]->getBlock(0) ) );

    //     // Extrapolation fuer BDF2
    //     // update(): this = alpha*A + beta*this
    //     extrapolatedVector->update( -1., *previousSolutions[1]->getBlock(0), 2. );

    //     u_rep_->importFromVector(extrapolatedVector, true);
    // }
    // else if(previousSolutions.size() == 1)
    // {
    //     MultiVectorConstPtr_Type u = previousSolutions[0]->getBlock(0);
    //     u_rep_->importFromVector(u, true);
    // }
    // else if (previousSolutions.size() == 0)
    // {
    //     MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
    //     u_rep_->importFromVector(u, true);
    // }

    // // (u-w)
    // *u_minus_w_rep_ = *u_rep_;
    // // update(): this = alpha*A + beta*this
    // u_minus_w_rep_->update(-1.0, *w_rep_, 1.0);


    // // ###############
    // // Neu assemblieren
    // // ###############

    // MatrixPtr_Type APN = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

    // MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    // this->feFactory_->assemblyAdvectionVecField( this->dim_, this->domain_FEType_vec_.at(0), N, u_minus_w_rep_, true );

    // P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    
    // this->feFactory_->assemblyAdditionalConvection( this->dim_, this->domain_FEType_vec_.at(0), P_, w_rep_, true );
    // P_->resumeFill();
    // P_->scale(density);
    // P_->scale(-1.0);
    // P_->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());


    // N->resumeFill();
    // N->scale(density);
    // N->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());

    // this->problemFluid_->A_->addMatrix(1.0, APN, 0.0);
    // P_->addMatrix(1.0, APN, 1.0);
    // N->addMatrix(1.0, APN, 1.0);

    // APN->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );

    // this->system_->addBlock( APN, 0, 0 );
}

// Residual of FSI is defined as:
// GI (Geometry Implicit)
//   |b_f     |   | F(u,p,d_f)  0       C_1^T   0 |
//R= |b_s     | - | 0           S       C_3^T   0 |
//   |C_2*d_s |   | C_1         C_2     0       0 |   
//   |0       |   | 0           C_4     0       G |
//   
// GE (Geometry Implicit)
//   |b_f     |   | F(u,p,d_f)  0       C_1^T   0 |
//R= |b_s     | - | 0           S       C_3^T   0 |
//   |C_2*d_s |   | C_1         C_2     0       0 |   

template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const
{
    
    MultiVectorConstPtr_Type geometrySolution;
    
    if(geometryExplicit_)
        geometrySolution = this->problemGeometry_->getSolution()->getBlock(0);
    else {
        moveMesh();
        geometrySolution = this->solution_->getBlock(4);
    }
    
    meshDisplacementNew_rep_->importFromVector(geometrySolution, true);
    
    MultiVectorConstPtr_Type fluidSolution = this->solution_->getBlock(0);
    u_rep_->importFromVector(fluidSolution, true);
    u_minus_w_rep_->importFromVector(fluidSolution, true); //    u_minus_w_rep_ = *u_rep_;
    
    *w_rep_ = *meshDisplacementNew_rep_;
    w_rep_->update(-1.0, *meshDisplacementOld_rep_, 1.0);
    double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);
    w_rep_->scale(1.0/dt);
    
    u_minus_w_rep_->update(-1.0, *w_rep_, 1.0);
    
    if (!geometryExplicit_) {
        
        P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        double density = this->problemTimeFluid_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
        
        this->feFactory_->assemblyAdditionalConvection( this->dim_, this->domain_FEType_vec_.at(0), P_, w_rep_, true );
        P_->resumeFill();
        P_->scale(density);
        P_->scale(-1.0);
        P_->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
        
        problemFluid_->assembleConstantMatrices();
        
        this->system_->addBlock( this->problemFluid_->system_->getBlock(0,1), 0, 1 );
        this->system_->addBlock( this->problemFluid_->system_->getBlock(1,0), 1, 0 );
        if (this->problemFluid_->system_->blockExists(1,1))
                this->system_->addBlock( this->problemFluid_->system_->getBlock(1,1), 1, 1 ); 
        // TEUCHOS_TEST_FOR_EXCEPTION(this->problemFluid_->system_->blockExists(1,1) , std::runtime_error, "Stabilization is used. Account for it."); ???
    }
    if ( this->verbose_ )
        std::cout << "Warning: Wrong consideration of temporal discretization for multi-stage RK methods!" << std::endl;
    
    this->problemFluid_->calculateNonLinResidualVecWithMeshVelo( "reverse", time, u_minus_w_rep_, P_ );
    this->system_->addBlock( problemFluid_->getSystem()->getBlock( 0, 0 ), 0, 0 );
   
    // we need to account for the coupling in the residuals
    if (materialModel_!="linear"){
        this->problemStructureNonLin_->calculateNonLinResidualVec( "reverse", time );
        this->residualVec_->addBlock( this->problemStructureNonLin_->getResidualVector()->getBlockNonConst(0) , 2);
        // we need to add a possible source term
    }
    else{
        MultiVectorPtr_Type residualSolidFSI =
            Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(2) );
        this->problemStructure_->getSystem()->getBlock(0,0)->apply( *this->problemStructure_->getSolution()->getBlock(0), *residualSolidFSI, Teuchos::NO_TRANS, -1. ); // y= -Ax + 0*y
        MultiVectorPtr_Type resSolidNonConst = Teuchos::rcp_const_cast<MultiVector_Type> ( this->residualVec_->getBlock(2) );
        resSolidNonConst->update(1., *this->problemStructure_->getRhs()->getBlock(0), 1.);
        // we need to add a possible source term
    }
    
    MultiVectorPtr_Type residualFluidVelocityFSI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(0) ); // residual of fluid part

    MultiVectorPtr_Type residualSolidFSI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(2) ); // residual of solid part

    MultiVectorPtr_Type residualCouplingFSI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(3) ); // residual of interface / lambda
    residualCouplingFSI->update( 1. , *this->rhs_->getBlock(3), 0. ); // change to -1 for standard
    
    //Now we need to add the coupling blocks
    // adding C_1^T * \lambda to fluid residual
    this->system_->getBlock(0,3)->apply( *this->solution_->getBlock(3) , *residualFluidVelocityFSI, Teuchos::NO_TRANS, -1., 1. ); 
    
    // adding C_3 ^T * \lambda to solid residual
    this->system_->getBlock(2,3)->apply( *this->solution_->getBlock(3) , *residualSolidFSI, Teuchos::NO_TRANS, -1., 1. ); 
     
     // adding C_1 * u to coupling/lambda 
    this->system_->getBlock(3,0)->apply( *this->solution_->getBlock(0) , *residualCouplingFSI, Teuchos::NO_TRANS, -1., 1. ); 
    
    // adding C_2 * d_s to coupling/lambda 
    this->system_->getBlock(3,2)->apply( *this->solution_->getBlock(2) , *residualCouplingFSI, Teuchos::NO_TRANS, -1., 1. );

    if (!geometryExplicit_) {
        // If we have GI coupling we need to add that component to residual as well
        MultiVectorPtr_Type residualGeometryFSI =
            Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(4) );

        residualGeometryFSI->update( 1. , *this->rhs_->getBlock(4), 0. ); // change to -1 for standard
        // G * d_f
        this->system_->getBlock(4,4)->apply( *this->solution_->getBlock(4) , *residualGeometryFSI, Teuchos::NO_TRANS, -1., 1. );
        // C_4 * d_f
        this->system_->getBlock(4,2)->apply( *this->solution_->getBlock(2) , *residualGeometryFSI, Teuchos::NO_TRANS, -1., 1. );
        
    }
    // might also be called in the sub calculateNonLinResidualVec() methods which where used above
    if (type == "reverse")
        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    else if (type == "standard"){
        this->residualVec_->scale(-1.);
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
    }


}

template <class SC, class LO, class GO, class NO>
void FSI<SC, LO, GO, NO>::calculateNonLinResidualVec(SmallMatrix<double> &coeff, std::string type, double time, BlockMatrixPtr_Type systemMass)
{
    NonLinearProblem<SC,LO,GO,NO>::calculateNonLinResidualVec(coeff,type,time);

    // for FSI we need to reassemble the massmatrix system if the mesh was moved for geometry implicit computations
    bool geometryExplicit = this->parameterList_->sublist("Parameter").get("Geometry Explicit",true);
    if( !geometryExplicit ) {
        MatrixPtr_Type massmatrix;
        this->setFluidMassmatrix( massmatrix );
        systemMass->addBlock( massmatrix, 0, 0 );
    }
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
void FSI<SC,LO,GO,NO>::setFromPartialVectorsInit() const
{
    
    //Fluid velocity
    this->solution_->addBlock( this->problemFluid_->getSolution()->getBlockNonConst(0), 0 );
    this->residualVec_->addBlock( this->problemFluid_->getResidualVector()->getBlockNonConst(0), 0 );
    this->residualVec_->addBlock( this->problemFluid_->getResidualVector()->getBlockNonConst(0), 0 );
    this->rhs_->addBlock( this->problemFluid_->getRhs()->getBlockNonConst(0), 0 );
    this->sourceTerm_->addBlock( this->problemFluid_->getSourceTerm()->getBlockNonConst(0), 0 );
    
    //Fluid pressure
    this->solution_->addBlock( this->problemFluid_->getSolution()->getBlockNonConst(1), 1 );
    this->residualVec_->addBlock( this->problemFluid_->getResidualVector()->getBlockNonConst(1), 1 );
    this->rhs_->addBlock( this->problemFluid_->getRhs()->getBlockNonConst(1), 1 );
    this->previousSolution_->addBlock( this->problemFluid_->getPreviousSolution()->getBlockNonConst(1), 1 );
    this->sourceTerm_->addBlock( this->problemFluid_->getSourceTerm()->getBlockNonConst(1), 1 );
    
    if (materialModel_=="linear"){
        this->solution_->addBlock( this->problemStructure_->getSolution()->getBlockNonConst(0), 2 );
        // we dont have a residual vector for linear problems
        this->rhs_->addBlock( this->problemStructure_->getRhs()->getBlockNonConst(0), 2 );
        this->sourceTerm_->addBlock( this->problemStructure_->getSourceTerm()->getBlockNonConst(0), 2 );
    }
    else{
        this->solution_->addBlock( this->problemStructureNonLin_->getSolution()->getBlockNonConst(0), 2 );
        this->residualVec_->addBlock( this->problemStructureNonLin_->getResidualVector()->getBlockNonConst(0), 2 );
        this->rhs_->addBlock( this->problemStructureNonLin_->getRhs()->getBlockNonConst(0), 2 );
        this->previousSolution_->addBlock( this->problemStructureNonLin_->getPreviousSolution()->getBlockNonConst(0), 2 );
        this->sourceTerm_->addBlock( this->problemStructureNonLin_->getSourceTerm()->getBlockNonConst(0), 2 );
    }
    
    if(!geometryExplicit_){
        this->solution_->addBlock( this->problemGeometry_->getSolution()->getBlockNonConst(0), 4 );
        // we dont have a previous solution for linear problems
        this->rhs_->addBlock( this->problemGeometry_->getRhs()->getBlockNonConst(0), 4 );
        this->sourceTerm_->addBlock( this->problemGeometry_->getSourceTerm()->getBlockNonConst(0), 4 );
    }
}
    
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::updateMeshDisplacement() const
{

     *meshDisplacementOld_rep_ = *meshDisplacementNew_rep_;

}


template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::solveGeometryProblem() const
{
    DomainPtr_Type domainFluid = Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(0));
    DomainPtr_Type domainStruct = Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(2));

    // Hole die dof-Map
    MapConstPtr_Type interfaceMapFluidVecField = domainFluid->getInterfaceMapVecFieldUnique();
    MapConstPtr_Type interfaceMapStructureVecField = domainStruct->getInterfaceMapVecFieldUnique();
    
    MapConstPtr_Type interfaceMapGlobalFluidVecField = domainFluid->getGlobalInterfaceMapVecFieldUnique();
    MapConstPtr_Type interfaceMapGlobalStructureVecField = domainStruct->getGlobalInterfaceMapVecFieldUnique();

    MeshUnstrPtr_Type meshUnstructuredFluid = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainFluid->getMesh() );
    vec3D_GO_ptr_Type indicesGlobalMatchedOriginFluid = meshUnstructuredFluid->getMeshInterface()->getIndicesGlobalMatchedOrigin();

    // Strukturloesung holen
    MultiVectorConstPtr_Type struc_sol_unique = this->solution_->getBlock(2);

    // Extrahiere nun aus der globalen Strukturloesung die Loesung auf dem Interface (in parallel).
    MultiVectorPtr_Type interfaceSolutionStruct = Teuchos::rcp( new MultiVector_Type( interfaceMapStructureVecField, 1 ) );

    {
        int flagCounter = 0;
        Teuchos::ArrayRCP< SC > valuesInterface = interfaceSolutionStruct->getDataNonConst(0); //single MultiVector
        Teuchos::ArrayRCP< SC > valuesStructure = struc_sol_unique->getDataNonConst(0); //single MultiVector

        for(UN i = 0; i < valuesInterface.size(); i++)
        {
            GO interfaceIDGlobal = interfaceMapGlobalStructureVecField->getGlobalElement( i ); // ID (vektorwertig, also dofID) in der Interface-Nummerierung
//            GO nodeID; // nodeID der interfaceID in der Interface-Nummerierung
//            LO localDofNumber; // Ranges from 0 to dofs-1
//            toNodeID(this->dim_, interfaceID, nodeID, localDofNumber); //This function assumes NodeWise ordering.
//
//            // Ggf. ist die ID auf einer anderen Flag.
//            // Bei nur einer Interface-Flag kommt man nicht hier hinein.
//            while( nodeID > indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size()-1 )
//            {
//                nodeID = nodeID - indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size();
//                flagCounter = flagCounter + 1; // hier steht dann die korrekte Flag
//            }
//
//            // Beachte: at(1) ist Struktur!!!
//            // GlobaleID des Interface-Knotens in der Struktur-Nummerierung
//            GO globalInterfaceIDNode = indicesGlobalMatchedOriginFluid->at(flagCounter).at(1).at(nodeID);
//            GO globalInterfaceIDinStructure; // dofID
//            toDofID(this->dim_, globalInterfaceIDNode, localDofNumber, globalInterfaceIDinStructure);

            // LokaleID auf dem Prozessor vom Interface-Knoten.
            LO localInterfaceIDinStructure = domainStruct->getMapVecFieldUnique()->getLocalElement(interfaceIDGlobal);

            valuesInterface[i] = valuesStructure[localInterfaceIDinStructure];
        }
    }

    MultiVectorPtr_Type interfaceSolutionFluid = Teuchos::rcp( new MultiVector_Type( interfaceMapFluidVecField, 1 ) );
    interfaceSolutionFluid->importFromVector( interfaceSolutionStruct );

    {
        // Strukturlsg als RW setzen (per Hand)
        this->problemGeometry_->setBoundariesSystem(); // RW im System setzen (mit den Flags von oben)
        this->problemGeometry_->getRhs()->putScalar(0.0);

        Teuchos::ArrayRCP< SC > valuesInterface = interfaceSolutionFluid->getDataNonConst(0); //single MultiVector
        Teuchos::ArrayRCP< SC > valuesFluidRhs = this->problemGeometry_->getRhs()->getBlock(0)->getDataNonConst(0); //single MultiVector

        int flagCounter = 0;
        for(UN i = 0; i < valuesInterface.size(); i++)
        {
            GO interfaceIDGlobal = interfaceMapGlobalFluidVecField->getGlobalElement( i );
//            GO interfaceID = interfaceMapFluidVecField->getGlobalElement( i ); // dofID in der Interface-Nummerierung
//            GO nodeID; // dofID der interfaceID in der Interface-Nummerierung
//            LO localDofNumber; // Ranges from 0 to dofs-1
//            toNodeID(this->dim_, interfaceID, nodeID, localDofNumber);//This function assumes NodeWise ordering.
//
//            // Ggf. ist die ID auf einer anderen Flag.
//            // Bei nur einer Interface-Flag kommt man nicht hier hinein.
//            while(nodeID > indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size()-1)
//            {
//                nodeID = nodeID - indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size();
//                flagCounter = flagCounter + 1; // hier steht dann die korrekte Flag
//            }
//
//            // Beachte: at(0) ist Fluid!!!
//            // GlobaleID des Interface-Knotens in der Fluid-Nummerierung
//            GO globalInterfaceIDNode = indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).at(nodeID);
//            GO globalInterfaceIDinFluid; // dofID
//            toDofID(this->dim_, globalInterfaceIDNode, localDofNumber, globalInterfaceIDinFluid);

            // LokaleID auf dem Prozessor des Interface-Knotens.
            LO localInterfaceIDinFluid = domainFluid->getMapVecFieldUnique()->getLocalElement(interfaceIDGlobal);

            valuesFluidRhs[localInterfaceIDinFluid] = valuesInterface[i];

            // valuesFluidRhs[localInterfaceIDinFluid.at(i)] = valuesInterface[i];
        }

        this->problemGeometry_->solve();
        
        if (!exporterGeo_.is_null())
            this->exporterGeo_->save( this->timeSteppingTool_->currentTime() );
        

    }
    
    

}


template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::setupSubTimeProblems(ParameterListPtr_Type parameterListFluid, ParameterListPtr_Type parameterListStructure) const
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
void FSI<SC,LO,GO,NO>::setFluidMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix fuer FSI combineSystems(), ggf nichtlinear.
    //######################
    double density = this->problemTimeFluid_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
    int size = this->problemTimeFluid_->getSystem()->size();

    this->problemTimeFluid_->systemMass_.reset(new BlockMatrix_Type(size));
    {
        massmatrix = Teuchos::rcp(new Matrix_Type( this->problemTimeFluid_->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
        // 0 = Fluid
        this->feFactory_->assemblyMass( this->dim_, this->problemTimeFluid_->getFEType(0), "Vector",  massmatrix, 0, true );
        massmatrix->resumeFill();
        massmatrix->scale(density);
        massmatrix->fillComplete( this->problemTimeFluid_->getDomain(0)->getMapVecFieldUnique(), this->problemTimeFluid_->getDomain(0)->getMapVecFieldUnique() );

        this->problemTimeFluid_->systemMass_->addBlock(massmatrix, 0, 0);
    }
}


// Function to compute the pressure boundary conditions for the fluid component
// Note it already contains some code connected to restarts
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::computePressureRHSInTime() const{

    // Type of pressure boundary condition
    std::string pressureRB = this->parameterList_->sublist("Parameter Fluid").get("Pressure Boundary Condition","None");

    // bool restart = this->parameterList_->sublist("Timestepping Parameter").get("Restart", false);
    // double timeStepRestart = this->parameterList_->sublist("Timestepping Parameter").get("Time step", 0.0); 
    
    // Resistance boundary condition based on 'A parallel two-level method for simulating blood ﬂows in branching arteries
    // with the resistive boundary condition -- Wu, Cai 2011'
    // We assemble
    // \int_{\Gamma_O} R flowrateOutlet \phi_f n ds + nu_f \int_{\Omega_O} \phi_f \cdot (\nabla u_f) \cdot n ds
    if (pressureRB == "Resistance")
    {
        if(this->verbose_)
            std::cout << " --- Computing resistance boundary condition .. " << std::endl;

        // Value added to the RHS of the fluid component
        MultiVectorPtr_Type FERhs = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ));

        // Parameters for pressure ramp
        vec_dbl_Type funcParameter(1,0.);
        funcParameter[0] = timeSteppingTool_->t_;            
       
        // We need the outlet are thus we need to know the flags, as they are not always the same
        int flagInlet =this->parameterList_->sublist("General").get("Flag Inlet Fluid", 4);
        int flagOutlet = this->parameterList_->sublist("General").get("Flag Outlet Fluid", 5);

        // We compute the initial flow rate in the first timestep
        if (timeSteppingTool_->currentTime()==0.) { 
            double flowRateInlet_n_1 = 0.;
            this->feFactory_->assemblyFlowRate(this->dim_, flowRateInlet_n_1, this->getDomain(0)->getFEType() , this->dim_, flagOutlet , u_rep_);  
            flowRateOutlet_n_1_ = flowRateInlet_n_1; 
        }  
        // else if(restart && timeStepRestart +1e-8 > timeSteppingTool_->currentTime() )
        // {
        //     if(this->verbose_)
        //         cout << " WARNING: Absorbing boundary condition is computed but the initial values usally corresponding to T=0 now correspond to the restart time " << endl;
            
        //     double flowRateInlet_n_1 = 0.;
        //     this->feFactory_->assemblyFlowRate(this->dim_, flowRateInlet_n_1, this->getDomain(0)->getFEType() , this->dim_, flagOutlet , u_rep_);  
        //     flowRateOutlet_n_1_ = flowRateInlet_n_1; 

        // }

        // Then we compute the flowrate in the current time step
        double flowRateInlet_n = 0.;
        this->feFactory_->assemblyFlowRate(this->dim_, flowRateInlet_n, this->getDomain(0)->getFEType() , this->dim_, flagOutlet , u_rep_);  
        flowRateOutlet_n_ = flowRateInlet_n; 

        vec_dbl_Type flowRateOutlet_timesteps(2);
        flowRateOutlet_timesteps[0] = flowRateOutlet_n_1_;
        flowRateOutlet_timesteps[1] = flowRateOutlet_n_;        

        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true); 
         
        // We compute the pressure value set on the outlet of the domain based on the resitive boundary
        pressureOutlet_ = this->feFactory_->assemblyResistanceBoundary(this->dim_, this->getDomain(0)->getFEType(),FERhs, u_rep_,flowRateOutlet_timesteps, funcParameter, pressureRamp,this->parameterList_,0);
                  
        flowRateOutlet_n_1_ = flowRateOutlet_n_;

        this->sourceTerm_->getBlockNonConst(0)->exportFromVector( FERhs, false, "Add" );

        //double density = this->parameterList_->sublist("Parameter").get("Density",1.);
        //this->problemTimeFluid_->getSourceTerm()->scale(density);
        // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
        // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!
        
        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
       // BlockMultiVectorPtr_Type tmpSourceterm = Teuchos::rcp(new BlockMultiVector_Type(1)) ;
       // tmpSourceterm->addBlock(this->sourceTerm_->getBlockNonConst(1),0);

            
        this->problemTimeFluid_->getRhs()->getBlockNonConst(0)->update(coeffSourceTermStructure, *this->sourceTerm_->getBlockNonConst(0), 1.);
        this->rhs_->addBlock( this->problemTimeFluid_->getRhs()->getBlockNonConst(0), 0 );

        if(this->verbose_)
            std::cout << "  .. done " << std::endl;

    
    }
    // Absorbing boundary condition implemented as introduced in 
    // 'AN EFFECTIVE FLUID-STRUCTURE INTERACTION FORMULATION FOR VASCULAR DYNAMICS BY GENERALIZED
    // ROBIN CONDITIONS, Nobile, Vergara, 2008'
    if (pressureRB == "Absorbing" || pressureRB == "Absorbing Paper")
    {
        if(this->verbose_)
            std::cout << " Computing absorbing boundary condition .. " << std::endl;

        MultiVectorPtr_Type FERhs = Teuchos::rcp(new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ));

        // Parameters for pressure ramp
        vec_dbl_Type funcParameter(1,0.);
        funcParameter[0] = timeSteppingTool_->t_;            
        
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
        u_rep_->importFromVector(u, true); 
         
        int flagInlet =this->parameterList_->sublist("General").get("Flag Inlet Fluid", 4);
        int flagOutlet = this->parameterList_->sublist("General").get("Flag Outlet Fluid", 5);

        if (timeSteppingTool_->currentTime()==0.) { 
            double areaInlet_init = 0.;
            double areaOutlet_init = 0.;

            this->feFactory_->assemblyArea(this->dim_,areaInlet_init, flagInlet);
            this->feFactory_->assemblyArea(this->dim_, areaOutlet_init, flagOutlet);

            areaInlet_init_ = areaInlet_init;
            areaOutlet_init_ = areaOutlet_init;

            double flowRateInlet_n_1 = 0.;
            this->feFactory_->assemblyFlowRate(this->dim_, flowRateInlet_n_1, this->getDomain(0)->getFEType() , this->dim_, flagOutlet , u_rep_);  
            flowRateOutlet_n_1_ = flowRateInlet_n_1; 

            if ( this->parameterList_->sublist("Timestepping Parameter").get("Checkpointing", false)){
                exporterBoundaryCondition_->exportData( "Area_Inlet ", areaInlet_init );
                exporterBoundaryCondition_->exportData( "Area_Outlet ", areaInlet_init );
            }
        } 
        // else if(restart && timeStepRestart +1e-8 > timeSteppingTool_->currentTime() )
        // {
        //     if(this->verbose_)
        //         cout << " WARNING: Absorbing boundary condition is computed but the initial values usally corresponding to T=0 now correspond to the restart time " << endl;
        //     double areaInlet_init = 0.;
        //     double areaOutlet_init = 0.;

        //     this->feFactory_->assemblyArea(this->dim_,areaInlet_init, flagInlet);
        //     this->feFactory_->assemblyArea(this->dim_, areaOutlet_init, flagOutlet);

        //     areaInlet_init_ = 0.0253605;//areaInlet_init;
        //     areaOutlet_init_ = 0.025605; //areaOutlet_init;

        //     double flowRateInlet_n_1 = 0.;
        //     this->feFactory_->assemblyFlowRate(this->dim_, flowRateInlet_n_1, this->getDomain(0)->getFEType() , this->dim_, flagOutlet , u_rep_);  
        //     flowRateOutlet_n_1_ = flowRateInlet_n_1;  
        // }
        double flowRateInlet_n = 0.;
        this->feFactory_->assemblyFlowRate(this->dim_, flowRateInlet_n, this->getDomain(0)->getFEType() , this->dim_, flagOutlet , u_rep_);  
        flowRateOutlet_n_ = flowRateInlet_n; 

        vec_dbl_Type flowRateOutlet_timesteps(2);
        flowRateOutlet_timesteps[0] = flowRateOutlet_n_1_;
        flowRateOutlet_timesteps[1] = flowRateOutlet_n_;

        // The traditional approach differs slightly from the approach used in Comparison of arterial wall models in fluid–structure interaction
        // simulations D. Balzani, A. Heinlein, A. Klawonn, O. Rheinbach, J. Schröder, 2023, and its predecessor. Thus, absorbing paper refers to 
        // the aforementioned paper.
        if(pressureRB == "Absorbing Paper"){

            double unsteadyStart = this->parameterList_->sublist("Parameter Fluid").get("Unsteady Start",0.1); 
            if( unsteadyStart +1e-10 > timeSteppingTool_->currentTime() &&  unsteadyStart -1e-10 < timeSteppingTool_->currentTime() )
            {
                double areaOutlet_T = 0.;
                this->feFactory_->assemblyArea(this->dim_, areaOutlet_T, flagOutlet);
                areaOutlet_T_ = areaOutlet_T;
                if(this->verbose_)
                    std::cout << " ---- Absorbing boundary condition: Start of unsteady Phase with areaOutlet_T=" << areaOutlet_T_<< " ---- " << std::endl;
            }

            pressureOutlet_ = this->feFactory_->assemblyAbsorbingBoundaryPaper(this->dim_, this->getDomain(0)->getFEType(),FERhs, u_rep_,flowRateOutlet_timesteps, funcParameter, pressureRamp,areaOutlet_init_, areaOutlet_T_,this->parameterList_,0);

        }
        else{
            double rampTime = this->parameterList_->sublist("Parameter").get("Max Ramp Time",0.1); 
            if( rampTime < timeSteppingTool_->currentTime() &&  rampTime + 0.05 > timeSteppingTool_->currentTime() )
            {
                double areaOutlet_T = 0.;
                this->feFactory_->assemblyArea(this->dim_, areaOutlet_T, flagOutlet);
                areaOutlet_T_ = areaOutlet_T;
                if(this->verbose_)
                    std::cout << " ---- Absorbing boundary condition: Start of unsteady Phase with areaOutlet_T=" << areaOutlet_T_<< " ---- " << std::endl;
            }
            pressureOutlet_ = this->feFactory_->assemblyAbsorbingBoundary(this->dim_, this->getDomain(0)->getFEType(),FERhs, u_rep_,flowRateOutlet_timesteps, funcParameter, pressureRamp,areaOutlet_init_,areaOutlet_T_,this->parameterList_,0);
       
        }

        // Adding the assembled RHS on the Fluid component to the fluid RHS
        this->sourceTerm_->getBlockNonConst(0)->exportFromVector( FERhs, false, "Add" );
        
        flowRateOutlet_n_1_ = flowRateOutlet_n_;
        if ( this->parameterList_->sublist("Timestepping Parameter").get("Checkpointing", false)){
            exporterBoundaryCondition_->exportData( "FlowrateOutlet_Previous_Timestep", flowRateOutlet_n_1_ );
        }
        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
       
            
        this->problemTimeFluid_->getRhs()->getBlockNonConst(0)->update(coeffSourceTermStructure, *this->sourceTerm_->getBlockNonConst(0), 1.);
        this->rhs_->addBlock( this->problemTimeFluid_->getRhs()->getBlockNonConst(0), 0 );

        if(this->verbose_)
            std::cout << "  .. done " << std::endl;

    
    }    
  
}

// TODO: updateMultistepRhsFSI() einbauen!
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::computeFluidRHSInTime( ) const
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


// Am Anfang der Zeititeration erst updateSolutionMultiPreviousStep() aufrufen und dann erst updateMultistepRhs(),
// damit die previousSolution_ initialisiert sind. Genauso fuer SystemMass
// TODO: updateSystemMassMultiPreviousStep() fertig programmieren
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::updateFluidInTime() const
{
    int nmbBDF = timeSteppingTool_->getBDFNumber();

    if(nmbBDF<2 && !this->parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation")) {
        if (timeSteppingTool_->currentTime()!=0.){
            this->problemTimeFluid_->updateSolutionMultiPreviousStep(2);
            this->problemTimeFluid_->updateSystemMassMultiPreviousStep(2);
        }
        else{
            this->problemTimeFluid_->updateSolutionMultiPreviousStep(1);
            this->problemTimeFluid_->updateSystemMassMultiPreviousStep(1);
        }
    }
    else{
        this->problemTimeFluid_->updateSolutionMultiPreviousStep(nmbBDF);
        this->problemTimeFluid_->updateSystemMassMultiPreviousStep(nmbBDF);
    }
}

template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::computeSolidRHSInTime() const {
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
void FSI<SC,LO,GO,NO>::setSolidMassmatrix( MatrixPtr_Type& massmatrix ) const
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
void FSI<SC,LO,GO,NO>::updateTime() const
{
    timeSteppingTool_->t_ = timeSteppingTool_->t_ + timeSteppingTool_->dt_prev_;
}

template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::moveMesh() const
{

    MultiVectorConstPtr_Type displacementUniqueConst;
    if(geometryExplicit_)
    {
        displacementUniqueConst = this->problemGeometry_->getSolution()->getBlock(0);
    }
    else
    {
        displacementUniqueConst = this->solution_->getBlock(4);
    }
    MultiVectorPtr_Type displacementRepeated = Teuchos::rcp( new MultiVector_Type( this->problemGeometry_->getDomain(0)->getMapVecFieldRepeated() ) );

    displacementRepeated->importFromVector( displacementUniqueConst );
    MultiVectorPtr_Type displacementUnique = Teuchos::rcp_const_cast<MultiVector_Type>(displacementUniqueConst);


    // Verschiebe die Gitter fuer Fluid-Velocity und Fluid-Druck
    // ACHTUNG: Klappt nur, weil die P2-Knoten hinter den P1-Knoten kommen.
    // Sonst muessen fuer den Druck die P1-Knoten extrahiert werden.
    // TODO: Wahrscheinlich reicht nur FSI-Domain, da selbes Objekt bei problemFluid_ und problemTimeFluid_.
    ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemFluid_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemFluid_->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeFluid_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeFluid_->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
}


template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::addInterfaceBlockRHS() const
{
    MultiVectorPtr_Type vectorToAdd = Teuchos::rcp( new MultiVector_Type( this->rhs_->getBlock(3) ) );

    C2_->apply(*(this->solution_->getBlock(2)), *vectorToAdd);
    this->rhs_->addBlock(vectorToAdd, 3);
}


// template<class SC,class LO,class GO,class NO>
// void FSI<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
//                                      const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
//                                     ) const
// {
//     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "implement NOX for steady FSI.");
//     std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
// //    if ( !type.compare("Monolithic"))
// //        evalModelImplMonolithic( inArgs, outArgs );
// //    else if ( !type.compare("FaCSI")){
// //        evalModelImplBlock( inArgs, outArgs );
// //    }
// //    else
// //        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Unkown preconditioner/solver type.");
// }

template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::getValuesOfInterest( vec_dbl_Type& values ){
    if (this->dim_==2)
        getValuesOfInterest2DBenchmark( values );
    else if(this->dim_==3)
        getValuesOfInterest3DBenchmark( values);
    
}
    
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::getValuesOfInterest2DBenchmark( vec_dbl_Type& values ){

    
    if ( valuesForExport_[0] >= 0.  ) { // we set the displacement of the 2D Turek Benchmark in x and y direction
        LO loc = valuesForExport_[0] + 10*Teuchos::ScalarTraits<SC>::eps();
        values[0] = this->getSolution()->getBlock(2)->getData(0)[2*loc];
        values[1] = this->getSolution()->getBlock(2)->getData(0)[2*loc+1];
    }
}
    
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::getValuesOfInterest3DBenchmark( vec_dbl_Type& values ){
    
    if ( valuesForExport_[0] >= 0.   ) { // we set the displacement of the 3D Richter Benchmark in x, y, and z direction
        LO loc = valuesForExport_[0] + 10*Teuchos::ScalarTraits<SC>::eps();
        values[0] = this->getSolution()->getBlock(2)->getData(0)[3*loc];
        values[1] = this->getSolution()->getBlock(2)->getData(0)[3*loc+1];
        values[2] = this->getSolution()->getBlock(2)->getData(0)[3*loc+2];
    }
}
    
template<class SC,class LO,class GO,class NO>
void FSI<SC,LO,GO,NO>::computeValuesOfInterestAndExport(){
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
void FSI<SC,LO,GO,NO>::initializeGE(){
/*All vectors (solution, rhs, previousSolution,...) should be initialized at this point (initializeProblem() was called)
 Therefore, all these BlockMVectors should have a length of 5 since the geometry problem is included in the general setup. We need to resize these BlockMVs here if the Geometry Explicit (GE) system is used.
*/
    if (geometryExplicit_) {
        this->solution_->resize( 4 );
        this->rhs_->resize( 4 );
        this->sourceTerm_->resize( 4 );
        this->rhsFuncVec_.resize( 4 );
        this->previousSolution_->resize( 4 );
        this->residualVec_->resize( 4 );
        this->initVectorSpaces();  //reinitialize NOX vector spaces
    }
}

}
#endif
