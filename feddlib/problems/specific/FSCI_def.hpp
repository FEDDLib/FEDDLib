#ifndef FSCI_def_hpp
#define FSCI_def_hpp
#include "FSCI_decl.hpp"

/*double OneFunc(double* x, int* parameter)
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
*/
namespace FEDD {
// Funktionen fuer die rechte Seite der Struktur/ Fluid/ Geometrie sind im jeweiligen Problem

template<class SC,class LO,class GO,class NO>
FSCI<SC,LO,GO,NO>::FSCI(const DomainConstPtr_Type &domainVelocity, std::string FETypeVelocity,
                    const DomainConstPtr_Type &domainPressure, std::string FETypePressure,
                    const DomainConstPtr_Type &domainStructure, std::string FETypeStructure,
                    const DomainConstPtr_Type &domainChem, std::string FETypeChem,
                    const DomainConstPtr_Type &domainInterface, std::string FETypeInterface,
                    const DomainConstPtr_Type &domainGeometry, std::string FETypeGeometry,
                    vec2D_dbl_Type diffusionTensor, RhsFunc_Type reactionFunc,
                    ParameterListPtr_Type parameterListFluid, ParameterListPtr_Type parameterListStructure,
                    ParameterListPtr_Type parameterListChem,
                    ParameterListPtr_Type parameterListFSCI, ParameterListPtr_Type parameterListGeometry,
                    Teuchos::RCP<SmallMatrix<int> > &defTS):
FSI<SC,LO,GO,NO>(domainVelocity,  FETypeVelocity,
                domainPressure,FETypePressure,
                domainStructure, FETypeStructure,
                domainInterface, FETypeInterface,
                domainGeometry, FETypeGeometry,
                parameterListFluid, parameterListStructure,
                parameterListFSCI, parameterListGeometry,defTS),
//NonLinearProblem<SC,LO,GO,NO>( parameterListFSCI, domainVelocity->getComm() ),
// hasSourceTerm = drittes Arguement. assembleSourceTerm() fuer NS nicht programmiert.
// Deswegen hier erstmal false (default Parameter).
// Fuer Struktur hingegen ist default Parameter true, da programmiert.
//P_(),
//problemFluid_(),
problemSCI_()
/*problemGeometry_(),
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
exporterGeo_()*/
{
    /*this->nonLinearTolerance_ = this->parameterList_->sublist("Parameter").get("relNonLinTol",1.0e-6);
    geometryExplicit_ = this->parameterList_->sublist("Parameter").get("Geometry Explicit",true);

    this->initNOXParameters();
    counterP = 0;
    
    std::string linearization = parameterListFSCI->sublist("General").get("Linearization","FixedPoint");
    
    TEUCHOS_TEST_FOR_EXCEPTION( !(linearization == "Newton" || linearization == "NOX")  && materialModel_ != "linear", std::runtime_error, "Nonlinear material models can only be used with Newton's method or FixedPoint (nonlinear material Jacobian will still be used).");
    
    this->addVariable( domainVelocity, FETypeVelocity, "u_f", domainVelocity->getDimension() ); // Fluid-Geschw.
    this->addVariable( domainPressure, FETypePressure, "p", 1); // Fluid-Druck
    this->addVariable( domainStructure, FETypeStructure, "d_s", domainStructure->getDimension() ); */// Struktur
    this->addVariable( domainChem, FETypeChem, "c", 1 ); // Chem
    /*this->addVariable( domainInterface, FETypeInterface, "lambda", domainInterface->getDimension() ); // Interface
    this->addVariable( domainGeometry, FETypeGeometry, "d_f", domainGeometry->getDimension() ); // Geometrie

    this->dim_ = this->getDomain(0)->getDimension();
    
    problemFluid_ = Teuchos::rcp( new FluidProblem_Type( domainVelocity, FETypeVelocity, domainPressure, FETypePressure, parameterListFluid ) );
    problemFluid_->initializeProblem();*/
    
    Teuchos::RCP<SmallMatrix<int>> defTSSCI; // Seperate Timestepping Matrix for SCI
    defTSSCI.reset( new SmallMatrix<int> (2) );
    (*defTSSCI)[0][0] = (*defTS)[2][2];
    (*defTSSCI)[1][1] = (*defTS)[3][3];
    this->problemSCI_ = Teuchos::rcp( new SCIProblem_Type( domainStructure, FETypeStructure, domainChem, FETypeChem, diffusionTensor,reactionFunc, parameterListStructure, parameterListChem, parameterListFSCI, defTSSCI ) );
    this->problemSCI_->initializeProblem();

    /*if (materialModel_=="linear"){
        problemStructure_ = Teuchos::rcp( new StructureProblem_Type( domainStructure, FETypeStructure, parameterListStructure ) );
        problemStructure_->initializeProblem();
    }
    else{
        problemStructureNonLin_ = Teuchos::rcp( new StructureNonLinProblem_Type( domainStructure, FETypeStructure, parameterListStructure) );
        problemStructureNonLin_->initializeProblem();
    }*/
    
    /*problemGeometry_ = Teuchos::rcp( new GeometryProblem_Type( domainGeometry, FETypeGeometry, parameterListGeometry ) );
    problemGeometry_->initializeProblem();
    //We initialize the subproblems. In the main routine, we need to call initializeFSCI(). There, we first initialize the vectors of the FSCI problem and then we set the pointers of the subproblems to the vectors of the full monolithic FSCI system. This way all values are only saved once in the subproblems and can be used by the monolithic FSCI system.
    
    meshDisplacementNew_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(5)->getMapVecFieldRepeated() ) );
    meshDisplacementOld_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(5)->getMapVecFieldRepeated() ) );
    u_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    u_minus_w_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(0)->getMapVecFieldRepeated() ) );
    
    // Exporting benchmark values
    if ( parameterListFSCI->sublist("General").get("Export Extra Data",false) ){
        if (this->dim_==2)
            findDisplacementTurek2DBenchmark();
        else if (this->dim_==3)
            findDisplacementRichter3DBenchmark();
    }
    if ( parameterListFSCI->sublist("General").get("Export drag and lift",false) ){
        exporterTxtDrag_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtDrag_->setup( "drag_force", this->comm_ );
        exporterTxtLift_ = Teuchos::rcp(new ExporterTxt () );
        exporterTxtLift_->setup( "lift_force", this->comm_ );
    }
    p_rep_ = Teuchos::rcp( new MultiVector_Type( this->getDomain(1)->getMapRepeated() ) );*/
    
}

    
template<class SC,class LO,class GO,class NO>
FSCI<SC,LO,GO,NO>::~FSCI()
{
    if (!this->exporterGeo_.is_null()) {
       this->exporterGeo_->closeExporter();
    }
}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::info()
{
    this->infoProblem();
    this->infoNonlinProblem();
}

// Alle linearen Probleme
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::assemble( std::string type ) const
{
    if (type == "") {
        if (this->verbose_)
        {
            std::cout << "-- Assembly FSCI ... " << std::endl;
        }

    //    P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), 10 ) );
        
        this->problemFluid_->assemble();
        
        // steady rhs wird hier assembliert.
        // rhsFunc auf 0 (=x) und 0 (=y) abaendern bei LinElas!
       /* if (materialModel_=="linear")
            this->problemStructure_->assemble();
        else
            this->problemStructureNonLin_->assemble();
        */
        this->problemSCI_->assemble();
        
        this->problemGeometry_->assemble();
        if ( this->geometryExplicit_ && this->parameterList_->sublist("Exporter").get("Export GE geometry solution",false)){
            exporterGeo_ = Teuchos::rcp(new Exporter_Type());
            
            DomainConstPtr_Type dom = this->getDomain(5);

            int exportEveryXTimesteps = this->parameterList_->sublist("Exporter").get( "Export every X timesteps", 1 );
            std::string suffix = this->parameterList_->sublist("Exporter").get("Geometry Suffix", "" );
            std::string varName = "d_f" + suffix;
            
            MeshPtr_Type meshNonConst = Teuchos::rcp_const_cast<Mesh_Type>( dom->getMesh() );
            this->exporterGeo_->setup(varName, meshNonConst, dom->getFEType(), exportEveryXTimesteps, this->parameterList_);
            
            MultiVectorConstPtr_Type exportVector = this->problemGeometry_->getSolution()->getBlock(0);
            
            this->exporterGeo_->addVariable( exportVector, varName, "Vector", this->dim_, dom->getMapUnique() );
        }
        
        
        if (!this->geometryExplicit_) {
            // RW fuer die Geometrie-Matrix H setzen, weil wenn wir das spaeter im FSCI-System
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
        MatrixPtr_Type C4(new Matrix_Type( this->getDomain(5)->getMapVecFieldUnique(), 1 ) ); // Geometrie (=Fluid)

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
                                                       this->getDomain(5)->getMapVecFieldUnique(), true);
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
        this->C2_ = C2;

        if(!geometryExplicit_)
        {
            C4->resumeFill();
            C4->scale(-1.0);
            // Domain = Spalten = Struktur; Range = Zeilen = Geometrie
            C4->fillComplete(this->getDomain(2)->getMapVecFieldUnique(), this->getDomain(5)->getMapVecFieldUnique());
        }
        

        
        // ###########################
        // Bloecke hinzufuegen
        // ###########################
        if(geometryExplicit_)
            this->system_.reset(new BlockMatrix_Type(5));
        else
            this->system_.reset(new BlockMatrix_Type(6));
        
        // Fluid
        this->system_->addBlock( this->problemFluid_->system_->getBlock(0,0), 0, 0 );
        this->system_->addBlock( this->problemFluid_->system_->getBlock(0,1), 0, 1 );
        this->system_->addBlock( this->problemFluid_->system_->getBlock(1,0), 1, 0 );
        if (this->getDomain(0)->getFEType()=="P1")
            this->system_->addBlock( this->problemFluid_->system_->getBlock(1,1), 1, 1 );
        
        // SCI
        /*if (materialModel_=="linear")
            this->system_->addBlock( this->problemStructure_->system_->getBlock(0,0), 2, 2 );
        else
            this->system_->addBlock( this->problemStructureNonLin_->system_->getBlock(0,0), 2, 2 );*/

        this->system_->addBlock( this->problemSCI_->system_->getBlock(0,0), 2, 2 ); // Structure
        this->system_->addBlock( this->problemSCI_->system_->getBlock(1,1), 3, 3); // Chem
        this->system_->addBlock( this->problemSCI_->system_->getBlock(0,1), 2, 3 ); // Coupling of chem
        this->system_->addBlock( this->problemSCI_->system_->getBlock(1,0), 3, 2 ); // Coupling of structure


        // Kopplung
        this->system_->addBlock( C1_T, 0, 4 );
        this->system_->addBlock( C3_T, 2, 4 );
        this->system_->addBlock( C1, 4, 0 );
        this->system_->addBlock( C2, 4, 2 );

        if (!dummyC.is_null())
            this->system_->addBlock( dummyC, 4, 4 );
        
        if(!geometryExplicit_)
        {
            // Geometrie
            this->system_->addBlock( this->problemGeometry_->system_->getBlock(0,0), 5, 5 );
            // Kopplung
            this->system_->addBlock( C4, 5, 2 );
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
        /*if (materialModel_=="linear")
            plStructure = this->problemStructure_->getParameterList();
        else
            plStructure = this->problemStructureNonLin_->getParameterList();*/

        this->setupSubTimeProblems(this->problemFluid_->getParameterList(), this->problemSCI_->getParameterList(),this->problemSCI_->getParameterList());
        
        if (this->verbose_)
        {
            std::cout << "Assembly done -- " << std::endl;
        }
    }
    else
        reAssemble(type);
}


template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::reAssemble(std::string type) const
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

    if(type == "UpdateChemInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateChemInTime)" << '\n';

        this->problemSCI_->updateChemInTime();
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

    if(type == "ComputeChemRHSInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (ComputeChemRHSInTime)" << '\n';
        
        this->problemSCI_->computeChemRHSInTime( );
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
        
        this->problemSCI_->computeSolidRHSInTime( );
        return;
    }

    
    // ###############
    // Berechne die Gittergeschwindigkeit w und FSCI-Konvektion-Anteil (u-w)
    // ###############
    MultiVectorConstPtr_Type fluidSolution = this->solution_->getBlock(0);
    this->u_rep_->importFromVector(fluidSolution, true);
    this->u_minus_w_rep_->importFromVector(fluidSolution, true); //    u_minus_w_rep_ = *u_rep_;

    MultiVectorConstPtr_Type geometrySolution;
    if(this->geometryExplicit_)
    {
        geometrySolution = this->problemGeometry_->getSolution()->getBlock(0);
    }
    else
    {
        geometrySolution = this->solution_->getBlock(5);
    }
    this->meshDisplacementNew_rep_->importFromVector(geometrySolution, true);

    *this->w_rep_ = *meshDisplacementNew_rep_;
    this->w_rep_->update(-1.0, *meshDisplacementOld_rep_, 1.0);
    this->w_rep_->scale( 1.0/dt );

    this->u_minus_w_rep_->update( -1.0, *w_rep_, 1.0 );

    // Selbiges fuer den Druck
    MultiVectorConstPtr_Type pressureSolution = this->solution_->getBlock(1);
    this->p_rep_->importFromVector(pressureSolution, true);


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
        if(this->geometryExplicit_)
        {
            // ACHTUNG: Fluid-Loesung wird hier auf Null gesetzt, wegen initializeVectors().
            // Somit dann auch die problemTimeFluid_ wodurch eine falsche BDF2-RHS entsteht.
            // Rufe im DAESolverInTime deswegen erneut setPartialSolutions() auf.
            
            
            this->problemFluid_->assembleConstantMatrices(); // Die Steifikeitsmatrix wird weiter unten erst genutzt
            
            // Es ist P = P_
            this->P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
//            counterP++;
//            std::string outNameMeshVelo = "meshVelo" + to_string(counterP) + ".mm";
//            w_rep_->writeMM(outNameMeshVelo);
//            Teuchos::ArrayRCP<SC> values = w_rep_->getDataNonConst(0);
//            for (int i=0; i<values.size()/2; i++) {
//                values[2*i] = i;
//            }
            this->feFactory_->assemblyAdditionalConvection( this->dim_, this->domain_FEType_vec_.at(0), P_, w_rep_, true );
            this->P_->resumeFill();
            this->P_->scale(density);
            this->P_->scale(-1.0);
//            P_->scale(.0);
            this->P_->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
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
            cout << "Assemble Newton for Fluid done" <<endl;
            //if (materialModel_ != "linear")
            this->problemSCI_->reAssemble("Newton");
            
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

            this->feFactory_->assemblyShapeDerivativeVelocity(this->dim_, this->domain_FEType_vec_.at(5), this->domain_FEType_vec_.at(1),
                            shapeVelocity, 5, u_rep_, w_rep_, p_rep_, dt, density, viscosity, true);
            this->feFactory_->assemblyShapeDerivativeDivergence(this->dim_, this->domain_FEType_vec_.at(5), this->domain_FEType_vec_.at(1),
                            shapeDiv, 1, 5, this->getDomain(1)->getMapUnique(), this->getDomain(5)->getMapVecFieldUnique(), u_rep_, true);
            shapeDiv->resumeFill();
            shapeDiv->scale(-1.0);
            shapeDiv->fillComplete(this->getDomain(4)->getMapVecFieldUnique(), this->getDomain(1)->getMapUnique());

            // Shape Reinschreiben
            this->system_->addBlock(shapeVelocity, 0, 5);
            this->system_->addBlock(shapeDiv, 1, 5);
            
            //if (materialModel_ != "linear")
            this->problemSCI_->reAssemble("Newton");

        }
    }

    this->system_->addBlock( problemFluid_->getSystem()->getBlock( 0, 0 ), 0, 0 );
    
    //if (materialModel_ != "linear")
    this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(0,0), 2, 2 );
    this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(1,1), 3, 3 );
    this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(0,1), 2, 3 );
    this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(1,0), 3, 2 );

}

/*template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::reAssembleExtrapolation(BlockMultiVectorPtrArray_Type previousSolutions)
{
    double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);
    // Fluid-Dichte
    double density = this->problemFluid_->getParameterList()->sublist("Parameter").get("Density",1.);

    // ###############
    // w bestimmen
    // ###############
    MultiVectorConstPtr_Type geometrySolution;
    if(geometryExplicit_)
    {
        geometrySolution = this->problemGeometry_->getSolution()->getBlock(0);
    }
    else
    {
        geometrySolution = this->solution_->getBlock(4);
    }
    this->meshDisplacementNew_rep_->importFromVector(geometrySolution, true);


    *this->w_rep_ = *meshDisplacementNew_rep_;
    // update(): this = alpha*A + beta*this
    this->w_rep_->update(-1.0, *meshDisplacementOld_rep_, 1.0);
    this->w_rep_->scale(1.0/dt);


    // ###############
    // u extrapolieren
    // ###############
    // Beachte: getBlock(0) ist hier der Richtige, da dies der u-Variable in FSCI entspricht.
    if (previousSolutions.size() >= 2)
    {
        MultiVectorPtr_Type extrapolatedVector = Teuchos::rcp( new MultiVector_Type( previousSolutions[0]->getBlock(0) ) );

        // Extrapolation fuer BDF2
        // update(): this = alpha*A + beta*this
        extrapolatedVector->update( -1., *previousSolutions[1]->getBlock(0), 2. );

        this->u_rep_->importFromVector(extrapolatedVector, true);
    }
    else if(previousSolutions.size() == 1)
    {
        MultiVectorConstPtr_Type u = previousSolutions[0]->getBlock(0);
       this->u_rep_->importFromVector(u, true);
    }
    else if (previousSolutions.size() == 0)
    {
        MultiVectorConstPtr_Type u = this->solution_->getBlock(0);
       this->u_rep_->importFromVector(u, true);
    }

    // (u-w)
    *u_minus_w_rep_ = *u_rep_;
    // update(): this = alpha*A + beta*this
    u_minus_w_rep_->update(-1.0, *w_rep_, 1.0);


    // ###############
    // Neu assemblieren
    // ###############

    MatrixPtr_Type APN = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );

    MatrixPtr_Type N = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    this->feFactory_->assemblyAdvectionVecField( this->dim_, this->domain_FEType_vec_.at(0), N, this->u_minus_w_rep_, true );

    this->P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
    
    this->feFactory_->assemblyAdditionalConvection( this->dim_, this->domain_FEType_vec_.at(0), this->P_, this->w_rep_, true );
    this->P_->resumeFill();
    this->P_->scale(density);
    this->P_->scale(-1.0);
    this->P_->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());


    N->resumeFill();
    N->scale(density);
    N->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());

    this->problemFluid_->A_->addMatrix(1.0, APN, 0.0);
    this->P_->addMatrix(1.0, APN, 1.0);
    N->addMatrix(1.0, APN, 1.0);

    APN->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique() );

    this->system_->addBlock( APN, 0, 0 );
}*/

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const
{
    
    MultiVectorConstPtr_Type geometrySolution;
    
    if(this->geometryExplicit_)
        geometrySolution = this->problemGeometry_->getSolution()->getBlock(0);
    else {
        moveMesh();
        geometrySolution = this->solution_->getBlock(4);
    }
    
    this->meshDisplacementNew_rep_->importFromVector(geometrySolution, true);
    
    MultiVectorConstPtr_Type fluidSolution = this->solution_->getBlock(0);
    this->u_rep_->importFromVector(fluidSolution, true);
    this->u_minus_w_rep_->importFromVector(fluidSolution, true); //    u_minus_w_rep_ = *u_rep_;
    
    *this->w_rep_ = *this->meshDisplacementNew_rep_;
    this->w_rep_->update(-1.0, *this->meshDisplacementOld_rep_, 1.0);
    double dt = this->parameterList_->sublist("Timestepping Parameter").get("dt",0.02);
    this->w_rep_->scale(1.0/dt);
    
    this->u_minus_w_rep_->update(-1.0, *w_rep_, 1.0);
    
    if (!this->geometryExplicit_) {
        
        P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        double density = this->problemTimeFluid_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
        
        this->feFactory_->assemblyAdditionalConvection( this->dim_, this->domain_FEType_vec_.at(0), this->P_, this->w_rep_, true );
        this->P_->resumeFill();
        this->P_->scale(density);
        this->P_->scale(-1.0);
        this->P_->fillComplete( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique());
        
        this->problemFluid_->assembleConstantMatrices();
        
        this->system_->addBlock( this->problemFluid_->system_->getBlock(0,1), 0, 1 );
        this->system_->addBlock( this->problemFluid_->system_->getBlock(1,0), 1, 0 );
        TEUCHOS_TEST_FOR_EXCEPTION(this->problemFluid_->system_->blockExists(1,1) , std::runtime_error, "Stabilization is used. Account for it.");
    }
    if ( this->verbose_ )
        std::cout << "Warning: Wrong consideration of temporal discretization for multi-stage RK methods!" << std::endl;
    
    this->problemFluid_->calculateNonLinResidualVecWithMeshVelo( "reverse", time, u_minus_w_rep_, P_ );
    this->system_->addBlock( problemFluid_->getSystem()->getBlock( 0, 0 ), 0, 0 );
    
    // we need to account for the coupling in the residuals
    this->problemSCI_->calculateNonLinResidualVec( "reverse", time );

    this->residualVec_->addBlock(  this->problemSCI_->getResidualVector()->getBlockNonConst(0) , 2);
    this->residualVec_->addBlock(  this->problemSCI_->getResidualVector()->getBlockNonConst(1) , 3);

    MultiVectorPtr_Type residualFluidVelocityFSCI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(0) );
    MultiVectorPtr_Type residualSolidFSCI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(2) );
    MultiVectorPtr_Type residualChemFSCI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(3) );


    MultiVectorPtr_Type residualCouplingFSCI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(4) );
    residualCouplingFSCI->update( 1. , *this->rhs_->getBlock(4), 0. ); // change to -1 for standard
    
    //Now we need to add the coupling blocks
    this->system_->getBlock(0,4)->apply( *this->solution_->getBlock(4) , *residualFluidVelocityFSCI, Teuchos::NO_TRANS, -1., 1. );
    
    this->system_->getBlock(2,4)->apply( *this->solution_->getBlock(4) , *residualSolidFSCI, Teuchos::NO_TRANS, -1., 1. );
    
    this->system_->getBlock(4,0)->apply( *this->solution_->getBlock(0) , *residualCouplingFSCI, Teuchos::NO_TRANS, -1., 1. );
    
    this->system_->getBlock(4,2)->apply( *this->solution_->getBlock(2) , *residualCouplingFSCI, Teuchos::NO_TRANS, -1., 1. );

    if (!geometryExplicit_) {
        
        MultiVectorPtr_Type residualGeometryFSCI =
            Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(5) );
        residualGeometryFSCI->update( 1. , *this->rhs_->getBlock(5), 0. ); // change to -1 for standard

        this->system_->getBlock(5,5)->apply( *this->solution_->getBlock(5) , *residualGeometryFSCI, Teuchos::NO_TRANS, -1., 1. );
        
        this->system_->getBlock(5,2)->apply( *this->solution_->getBlock(2) , *residualGeometryFSCI, Teuchos::NO_TRANS, -1., 1. );
        
    }
    // might also be called in the sub calculateNonLinResidualVec() methods which where used above
    if (type == "reverse")
        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    else if (type == "standard"){
        this->residualVec_->scale(-1.);
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
    }


}

// Muss derzeit nur am Anfang jeder Zeititeration aufgerufen werden, damit
// problemTimeFluid_ und problemTimeStructure_ die aktuelle Loesung haben.
// ACHTUNG: Wenn wir irgendwann einmal anfangen reAssemble() auf problemFluid_ und
// problemStructure_ aufzurufen, dann muessen wir in jeder nichtlinearen Iteration
// diese setPartialSolutions() aufrufen, damit problemFluid_ und problemStructure_
// den korrekten nichtlinearen Term ausrechnen koennen.
// CH: Ist das noch relevant?
// We need to build FSCI so this method is not needed anymore
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setFromPartialVectorsInit() const
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
    
    /*if (materialModel_=="linear"){
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
    }*/
    // Structure 
     this->solution_->addBlock( this->problemSCI_->getSolution()->getBlockNonConst(0), 2 );
    this->residualVec_->addBlock( this->problemSCI_->getResidualVector()->getBlockNonConst(0), 2 );
    this->rhs_->addBlock( this->problemSCI_->getRhs()->getBlockNonConst(0), 2 );
    this->previousSolution_->addBlock( this->problemSCI_->getPreviousSolution()->getBlockNonConst(0), 2 );
    this->sourceTerm_->addBlock( this->problemSCI_->getSourceTerm()->getBlockNonConst(0), 2 );

    // Diffusion 
    this->solution_->addBlock( this->problemSCI_->getSolution()->getBlockNonConst(1), 3 );
    this->residualVec_->addBlock( this->problemSCI_->getResidualVector()->getBlockNonConst(1), 3 );
    this->rhs_->addBlock( this->problemSCI_->getRhs()->getBlockNonConst(1), 3 );
    this->previousSolution_->addBlock( this->problemSCI_->getPreviousSolution()->getBlockNonConst(1), 3 );
    this->sourceTerm_->addBlock( this->problemSCI_->getSourceTerm()->getBlockNonConst(1), 3 );

    if(!this->geometryExplicit_){
        this->solution_->addBlock( this->problemGeometry_->getSolution()->getBlockNonConst(0), 5 );
        // we dont have a previous solution for linear problems
        this->rhs_->addBlock( this->problemGeometry_->getRhs()->getBlockNonConst(0), 5 );
        this->sourceTerm_->addBlock( this->problemGeometry_->getSourceTerm()->getBlockNonConst(0), 5 );
    }
}
    
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::updateMeshDisplacement() const
{

     *meshDisplacementOld_rep_ = *meshDisplacementNew_rep_;

}


template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::solveGeometryProblem() const
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
void FSCI<SC,LO,GO,NO>::setupSubTimeProblems(ParameterListPtr_Type parameterListFluid, ParameterListPtr_Type parameterListStructure,ParameterListPtr_Type parameterListChem ) const
{
    if(this->verbose_)
        std::cout << "-- Setup FSCI Sub-TimeProblems \n" << std::flush;

    double dt = timeSteppingTool_->get_dt();
    double beta = timeSteppingTool_->get_beta();

    int sizeFluid = this->problemFluid_->getSystem()->size();
    /*int sizeStructure;
    if (materialModel_=="linear")
        sizeStructure = this->problemStructure_->getSystem()->size();
    else
        sizeStructure = this->problemStructureNonLin_->getSystem()->size();
    */
    problemTimeFluid_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemFluid_, this->comm_));
   
   /* if (materialModel_=="linear")
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructure_, this->comm_));
    else
        problemTimeStructure_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemStructureNonLin_, this->comm_));
    */
    problemTimeSCI_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemSCI_, this->comm_));
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
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Implement other FSCI fluid time stepping than BDF.");
    }
    // ######################
    // Struktur: Mass-, Problem, SourceTerm Koeffizienten
    // ######################
    // Koeffizienten vor der Massematrix und vor der Systemmatrix des steady-Problems
    /*SmallMatrix<double> massCoeffStructure(sizeStructure);
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
            if((*defTS_)[i + sizeFluid][j + sizeFluid] == 1 && i == j) // Weil: (u_f, p, d_s,...) und timeStepDef_ von FSCI
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
    this->problemTimeStructure_->setTimeParameters(massCoeffStructure,problemCoeffStructure);*/

   this->problemTimeFluid_->assemble( "MassSystem" );
   // this->problemTimeStructure_->assemble( "MassSystem" );
  // this->problemSCI_->setupSubTimeProblems(parameterListStructure,parameterListChem); // already called in SCI
    if(this->verbose_)
        std::cout << " Setup Sub-Timeproblems done-- \n" << endl;

}


template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setFluidMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix fuer FSCI combineSystems(), ggf nichtlinear.
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


// TODO: updateMultistepRhsFSCI() einbauen!
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::computeFluidRHSInTime( ) const
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

    if (this->timeSteppingTool_->currentTime()==0.) {
        SmallMatrix<double> tmpmassCoeff(sizeFluid);
        SmallMatrix<double> tmpproblemCoeff(sizeFluid);
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++) {
                if ((*this->defTS_)[i][j]==1 && i==j) {
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
    if (this->timeSteppingTool_->currentTime()==0.) {
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
    if (this->timeSteppingTool_->currentTime()==0.) {
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
                if ((*this->defTS_)[i][j]==1){
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
void FSCI<SC,LO,GO,NO>::updateFluidInTime() const
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
void FSCI<SC,LO,GO,NO>::computeSolidRHSInTime() const {
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
    this->problemSCI_->problemTimeStructure_->updateSolutionNewmarkPreviousStep(dt, beta, gamma);
    
    // Stelle die rechte Seite des zeitdiskretisierten Systems auf (ohne f_{n+1}).
    // Bei Newmark lautet dies:
    // M*[\frac{1}{dt^2*beta}*u_n + \frac{1}{dt*beta}*u'_n + \frac{0.5 - beta}{beta}*u''_n],
    // wobei u' = v (velocity) und u'' = w (acceleration).
    this->problemSCI_->problemTimeStructure_->updateNewmarkRhs(dt, beta, gamma, coeffTemp);
    
    //can we get rid of this?
    double time = timeSteppingTool_->currentTime() + dt;
    
    // TODO: SourceTerm wird in jedem Zeitschritt neu berechnet; auch wenn konstant!!!
    // if(time == 0){nur dann konstanten SourceTerm berechnen}
    if (this->problemSCI_->problemTimeStructure_->hasSourceTerm())
    {
        this->problemSCI_->problemTimeStructure_->assembleSourceTerm( time );
        
        // Fuege die rechte Seite der DGL (f bzw. f_{n+1}) der rechten Seite hinzu (skaliert mit coeffSourceTerm)
        // Die Skalierung mit der Dichte erfolgt schon in der Assemblierungsfunktion!
        
        // addSourceTermToRHS() aus DAESolverInTime
        double coeffSourceTermStructure = 1.0;
        BlockMultiVectorPtr_Type tmpPtr =  this->problemSCI_->problemTimeStructure_->getSourceTerm();
         this->problemSCI_->problemTimeStructure_->getRhs()->update(coeffSourceTermStructure, *tmpPtr, 1.);
    }

}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setSolidMassmatrix( MatrixPtr_Type& massmatrix ) const
{

    //######################
    // Massematrix
    //######################
    this->problemSCI_->setSolidMassmatrix(massmatrix);


    /*double density = this->problemSCI_->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
    int size = this->problemSCI_->problemTimeStructure_->getSystem()->size();

    if(timeSteppingTool_->currentTime() == 0.0)
    {
        this->problemSCI_->problemTimeStructure_->systemMass_.reset(new BlockMatrix_Type(size));
        {

            massmatrix = Teuchos::rcp(new Matrix_Type( this->problemSCI_->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // 2 = Struktur
            this->feFactory_->assemblyMass(this->dim_, this->problemSCI_->problemTimeStructure_->getFEType(0), "Vector", massmatrix, 2, true);
            massmatrix->resumeFill();
            massmatrix->scale(density);
            massmatrix->fillComplete( this->problemSCI_->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->problemSCI_->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique());

            this->problemSCI_->problemTimeStructure_->systemMass_->addBlock( massmatrix, 0, 0 );
        }
    }*/
}

// --------------
// Set chem mass matrix
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setChemMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    this->problemSCI_->setChemMassmatrix(massmatrix);
    //######################
    // Massematrix
    //######################

    /*double density = this->problemSCI_->problemTimeStructure_->getParameterList()->sublist("Parameter").get("Density",1000.e-0);
    int size = this->problemSCI_->problemTimeStructure_->getSystem()->size();

    if(timeSteppingTool_->currentTime() == 0.0)
    {
        this->problemSCI_->problemTimeStructure_->systemMass_.reset(new BlockMatrix_Type(size));
        {

            massmatrix = Teuchos::rcp(new Matrix_Type( this->problemSCI_->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getApproxEntriesPerRow() ) );
            // 2 = Struktur
            this->feFactory_->assemblyMass(this->dim_, this->problemSCI_->problemTimeStructure_->getFEType(0), "Vector", massmatrix, 2, true);
            massmatrix->resumeFill();
            massmatrix->scale(density);
            massmatrix->fillComplete( this->problemSCI_->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique(), this->problemSCI_->problemTimeStructure_->getDomain(0)->getMapVecFieldUnique());

            this->problemSCI_->problemTimeStructure_->systemMass_->addBlock( massmatrix, 0, 0 );
        }
    }*/
}


// Damit die richtige timeSteppingTool_->currentTime() genommen wird.
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::updateTime() const
{
    this->timeSteppingTool_->t_ = this->timeSteppingTool_->t_ + this->timeSteppingTool_->dt_prev_;
}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::moveMesh() const
{

    MultiVectorConstPtr_Type displacementUniqueConst;
    if(this->geometryExplicit_)
    {
        displacementUniqueConst = this->problemGeometry_->getSolution()->getBlock(0);
    }
    else
    {
        displacementUniqueConst = this->solution_->getBlock(5);
    }
    MultiVectorPtr_Type displacementRepeated = Teuchos::rcp( new MultiVector_Type( this->problemGeometry_->getDomain(0)->getMapVecFieldRepeated() ) );

    displacementRepeated->importFromVector( displacementUniqueConst );
    MultiVectorPtr_Type displacementUnique = Teuchos::rcp_const_cast<MultiVector_Type>(displacementUniqueConst);


    // Verschiebe die Gitter fuer Fluid-Velocity und Fluid-Druck
    // ACHTUNG: Klappt nur, weil die P2-Knoten hinter den P1-Knoten kommen.
    // Sonst muessen fuer den Druck die P1-Knoten extrahiert werden.
    // TODO: Wahrscheinlich reicht nur FSCI-Domain, da selbes Objekt bei problemFluid_ und problemTimeFluid_.
    ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemFluid_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemFluid_->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeFluid_->getDomain(0)) )->moveMesh(displacementUnique, displacementRepeated);
    ( Teuchos::rcp_const_cast<Domain_Type>(this->problemTimeFluid_->getDomain(1)) )->moveMesh(displacementUnique, displacementRepeated);
}


template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::addInterfaceBlockRHS() const
{
    MultiVectorPtr_Type vectorToAdd = Teuchos::rcp( new MultiVector_Type( this->rhs_->getBlock(4) ) );

    this->C2_->apply(*(this->solution_->getBlock(2)), *vectorToAdd);
    this->rhs_->addBlock(vectorToAdd, 4);
}


template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                     const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                    ) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "implement NOX for steady FSCI.");
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
void FSCI<SC,LO,GO,NO>::getValuesOfInterest( vec_dbl_Type& values ){
    if (this->dim_==2)
        getValuesOfInterest2DBenchmark( values );
    else if(this->dim_==3)
        getValuesOfInterest3DBenchmark( values);
    
}
    
    
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::computeValuesOfInterestAndExport(){
    /*if ( this->getParameterList()->sublist("General").get("Export drag and lift",false) ) {
        
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
            bcFactoryDrag->addBC(FSI::drag2D, 4, 0, domainVelocity, "Dirichlet", dim); // obstacle
            bcFactoryDrag->addBC(drag2D, 5, 0, domainVelocity, "Dirichlet", dim); // interface; check main FSCI for matching flags at the obstacle and interface
            bcFactoryLift->addBC(lift2D, 4, 0, domainVelocity, "Dirichlet", dim);
            bcFactoryLift->addBC(lift2D, 5, 0, domainVelocity, "Dirichlet", dim);
        }
        else if( dim == 3 ){
            bcFactoryDrag->addBC(drag3D, 3, 0, domainVelocity, "Dirichlet", dim); // check main FSCI for matching
            bcFactoryDrag->addBC(drag3D, 6, 0, domainVelocity, "Dirichlet", dim); // check main FSCI for matching flags at the obstacle and interface
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
    }*/
}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::initializeGE(){
/*All vectors (solution, rhs, previousSolution,...) should be initialized at this point (initializeProblem() was called)
 Therefore, all these BlockMVectors should have a length of 5 since the geometry problem is included in the general setup. We need to resize these BlockMVs here if the Geometry Explicit (GE) system is used.
*/
    if (this->geometryExplicit_) {
        this->solution_->resize( 5 );
        this->rhs_->resize( 5 );
        this->sourceTerm_->resize( 5 );
        this->rhsFuncVec_.resize( 5 );
        this->previousSolution_->resize( 5 );
        this->residualVec_->resize( 5 );
        this->initVectorSpaces();  //reinitialize NOX vector spaces
    }
}

}
#endif
