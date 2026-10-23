#ifndef FSCI_def_hpp
#define FSCI_def_hpp
#include "FSCI_decl.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"
#include "feddlib/core/General/BCBuilder.hpp"

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
                    ParameterListPtr_Type parameterListChem, ParameterListPtr_Type parameterListSCI,
                    ParameterListPtr_Type parameterListFSCI, ParameterListPtr_Type parameterListGeometry,
                    Teuchos::RCP<SmallMatrix<int> > &defTS):
FSI<SC,LO,GO,NO>(domainVelocity,  FETypeVelocity,
                domainPressure,FETypePressure,
                domainStructure, FETypeStructure,
                domainInterface, FETypeInterface,
                domainGeometry, FETypeGeometry,
                parameterListFluid, parameterListStructure,
                parameterListFSCI, parameterListGeometry,defTS),
problemSCI_()
{
    this->addVariable( domainChem, FETypeChem, "c", 1 ); // Chem last added component!
    this->initNOXParameters();
    
    chemistryExplicit_ =    parameterListSCI->sublist("Parameter").get("Chemistry Explicit",false);

    
    Teuchos::RCP<SmallMatrix<int>> defTSSCI; // Seperate Timestepping Matrix for SCI
    defTSSCI.reset( new SmallMatrix<int> (1) );
    if(!chemistryExplicit_){
        defTSSCI.reset( new SmallMatrix<int> (2) );
        (*defTSSCI)[1][1] = (*defTS)[4][4];

    }
    (*defTSSCI)[0][0] = (*defTS)[2][2];

    this->problemSCI_ = Teuchos::rcp( new SCIProblem_Type( domainStructure, FETypeStructure, domainChem, FETypeChem, diffusionTensor,reactionFunc, parameterListStructure, parameterListChem, parameterListSCI, defTSSCI ) );
    this->problemSCI_->initializeProblem();

   // this->initVectorSpaces();  //reinitialize NOX vector spaces

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

        this->problemFluid_->assemble();
        
        this->problemSCI_->assemble();
        
        this->problemGeometry_->assemble();
        if ( this->geometryExplicit_ && this->parameterList_->sublist("Exporter").get("Export GE geometry solution",false)){
            this->exporterGeo_ = Teuchos::rcp(new Exporter_Type());
            
            DomainConstPtr_Type dom = this->getDomain(4);

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
        this->feFactory_->assemblyFSICoupling(this->dim_, this->domain_FEType_vec_.at(0), C1, C1_T, 0, 0,this->getDomain(0)->getInterfaceMapVecFieldUnique(), this->getDomain(0)->getMapVecFieldUnique(), true);
        
        // Struktur-Bloecke
        // ACHTUNG: Die Interface-Variable \lambda wird eindeutig von der Fluid-Seite gehalten.
        // Deswegen auch getDomain(0) fuer die Spalten der Kopplungsbloecke.
        this->feFactory_->assemblyFSICoupling(this->dim_, this->domain_FEType_vec_.at(2), C2, C3_T, 0, 2,this->getDomain(0)->getInterfaceMapVecFieldUnique(), this->getDomain(2)->getMapVecFieldUnique(), true);

        
        // Falls geometrisch implizit
        /*if(!this->geometryExplicit_)
        {
            // TODO: Wegen IndicesGlobalMatched_ vlt .sicherheitshalber FEloc = 0 nehmen.
            // TODO: Check C4
            this->feFactory_->assemblyGeometryCoupling(this->dim_, this->domain_FEType_vec_.at(4), C4, 4,
                                                       this->getDomain(0)->getGlobalInterfaceMapUnique(),
                                                       this->getDomain(2)->getMapVecFieldUnique(),
                                                       this->getDomain(5)->getMapVecFieldUnique(), true);
        }*/
        
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

        if(!this->geometryExplicit_)
        {
            C4->resumeFill();
            C4->scale(-1.0);
            // Domain = Spalten = Struktur; Range = Zeilen = Geometrie
            C4->fillComplete(this->getDomain(2)->getMapVecFieldUnique(), this->getDomain(5)->getMapVecFieldUnique());
        }
        

        
        // ###########################
        // Bloecke hinzufuegen
        // ###########################
        if(this->geometryExplicit_){
            this->system_.reset(new BlockMatrix_Type(5));
            if(chemistryExplicit_)
                this->system_.reset(new BlockMatrix_Type(4));
        }
        else{
            this->system_.reset(new BlockMatrix_Type(6));
            if(chemistryExplicit_)
                this->system_.reset(new BlockMatrix_Type(5));
        }
        
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

        // Chemistry goes to the last block
        this->system_->addBlock( this->problemSCI_->system_->getBlock(0,0), 2, 2 ); // Structure
        if(!chemistryExplicit_){
            this->system_->addBlock( this->problemSCI_->system_->getBlock(1,1), 4, 4); // Chem
            this->system_->addBlock( this->problemSCI_->system_->getBlock(0,1), 2, 4 ); // Coupling of chem
            this->system_->addBlock( this->problemSCI_->system_->getBlock(1,0), 4, 2 ); // Coupling of structure
        }

        // Kopplung
        this->system_->addBlock( C1_T, 0, 3 );
        this->system_->addBlock( C3_T, 2, 3 );
        this->system_->addBlock( C1, 3, 0 );
        this->system_->addBlock( C2, 3, 2 );

        if (!dummyC.is_null())
            this->system_->addBlock( dummyC, 3, 3 );
        
        if(!this->geometryExplicit_)
        {
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Only Geometry explicit available here");

        }

        // Sollte (bzw. muss) erst aufgerufen werden, wenn alle Bloecke aus assemble()
        // in das ganze System hineingeschrieben worden ist. Ansonsten findet
        // blockExists() keinen Block
    //    this->initializeVectors();
        //this->initializeVectorsNonLinear();
        //NOX

        // We set the vector from the partial problems
        this->setFromPartialVectorsInit();
        
        // Fuer die Zeitprobleme
        this->timeSteppingTool_ = Teuchos::rcp(new TimeSteppingTools(sublist(this->parameterList_,"Timestepping Parameter") , this->comm_));
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

    double dt = this->timeSteppingTool_->dt_;

    // Fluid-Dichte
    double density = this->problemFluid_->getParameterList()->sublist("Parameter").get("Density",1.);
    double viscosity = this->problemFluid_->getParameterList()->sublist("Parameter").get("Viscosity",1.);

    if(type == "UpdateMeshDisplacement")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (UpdateMeshDisplacement)" << '\n';

        // Da dieser Abschnitt zu Beginn der neuen Zeititeration aufgerufen wird, muss
        // auch die alte Gitterverschiebung durch die neue Loesung aktualisiert werden.
        this->updateMeshDisplacement();
        this->problemSCI_->reAssemble("UpdateMeshDisplacement");
        return;
    }

    if(type == "SolveGeometryProblem")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (SolveGeometryProblem)" << '\n';

        this->solveGeometryProblem();
        return;
    }

    if(type == "SolveChemistryProblem")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (SolveChemistryProblem)" << '\n';

        this->problemSCI_->solveChemistryProblem();
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

        this->updateFluidInTime();
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

        this->moveMesh();
        return;
    }

    if(type == "AddInterfaceBlockRHS")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (AddInterfaceBlockRHS)" << '\n';

        this->addInterfaceBlockRHS();
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
        
        this->computeFluidRHSInTime( );
        return;
    }

    if(type == "ComputeSolidRHSInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (ComputeSolidRHSInTime)" << '\n';
        
        this->computeSolidRHSInTime( );
        return;
    }
    if(type == "ComputePressureRHSInTime")
    {
        if(this->verbose_)
            std::cout << "-- Reassembly (ComputePressureRHSInTime)" << '\n';
        
        this->computePressureRHSInTime();
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

    *this->w_rep_ = *this->meshDisplacementNew_rep_;
    this->w_rep_->update(-1.0, *this->meshDisplacementOld_rep_, 1.0);
    this->w_rep_->scale( 1.0/dt );

    this->u_minus_w_rep_->update( -1.0, *this->w_rep_, 1.0 );

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
            this->feFactory_->assemblyAdditionalConvection( this->dim_, this->domain_FEType_vec_.at(0), this->P_, this->w_rep_, true );
            this->P_->resumeFill();
            this->P_->scale(density);
            this->P_->scale(-1.0);
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

    if(this->geometryExplicit_)
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
            
            this->problemFluid_->reAssemble( "Newton" );
            this->problemSCI_->reAssemble("Newton");
            
        }
    }
    /*else
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

            this->problemFluid_->reAssemble( "Newton" );
            
            // TODO: Shape
            // domain(0) = domain(4)
            MatrixPtr_Type shapeVelocity = Teuchos::rcp(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
            MatrixPtr_Type shapeDiv = Teuchos::rcp(new Matrix_Type( this->getDomain(1)->getMapUnique(), this->getDomain(1)->getDimension() * this->getDomain(1)->getApproxEntriesPerRow() ) ); // shape fuer div-Nebenbedingung

            this->feFactory_->assemblyShapeDerivativeVelocity(this->dim_, this->domain_FEType_vec_.at(5), this->domain_FEType_vec_.at(1),
                            shapeVelocity, 5, this->u_rep_, this->w_rep_, this->p_rep_, dt, density, viscosity, true);
            this->feFactory_->assemblyShapeDerivativeDivergence(this->dim_, this->domain_FEType_vec_.at(5), this->domain_FEType_vec_.at(1),
                            shapeDiv, 1, 5, this->getDomain(1)->getMapUnique(), this->getDomain(5)->getMapVecFieldUnique(), this->u_rep_, true);
            shapeDiv->resumeFill();
            shapeDiv->scale(-1.0);
            shapeDiv->fillComplete(this->getDomain(4)->getMapVecFieldUnique(), this->getDomain(1)->getMapUnique());

            // Shape Reinschreiben
            this->system_->addBlock(shapeVelocity, 0, 5);
            this->system_->addBlock(shapeDiv, 1, 5);
            
            this->problemSCI_->reAssemble("Newton");

        }
    }*/

    this->system_->addBlock( this->problemFluid_->getSystem()->getBlock( 0, 0 ), 0, 0 );
    
    this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(0,0), 2, 2 );
    if(!chemistryExplicit_){
        this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(1,1), 4, 4 );
        this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(0,1), 2, 4 );
        this->system_->addBlock(  this->problemSCI_->getSystem()->getBlock(1,0), 4, 2 );
    }
}
// Damit die richtige timeSteppingTool_->currentTime() genommen wird.
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::updateTime() const
{
    this->timeSteppingTool_->t_ = this->timeSteppingTool_->t_ + this->timeSteppingTool_->dt_prev_;
    this->problemTimeFluid_->updateTime(this->timeSteppingTool_->t_);
    this->problemSCI_->reAssemble("UpdateTime");

}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::calculateNonLinResidualVec(std::string type, double time) const
{
    
    MultiVectorConstPtr_Type geometrySolution;
    
    if(this->geometryExplicit_)
        geometrySolution = this->problemGeometry_->getSolution()->getBlock(0);
    /*else {
        this->moveMesh();
        geometrySolution = this->solution_->getBlock(4);
    }*/
    
     
    this->meshDisplacementNew_rep_->importFromVector(geometrySolution, true);
    
    MultiVectorConstPtr_Type fluidSolution = this->solution_->getBlock(0);
    this->u_rep_->importFromVector(fluidSolution, true);
    this->u_minus_w_rep_->importFromVector(fluidSolution, true); //    u_minus_w_rep_ = *u_rep_;
    
    *this->w_rep_ = *this->meshDisplacementNew_rep_;
    this->w_rep_->update(-1.0, *this->meshDisplacementOld_rep_, 1.0);
    double dt = this->timeSteppingTool_->dt_;
    this->w_rep_->scale(1.0/dt);
    
    this->u_minus_w_rep_->update(-1.0, *this->w_rep_, 1.0);
    
    if (!this->geometryExplicit_) {
        
        this->P_.reset(new Matrix_Type( this->getDomain(0)->getMapVecFieldUnique(), this->getDomain(0)->getDimension() * this->getDomain(0)->getApproxEntriesPerRow() ) );
        double density = this->problemTimeFluid_->getParameterList()->sublist("Parameter").get("Density",1.e-0);
        
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
    

    this->problemFluid_->calculateNonLinResidualVecWithMeshVelo( "reverse", time, this->u_minus_w_rep_, this->P_ );
    this->system_->addBlock( this->problemFluid_->getSystem()->getBlock( 0, 0 ), 0, 0 );
    
    // we need to account for the coupling in the residuals
    this->problemSCI_->calculateNonLinResidualVec( "reverse", time );

    this->residualVec_->addBlock(  this->problemSCI_->getResidualVector()->getBlockNonConst(0) , 2);
    if(!chemistryExplicit_) {
        this->residualVec_->addBlock(this->problemSCI_->getResidualVector()->getBlockNonConst(1), 4);
    }
    MultiVectorPtr_Type residualFluidVelocityFSCI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(0) );
    MultiVectorPtr_Type residualSolidFSCI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(2) );
    
    if(!chemistryExplicit_) {
        MultiVectorPtr_Type residualChemFSCI = Teuchos::rcp_const_cast<MultiVector_Type>(this->residualVec_->getBlock(4));
    }

    MultiVectorPtr_Type residualCouplingFSCI =
        Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(3) );
    residualCouplingFSCI->update( 1. , *this->rhs_->getBlock(3), 0. ); // change to -1 for standard
    
  //Now we need to add the coupling blocks
    this->system_->getBlock(0,3)->apply( *this->solution_->getBlock(3) , *residualFluidVelocityFSCI, Teuchos::NO_TRANS, -1., 1. );
    
    this->system_->getBlock(2,3)->apply( *this->solution_->getBlock(3) , *residualSolidFSCI, Teuchos::NO_TRANS, -1., 1. );
    
    this->system_->getBlock(3,0)->apply( *this->solution_->getBlock(0) , *residualCouplingFSCI, Teuchos::NO_TRANS, -1., 1. );
    
    this->system_->getBlock(3,2)->apply( *this->solution_->getBlock(2) , *residualCouplingFSCI, Teuchos::NO_TRANS, -1., 1. );

   /* if (!this->geometryExplicit_) {
        
        MultiVectorPtr_Type residualGeometryFSCI =
            Teuchos::rcp_const_cast<MultiVector_Type>( this->residualVec_->getBlock(4) );
        residualGeometryFSCI->update( 1. , *this->rhs_->getBlock(4), 0. ); // change to -1 for standard

        this->system_->getBlock(4,4)->apply( *this->solution_->getBlock(4) , *residualGeometryFSCI, Teuchos::NO_TRANS, -1., 1. );
        
        this->system_->getBlock(4,2)->apply( *this->solution_->getBlock(2) , *residualGeometryFSCI, Teuchos::NO_TRANS, -1., 1. );
        
    }*/
    // might also be called in the sub calculateNonLinResidualVec() methods which where used above
    if (type == "reverse"){
        this->bcFactory_->setBCMinusVector( this->residualVec_, this->solution_, time );
    }
    else if (type == "standard"){
        this->residualVec_->scale(-1.);
        this->bcFactory_->setVectorMinusBC( this->residualVec_, this->solution_, time );
    } 
    this->setBoundariesRHS(this->timeSteppingTool_->currentTime());
    /*bool plotResVector = this->getParameterList()->sublist("General").get("Plot Residual Vector",true);
    double range1 = this->getParameterList()->sublist("General").get("Plot Residual Vector Start",0.0);
    double range2 = this->getParameterList()->sublist("General").get("Plot Residual Vector End",1.0);
    if(plotResVector && this->timeSteppingTool_->currentTime() >= range1 && this->timeSteppingTool_->currentTime() <= range2)
        this->plotResidualVec(this->timeSteppingTool_->currentTime());*/



}

// This connects the here computed 'solution_' to the solutions of the subproblems. This way the subproblems have
// the correct solution when computing residual and co
// We need to build FSCI so this method is not needed anymore
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setFromPartialVectorsInit() const
{
    
    //Fluid velocity
    this->solution_->addBlock( this->problemFluid_->getSolution()->getBlockNonConst(0), 0 );
    this->residualVec_->addBlock( this->problemFluid_->getResidualVector()->getBlockNonConst(0), 0 );
    this->previousSolution_->addBlock( this->problemFluid_->getPreviousSolution()->getBlockNonConst(0), 0 );
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
    if(!chemistryExplicit_){
        this->solution_->addBlock( this->problemSCI_->getSolution()->getBlockNonConst(1), 4 );
        this->residualVec_->addBlock( this->problemSCI_->getResidualVector()->getBlockNonConst(1), 4 );
        this->rhs_->addBlock( this->problemSCI_->getRhs()->getBlockNonConst(1), 4 );
        this->previousSolution_->addBlock( this->problemSCI_->getPreviousSolution()->getBlockNonConst(1), 4 );
        this->sourceTerm_->addBlock( this->problemSCI_->getSourceTerm()->getBlockNonConst(1), 4 );
    }
   /* if(!this->geometryExplicit_){
        this->solution_->addBlock( this->problemGeometry_->getSolution()->getBlockNonConst(0), 5 );
        // we dont have a previous solution for linear problems
        this->rhs_->addBlock( this->problemGeometry_->getRhs()->getBlockNonConst(0), 5 );
        this->sourceTerm_->addBlock( this->problemGeometry_->getSourceTerm()->getBlockNonConst(0), 5 );
    }*/

}
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setupSubTimeProblems(ParameterListPtr_Type parameterListFluid, ParameterListPtr_Type parameterListStructure,ParameterListPtr_Type parameterListChem ) const
{
    if(this->verbose_)
        std::cout << "-- Setup FSCI Sub-TimeProblems \n" << std::flush;

    double dt = this->timeSteppingTool_->get_dt();
    double beta = this->timeSteppingTool_->get_beta();

    int sizeFluid = this->problemFluid_->getSystem()->size();
    /*int sizeStructure;
    if (materialModel_=="linear")
        sizeStructure = this->problemStructure_->getSystem()->size();
    else
        sizeStructure = this->problemStructureNonLin_->getSystem()->size();
    */
    this->problemTimeFluid_.reset(new TimeProblem<SC,LO,GO,NO>(*this->problemFluid_, this->comm_));
   
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
                if ((*this->defTS_)[i][j]==1 && i==j) {
                    defFluid[i][j] = 1;
                    massCoeffFluid[i][j] = this->timeSteppingTool_->getInformationBDF(0) / dt;
                }
                else{
                    massCoeffFluid[i][j] = 0.0;
                }
            }
        }
        for (int i=0; i<sizeFluid; i++) {
            for (int j=0; j<sizeFluid; j++){
                if ((*this->defTS_)[i][j]==1){
                    problemCoeffFluid[i][j] = this->timeSteppingTool_->getInformationBDF(1);
                    coeffSourceTermFluid = this->timeSteppingTool_->getInformationBDF(1);
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
   this->problemTimeFluid_->assemble( "MassSystem" );
   // this->problemTimeStructure_->assemble( "MassSystem" );
  // this->problemSCI_->setupSubTimeProblems(parameterListStructure,parameterListChem); // already called in SCI
    if(this->verbose_)
        std::cout << " Setup Sub-Timeproblems done-- \n" << std::endl;

}


template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setFluidMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix fuer FSCI combineSystems(), ggf nichtlinear.
    //######################
    double density = this->problemTimeFluid_->getParameterList()->sublist("Parameter").get("Density",1.e-0);
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

// Am Anfang der Zeititeration erst updateSolutionMultiPreviousStep() aufrufen und dann erst updateMultistepRhs(),
// damit die previousSolution_ initialisiert sind. Genauso fuer SystemMass
// TODO: updateSystemMassMultiPreviousStep() fertig programmieren
/*template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::updateFluidInTime() const
{
    int nmbBDF = this->timeSteppingTool_->getBDFNumber();

    if(nmbBDF<2 && !this->parameterList_->sublist("General").get("Linearization","FixedPoint").compare("Extrapolation")) {
        if (this->timeSteppingTool_->currentTime()!=0.){
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
}*/

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::computeSolidRHSInTime() const {
    //######################
    // RHS nach Newmark
    //######################
    double dt = this->timeSteppingTool_->get_dt();
    double beta = this->timeSteppingTool_->get_beta();
    double gamma = this->timeSteppingTool_->get_gamma();
    
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
    double time = this->timeSteppingTool_->currentTime() + dt;
    
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
}

// --------------
// Set chem mass matrix
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::setChemMassmatrix( MatrixPtr_Type& massmatrix ) const
{
    //######################
    // Massematrix
    //######################
    this->problemSCI_->setChemMassmatrix(massmatrix);
    

}


template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
                                     const Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs
                                    ) const
{
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "implement NOX for steady FSCI.");
    std::string type = this->parameterList_->sublist("General").get("Preconditioner Method","Monolithic");
}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::getValuesOfInterest( vec_dbl_Type& values ){    
}
    
    
template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::computeValuesOfInterestAndExport(){
}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::initializeGE(){
/*All vectors (solution, rhs, previousSolution,...) should be initialized at this point (initializeProblem() was called)
 Therefore, all these BlockMVectors should have a length of 5 since the geometry problem is included in the general setup. We need to resize these BlockMVs here if the Geometry Explicit (GE) system is used.

 In first initialization all variables were added including geometry. Thus, the general vectors are on entry too long. But the geometry is no longer the last component
 and we need to account for the diffusion

*/
    if (this->geometryExplicit_) {
        if(chemistryExplicit_){
            this->solution_->resize( 4 );
            this->rhs_->resize( 4 );
            this->sourceTerm_->resize( 4 );
            this->rhsFuncVec_.resize( 4 );
            this->previousSolution_->resize( 4 );
            this->residualVec_->resize( 4 );
            // Add Diffusion again 
        
        }
        else{
            this->solution_->resize( 5 );
            this->rhs_->resize( 5 );
            this->sourceTerm_->resize( 5 );
            this->rhsFuncVec_.resize( 5 );
            this->previousSolution_->resize( 5 );
            this->residualVec_->resize( 5 );
            // Add Diffusion again 
            this->solution_->addBlock( this->problemSCI_->getSolution()->getBlockNonConst(1), 4 );
            this->residualVec_->addBlock( this->problemSCI_->getResidualVector()->getBlockNonConst(1), 4 );
            this->rhs_->addBlock( this->problemSCI_->getRhs()->getBlockNonConst(1), 4 );
            this->previousSolution_->addBlock( this->problemSCI_->getPreviousSolution()->getBlockNonConst(1), 4 );
            this->sourceTerm_->addBlock( this->problemSCI_->getSourceTerm()->getBlockNonConst(1), 4 );
        }
        this->initVectorSpaces();  //reinitialize NOX vector spaces
    }
}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::exportValuesOfInterest(double time)
{
    problemSCI_->exportValuesOfInterest(time);
    if(this->geometryExplicit_)
    {
        std::string varName = std::to_string(time);
        this->exporterGeometry_->writeVariablesHDF5(varName,this->problemGeometry_->getSolution()->getBlock(0));
    }
}

template<class SC,class LO,class GO,class NO>
void FSCI<SC,LO,GO,NO>::importValuesOfInterest(double time)
{
    problemSCI_->importValuesOfInterest(time);
}

}
#endif
