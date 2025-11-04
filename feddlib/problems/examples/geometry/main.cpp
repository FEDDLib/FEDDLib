#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/problems/specific/LinElas.hpp"

void geometryBC2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.01;
    res[1] = 0.02;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroBC(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;

    return;
}

void x1(double* x, double* res, double t, const double* parameters){

    res[0] = 10.;

    return;
}
void x2(double* x, double* res, double t, const double* parameters){

    res[0] = 20.;

    return;
}
void x3(double* x, double* res, double t, const double* parameters){

    res[0] = 30.;

    return;
}
void x4(double* x, double* res, double t, const double* parameters){

    res[0] = 40.;

    return;
}
void x5(double* x, double* res, double t, const double* parameters){

    res[0] = 50.;

    return;
}
void rhs2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    
    return;
}

void rhs(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = parameters[1];
    res[2] = 0.;
    return;
}


typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

// Berechnet von einer dofID, d.h. dim*nodeID+(0,1,2), die entsprechende nodeID.
// IN localDofNumber steht dann, ob es die x- (=0), y- (=1) oder z-Komponente (=2) ist.
void toNodeID(UN dim, GO dofID, GO& nodeID, LO& localDofNumber )
{
    nodeID = (GO) (dofID/dim);
    localDofNumber = (LO) (dofID%dim);
}

// Diese Funktion berechnet genau das umgekehrte. Also von einer nodeID die entsprechende dofID
void toDofID(UN dim, GO nodeID, LO localDofNumber, GO& dofID )
{
    dofID = (GO) ( dim * nodeID + localDofNumber);
}

using namespace FEDD;
using namespace Teuchos;
using namespace std;

int main(int argc, char *argv[])
{


    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef ExporterParaView<SC,LO,GO,NO> ExporterPV_Type;
    typedef RCP<ExporterPV_Type> ExporterPVPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    
    typedef Map<LO,GO,NO> Map_Type;
    typedef RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
   
    
    typedef std::vector<GO> vec_GO_Type;
    typedef std::vector<vec_GO_Type> vec2D_GO_Type;
    typedef std::vector<vec2D_GO_Type> vec3D_GO_Type;
    typedef Teuchos::RCP<vec3D_GO_Type> vec3D_GO_ptr_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        return EXIT_SUCCESS;
    }

    bool verbose (comm->getRank() == 0);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        int 		volumeID        = parameterListProblem->sublist("Parameter").get("Volume ID",0);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        string		meshName_fluid    	= parameterListProblem->sublist("Parameter").get("Mesh Name Fluid","dfg_fsi_fluid_h002.mesh");
        string		meshName_struc    	= parameterListProblem->sublist("Parameter").get("Mesh Name Structure","dfg_fsi_solid_h002.mesh");
        string		meshDelimiter   = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
        string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        int         n;

        Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
        Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
        Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        bool serialP2Mesh = parameterListProblem->sublist("General").get("Build serial P2 Mesh",true);

        {
            if (verbose)
            {
                cout << "#################################################" <<endl;
                cout << "############ Starting geometry test ... #########" <<endl;
                cout << "#################################################" <<endl;
            }

            DomainPtr_Type domainP1fluid;
            DomainPtr_Type domainP1struct;
            DomainPtr_Type domainP2fluid;
            DomainPtr_Type domainP2struct;

            {
                Teuchos::TimeMonitor totalTimeMonitor(*totalTime);
                {
                    Teuchos::TimeMonitor buildMeshMonitor(*buildMesh);
                    if (verbose)
                    {
                        cout << "-- Building Mesh ... " << flush;
                    }

                    domainP1fluid.reset( new Domain_Type( comm, dim ) );
                    domainP1struct.reset( new Domain_Type( comm, dim ) );
                    domainP2fluid.reset( new Domain_Type( comm, dim ) );
                    domainP2struct.reset( new Domain_Type( comm, dim ) );

                    if (!meshType.compare("unstructured")) {
                        
                        vec_int_Type idsInterface(1,-1);
                        idsInterface[0] = 5;
                        
                        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
                        MeshPartitioner_Type::DomainPtrArray_Type domainP2Array(2);
                        domainP1Array[0] = domainP1fluid;
                        domainP1Array[1] = domainP1struct;
                        
                        ParameterListPtr_Type pListPartitioner = sublist( parameterListProblem, "Mesh Partitioner" );
                        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                        
                        partitionerP1.readAndPartition();
                        domainP2fluid->buildP2ofP1Domain( domainP1fluid );
                        domainP2struct->buildP2ofP1Domain( domainP1struct );
                        
                        // Calculate distances is done in: identifyInterfaceParallelAndDistance
                        domainP1fluid->identifyInterfaceParallelAndDistance(domainP1struct, idsInterface);
                        domainP2fluid->identifyInterfaceParallelAndDistance(domainP2struct, idsInterface);

                        
                    }
                    else{
                        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test for unstructured meshes read from .mesh-file. Change mesh type in setup file to 'unstructured'.");
                    }
                    if (verbose){
                        cout << "done! -- " << endl;
                    }
                }
            }


            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluid( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStruct( new BCBuilder<SC,LO,GO,NO>( ) );

            DomainPtr_Type domainFluid;
            DomainPtr_Type domainStruct;
            if (!discType.compare("P2")){
                domainFluid = domainP2fluid;
                domainStruct = domainP2struct;
            }
            else{
                domainFluid = domainP1fluid;
                domainStruct = domainP1struct;
            }

            bcFactoryFluid->addBC(x1, 0, 0, domainFluid, "Dirichlet", 1);
            bcFactoryFluid->addBC(x2, 1, 0, domainFluid, "Dirichlet", 1);
            bcFactoryFluid->addBC(x3, 2, 0, domainFluid, "Dirichlet", 1);
            bcFactoryFluid->addBC(x4, 3, 0, domainFluid, "Dirichlet", 1);
            bcFactoryFluid->addBC(x5, 4, 0, domainFluid, "Dirichlet", 1);

            bcFactoryStruct->addBC(x1, 0, 0, domainStruct, "Dirichlet", 1);
            bcFactoryStruct->addBC(x2, 1, 0, domainStruct, "Dirichlet", 1);
            bcFactoryStruct->addBC(x3, 2, 0, domainStruct, "Dirichlet", 1);
            bcFactoryStruct->addBC(x4, 3, 0, domainStruct, "Dirichlet", 1);
            bcFactoryStruct->addBC(x5, 4, 0, domainStruct, "Dirichlet", 1);


            MultiVectorPtr_Type valuesFluid = rcp(new MultiVector_Type( domainFluid->getMapUnique() ) );
            MultiVectorPtr_Type valuesStruct = rcp(new MultiVector_Type( domainStruct->getMapUnique() ) );
            BlockMultiVectorPtr_Type valuesFluidBlock = rcp(new BlockMultiVector_Type( 1 ) );
            BlockMultiVectorPtr_Type valuesStructBlock = rcp(new BlockMultiVector_Type( 1 ) );

            valuesFluidBlock->addBlock( valuesFluid, 0 );
            valuesStructBlock->addBlock( valuesStruct, 0 );

            bcFactoryFluid->setRHS( valuesFluidBlock );
            bcFactoryStruct->setRHS( valuesStructBlock );

            bool exportParaView = parameterListAll->sublist("General").get("ParaViewExport",true);
            if (exportParaView) {
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

                DomainPtr_Type domain = domainFluid;
                                
                exPara->setup("fluid", domain->getMesh(), discType);

                MultiVectorConstPtr_Type valuesFluidConst = valuesFluidBlock->getBlock( 0 );
                exPara->addVariable( valuesFluidConst, "values", "Scalar", 1, domain->getMapUnique());

                exPara->save(0.0);
                exPara->closeExporter();
            }

            bool exportParaViewSubomains = parameterListAll->sublist("General").get("ParaViewExport Subdomains",false);
            if (exportParaViewSubomains) {

                DomainPtr_Type domain = domainP1fluid;
                MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domain->getElementMap() ) );
                MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                vecDecomposition->putScalar(comm->getRank()+1.);

                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

                exPara->setup("subdomainsFluid", domain->getMesh(), "P0");
                
                exPara->addVariable( vecDecompositionConst, "subdomain", "Scalar", 1, domain->getElementMap());

                exPara->save(0.0);
                exPara->closeExporter();
            }


            if (exportParaView) {
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

                DomainPtr_Type domain = domainStruct;
                
                exPara->setup("structure", domain->getMesh(), discType);

                MultiVectorConstPtr_Type valuesStructConst = valuesStructBlock->getBlock( 0 );
                exPara->addVariable( valuesStructConst, "values", "Scalar", 1, domain->getMapUnique());

                exPara->save(0.0);
            }


            if (exportParaViewSubomains) {
                DomainPtr_Type domain = domainP1struct;

                MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domain->getElementMap() ) );
                MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;

                vecDecomposition->putScalar(comm->getRank()+1.);

                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup("subdomainsStruct", domain->getMesh(), "P0");

                exPara->addVariable( vecDecompositionConst, "subdomain", "Scalar", 1, domain->getElementMap());

                exPara->save(0.0);
                exPara->closeExporter();
            }


            // ########################
            // Export fuer den Abstand zum Interface
            // ########################
            if (exportParaView)
            {
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

                DomainPtr_Type domain = domainFluid;
                exPara->setup("distToInterface", domain->getMesh(), discType);
                
                Teuchos::ArrayRCP< SC > values = valuesFluid->getDataNonConst(0); // only a single MultiVector

                vec_dbl_ptr_Type distancesRepeated = domain->getDistancesToInterface();
                for(UN i = 0; i < values.size(); i++)
                {
                    GO indexGlob = domain->getMapUnique()->getGlobalElement(i);
                    LO index = domain->getMapRepeated()->getLocalElement( indexGlob );
                    values[i] = distancesRepeated->at(index);
                }

                MultiVectorConstPtr_Type valuesFluidConst = valuesFluid;
                exPara->addVariable( valuesFluidConst, "distances", "Scalar", 1, domain->getMapUnique());

                exPara->save(0.0);
                exPara->closeExporter();
            }

            // ############################################################
            // ########### Test fuer das Geometry-Problem #################
            // ############################################################
            // Strukturproblem loesen
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );
            bcFactoryStructure->addBC(zeroDirichlet2D, 1, 0, domainStruct, "Dirichlet", dim); // 1 ist homogener Dirichletrand der Struktur

            LinElas<SC,LO,GO,NO> LinElas(domainStruct,discType,parameterListAll);

            {
                Teuchos::TimeMonitor solveTimeMonitor(*solveTime);

                LinElas.addBoundaries(bcFactoryStructure); // Dem Problem RW hinzufuegen

                if (dim==2)
                    LinElas.addRhsFunction( rhs2D );
                else if(dim==3)
                    LinElas.addRhsFunction( rhs );
                
                double force = parameterListAll->sublist("Parameter").get("Volume force",0.);
                double degree = 0;
                
                LinElas.addParemeterRhs( force );
                LinElas.addParemeterRhs( degree );
                
                // ######################
                // Matrix assemblieren, RW setzen und System loesen
                // ######################
                LinElas.initializeProblem();
                LinElas.assemble();
                LinElas.addToRhs( LinElas.getSourceTerm() );
                LinElas.setBoundaries(); // In der Klasse Problem
                LinElas.solve();

            }

            // // ######################
            // // Exporter fuer die Strukturloesung
            // // ######################
            if(exportParaView)
            {
                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

                DomainPtr_Type domain =  domainStruct;

                Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolution = LinElas.getSolution()->getBlock(0);
                
                exPara->setup("linearElasticity", domain->getMesh(), discType);

                exPara->addVariable(exportSolution, "d_s", "Vector", dim, domain->getMapUnique());

                exPara->save(0.0);
                exPara->closeExporter();
            }

            // ########################
            // RW Geometry-Problem
            // ########################
            // Temporaer einfach RW setzen. Spaeter die Strukturloesung als RW setzen
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryGeometry( new BCBuilder<SC,LO,GO,NO>( ) );

            bcFactoryGeometry->addBC(zeroDirichlet2D, 1, 0, domainFluid, "Dirichlet", dim); // wall
            bcFactoryGeometry->addBC(zeroDirichlet2D, 2, 0, domainFluid, "Dirichlet", dim); // inflow
            bcFactoryGeometry->addBC(zeroDirichlet2D, 3, 0, domainFluid, "Dirichlet", dim); // outflow
            bcFactoryGeometry->addBC(zeroDirichlet2D, 4, 0, domainFluid, "Dirichlet", dim); // fluid-obstacle
            bcFactoryGeometry->addBC(geometryBC2D, 5, 0, domainFluid, "Dirichlet", dim); // 5 ist die Interface-Flag


            // Wir muessen aus der kompletten Strukturloesung die Loesung auf dem Interface
            // (Struktur) extrahieren und diesen dann an die korrekten Prozessoren,
            // die die entsprechenden, passenden Knoten auf dem Fluid enthalten, senden,
            // damit wir RW setzen koennen.

            // Mesh_ umcasten  in unstructured und dann meshUnstructured nutzen,
            // da nur unstructured Attribut MeshInterface_ besitzt.
            // Jeder Prozessor kennt also das komplette matched Interface
            // MeshUnstrPtr_Type meshUnstructuredFluid = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainFluid->getMesh() );
            //
            // vec3D_long_ptr_Type indicesGlobalMatchedOriginFluid = meshUnstructuredFluid->getMeshInterface()->getIndicesGlobalMatchedOrigin();
            //
            // // lokale Interface ID, die ich dem Knoten zuschreibe; Wir zaehlen von 0 bis #\Gamma-1
            // GO localInterfaceID = 0; // long long
            //
            // // In diesen beiden Vektoren stehen die Interface IDs in der Interface-Nummerierung,
            // // die der Prozessor haelt, und dies einmal von der Struktur und dem Fluid aus. Unique!!!
            // vec_long_Type vecInterfaceMapStructure; // vec_long ist long long wg. 64
            // vec_long_Type vecInterfaceMapFluid;
            //
            // // indicesGlobalMatchedOriginFluid sind die IndicesGlobalMatched von Fluid-Sicht aus.
            // // D.h. in indicesGlobalMatchedOriginFluid.at(0).at(0) ist Fluid GID und
            // // indicesGlobalMatchedOriginFluid.at(0).at(1) ist Struktur GID.
            //
            // // Fluid und Struktur haben gleich viele Interfaceknoten, deswegen
            // // muessen wir hier nicht zwischen Fluid und Struktur differenzieren
            // // ACHTUNG: indicesGlobalMatchedOriginFluid ist hier nicht partitioniert!!!!
            // for(int i = 0; i < indicesGlobalMatchedOriginFluid->size(); i++) // Schleife ueber jede flag
            // {
            //     for(int j = 0; j < indicesGlobalMatchedOriginFluid->at(i).at(0).size(); j++) // GIDs innerhalb der flag (vom Fluid)
            //     {
            //         // Wir muessen long long anstatt int nutzen, da MyGID auf 64 gestellt ist/ genutzt wird
            //         GO globalIDOfInterfaceNodeFluid = indicesGlobalMatchedOriginFluid->at(i).at(0).at(j);
            //         GO globalIDOfInterfaceNodeStructure = indicesGlobalMatchedOriginFluid->at(i).at(1).at(j);
            //
            //         // liefert true, falls Proz. GID besitzt
            //         if( domainFluid->getMapUnique()->getLocalElement(globalIDOfInterfaceNodeFluid) != OrdinalTraits<LO>::invalid())
            //         {
            //             // Schreibe die lokale Interface ID hinein
            //             vecInterfaceMapFluid.push_back(localInterfaceID);
            //         }
            //
            //         // liefert true, falls Proz. GID besitzt
            //         if( domainStruct->getMapUnique()->getLocalElement(globalIDOfInterfaceNodeStructure) != OrdinalTraits<LO>::invalid())
            //         {
            //             // Schreibe die lokale Interface ID hinein
            //             vecInterfaceMapStructure.push_back(localInterfaceID);
            //         }
            //
            //         localInterfaceID = localInterfaceID + 1;
            //
            //     }
            // }
            //
            // // Am Ende steht in localInterfaceID wie viele Interface-Knoten es insgesamt gibt
            // GO numberInterfaceNodes = localInterfaceID; // long long wg. 64
            //
            // // Baue nun die InterfaceMap fuer Fluid und Struktur
            //
            // Teuchos::ArrayView<GO> vecInterfaceMapFluidArray =  Teuchos::arrayViewFromVector( vecInterfaceMapFluid );
            // MapPtr_Type interfaceMapFluid = rcp(new Map_Type( numberInterfaceNodes, vecInterfaceMapFluidArray, 0, comm ) ); //maybe numberInterfaceNodes instead of -1
            //
            // Teuchos::ArrayView<GO> vecInterfaceMapStructureArray =  Teuchos::arrayViewFromVector( vecInterfaceMapStructure );
            // MapPtr_Type interfaceMapStructure = rcp(new Map_Type( numberInterfaceNodes, vecInterfaceMapStructureArray, 0, comm ) ); //maybe numberInterfaceNodes

            // Baue die Interface-Maps in der Interface-Nummerierung
            domainFluid->buildUniqueInterfaceMaps();
            domainStruct->buildUniqueInterfaceMaps();

            // Hole die dof-Map
            MapConstPtr_Type interfaceMapFluidVecField = domainFluid->getInterfaceMapVecFieldUnique();
            MapConstPtr_Type interfaceMapStructureVecField = domainStruct->getInterfaceMapVecFieldUnique();

            // TODO: Noch fehlerhaft domainFluid, wie es scheint.

            // vec_long_Type localInterfaceIDinFluid = domainFluid->getLocalInterfaceIDInGlobal();
            // vec_long_Type localInterfaceIDinStructure = domainStruct->getLocalInterfaceIDInGlobal();

            MeshUnstrPtr_Type meshUnstructuredFluid = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainFluid->getMesh() );
            vec3D_GO_ptr_Type indicesGlobalMatchedOriginFluid = meshUnstructuredFluid->getMeshInterface()->getIndicesGlobalMatchedOrigin();

            // Strukturloesung holen
            MultiVectorConstPtr_Type struc_sol_unique = LinElas.getSolution()->getBlock(0);

            // Extrahiere nun aus der globalen Strukturloesung die Loesung auf dem Interface (in parallel).
            MultiVectorPtr_Type interfaceSolutionStruct = rcp( new MultiVector_Type( interfaceMapStructureVecField, 1 ) );

            {
                int flagCounter = 0;
                Teuchos::ArrayRCP< SC > valuesInterface = interfaceSolutionStruct->getDataNonConst(0); //single MultiVector
                Teuchos::ArrayRCP< SC > valuesStructure = struc_sol_unique->getDataNonConst(0); //single MultiVector

                for(UN i = 0; i < valuesInterface.size(); i++)
                {
                    GO interfaceID = interfaceMapStructureVecField->getGlobalElement( i ); // ID (vektorwertig, also dofID) in der Interface-Nummerierung
                    GO nodeID; // nodeID der interfaceID in der Interface-Nummerierung
                    LO localDofNumber; // Ranges from 0 to dofs-1
                    toNodeID(dim, interfaceID, nodeID, localDofNumber); //This function assumes NodeWise ordering.

                    // Ggf. ist die ID auf einer anderen Flag.
                    // Bei nur einer Interface-Flag kommt man nicht hier hinein.
                    while( nodeID > indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size()-1 )
                    {
                        nodeID = nodeID - indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size();
                        flagCounter = flagCounter + 1; // hier steht dann die korrekte Flag
                    }

                    // Beachte: at(1) ist Struktur!!!
                    // GlobaleID des Interface-Knotens in der Struktur-Nummerierung
                    GO globalInterfaceIDNode = indicesGlobalMatchedOriginFluid->at(flagCounter).at(1).at(nodeID);
                    GO globalInterfaceIDinStructure; // dofID
                    toDofID(dim, globalInterfaceIDNode, localDofNumber, globalInterfaceIDinStructure);

                    // LokaleID auf dem Prozessor vom Interface-Knoten.
                    LO localInterfaceIDinStructure = domainStruct->getMapVecFieldUnique()->getLocalElement(globalInterfaceIDinStructure);

                    valuesInterface[i] = valuesStructure[localInterfaceIDinStructure];
                }
            }



            // {
            //     Teuchos::ArrayRCP< SC > valuesInterface = interfaceSolutionStruct->getDataNonConst(0); //single MultiVector
            //     Teuchos::ArrayRCP< SC > valuesStructure = struc_sol_unique->getDataNonConst(0); //single MultiVector
            //
            //     for(UN i = 0; i < valuesInterface.size(); i++)
            //     {
            //         valuesInterface[i] = valuesStructure[localInterfaceIDinStructure.at(i)];
            //     }
            // }



            // Baue Import indem wir eine TargetMap und eine SourceMap bauen/ angeben
            // SourceMap: Welche Indizes besitze ich (= der derzeitige Prozzesor)
            // TargetMap: Welche Indizes benoetige ich (= der derzeitge Prozessor)
            // Auf gut Deutsch also: Welche Indizes will ich (= der Prozessor) importiert haben/ am Ende besitzen.
            // Wir umgehen hiermit das explizite Aufstellen des Imports.
            MultiVectorPtr_Type interfaceSolutionFluid = rcp( new MultiVector_Type( interfaceMapFluidVecField, 1 ) );
            interfaceSolutionFluid->importFromVector( interfaceSolutionStruct );

            // ########################
            // Geometry-Problem loesen
            // ########################
            // Hole den Extension-Operator auf dem Fluid-Operator
            Geometry<SC,LO,GO,NO> Geometry( domainFluid, discType, parameterListAll );

            {
                Teuchos::TimeMonitor solveTimeMonitor(*solveTime);

                Geometry.addBoundaries(bcFactoryGeometry); // Dem Problem RW hinzufuegen

                Geometry.initializeProblem();
                // ######################
                // Matrix assemblieren, RW setzen und System loesen
                // ######################
                Geometry.assemble();

                // Geometry.SetBoundaries(); // In der Klasse Problem

                // Strukturlsg als RW setzen (per Hand)
                Geometry.setBoundariesSystem(); // RW im System setzen (mit den Flags von oben)
                Geometry.getRhs()->putScalar(0.0);

                Teuchos::ArrayRCP< SC > valuesInterface = interfaceSolutionFluid->getDataNonConst(0); //single MultiVector
                Teuchos::ArrayRCP< SC > valuesFluidRhs = Geometry.getRhs()->getBlock(0)->getDataNonConst(0); //single MultiVector

                int flagCounter = 0;
                for(UN i = 0; i < valuesInterface.size(); i++)
                {
                    GO interfaceID = interfaceMapFluidVecField->getGlobalElement( i ); // dofID in der Interface-Nummerierung
                    GO nodeID; // dofID der interfaceID in der Interface-Nummerierung
                    LO localDofNumber; // Ranges from 0 to dofs-1
                    toNodeID(dim, interfaceID, nodeID, localDofNumber);//This function assumes NodeWise ordering.

                    // Ggf. ist die ID auf einer anderen Flag.
                    // Bei nur einer Interface-Flag kommt man nicht hier hinein.
                    while(nodeID > indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size()-1)
                    {
                        nodeID = nodeID - indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).size();
                        flagCounter = flagCounter + 1; // hier steht dann die korrekte Flag
                    }

                    // Beachte: at(0) ist Fluid!!!
                    // GlobaleID des Interface-Knotens in der Fluid-Nummerierung
                    GO globalInterfaceIDNode = indicesGlobalMatchedOriginFluid->at(flagCounter).at(0).at(nodeID);
                    GO globalInterfaceIDinFluid; // dofID
                    toDofID(dim, globalInterfaceIDNode, localDofNumber, globalInterfaceIDinFluid);

                    // LokaleID auf dem Prozessor des Interface-Knotens.
                    GO localInterfaceIDinFluid = domainFluid->getMapVecFieldUnique()->getLocalElement(globalInterfaceIDinFluid);

                    valuesFluidRhs[localInterfaceIDinFluid] = valuesInterface[i];

                    // valuesFluidRhs[localInterfaceIDinFluid.at(i)] = valuesInterface[i];
                }

                Geometry.solve();
            }

            // ########################
            // Export fuer die Fortsetzungs-Loesung
            // ########################
            if ( parameterListAll->sublist("General").get("ParaViewExport",false) ) {

                ExporterPVPtr_Type exPara(new ExporterPV_Type());
                DomainPtr_Type domain = domainFluid;
                MultiVectorConstPtr_Type exportSolution = Geometry.getSolution()->getBlock(0);
                
                exPara->setup("Extension", domain->getMesh(), discType);

                exPara->addVariable(exportSolution, "Extension", "Vector", dim, domain->getMapUnique());

                exPara->save(0.0);
            }

            // ########################
            // Fluid-Gitter bewegen und exportieren
            // ########################
            // Setze die aktuelle (nicht-deformierte) Konfiguration als Referenzkonfiguration
            domainFluid->setReferenceConfiguration();

            // Extrahiere die Loesung; hier unique, da GLS eindeutig geloest werden muss

            // Geometrieloesung holen
            MultiVectorConstPtr_Type displacementUniqueConst = Geometry.getSolution()->getBlock(0);
            MultiVectorPtr_Type displacementRepeated = rcp( new MultiVector_Type( Geometry.getDomain(0)->getMapVecFieldRepeated() ) );

            displacementRepeated->importFromVector( displacementUniqueConst );
            MultiVectorPtr_Type displacementUnique = rcp_const_cast<MultiVector_Type>(displacementUniqueConst);

            // Verschiebe das Gitter
            domainFluid->moveMesh(displacementUnique, displacementRepeated);

            // Exportiere die berechnete Loesung, jedoch auf dem schon deformierten Gitter (die Loesung ist egal. Wir wollen nur das Gitter haben)
            if ( parameterListAll->sublist("General").get("ParaViewExport",false) ) {

                ExporterPVPtr_Type exPara(new ExporterPV_Type());
                DomainPtr_Type domain = domainFluid;
                MultiVectorConstPtr_Type exportSolution = Geometry.getSolution()->getBlock(0);
                
                exPara->setup("MovedFluidMesh", domain->getMesh(), discType);
                                
                exPara->addVariable(exportSolution, "MovedFluidMesh", "Vector", dim, domain->getMapUnique());

                exPara->save(0.0);

            }
        }
    }

    Teuchos::TimeMonitor::report(std::cout);

    return EXIT_SUCCESS;
}
