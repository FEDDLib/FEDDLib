#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 MeshInterface test
 
 @brief  MeshInterface test
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */



void zeroBC(double* x, double* res, double t, const double* parameters){

    res[0] = 0.;

    return;
}

void x1(double* x, double* res, double t, const double* parameters){

    res[0] = 1.;

    return;
}
void x2(double* x, double* res, double t, const double* parameters){

    res[0] = 2.;

    return;
}
void x3(double* x, double* res, double t, const double* parameters){

    res[0] = 3.;

    return;
}
void x4(double* x, double* res, double t, const double* parameters){

    res[0] = 4.;

    return;
}
void x5(double* x, double* res, double t, const double* parameters){

    res[0] = 5.;

    return;
}
void x6(double* x, double* res, double t, const double* parameters){
    
    res[0] = 6.;
    
    return;
}
void x7(double* x, double* res, double t, const double* parameters){
    
    res[0] = 7.;
    
    return;
}
void x10(double* x, double* res, double t, const double* parameters){
    
    res[0] = 10.;
    
    return;
}

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;

using namespace FEDD;
using namespace Teuchos;
using namespace std;
int main(int argc, char *argv[]) {
    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    Teuchos::CommandLineProcessor My_CLP;


    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string xmlFile = "meshes_interface.xml";
    myCLP.setOption("file",&xmlFile,"xmlFile");

    int findInterfaceP2 = 1;
    myCLP.setOption("interface",&findInterfaceP2,"interface");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    bool verbose (comm->getRank() == 0);

    int 		dim				= Teuchos::getParametersFromXmlFile(xmlFile)->sublist("Parameter").get("Dimension",2);
    string		meshType    	= Teuchos::getParametersFromXmlFile(xmlFile)->sublist("Parameter").get("Mesh Type","unstructured");
    string		meshName_fluid    	= Teuchos::getParametersFromXmlFile(xmlFile)->sublist("Parameter").get("Mesh Name Fluid","dfg_fsi_fluid_h004.mesh");
    string		meshName_struc    	= Teuchos::getParametersFromXmlFile(xmlFile)->sublist("Parameter").get("Mesh Name Structure","dfg_fsi_solid_h004.mesh");
    string		meshDelimiter   = Teuchos::getParametersFromXmlFile(xmlFile)->sublist("Parameter").get("Mesh Delimiter"," ");
    string      discType = Teuchos::getParametersFromXmlFile(xmlFile)->sublist("Parameter").get("Discretization","P1");

    ParameterListPtr_Type pList = Teuchos::getParametersFromXmlFile(xmlFile);
    
    Teuchos::RCP<Teuchos::Time> totalTime(Teuchos::TimeMonitor::getNewCounter("main: Total Time"));
    Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Build Mesh"));
    Teuchos::RCP<Teuchos::Time> solveTime(Teuchos::TimeMonitor::getNewCounter("main: Solve problem time"));

    int numProcsCoarseSolve = Teuchos::getParametersFromXmlFile(xmlFile)->sublist("General").get("Mpi Ranks Coarse",0);

    int size = comm->getSize() - numProcsCoarseSolve;

    int nmbInterfaceFlags = Teuchos::getParametersFromXmlFile(xmlFile)->sublist("Parameter").get("Interface Flags",1);
    int partition = 0;
    {
        if (verbose) {
            cout << "#################################################" <<endl;
            cout << "############ Starting meshes interfaces ... #####" <<endl;
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
                if (verbose) {
                    cout << "-- Building Mesh ... " << flush;
                }

                domainP1fluid.reset( new Domain_Type( comm, dim ) );
                domainP1struct.reset( new Domain_Type( comm, dim ) );
                domainP2fluid.reset( new Domain_Type( comm, dim ) );
                domainP2struct.reset( new Domain_Type( comm, dim ) );

                if (!meshType.compare("unstructured")) {

                    int volumeID = 10; // default value which would be used anyway

                    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
                    MeshPartitioner_Type::DomainPtrArray_Type domainP2Array(2);
                    domainP1Array[0] = domainP1fluid;
                    domainP1Array[1] = domainP1struct;
                    
                    ParameterListPtr_Type pListPartitioner = sublist( pList, "Mesh Partitioner" );
                    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                    
                    partitionerP1.readAndPartition();

                    domainP2fluid->buildP2ofP1Domain( domainP1fluid );
                    domainP2struct->buildP2ofP1Domain( domainP1struct );

                    vec_int_Type idsInterface(1,5);
                    if (dim==3) {
                        idsInterface[0] = 6;
                        if (nmbInterfaceFlags==2)
                            idsInterface.push_back(3);
                    }
                    

                    if (verbose){
                        cout << "\n";
                        cout << "-- Determine P1 interface ... " << flush;
                    }
                    domainP1fluid->identifyInterfaceParallelAndDistance(domainP1struct, idsInterface);
                    if (verbose){
                        cout << " done --" << endl;
                        cout << "-- Determine P2 interface ... " << flush;
                    }
                    if (findInterfaceP2)
                        domainP2fluid->identifyInterfaceParallelAndDistance(domainP2struct, idsInterface);

                    
                    if (verbose)
                        cout << " done --" << endl;

                }
                else
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Test for unstructured meshes read from .mesh-file. Change mesh type in setup file to 'unstructured'.");
                if (verbose)
                    cout << "done! -- " << endl;
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
        bcFactoryFluid->addBC(x6, 5, 0, domainFluid, "Dirichlet", 1);
        bcFactoryFluid->addBC(x7, 6, 0, domainFluid, "Dirichlet", 1);
        bcFactoryFluid->addBC(x10, 10, 0, domainFluid, "Dirichlet", 1);
        
        bcFactoryStruct->addBC(x1, 0, 0, domainStruct, "Dirichlet", 1);
        bcFactoryStruct->addBC(x2, 1, 0, domainStruct, "Dirichlet", 1);
        bcFactoryStruct->addBC(x3, 2, 0, domainStruct, "Dirichlet", 1);
        bcFactoryStruct->addBC(x4, 3, 0, domainStruct, "Dirichlet", 1);
        bcFactoryStruct->addBC(x5, 4, 0, domainStruct, "Dirichlet", 1);
        bcFactoryStruct->addBC(x6, 5, 0, domainStruct, "Dirichlet", 1);
        bcFactoryStruct->addBC(x7, 6, 0, domainStruct, "Dirichlet", 1);
        bcFactoryStruct->addBC(x10, 10, 0, domainStruct, "Dirichlet", 1);

        MultiVectorPtr_Type valuesFluid = rcp(new MultiVector_Type( domainFluid->getMapUnique() ) );
        MultiVectorPtr_Type valuesStruct = rcp(new MultiVector_Type( domainStruct->getMapUnique() ) );
        BlockMultiVectorPtr_Type valuesFluidBlock = rcp(new BlockMultiVector_Type( 1 ) );
        BlockMultiVectorPtr_Type valuesStructBlock = rcp(new BlockMultiVector_Type( 1 ) );

        valuesFluidBlock->addBlock( valuesFluid, 0 );
        valuesStructBlock->addBlock( valuesStruct, 0 );

        bcFactoryFluid->setRHS( valuesFluidBlock );
        bcFactoryStruct->setRHS( valuesStructBlock );

        bool exportParaView = Teuchos::getParametersFromXmlFile(xmlFile)->sublist("General").get("ParaViewExport",true);
        if (exportParaView) {
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            DomainPtr_Type domain = domainFluid;
            
            exPara->setup("fluid", domain->getMesh(), discType);
            
            MultiVectorConstPtr_Type valuesFluidConst = valuesFluidBlock->getBlock( 0 );
            exPara->addVariable( valuesFluidConst, "values", "Scalar", 1, domain->getMapUnique());

            exPara->save(0.0);
            exPara->closeExporter();
        }

        bool exportParaViewSubomains = Teuchos::getParametersFromXmlFile(xmlFile)->sublist("General").get("ParaViewExport Subdomains",false);
        if (exportParaViewSubomains) {

            DomainPtr_Type domain = domainP1fluid;
            MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domain->getElementMap() ) );
            MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
            vecDecomposition->putScalar(comm->getRank()+1.);

            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            exPara->setup( "subdomains_fluid", domain->getMesh(), "P0" );
            
            exPara->addVariable( vecDecompositionConst, "subdomain", "Scalar", 1, domain->getElementMap());

            exPara->save(0.0);
            exPara->closeExporter();
        }


        if (exportParaView) {
            Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            DomainPtr_Type domain = domainStruct;

            exPara->setup("solid", domain->getMesh(), discType);
            
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

            exPara->setup( "subdomains_solid", domain->getMesh(), "P0" );
            
            exPara->addVariable( vecDecompositionConst, "subdomain", "Scalar", 1, domain->getElementMap() );


            exPara->save(0.0);
            exPara->closeExporter();
        }



    }

    Teuchos::TimeMonitor::report(cout);

    return(EXIT_SUCCESS);
}
