#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

void rhsDummy2D(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    return;
}

void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet2D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflow2D(double* x, double* res, double t, const double* parameters)
{
    double H = parameters[1];

    if(t < .5)
    {
        res[0] = (4.0*1.5*parameters[0]*x[1]*(H-x[1])/(H*H)) * 0.5 * ( ( 1 - cos(2.*M_PI*t) ) );
        res[1] = 0.;
    }
    else
    {
        res[0] = (4.0*1.5*parameters[0]*x[1]*(H-x[1])/(H*H));
        res[1] = 0.;
    }

    return;
}

void inflow3DRichter(double* x, double* res, double t, const double* parameters)
{
    double H = parameters[1];

    if(t < 2.)
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) ) * 0.5 * ( ( 1 - cos( M_PI*t/2.0)  ));
        res[1] = 0.;
        res[2] = 0.;
    }
    else
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) );
        res[1] = 0.;
        res[2] = 0.;
    }

    return;
}

void inflow3DRichterFaster(double* x, double* res, double t, const double* parameters)
{
    double H = parameters[1];
    
    if(t < .5)
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) ) * 0.5 * ( ( 1 - cos( 2.*M_PI*t )  ));
        res[1] = 0.;
        res[2] = 0.;
    }
    else
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) );
        res[1] = 0.;
        res[2] = 0.;
    }
    
    return;
}

void inflow3DRichterSuperFast(double* x, double* res, double t, const double* parameters)
{
    double H = parameters[1];
    
    if(t < .1)
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) ) * 0.5 * ( ( 1 - cos( 10.*M_PI*t )  ));
        res[1] = 0.;
        res[2] = 0.;
    }
    else
    {
        res[0] = 9./8 * parameters[0] *x[1]*(H-x[1])*(H*H-x[2]*x[2])/( H*H*(H/2.)*(H/2.) );
        res[1] = 0.;
        res[2] = 0.;
    }
    
    return;
}

void dummyFunc(double* x, double* res, double t, const double* parameters)
{
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

    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
    string xmlProblemFile = "parametersProblemFSI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFileGE = "parametersPrecGE.xml"; // GE
    string xmlPrecFileGI = "parametersPrecGI.xml"; // GI
    myCLP.setOption("precfileGE",&xmlPrecFileGE,".xml file with Inputparameters.");
    myCLP.setOption("precfileGI",&xmlPrecFileGI,".xml file with Inputparameters.");
    string xmlSolverFileFSI = "parametersSolverFSI.xml"; // GI
    myCLP.setOption("solverfileFSI",&xmlSolverFileFSI,".xml file with Inputparameters.");
    string xmlSolverFileGeometry = "parametersSolverGeometry.xml"; // GE
    myCLP.setOption("solverfileGeometry",&xmlSolverFileGeometry,".xml file with Inputparameters.");

    string xmlPrecFileFluidMono = "parametersPrecFluidMono.xml";
    string xmlPrecFileFluidTeko = "parametersPrecFluidTeko.xml";
    myCLP.setOption("precfileFluidMono",&xmlPrecFileFluidMono,".xml file with Inputparameters.");
    myCLP.setOption("precfileFluidTeko",&xmlPrecFileFluidTeko,".xml file with Inputparameters.");
    string xmlProblemFileFluid = "parametersProblemFluid.xml";
    myCLP.setOption("problemFileFluid",&xmlProblemFileFluid,".xml file with Inputparameters.");
    string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileGeometry = "parametersPrecGeometry.xml";
    myCLP.setOption("precfileGeometry",&xmlPrecFileGeometry,".xml file with Inputparameters.");
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }

    bool verbose (comm->getRank() == 0);

    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrecGE = Teuchos::getParametersFromXmlFile(xmlPrecFileGE);
        ParameterListPtr_Type parameterListPrecGI = Teuchos::getParametersFromXmlFile(xmlPrecFileGI);
        ParameterListPtr_Type parameterListSolverFSI = Teuchos::getParametersFromXmlFile(xmlSolverFileFSI);
        ParameterListPtr_Type parameterListSolverGeometry = Teuchos::getParametersFromXmlFile(xmlSolverFileGeometry);
//        ParameterListPtr_Type parameterListProblemFluid = Teuchos::getParametersFromXmlFile(xmlProblemFileFluid);
//        ParameterListPtr_Type parameterListProblemStructure = Teuchos::getParametersFromXmlFile("parametersProblemStructure.xml");
        ParameterListPtr_Type parameterListPrecGeometry = Teuchos::getParametersFromXmlFile(xmlPrecFileGeometry);

        ParameterListPtr_Type parameterListPrecFluidMono = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidMono);
        ParameterListPtr_Type parameterListPrecFluidTeko = Teuchos::getParametersFromXmlFile(xmlPrecFileFluidTeko);

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        
        bool geometryExplicit = parameterListProblem->sublist("Parameter").get("Geometry Explicit",true);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        if(geometryExplicit)
            parameterListAll->setParameters(*parameterListPrecGE);
        
        else
            parameterListAll->setParameters(*parameterListPrecGI);
        
        parameterListAll->setParameters(*parameterListSolverFSI);

        
        ParameterListPtr_Type parameterListFluidAll(new Teuchos::ParameterList(*parameterListPrecFluidMono)) ;
        sublist(parameterListFluidAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Fluid") );
        parameterListFluidAll->setParameters(*parameterListPrecFluidTeko);

        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );

        parameterListStructureAll->setParameters(*parameterListPrecStructure);


        
        // Fuer das Geometrieproblem, falls GE
        // CH: We might want to add a paramterlist, which defines the Geometry problem
        ParameterListPtr_Type parameterListGeometry(new Teuchos::ParameterList(*parameterListPrecGeometry));
        parameterListGeometry->setParameters(*parameterListSolverGeometry);
        // we only compute the preconditioner for the geometry problem once
        sublist( parameterListGeometry, "General" )->set( "Preconditioner Method", "MonolithicConstPrec" );
        sublist( parameterListGeometry, "Parameter" )->set( "Model", parameterListProblem->sublist("Parameter").get("Model Geometry","Laplace") );
        
        sublist( parameterListGeometry, "Parameter" )->set( "Poisson Ratio", 0.4 );
        sublist( parameterListGeometry, "Parameter" )->set( "Mu", 2.0e+6 );
            
        int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        
        string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        string preconditionerMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;

        TimePtr_Type totalTime(TimeMonitor_Type::getNewCounter("FEDD - main - Total Time"));
        TimePtr_Type buildMesh(TimeMonitor_Type::getNewCounter("FEDD - main - Build Mesh"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        // #####################
        // Mesh bauen und wahlen
        // #####################
        {
            if (verbose)
            {
                cout << "###############################################" <<endl;
                cout << "############ Starting FSI  ... ################" <<endl;
                cout << "###############################################" <<endl;
            }

            DomainPtr_Type domainP1fluid;
            DomainPtr_Type domainP1struct;
            DomainPtr_Type domainP2fluid;
            DomainPtr_Type domainP2struct;
          
            std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","parabolic");
            
            {
                TimeMonitor_Type totalTimeMonitor(*totalTime);
                {
                    TimeMonitor_Type buildMeshMonitor(*buildMesh);
                    if (verbose)
                    {
                        cout << " -- Building Mesh ... " << flush;
                    }

                    domainP1fluid.reset( new Domain_Type( comm, dim ) );
                    domainP1struct.reset( new Domain_Type( comm, dim ) );
                    domainP2fluid.reset( new Domain_Type( comm, dim ) );
                    domainP2struct.reset( new Domain_Type( comm, dim ) );
                    //                    
                    if (!meshType.compare("unstructured")) {

                        vec_int_Type idsInterface(1,-1);
                        if (bcType == "partialCFD")
                            idsInterface[0] = 5;
                        else if ( bcType == "Richter3D" || bcType == "Richter3DFull" || bcType == "Richter3DFullFaster" || bcType == "Richter3DFullSuperFast"){
                            idsInterface[0] = 6;
                        }
                        else if (bcType == "Tube2D"){
                            idsInterface[0] = 4;
                        }
                        else{
                            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"Interface ID not known for this problem.");
                        }

                        
                        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
                        domainP1Array[0] = domainP1fluid;
                        domainP1Array[1] = domainP1struct;
                        
                        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
                        if (!discType.compare("P2")){
                            pListPartitioner->set("Build Edge List",true);
                            pListPartitioner->set("Build Surface List",true);
                        }
                        else{
                            pListPartitioner->set("Build Edge List",false);
                            pListPartitioner->set("Build Surface List",false);
                        }
                        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                        
                        partitionerP1.readAndPartition();
                        
                        if (!discType.compare("P2")){
                            domainP2fluid->buildP2ofP1Domain( domainP1fluid );
                            domainP2struct->buildP2ofP1Domain( domainP1struct );
                        }
                        // Calculate distances is done in: identifyInterfaceParallelAndDistance
                        domainP1fluid->identifyInterfaceParallelAndDistance(domainP1struct, idsInterface);
                        if (!discType.compare("P2"))
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

            DomainPtr_Type domainFluidVelocity;
            DomainPtr_Type domainFluidPressure;
            DomainPtr_Type domainStructure;
            DomainPtr_Type domainGeometry;
            if (!discType.compare("P2"))
            {
                domainFluidVelocity = domainP2fluid;
                domainFluidPressure = domainP1fluid;
                domainStructure = domainP2struct;
                domainGeometry = domainP2fluid;
            }
            else
            {
                domainFluidVelocity = domainP1fluid;
                domainFluidPressure = domainP1fluid;
                domainStructure = domainP1struct;
                domainGeometry = domainP1fluid;
//                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"P1/P1 for FSI not implemented!");
            }


            if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
                
                if (verbose)
                    std::cout << "\t### Exporting fluid and solid subdomains ###\n";

                typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
                typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
                typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
                typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
                typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

                {
                    MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domainFluidVelocity->getElementMap() ) );
                    MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                    vecDecomposition->putScalar(comm->getRank()+1.);
                    
                    Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                    
                    exPara->setup( "subdomains_fluid", domainFluidVelocity->getMesh(), "P0" );
                    
                    exPara->addVariable( vecDecompositionConst, "subdomains", "Scalar", 1, domainFluidVelocity->getElementMap());
                    exPara->save(0.0);
                    exPara->closeExporter();
                }
                {
                    MultiVectorPtr_Type vecDecomposition = rcp(new MultiVector_Type( domainStructure->getElementMap() ) );
                    MultiVectorConstPtr_Type vecDecompositionConst = vecDecomposition;
                    vecDecomposition->putScalar(comm->getRank()+1.);
                    
                    Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                    
                    exPara->setup( "subdomains_solid", domainStructure->getMesh(), "P0" );
                    
                    exPara->addVariable( vecDecompositionConst, "subdomains", "Scalar", 1, domainStructure->getElementMap());
                    exPara->save(0.0);
                    exPara->closeExporter();
                }

            }
            
            // Baue die Interface-Maps in der Interface-Nummerierung
            domainFluidVelocity->buildInterfaceMaps();
            
            domainStructure->buildInterfaceMaps();

            // domainInterface als dummyDomain mit mapVecFieldRepeated_ als interfaceMapVecFieldUnique_.
            // Wird fuer den Vorkonditionierer und Export gebraucht.
            // mesh is needed for rankRanges
            DomainPtr_Type domainInterface;
            domainInterface.reset( new Domain_Type( comm ) );
            domainInterface->setDummyInterfaceDomain(domainFluidVelocity);
            
            domainFluidVelocity->setReferenceConfiguration();
            domainFluidPressure->setReferenceConfiguration();

            {
            ///////////////////////////
            // Mesh Export fuer Matlab
            ///////////////////////////
//            ofstream myFile;
//            FILE * pFile;
//            std::cout << "FluidP2Elements..." << '\n';
//            myFile.open ("FluidP2ElementsH00.txt");
//            for(int i = 0; i < domainFluidVelocity->getElements()->size(); i++)
//            {
//                for(int j= 0; j < domainFluidVelocity->getElements()->at(i).size(); j++)
//                {
//                    myFile << domainFluidVelocity->getElements()->at(i).at(j);
//                    myFile << " ";
//                }
//                myFile << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "FluidP1Elements..." << '\n';
//            myFile.open ("FluidP1ElementsH00.txt");
//            for(int i = 0; i < domainFluidPressure->getElements()->size(); i++)
//            {
//                for(int j= 0; j < domainFluidPressure->getElements()->at(i).size(); j++)
//                {
//                    myFile << domainFluidPressure->getElements()->at(i).at(j);
//                    myFile << " ";
//                }
//                myFile << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "StrucP2Elements..." << '\n';
//            myFile.open ("StrucP2ElementsH00.txt");
//            for(int i = 0; i < domainStructure->getElements()->size(); i++)
//            {
//                for(int j= 0; j < domainStructure->getElements()->at(i).size(); j++)
//                {
//                    myFile << domainStructure->getElements()->at(i).at(j);
//                    myFile << " ";
//                }
//                myFile << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "FluidP2Nodes..." << '\n';
////            myFile.open ("FluidP2NodesH00.txt");
//            pFile = fopen ("FluidP2NodesH00.txt","w");
//
//            for(int i = 0; i < domainFluidVelocity->getPointsUnique()->size(); i++)
//            {
//                for(int j= 0; j < domainFluidVelocity->getPointsUnique()->at(i).size(); j++)
//                {
//                    fprintf(pFile,"%4.10f ", domainFluidVelocity->getPointsUnique()->at(i).at(j) );
//
////                    printf("%4.10f ", domainFluidVelocity->getPointsUnique()->at(i).at(j) );
////                    myFile << domainFluidVelocity->getPointsUnique()->at(i).at(j);
////                    myFile << " ";
//                }
//                fprintf(pFile,"\n");
////                myFile << endl;
//            }
//            fclose(pFile);
////            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "FluidP1Nodes..." << '\n';
////            myFile.open ("FluidP1NodesH00.txt");
//            pFile = fopen ("FluidP1NodesH00.txt","w");
//
//            for(int i = 0; i < domainFluidPressure->getPointsUnique()->size(); i++)
//            {
//                for(int j= 0; j < domainFluidPressure->getPointsUnique()->at(i).size(); j++)
//                {
//                    fprintf(pFile,"%4.10f ", domainFluidPressure->getPointsUnique()->at(i).at(j) );
////                    myFile << domainFluidPressure->getPointsUnique()->at(i).at(j);
////                    myFile << " ";
//                }
//                fprintf(pFile,"\n");
////                myFile << endl;
//            }
//            fclose(pFile);
////            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "StrucP2Nodes..." << '\n';
//            
//            pFile = fopen ("StrucP2NodesH00.txt","w");
//            for(int i = 0; i < domainStructure->getPointsUnique()->size(); i++)
//            {
//                for(int j= 0; j < domainStructure->getPointsUnique()->at(i).size(); j++)
//                {
//                    fprintf(pFile,"%4.10f ", domainStructure->getPointsUnique()->at(i).at(j) );
////                    myFile << domainStructure->getPointsUnique()->at(i).at(j);
////                    myFile << " ";
//                }
//                fprintf(pFile,"\n");
//            }
//            fclose(pFile);
////            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "FluidP2Flags..." << '\n';
//            myFile.open ("FluidP2FlagsH00.txt");
//            for(int i = 0; i < domainFluidVelocity->getBCFlagUnique()->size(); i++)
//            {
//                myFile << domainFluidVelocity->getBCFlagUnique()->at(i) << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "FluidPFlags..." << '\n';
//            myFile.open ("FluidP1FlagsH00.txt");
//            for(int i = 0; i < domainFluidPressure->getBCFlagUnique()->size(); i++)
//            {
//                myFile << domainFluidPressure->getBCFlagUnique()->at(i) << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
//            
//            std::cout << "StrucP2Flags..." << '\n';
//            myFile.open ("StrucP2FlagsH00.txt");
//            for(int i = 0; i < domainStructure->getBCFlagUnique()->size(); i++)
//            {
//                myFile << domainStructure->getBCFlagUnique()->at(i) << endl;
//            }
//            myFile.close();
//            std::cout << "done" << '\n';
        }
                            
            // #####################
            // Problem definieren
            // #####################
            Teuchos::RCP<SmallMatrix<int>> defTS;
            if(geometryExplicit)
            {
                // SmallMatrix<int> defTS(4);
                defTS.reset( new SmallMatrix<int> (4) );

                // Fluid
                (*defTS)[0][0] = 1;
                (*defTS)[0][1] = 1;

                // Struktur
                (*defTS)[2][2] = 1;
            }
            else
            {
                // SmallMatrix<int> defTS(5);
                defTS.reset( new SmallMatrix<int> (5) );

                // Fluid
                (*defTS)[0][0] = 1;
                (*defTS)[0][1] = 1;

                // TODO: [0][4] und [1][4] bei GI + Newton noetig?

                // Struktur
                (*defTS)[2][2] = 1;
            }

            FSI<SC,LO,GO,NO> fsi(domainFluidVelocity, discType,
                                 domainFluidPressure, "P1",
                                 domainStructure, discType,
                                 domainInterface, discType,
                                 domainGeometry, discType,
                                 parameterListFluidAll,
                                 parameterListStructureAll,
                                 parameterListAll,
                                 parameterListGeometry,
                                 defTS);


            domainFluidVelocity->info();
            domainFluidPressure->info();
            domainStructure->info();
            domainGeometry->info();
            
            fsi.info();
            
            // #####################
            // Randwerte
            // Zusammenfassung der Flags:
            // Fluid/ Geometrie: 1 = wall; 2 = inflow; 3 = outflow; 4 = Fluid-obstacle; 5 = Interface in 2d
            // Struktur: 1 = linke Seite (homogener Dirichletrand); 5 = Interface in 2d
            // #####################
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("MeanVelocity",2.0));
            
            if(!bcType.compare("partialCFD"))
            {
                parameter_vec.push_back(.41);//height of inflow region
            }
            else if( !bcType.compare("Richter3D") || !bcType.compare("Richter3DFull") || !bcType.compare("Richter3DFullFaster") || !bcType.compare("Richter3DFullSuperFast") )
            {
                parameter_vec[0] = 1.;
                parameter_vec.push_back(.4);
            }
            else if(!bcType.compare("Tube2D"))
            {
                parameter_vec.push_back(1.);
                parameter_vec[0] = 1.5;
            }
            else
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Select a valid boundary condition. Only partialCFD available.");
            }


            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );


            // TODO: Vermutlich braucht man keine bcFactoryFluid und bcFactoryStructure,
            // da die RW sowieso auf dem FSI-Problem gesetzt werden.

            // Fluid-RW
            {
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluid( new BCBuilder<SC,LO,GO,NO>( ) );
                
                if (dim==2)
                {
                    if(!bcType.compare("partialCFD"))
                    {
                        bcFactory->addBC(zeroDirichlet2D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                        bcFactory->addBC(inflow2D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow

                        bcFactoryFluid->addBC(zeroDirichlet2D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                        bcFactoryFluid->addBC(inflow2D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow

                        bcFactory->addBC(zeroDirichlet2D, 4, 0, domainFluidVelocity, "Dirichlet", dim); // obstacle
                        bcFactoryFluid->addBC(zeroDirichlet2D, 4, 0, domainFluidVelocity, "Dirichlet", dim); // obstacle                        
                        
                    }
                    else if(bcType == "Tube2D"){
                        bcFactory->addBC(zeroDirichlet2D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                        bcFactory->addBC(inflow2D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow
                        bcFactoryFluid->addBC(zeroDirichlet2D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                        bcFactoryFluid->addBC(inflow2D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow
                    }
                    
                }
                else if(dim==3)
                {
                    bcFactory->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                    if (bcType == "Richter3D" || bcType == "Richter3DFull")
                        bcFactory->addBC(inflow3DRichter, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow
                    else if (bcType == "Richter3DFullFaster")
                        bcFactory->addBC(inflow3DRichterFaster, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow
                    else if (bcType == "Richter3DFullSuperFast")
                        bcFactory->addBC(inflow3DRichterSuperFast, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow
                    
                    bcFactoryFluid->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                    if (bcType == "Richter3D" || bcType == "Richter3DFull")
                        bcFactoryFluid->addBC(inflow3DRichter, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);
                    else if (bcType == "Richter3DFullFaster")
                        bcFactoryFluid->addBC(inflow3DRichterFaster, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);
                    else if (bcType == "Richter3DFullSuperFast")
                        bcFactoryFluid->addBC(inflow3DRichterSuperFast, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);
                }

                // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
                // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
                fsi.problemFluid_->addBoundaries(bcFactoryFluid);
                

            }

            // Struktur-RW
            {
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );

                if(dim == 2)
                {
                    bcFactory->addBC(zeroDirichlet2D, 1, 2, domainStructure, "Dirichlet", dim); // linke Seite
                    bcFactoryStructure->addBC(zeroDirichlet2D, 1, 0, domainStructure, "Dirichlet", dim); // linke Seite
                }
                else if(dim == 3)
                {
                    bcFactory->addBC(zeroDirichlet3D, 1, 2, domainStructure, "Dirichlet", dim); // linke Seite
                    bcFactoryStructure->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet", dim); // linke Seite
                }
                
                // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
                // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
                if (!fsi.problemStructure_.is_null())
                    fsi.problemStructure_->addBoundaries(bcFactoryStructure);
                else
                    fsi.problemStructureNonLin_->addBoundaries(bcFactoryStructure);
            }
            // RHS dummy for structure
            if (dim==2) {
                if (!fsi.problemStructure_.is_null())
                    fsi.problemStructure_->addRhsFunction( rhsDummy2D );
                else
                    fsi.problemStructureNonLin_->addRhsFunction( rhsDummy2D );
                
            }
            else if (dim==3) {
             
                if (!fsi.problemStructure_.is_null())
                    fsi.problemStructure_->addRhsFunction( rhsDummy );
                else
                    fsi.problemStructureNonLin_->addRhsFunction( rhsDummy );
                
            }
            
            

            // Geometrie-RW separat, falls geometrisch explizit.
            // Bei Geometrisch implizit: Keine RW in die factoryFSI fuer das
            // Geometrie-Teilproblem, da sonst (wg. dem ZeroDirichlet auf dem Interface,
            // was wir brauchen wegen Kopplung der Struktur) der Kopplungsblock C4
            // in derselben Zeile, der nur Werte auf dem Interface haelt, mit eliminiert.
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryGeometry( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluidInterface;
            if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
                bcFactoryFluidInterface = Teuchos::rcp( new BCBuilder<SC,LO,GO,NO>( ) );

            if(dim == 2)
            {
                if(!bcType.compare("partialCFD")) {
                bcFactoryGeometry->addBC(zeroDirichlet2D, 1, 0, domainGeometry, "Dirichlet", dim); // wall
                bcFactoryGeometry->addBC(zeroDirichlet2D, 2, 0, domainGeometry, "Dirichlet", dim); // inflow
                bcFactoryGeometry->addBC(zeroDirichlet2D, 3, 0, domainGeometry, "Dirichlet", dim); // outflow
                bcFactoryGeometry->addBC(zeroDirichlet2D, 4, 0, domainGeometry, "Dirichlet", dim); // fluid-obstacle
                bcFactoryGeometry->addBC(zeroDirichlet2D, 5, 0, domainGeometry, "Dirichlet", dim); // 5 ist die Interface-Flag
                if ( preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko" )
                    bcFactoryFluidInterface->addBC(zeroDirichlet2D, 5, 0, domainFluidVelocity, "Dirichlet", dim);
                // Die RW, welche nicht Null sind in der rechten Seite (nur Interface) setzen wir spaeter per Hand
                }
                else if(!bcType.compare("Tube2D")){
                    bcFactoryGeometry->addBC(zeroDirichlet2D, 1, 0, domainGeometry, "Dirichlet", dim); // wall
                    bcFactoryGeometry->addBC(zeroDirichlet2D, 2, 0, domainGeometry, "Dirichlet", dim); // inflow
                    bcFactoryGeometry->addBC(zeroDirichlet2D, 3, 0, domainGeometry, "Dirichlet", dim); // outflow
                    bcFactoryGeometry->addBC(zeroDirichlet2D, 4, 0, domainGeometry, "Dirichlet", dim); // interface
                }
            }
            else if(dim == 3)
            {
//                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "fix Randwerte für Richter Benchmark");
                bcFactoryGeometry->addBC(zeroDirichlet3D, 1, 0, domainGeometry, "Dirichlet", dim); // wall
                bcFactoryGeometry->addBC(zeroDirichlet3D, 2, 0, domainGeometry, "Dirichlet", dim); // inflow
                if (bcType == "Richter3D"){
                    bcFactoryGeometry->addBC(zeroDirichlet3D, 3, 0, domainGeometry, "Dirichlet", dim); // interface line and symmetry wall
                    bcFactoryGeometry->addBC(zeroDirichlet3D, 4, 0, domainGeometry, "Dirichlet", dim); // symmetry-wall
                }
                bcFactoryGeometry->addBC(zeroDirichlet3D, 5, 0, domainGeometry, "Dirichlet", dim); // outflow
                bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // interface
                if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
                    bcFactoryFluidInterface->addBC(zeroDirichlet2D, 6, 0, domainFluidVelocity, "Dirichlet", dim);
                // Die RW, welche nicht Null sind in der rechten Seite (nur Interface) setzen wir spaeter per Hand
            }

            fsi.problemGeometry_->addBoundaries(bcFactoryGeometry);
            if ( preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
                fsi.getPreconditioner()->setFaCSIBCFactory( bcFactoryFluidInterface );


            // #####################
            // Zeitintegration
            // #####################
            fsi.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

            fsi.initializeProblem();
            fsi.initializeGE();
            // Matrizen assemblieren
            fsi.assemble();
                        
            DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

            // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
            // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
            daeTimeSolver.defineTimeStepping(*defTS);

            // Uebergebe das (nicht) lineare Problem
            daeTimeSolver.setProblem(fsi);

            // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
            // defTS definiert worden sind.
            daeTimeSolver.setupTimeStepping();

            daeTimeSolver.advanceInTime();
        }
    }

    TimeMonitor_Type::report(std::cout);

    return(EXIT_SUCCESS);
}
