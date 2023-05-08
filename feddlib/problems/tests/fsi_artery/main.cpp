#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

void zeroBC(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void parabolicInflow3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[0] = parameters[0] / parameters[2] * x[0];
        res[1] = 0.;
        res[2] = 0.;
    }

    return;
}

void parabolicInflow3DLin(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = parameters[0] / parameters[2] * x[0];
        res[1] = 0.;
        res[2] = 0.;
    }

    return;
}

void rhsDummy(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
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
    
    string xmlProbL = "plistProblemLaplace.xml";
    myCLP.setOption("probLaplace",&xmlProbL,".xml file with Inputparameters.");
    string xmlPrecL = "plistPrecLaplace.xml";
    myCLP.setOption("precLaplace",&xmlPrecL,".xml file with Inputparameters.");
    string xmlSolverL = "plistSolverLaplace.xml";
    myCLP.setOption("solverLaplace",&xmlSolverL,".xml file with Inputparameters.");
    
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
        ParameterListPtr_Type parameterListSolverFSI = Teuchos::getParametersFromXmlFile(xmlSolverFileFSI);
        ParameterListPtr_Type parameterListSolverGeometry = Teuchos::getParametersFromXmlFile(xmlSolverFileGeometry);
        ParameterListPtr_Type parameterListPrecGeometry = Teuchos::getParametersFromXmlFile(xmlPrecFileGeometry);

        ParameterListPtr_Type parameterListPrecGE = Teuchos::getParametersFromXmlFile(xmlPrecFileGE);
        ParameterListPtr_Type parameterListPrecGI = Teuchos::getParametersFromXmlFile(xmlPrecFileGI);
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
        
        double poissonRatio = parameterListProblem->sublist("Parameter Geometry").get("Poisson Ratio",0.3);
        double mu = parameterListProblem->sublist("Parameter Geometry").get("Mu",2.0e+6);
        double distanceLaplace = parameterListProblem->sublist("Parameter Geometry").get("Distance Laplace",0.1);
        double coefficientLaplace = parameterListProblem->sublist("Parameter Geometry").get("Coefficient Laplace",1000.);
        
        sublist( parameterListGeometry, "Parameter" )->set( "Poisson Ratio", poissonRatio );
        sublist( parameterListGeometry, "Parameter" )->set( "Mu", mu );
        sublist( parameterListGeometry, "Parameter" )->set( "Distance Laplace", distanceLaplace );
        sublist( parameterListGeometry, "Parameter" )->set( "Coefficient Laplace", coefficientLaplace );
            
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
          
            std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","Compute Inflow");
            
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

                        vec_int_Type idsInterface(1,6);
                                                
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
                if (verbose)
                    std::cout << "### Double check temporal discretization of Shape Derivatives! ###" << std::endl;
                
                (*defTS)[0][4] = 1;
                (*defTS)[1][4] = 1;
                
                // Struktur
                (*defTS)[2][2] = 1;
            }

            FSI<SC,LO,GO,NO> fsi(domainFluidVelocity, discType,
                                 domainFluidPressure, "P1",
                                 domainStructure, discType,
                                 domainInterface, discType,
                                 domainGeometry, discType,
                                 parameterListFluidAll, parameterListStructureAll, parameterListAll,
                                 parameterListGeometry, defTS);


            domainFluidVelocity->info();
            domainFluidPressure->info();
            domainStructure->info();
            domainGeometry->info();
            fsi.info();
                     
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("Max Velocity",1.));
            parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Max Ramp Time",2.) );
            
            TEUCHOS_TEST_FOR_EXCEPTION(bcType != "Compute Inflow", std::logic_error, "Select a valid boundary condition. Only Compute Inflow available.");

            //#############################################
            //#############################################
            //#### Compute parabolic inflow with laplacian
            //#############################################
            //#############################################
            MultiVectorConstPtr_Type solutionLaplace;
            {
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryLaplace(new BCBuilder<SC,LO,GO,NO>( ));
                
                bcFactoryLaplace->addBC(zeroBC, 2, 0, domainFluidVelocity, "Dirichlet", 1); //inflow ring
                bcFactoryLaplace->addBC(zeroBC, 3, 0, domainFluidVelocity, "Dirichlet", 1); //outflow ring
                bcFactoryLaplace->addBC(zeroBC, 6, 0, domainFluidVelocity, "Dirichlet", 1); //surface
                
                ParameterListPtr_Type parameterListProblemL = Teuchos::getParametersFromXmlFile(xmlProbL);
                ParameterListPtr_Type parameterListPrecL = Teuchos::getParametersFromXmlFile(xmlPrecL);
                ParameterListPtr_Type parameterListSolverL = Teuchos::getParametersFromXmlFile(xmlSolverL);

                ParameterListPtr_Type parameterListLaplace(new Teuchos::ParameterList(*parameterListProblemL)) ;
                parameterListLaplace->setParameters(*parameterListPrecL);
                parameterListLaplace->setParameters(*parameterListSolverL);
                
                Laplace<SC,LO,GO,NO> laplace( domainFluidVelocity, discType, parameterListLaplace, false );
                {
                    laplace.addRhsFunction(oneFunc);
                    laplace.addBoundaries(bcFactoryLaplace);
                    
                    laplace.initializeProblem();
                    laplace.assemble();
                    laplace.setBoundaries();
                    laplace.solve();
                }
                
                //We need the values in the inflow area. Therefore, we use the above bcFactory and the volume flag 10 and the outlet flag 5 and set zero Dirichlet boundary values
                bcFactoryLaplace->addBC(zeroBC, 5, 0, domainFluidVelocity, "Dirichlet", 1);
                bcFactoryLaplace->addBC(zeroBC, 10, 0, domainFluidVelocity, "Dirichlet", 1);
                bcFactoryLaplace->setRHS( laplace.getSolution(), 0./*time; does not matter here*/ );
                solutionLaplace = laplace.getSolution()->getBlock(0);
            
                SC maxValue = solutionLaplace->getMax();
                
                parameter_vec.push_back(maxValue);

                Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
                
                exPara->setup("parabolicInflow", domainFluidVelocity->getMesh(), discType);
                
//                exPara->setup(domainFluidVelocity->getDimension(), domainFluidVelocity->getNumElementsGlobal(), domainFluidVelocity->getElements(), domainFluidVelocity->getPointsUnique(), domainFluidVelocity->getMapUnique(), domainFluidVelocity->getMapRepeated(), discType, "parabolicInflow", 1, comm);

                MultiVectorConstPtr_Type valuesConst = laplace.getSolution()->getBlock(0);
                exPara->addVariable( valuesConst, "values", "Scalar", 1, domainFluidVelocity->getMapUnique() );

                exPara->save(0.0);
                exPara->closeExporter();

            }
            
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

            // TODO: Vermutlich braucht man keine bcFactoryFluid und bcFactoryStructure,
            // da die RW sowieso auf dem FSI-Problem gesetzt werden.

            // Fluid-RW
            {
                bool zeroPressure = parameterListProblem->sublist("Parameter Fluid").get("Set Outflow Pressure to Zero",false);
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluid( new BCBuilder<SC,LO,GO,NO>( ) );
                               
                //bcFactory->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                 string rampType = parameterListProblem->sublist("Parameter Fluid").get("Ramp type","cos");
                if (rampType == "cos") {
                    bcFactory->addBC(parabolicInflow3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow ring
                    bcFactory->addBC(parabolicInflow3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow
                    
                    
                    //bcFactoryFluid->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                    bcFactoryFluid->addBC(parabolicInflow3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow ring
                    bcFactoryFluid->addBC(parabolicInflow3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow
                    
                }
                else if(rampType == "linear"){
                    bcFactory->addBC(parabolicInflow3DLin, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow ring
                    bcFactory->addBC(parabolicInflow3DLin, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow
                    
                    //bcFactoryFluid->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                    bcFactoryFluid->addBC(parabolicInflow3DLin, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow ring
                    bcFactoryFluid->addBC(parabolicInflow3DLin, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow
                }
                
                bcFactory->addBC(zeroDirichlet3D, 3, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // outflow ring
                
                bcFactoryFluid->addBC(zeroDirichlet3D, 3, 0, domainFluidVelocity, "Dirichlet", dim); // outflow ring
                
                if (zeroPressure) {
                    bcFactory->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                    bcFactory->addBC(zeroBC, 5, 1, domainFluidPressure, "Dirichlet", 1); // outflow
                    
                    bcFactoryFluid->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                    bcFactoryFluid->addBC(zeroBC, 5, 1, domainFluidPressure, "Dirichlet", 1); // outflow
                }
                
                // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
                // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
                fsi.problemFluid_->addBoundaries(bcFactoryFluid);
            }

            // Struktur-RW
            {
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );
                bcFactory->addBC(zeroDirichlet3D, 2, 2, domainStructure, "Dirichlet", dim); // ring inflow
                bcFactory->addBC(zeroDirichlet3D, 3, 2, domainStructure, "Dirichlet", dim); // ring outflow
                bcFactory->addBC(zeroDirichlet3D, 4, 2, domainStructure, "Dirichlet", dim); // donut inflow
                bcFactory->addBC(zeroDirichlet3D, 5, 2, domainStructure, "Dirichlet", dim); // donut outflow
                
                bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet", dim); // ring inflow
                bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet", dim); // ring outflow
                bcFactoryStructure->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet", dim); // donut inflow
                bcFactoryStructure->addBC(zeroDirichlet3D, 5, 0, domainStructure, "Dirichlet", dim); // donut outflow
                // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
                // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
                if (!fsi.problemStructure_.is_null())
                    fsi.problemStructure_->addBoundaries(bcFactoryStructure);
                else
                    fsi.problemStructureNonLin_->addBoundaries(bcFactoryStructure);
            }
            // RHS dummy for structure
                    
            if (!fsi.problemStructure_.is_null())
                fsi.problemStructure_->addRhsFunction( rhsDummy );
            else
                fsi.problemStructureNonLin_->addRhsFunction( rhsDummy );
        
            // Geometrie-RW separat, falls geometrisch explizit.
            // Bei Geometrisch implizit: Keine RW in die factoryFSI fuer das
            // Geometrie-Teilproblem, da sonst (wg. dem ZeroDirichlet auf dem Interface,
            // was wir brauchen wegen Kopplung der Struktur) der Kopplungsblock C4
            // in derselben Zeile, der nur Werte auf dem Interface haelt, mit eliminiert.
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryGeometry( new BCBuilder<SC,LO,GO,NO>( ) );
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluidInterface;
            if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
                bcFactoryFluidInterface = Teuchos::rcp( new BCBuilder<SC,LO,GO,NO>( ) );

//                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "fix Randwerte für Richter Benchmark");
            bcFactoryGeometry->addBC(zeroDirichlet3D, 2, 0, domainGeometry, "Dirichlet", dim); // inflow ring
            bcFactoryGeometry->addBC(zeroDirichlet3D, 3, 0, domainGeometry, "Dirichlet", dim); // outflow ring
            bcFactoryGeometry->addBC(zeroDirichlet3D, 4, 0, domainGeometry, "Dirichlet", dim); // inflow
            bcFactoryGeometry->addBC(zeroDirichlet3D, 5, 0, domainGeometry, "Dirichlet", dim); // outflow
            // Die RW, welche nicht Null sind in der rechten Seite (nur Interface) setzen wir spaeter per Hand.
            // Hier erstmal Dirichlet Nullrand, wird spaeter von der Sturkturloesung vorgegeben
            bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // interface
            if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 6, 0, domainFluidVelocity, "Dirichlet", dim);

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
