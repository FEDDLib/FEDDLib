#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/FSCI.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*! Test case for specific artery geometrie or straight tube geometry. Inflow depends on inflow region
	-> artery: Inflow scaled with normal vector on inflow (x,y,z) * laplaceInflow	

*/



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
        res[2] = parameters[0] / 1. * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = parameters[0] / 1. * x[0];

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
        res[2] = -parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = -parameters[0] / parameters[2] * x[0];
    }

    return;
}

void parabolicInflow3DArtery(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = -parameters[0] / parameters[2] * x[0] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else
    {
        res[1] = 0.;
        res[0] = 0.;
        res[2] = -parameters[0] / parameters[2] * x[0];
    }

    return;
}

void parabolicInflow3DLinArtery(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // parameters[2] is the maxium solution value of the laplacian parabolic inflow problme
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    
    if(t < parameters[1])
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = -parameters[0] / parameters[2] * x[0] * t / parameters[1];
    }
    else
    {
        res[0] = 0.;
        res[1] = 0.;
        res[2] = -parameters[0] / parameters[2] * x[0];
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

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 0.0;	
    res[0] = m * x[0];

}
void inflowChem(double* x, double* res, double t, const double* parameters)
{
    res[0] = 1.;
    
    return;
}

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType NO;

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
    string xmlProblemFile = "parametersProblemFSCI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFileGE = "parametersPrecGE.xml"; // GE
    string xmlPrecFileGI = "parametersPrecGI.xml"; // GI
    myCLP.setOption("precfileGE",&xmlPrecFileGE,".xml file with Inputparameters.");
    myCLP.setOption("precfileGI",&xmlPrecFileGI,".xml file with Inputparameters.");
    string xmlSolverFileFSI = "parametersSolverFSCI.xml"; // GI
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
    
    string xmlPrecFileChem = "parametersPrecChem.xml";
    myCLP.setOption("precfileChem",&xmlPrecFileChem,".xml file with Inputparameters.");
    
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
        
        ParameterListPtr_Type parameterListPrecChem = Teuchos::getParametersFromXmlFile(xmlPrecFileChem);
        
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
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );

        ParameterListPtr_Type parameterListChemAll(new Teuchos::ParameterList(*parameterListPrecChem));
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Chem") );
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );

        
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
        if (verbose)
        {
            cout << "###############################################" <<endl;
            cout << "############ Starting FSCI  ... ################" <<endl;
            cout << "###############################################" <<endl;
        }

        DomainPtr_Type domainP1fluid;
        DomainPtr_Type domainP1struct;
        DomainPtr_Type domainP1chem;
        DomainPtr_Type domainP2fluid;
        DomainPtr_Type domainP2struct;
        DomainPtr_Type domainP2chem;
        
        DomainPtr_Type domainFluidVelocity;
        DomainPtr_Type domainFluidPressure;
        DomainPtr_Type domainChem;
        DomainPtr_Type domainStructure;
        DomainPtr_Type domainGeometry;
        
        std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","Compute Inflow");
        std::string geometryType = parameterListAll->sublist("Parameter").get("Geometry Type","Artery");
        
        
        TimeMonitor_Type totalTimeMonitor(*totalTime);
    
        TimeMonitor_Type buildMeshMonitor(*buildMesh);
        if (verbose)
        {
            cout << " -- Building Mesh ... " << flush;
        }

        domainP1fluid.reset( new Domain_Type( comm, dim ) );
        domainP1struct.reset( new Domain_Type( comm, dim ) );
        domainP1chem.reset(new Domain_Type(comm,dim));
        
        domainP2fluid.reset( new Domain_Type( comm, dim ) );
        domainP2struct.reset( new Domain_Type( comm, dim ) );
        domainP2chem.reset( new Domain_Type(comm,dim));
        
        //                    

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
        
        
        if (!discType.compare("P2"))
        {
            domainFluidVelocity = domainP2fluid;
            domainFluidPressure = domainP1fluid;
            domainChem = domainP2struct;
            domainStructure = domainP2struct;
            domainGeometry = domainP2fluid;
        }
        else
        {
            domainFluidVelocity = domainP1fluid;
            domainFluidPressure = domainP1fluid;
            domainStructure = domainP1struct;
            domainChem = domainP1struct;
            domainGeometry = domainP1fluid;
            //                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,"P1/P1 for FSI not implemented!");
        }

        /*Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domainFluidVelocity->getMapUnique()));
        vec_int_ptr_Type BCFlags = domainFluidVelocity->getBCFlagUnique();

        Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
        for(int i=0; i< entries.size(); i++){
            entries[i] = BCFlags->at(i);
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

        exPara->setup("FlagsFluid",domainFluidVelocity->getMesh(), discType);

        exPara->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domainFluidVelocity->getMapUnique(), domainFluidVelocity->getMapUniqueP2());

        exPara->save(0.0);

        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara2(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution2(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapUnique()));
        vec_int_ptr_Type BCFlags2 = domainStructure->getBCFlagUnique();

        Teuchos::ArrayRCP< SC > entries2  = exportSolution2->getDataNonConst(0);
        for(int i=0; i< entries2.size(); i++){
            entries2[i] = BCFlags2->at(i);
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst2 = exportSolution2;

        exPara2->setup("FlagsStructure", domainStructure->getMesh(), discType);

        exPara2->addVariable(exportSolutionConst2, "Flags", "Scalar", 1,domainStructure->getMapUnique(), domainStructure->getMapUniqueP2());

        exPara2->save(0.0);


        if (verbose){
            cout << "done! -- " << endl;
        }*/
            
                        
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
        

        
        /*if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
            
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
        }*/



            // #####################
        // Problem definieren
        // #####################
        Teuchos::RCP<SmallMatrix<int>> defTS;
        if(geometryExplicit)
        {
            // SmallMatrix<int> defTS(4);
            defTS.reset( new SmallMatrix<int> (5) );

            // Fluid
            (*defTS)[0][0] = 1;
            (*defTS)[0][1] = 1;

            // Struktur
            (*defTS)[2][2] = 1;
            
            //Chem
            (*defTS)[3][3] = 1;
        }
        else
        {
            // SmallMatrix<int> defTS(5);
            defTS.reset( new SmallMatrix<int> (6) );

            // Fluid
            (*defTS)[0][0] = 1;
            (*defTS)[0][1] = 1;
            // TODO: [0][4] und [1][4] bei GI + Newton noetig?
            /* if (verbose)
                std::cout << "### Double check temporal discretization of Shape Derivatives! ###" << std::endl;
            
            (*defTS)[0][5] = 1;
            (*defTS)[1][5] = 1;*/
            
            // Struktur
            (*defTS)[2][2] = 1;
            
            (*defTS)[3][3] = 1;
            
        }

        vec2D_dbl_Type diffusionTensor(dim,vec_dbl_Type(3));
        double D0 = parameterListAll->sublist("Parameter Diffusion").get("D0",1.);
        for(int i=0; i<dim; i++){
            diffusionTensor[0][0] =D0;
            diffusionTensor[1][1] =D0;
            diffusionTensor[2][2] =D0;

            if(i>0){
            diffusionTensor[i][i-1] = 0;
            diffusionTensor[i-1][i] = 0;
            }
            else
            diffusionTensor[i][i+1] = 0;				
        }
 
        FSCI<SC,LO,GO,NO> fsci(domainFluidVelocity, discType,
                                domainFluidPressure, "P1",
                                domainStructure, discType,
                                domainChem, discType,
                                domainInterface, discType,
                                domainGeometry, discType,
                                diffusionTensor, reactionFunc,
                                parameterListFluidAll, parameterListStructureAll, parameterListChemAll, parameterListAll,
                                parameterListGeometry, defTS);


        domainFluidVelocity->info();
        domainFluidPressure->info();
        domainStructure->info();
        domainGeometry->info();

        fsci.info();
                    
        std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter").get("Max Velocity",1.));
        parameter_vec.push_back( parameterListProblem->sublist("Parameter").get("Max Ramp Time",2.) );
        
        TEUCHOS_TEST_FOR_EXCEPTION(bcType != "Compute Inflow", std::logic_error, "Select a valid boundary condition. Only Compute Inflow available.");

        //#############################################
        //#############################################
        //#### Compute parabolic inflow with laplacian
        //#############################################
        //#############################################
        /*MultiVectorConstPtr_Type solutionLaplace;
        {
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryLaplace(new BCBuilder<SC,LO,GO,NO>( ));
            
                bcFactoryLaplace->addBC(zeroBC, 4, 0, domainFluidVelocity, "Dirichlet", 1); //inflow ring
                bcFactoryLaplace->addBC(zeroBC, 4, 0, domainFluidVelocity, "Dirichlet", 1); //outflow ring
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
                bcFactoryLaplace->addBC(zeroBC, 3, 0, domainFluidVelocity, "Dirichlet", 1);
            bcFactoryLaplace->addBC(zeroBC, 10, 0, domainFluidVelocity, "Dirichlet", 1);
            bcFactoryLaplace->setRHS( laplace.getSolution(), 0.);
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
        */
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

        // TODO: Vermutlich braucht man keine bcFactoryFluid und bcFactoryStructure,
        // da die RW sowieso auf dem FSI-Problem gesetzt werden.

        // Fluid-RW
        {
            bool zeroPressure = parameterListProblem->sublist("Parameter Fluid").get("Set Outflow Pressure to Zero",false);
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluid( new BCBuilder<SC,LO,GO,NO>( ) );
                            
            bcFactory->addBC(parabolicInflow3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);//, solutionLaplace); // inflow 
            //bcFactory->addBC(parabolicInflow3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow
            bcFactoryFluid->addBC(parabolicInflow3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);//, solutionLaplace); // inflow 
            //bcFactoryFluid->addBC(parabolicInflow3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow

            bcFactoryFluid->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall

            
            bcFactory->addBC(zeroDirichlet3D, 2, 0, domainFluidVelocity, "Dirichlet", dim); // inflow ring                
            bcFactoryFluid->addBC(zeroDirichlet3D, 2, 0, domainFluidVelocity, "Dirichlet", dim); // inflow ring

                bcFactory->addBC(zeroDirichlet3D, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec);//, solutionLaplace); // inflow ring                
                bcFactoryFluid->addBC(zeroDirichlet3D, 4, 0, domainFluidVelocity, "Dirichlet", dim); // inflow ring
            
            if (zeroPressure) {
                //bcFactory->addBC(zeroBC, 4, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                    bcFactory->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow
                
                //bcFactoryFluid->addBC(zeroBC, 4, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                    bcFactoryFluid->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow
            }
            
            // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
            // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
            fsci.problemFluid_->addBoundaries(bcFactoryFluid);
        }

        // Struktur-RW
        {
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );
                bcFactory->addBC(zeroDirichlet3D, 0, 2, domainStructure, "Dirichlet_Y_Z", dim); // inflow/outflow strip fixed in y direction
                bcFactory->addBC(zeroDirichlet3D, 1, 2, domainStructure, "Dirichlet_X_Y", dim); // inflow/outflow strip fixed in y direction
                bcFactory->addBC(zeroDirichlet3D, 2, 2, domainStructure, "Dirichlet_Z", dim); // inlet fixed in Z direction
                bcFactory->addBC(zeroDirichlet3D, 3, 2, domainStructure, "Dirichlet_X", dim); // outlet fixed in X direction

                bcFactoryStructure->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet_Y_Z", dim); 
                bcFactoryStructure->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_X_Y", dim); 
                bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Z", dim);           
                bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_X", dim); 
            // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
            // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
            if (!fsci.problemSCI_->problemStructure_.is_null())
                fsci.problemSCI_->problemStructure_->addBoundaries(bcFactoryStructure);
            else
                fsci.problemSCI_->problemStructureNonLin_->addBoundaries(bcFactoryStructure);
        }
        // RHS dummy for structure
                
        if (!fsci.problemSCI_->problemStructure_.is_null())
            fsci.problemSCI_->problemStructure_->addRhsFunction( rhsDummy );
        else
            fsci.problemSCI_->problemStructureNonLin_->addRhsFunction( rhsDummy );
    
        // Geometrie-RW separat, falls geometrisch explizit.
        // Bei Geometrisch implizit: Keine RW in die factoryFSI fuer das
        // Geometrie-Teilproblem, da sonst (wg. dem ZeroDirichlet auf dem Interface,
        // was wir brauchen wegen Kopplung der Struktur) der Kopplungsblock C4
        // in derselben Zeile, der nur Werte auf dem Interface haelt, mit eliminiert.
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryGeometry( new BCBuilder<SC,LO,GO,NO>( ) );
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluidInterface;
        if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
            bcFactoryFluidInterface = Teuchos::rcp( new BCBuilder<SC,LO,GO,NO>( ) );

            bcFactoryGeometry->addBC(zeroDirichlet3D, 0, 0, domainGeometry, "Dirichlet", dim); // inflow/outflow strip fixed in y direction
        bcFactoryGeometry->addBC(zeroDirichlet3D, 1, 0, domainGeometry, "Dirichlet", dim); // inflow/outflow strip fixed in y direction
        bcFactoryGeometry->addBC(zeroDirichlet3D, 2, 0, domainGeometry, "Dirichlet", dim); // inlet fixed in Z direction
        bcFactoryGeometry->addBC(zeroDirichlet3D, 3, 0, domainGeometry, "Dirichlet", dim); // inlet fixed in X direction
        bcFactoryGeometry->addBC(zeroDirichlet3D, 4, 0, domainGeometry, "Dirichlet", dim); // inlet/outlet Ring
        bcFactoryGeometry->addBC(zeroDirichlet3D, 5, 0, domainGeometry, "Dirichlet", dim); // ?
            bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // Interface
        
        // Die RW, welche nicht Null sind in der rechten Seite (nur Interface) setzen wir spaeter per Hand.
        // Hier erstmal Dirichlet Nullrand, wird spaeter von der Sturkturloesung vorgegeben
            bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // interface
        if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 6, 0, domainFluidVelocity, "Dirichlet", dim);

        fsci.problemGeometry_->addBoundaries(bcFactoryGeometry);
        if ( preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko")
            fsci.getPreconditioner()->setFaCSIBCFactory( bcFactoryFluidInterface );

        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryChem( new BCBuilder<SC,LO,GO,NO>( ) ); 
        if (dim==2)
        {
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
                            
        }
        else if(dim==3)
        {

            bcFactory->addBC(inflowChem, 6, 3, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactoryChem->addBC(inflowChem, 6, 0, domainChem, "Dirichlet", 1); // inflow of Chem
            // bcFactoryChem->addBC(inflowChem, 4, 1, domainChem, "Dirichlet", 1); // inflow of Chem
        }

    // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
    // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()

        fsci.problemSCI_->problemChem_->addBoundaries(bcFactoryChem);

        // #####################
        // Zeitintegration
        // #####################
        fsci.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen
        fsci.problemSCI_->addBoundaries(bcFactory);
        fsci.initializeProblem();
        
        fsci.initializeGE();
        // Matrizen assemblieren

        fsci.assemble();
    
        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

        // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
        // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
        daeTimeSolver.defineTimeStepping(*defTS);

        // Uebergebe das (nicht) lineare Problem
        daeTimeSolver.setProblem(fsci);

        // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
        // defTS definiert worden sind.
        daeTimeSolver.setupTimeStepping();

        daeTimeSolver.advanceInTime();
    }
    

    TimeMonitor_Type::report(std::cout);

    return(EXIT_SUCCESS);
}
