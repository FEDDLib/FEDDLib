#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/General/HDF5Import.hpp"
#include "feddlib/core/General/HDF5Export.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include "feddlib/problems/specific/LinElas.hpp"
#include "feddlib/problems/specific/NonLinElasticity.hpp"
#include "feddlib/problems/specific/NavierStokes.hpp"
#include "feddlib/problems/specific/Geometry.hpp"
#include "feddlib/problems/specific/FSI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include "feddlib/problems/Solver/Preconditioner.hpp"


/*! Test case for the curved artery gemeotry proposed in 
    [1] Numerical modeling of fluid–structure interaction in arteries with
      anisotropic polyconvex hyperelastic and anisotropic viscoelastic
      material models at finite strains
      Daniel Balzani, Simone Deparis, Simon Fausten, Davide Forti, Alexander Heinlein,
      Axel Klawonn, Alfio Quarteroni, Oliver Rheinbach,and Joerg Schröder
    [2] Comparison of arterial wall models in fluid–structure interaction simulations
      D. Balzani, A. Heinlein, A. Klawonn,O. Rheinbach, J. Schröder (2023)
    - The Parameters for fluid and structure are currently set to the NH3/NH4 model in [2]
    - The parabolic inflow profile was computed beforehand (e.g. in laplace example) and safed
      to be reloaded for specific mesh resolution and discretization
    - The uses meshes are #2, #3, #4 from table 4 in [1]. The other resolutions were to large for P2 disc.
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

// void flowrate3D(double* x, double* res, double t, const double* parameters)
// {
//     // parameters[0] is the maxium desired velocity
//     // parameters[1] rampTime
//     // parameters[2] flowrate

//     // The center point of the inlet is (0,0,0)   

//     // Distance from center
//     double Q = 0.;
//     if(t < parameters[1])
//     {
       
//         Q = parameters[2] * 0.5*( 1. - cos( M_PI*t/parameters[1] ));
//     }
//     else
//     {
//         Q = parameters[2];
//     }

//     res[0] = Q;
//     return;
// }


void flowrate3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] rampTime
    // parameters[2] flowrate
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing
    double heartBeatStart = parameters[3];

    if(t < parameters[1])
    {
        res[0] = parameters[2] * 0.5 * ( ( 1 - cos( M_PI*t/parameters[1]) ));
    }
    else if(t > heartBeatStart)
    {
    
        double a0    = 11.693284502463376;
        double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
            0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
            0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
        double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
            -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
            0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
                    
        double Q = 0.5*a0;
        

        double t_min = t - fmod(t,1.0)+heartBeatStart-std::floor(t); ; //FlowConditions::t_start_unsteady;
        double t_max = t_min + 1.0; // One heartbeat lasts 1.0 second    
        double y = M_PI * ( 2.0*( t-t_min ) / ( t_max - t_min ) -1.0)  ;
        
        for(int i=0; i< 20; i++)
            Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
        
        
        // Remove initial offset due to FFT
        Q -= 0.026039341343493;
        Q = (Q - 2.85489)/(7.96908-2.85489);

        res[0] =  parameters[2] + parameters[2]* Q *1.6563 - 0.1 ;
        
    }
    else
    {
        res[0] = parameters[2] ;

    }

    return;
}

void flowrate3DLinear(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] rampTime
    // parameters[3] flowrate

    // The center point of the inlet is (0,0,0)   

    // Distance from center
    double Q = 0.;
    if(t < parameters[1])
    {
       
        Q = parameters[2] *  t / parameters[1];
    }
    else
    {
        Q = parameters[2];
    }

    res[0] = Q;
    return;
}

void parabolicInflow3D(double* x, double* res, double t, const double* parameters)
{
    // parameters[0] is the maxium desired velocity
    // parameters[1] end of ramp
    // we use x[0] for the laplace solution in the considered point. Therefore, point coordinates are missing

    res[0] = 0.;
    res[1] = 0.;
    res[2] = -parameters[0]  * x[0];


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

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
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
     string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileGeometry = "parametersPrecGeometry.xml";
    myCLP.setOption("precfileGeometry",&xmlPrecFileGeometry,".xml file with Inputparameters.");
      
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
        sublist(parameterListGeometry, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Geometry") );
    // we only compute the preconditioner for the geometry problem once
        sublist( parameterListGeometry, "General" )->set( "Preconditioner Method", "MonolithicConstPrec" );
                   
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
                    
            DomainPtr_Type domainFluidVelocity;
            DomainPtr_Type domainFluidPressure;
            DomainPtr_Type domainStructure;
            DomainPtr_Type domainGeometry;
                        
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
                        idsInterface.push_back(4);
                        idsInterface.push_back(5);
               
                        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
                        domainP1Array[0] = domainP1fluid;
                        domainP1Array[1] = domainP1struct;
                        
                        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
                        if (!discType.compare("P2")){
                            pListPartitioner->set("Build Edge List",true);
                            pListPartitioner->set("Build Surface List",true);
                        }
                        else{
                            pListPartitioner->set("Build Edge List",true);
                            pListPartitioner->set("Build Surface List",true);
                        }
                        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
                        
                        partitionerP1.readAndPartition();
                        
                        if (!discType.compare("P2")){
                            domainP2fluid->buildP2ofP1Domain( domainP1fluid );
                            domainP2struct->buildP2ofP1Domain( domainP1struct );
                        }
                        
                        
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
                        // domainFluidVelocity->preProcessMesh(true,false);
                        domainFluidVelocity->exportNodeFlags("Fluid");
                        domainStructure->exportNodeFlags("Solid");
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
         
            if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
                
                if (verbose)
                    std::cout << "\t### Exporting fluid and solid subdomains ###\n";

               domainFluidVelocity->exportDistribution("Fluid");
               domainStructure->exportDistribution("Solid");

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
                     
            std::vector<double> parameter_vec(1, parameterListProblem->sublist("Parameter Fluid").get("Max Velocity",1.));
            parameter_vec.push_back( parameterListProblem->sublist("Parameter Fluid").get("Max Ramp Time",0.1) );
            parameter_vec.push_back(parameterListProblem->sublist("Parameter Fluid").get("Flowrate",3.0));
            parameter_vec.push_back(parameterListProblem->sublist("Parameter Fluid").get("Heart Beat Start",0.2));

        
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

            // TODO: Vermutlich braucht man keine bcFactoryFluid und bcFactoryStructure,
            // da die RW sowieso auf dem FSI-Problem gesetzt werden.

            // Fluid-RW
            {

                MultiVectorConstPtr_Type solutionLaplace;
                string meshNumber = parameterListProblem->sublist("Mesh Partitioner").get("Mesh Number","2");

                HDF5Import<SC,LO,GO,NO> importer(domainFluidVelocity->getMapUnique() ,"laplace_parabolic_fluidBenchmark"+ meshNumber+"_"+discType);
                Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > solutionImported = importer.readVariablesHDF5("solution");
                solutionLaplace = solutionImported; // This must me normalized to 1!!

                bool zeroPressure = parameterListProblem->sublist("Parameter Fluid").get("Set Outflow Pressure to Zero",false);
                Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryFluid( new BCBuilder<SC,LO,GO,NO>( ) );
                               
                //bcFactory->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                string rampType = parameterListProblem->sublist("Parameter Fluid").get("Ramp type","cos");
                if (rampType == "cos") {
                
                    bcFactory->addBC(parabolicInflow3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace,true, flowrate3D); // inflow ring
                    //bcFactory->addBC(parabolicInflow3DArtery, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow
                    bcFactoryFluid->addBC(parabolicInflow3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace,true, flowrate3D); // inflow ring
                    //bcFactoryFluid->addBC(parabolicInflow3DArtery, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow                  
                }
                else if(rampType == "linear"){
                
                    bcFactory->addBC(parabolicInflow3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace,true, flowrate3DLinear); // inflow 
                    //bcFactory->addBC(parabolicInflow3DLinArtery, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow ring
                    bcFactoryFluid->addBC(parabolicInflow3D, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace,true, flowrate3DLinear); // inflow
                    //bcFactoryFluid->addBC(parabolicInflow3DLinArtery, 4, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec, solutionLaplace); // inflow ring               
                }
                
                // bcFactory->addBC(zeroDirichlet3D, 4, 0, domainFluidVelocity, "Dirichlet_Z", dim); // inflow ring                
                // bcFactoryFluid->addBC(zeroDirichlet3D, 4, 0, domainFluidVelocity, "Dirichlet_Z", dim); // inflow ring
                
                // bcFactory->addBC(zeroDirichlet3D, 5, 0, domainFluidVelocity, "Dirichlet_X", dim); // outflow ring                
                // bcFactoryFluid->addBC(zeroDirichlet3D, 5, 0, domainFluidVelocity, "Dirichlet_X", dim); // outflow ring

                if (zeroPressure) {
                    //bcFactory->addBC(zeroBC, 4, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                    bcFactory->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow
                    
                    //bcFactoryFluid->addBC(zeroBC, 4, 1, domainFluidPressure, "Dirichlet", 1); // outflow ring
                    bcFactoryFluid->addBC(zeroBC, 3, 1, domainFluidPressure, "Dirichlet", 1); // outflow
                }
                
                // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
                // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
                fsi.problemFluid_->addBoundaries(bcFactoryFluid);
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

            // bcFactoryGeometry->addBC(zeroDirichlet3D, 0, 0, domainGeometry, "Dirichlet", dim); // inflow/outflow strip fixed in y direction
            // bcFactoryGeometry->addBC(zeroDirichlet3D, 1, 0, domainGeometry, "Dirichlet", dim); // inflow/outflow strip fixed in y direction
            // bcFactoryGeometry->addBC(zeroDirichlet3D, 2, 0, domainGeometry, "Dirichlet", dim); // inlet fixed in Z direction
            // bcFactoryGeometry->addBC(zeroDirichlet3D, 3, 0, domainGeometry, "Dirichlet", dim); // inlet fixed in X direction
            bcFactoryGeometry->addBC(zeroDirichlet3D, 4, 0, domainGeometry, "Dirichlet", dim); // inlet Ring
            bcFactoryGeometry->addBC(zeroDirichlet3D, 5, 0, domainGeometry, "Dirichlet", dim); // outlet Ring
            bcFactoryGeometry->addBC(zeroDirichlet3D, 6, 0, domainGeometry, "Dirichlet", dim); // Interface
          
            // Die RW, welche nicht Null sind in der rechten Seite (nur Interface) setzen wir spaeter per Hand.
            // Hier erstmal Dirichlet Nullrand, wird spaeter von der Sturkturloesung vorgegeben
            if (preconditionerMethod == "FaCSI" || preconditionerMethod == "FaCSI-Teko"){
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 4, 0, domainFluidVelocity, "Dirichlet", dim);
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 5, 0, domainFluidVelocity, "Dirichlet", dim);
                bcFactoryFluidInterface->addBC(zeroDirichlet3D, 6, 0, domainFluidVelocity, "Dirichlet", dim);
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
    return EXIT_SUCCESS;
}
