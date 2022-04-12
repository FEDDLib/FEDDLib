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




void dummyFunc(double* x, double* res, double t, const double* parameters)
{
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
    string xmlProblemFile = "parametersProblemFSI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");       
    string xmlSolverFileFSI = "parametersSolverFSI.xml"; // GI
    myCLP.setOption("solverfileFSI",&xmlSolverFileFSI,".xml file with Inputparameters.");
    

    
    string xmlProblemFileFluid = "parametersProblemChem.xml";
    myCLP.setOption("problemFileFluid",&xmlProblemFileFluid,".xml file with Inputparameters.");
    string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    
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

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        ParameterListPtr_Type parameterListPrecChem = Teuchos::getParametersFromXmlFile(xmlPrecFileChem);
        


        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
      
        
        parameterListAll->setParameters(*parameterListSolverFSI);

        
        ParameterListPtr_Type parameterListFluidAll(new Teuchos::ParameterList(*parameterListPrecFluidMono)) ;
        sublist(parameterListFluidAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Fluid") );
        parameterListFluidAll->setParameters(*parameterListPrecFluidTeko);

        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );

        parameterListStructureAll->setParameters(*parameterListPrecStructure);
                 
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
                cout << "############ Starting SCI  ... ################" <<endl;
                cout << "###############################################" <<endl;
            }

            DomainPtr_Type domainP1chem;
            DomainPtr_Type domainP1struct;
            DomainPtr_Type domainP2chem;
            DomainPtr_Type domainP2struct;
          
          
            DomainPtr_Type domainChem;
            DomainPtr_Type domainStructure;
            
            std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","parabolic");
            
     
            domainP1chem.reset( new Domain_Type( comm, dim ) );
            domainP1struct.reset( new Domain_Type( comm, dim ) );
            domainP2chem.reset( new Domain_Type( comm, dim ) );
            domainP2struct.reset( new Domain_Type( comm, dim ) );
                                    
			MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(2);
			domainP1Array[0] = domainP1chem;
			domainP1Array[1] = domainP1struct;
		

            pListPartitioner->set("Build Edge List",true);
            pListPartitioner->set("Build Surface List",true);
		                    
            MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
            
            partitionerP1.readAndPartition();
                        
            if (!discType.compare("P2")){
                domainP2chem->buildP2ofP1Domain( domainP1chem );
                domainP2struct->buildP2ofP1Domain( domainP1struct );
                
                domainChem = domainP2chem;
            	domainStructure = domainP2Structure;   
            }
         
           


            if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
                
                if (verbose)
                    std::cout << "\t### Exporting fluid and solid subdomains ###\n";

                typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
                typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
                typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
                typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
                typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
				// Same subdomain for solid and chemistry, as they have same domain
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
    
            
            domainChem->setReferenceConfiguration();


            SCI<SC,LO,GO,NO> sci(domainFluidVelocity, discType,
                                 domainFluidPressure, "P1",
                                 domainStructure, discType,
                                 domainInterface, discType,
                                 parameterListFluidAll,
                                 parameterListStructureAll,
                                 parameterListAll,
                                 defTS);
            
            sci.info();
            
            
            Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );


        
                
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
                
            }
            else if(dim==3)
            {
                bcFactory->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                
                    bcFactory->addBC(inflow3DRichterSuperFast, 2, 0, domainFluidVelocity, "Dirichlet", dim, parameter_vec); // inflow
                
                bcFactoryFluid->addBC(zeroDirichlet3D, 1, 0, domainFluidVelocity, "Dirichlet", dim); // wall
                
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
