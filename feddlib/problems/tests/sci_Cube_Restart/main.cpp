#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/problems/specific/SCI.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"
#include "feddlib/problems/Solver/NonLinearSolver.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include "feddlib/core/General/AceGenInterfaceCheck.hpp"
#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_StackedTimer.hpp>

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


void zeroDirichlet(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;

    return;
}

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 0.0;	
    res[0] = m * x[0];

}

void zeroDirichlet3D(double* x, double* res, double t, const double* parameters)
{
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;

    return;
}

void inflowChem(double* x, double* res, double t, const double* parameters)
{
	if(t>=parameters[0])
    	res[0] = 1.;
    else	
    	res[0] = 0.;
    return;
}


void rhsX(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    
    res[0] = parameters[1];
    res[1] = 0.;
    res[2] = 0.;
    return;
}

void rhsY(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] =  parameters[1];
    res[2] = 0.;
    return;
}

void rhsZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    res[1] = 0.;
    res[2] = parameters[1];
    return;
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStep (lambda)
// 3 : LoadStep end time
// 4 : Flag 

void rhsYZ(double* x, double* res, double* parameters){

    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];

  	res[0] =0.;
    res[1] =0.;
    res[2] =0.;

    if(parameters[0]+1.e-12 < TRamp)
        force = (parameters[0]+loadStepSize) * parameters[1] / TRamp ;
    else
        force = parameters[1];


    if(parameters[6] == 4  || parameters[6] == 5){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    
    return;
}


// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStep (lambda)
// 3 : LoadStep end time
// 4 : Flag 

void rhsArtery(double* x, double* res, double* parameters){

    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];

  	res[0] =0.;
    res[1] =0.;
    res[2] =0.;
    
    if(parameters[0]+1.e-12 < TRamp)
        force = (parameters[0]+loadStepSize) * parameters[1] / TRamp ;
    else
        force = parameters[1];


    if(parameters[5] == 5){
      	res[0] = force;
        res[1] = force;
        res[2] = force;
    }
    
    return;
}

// Parameter Structure
// 0 : time
// 1 : force
// 2 : loadStepSize
// 3 : LoadStep end time
// 4 : Flag 

void rhsHeartBeatCube(double* x, double* res, double* parameters){

    res[0] =0.;
    res[1] =0.;
    res[2] = 0.;
    double lambda=0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double heartBeatStart = parameters[4];
    
	double a0    = 11.693284502463376;
	double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
	double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		         
    double Q = 0.5*a0;
    

    double t_min = parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
    double t_max = t_min + 0.52; // One heartbeat lasts 1.0 second    
    double y = M_PI * ( 2.0*( parameters[0]-t_min ) / ( t_max - t_min ) -0.87  );
    
    for(int i=0; i< 20; i++)
        Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
    
    
    // Remove initial offset due to FFT
    Q -= 0.026039341343493;
    Q = (Q - 2.85489)/(7.96908-2.85489);
    
    bool Qtrue=false;
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0]+1.0e-10 < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))+1.e-10< 0.5)
    	lambda= 0.75;
    else{
        lambda = 0.75+0.25*Q;//*0.005329; // 0.775+0.125 * cos(4*M_PI*(parameters[0]));
        Qtrue = true; 
    } 
  
    double forceDirection = force/fabs(force);
    if(parameters[6]==5 || parameters[6]==4){
        res[0] =lambda*force;//+forceDirection*Q;
        res[1] =lambda*force;//+forceDirection*Q;
        res[2] =lambda*force;//+forceDirection*Q;        
    } 
    
   /* if(parameters[0]< heartBeatStart){
    	Q = 0.;
    }
    
    if(parameters[0]+1e-12 < TRamp)
        force = force * (parameters[0]+loadStepSize);
    
    if(parameters[6] == 5 || parameters[6] == 4){
     	res[0] = force+Q*0.005329;
        res[1] = force+Q*0.005329;
       	res[2] = force+Q*0.005329;
       	      	
    }*/
      
}


void rhsCubePaper(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[2] = 0.;
    double force = parameters[1];
    double loadStepSize = parameters[2];
    double TRamp = parameters[3];
    double lambda=0.;
    
    double heartBeatStart = parameters[4];
    
    if(parameters[0]+1e-12 < TRamp)
        lambda = 0.875*(parameters[0]+loadStepSize)/ TRamp;
    else if(parameters[0] <= TRamp+1.e-12)
    	lambda = 0.875;
    else if (parameters[0] < heartBeatStart)
    	lambda = 0.875;
    else if( parameters[0] < heartBeatStart + 0.5)
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= heartBeatStart + 0.5 && (parameters[0] - std::floor(parameters[0]))< 0.5)
    	lambda= 0.75;
    else
        lambda = 0.875 - 0.125 * cos(4*M_PI*(parameters[0]));
     
     
    if(parameters[6] == 5 || parameters[6] == 4){
        res[0] =lambda*force;
        res[1] =lambda*force;
        res[2] =lambda*force; 
    }
    
    

            
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
    if (!FEDD::aceGenInterfaceAvailable())
        return EXIT_FAILURE;

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

       // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    std::string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
   
    std::string xmlProblemFile = "parametersProblemSCI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");    
    
    std::string xmlProblemStructureFile = "parametersProblemStructure.xml";  
    myCLP.setOption("problemfileStructure",&xmlProblemStructureFile,".xml file with Inputparameters.");    
 
    std::string xmlSolverFileSCI = "parametersSolverSCI.xml"; 
    myCLP.setOption("solverfileSCI",&xmlSolverFileSCI,".xml file with Inputparameters.");
    
    std::string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    
    std::string xmlPrecFileChem = "parametersPrecChem.xml";
    myCLP.setOption("precfileChem",&xmlPrecFileChem,".xml file with Inputparameters.");
    
    std::string xmlPrecCEFile = "parametersPrecCE.xml";
    myCLP.setOption("precCEfile",&xmlPrecCEFile,".xml file with Inputparameters.");

   
    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    bool testPassed = true;
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    {
        mpiSession.~GlobalMPISession();
        return 0;
    }
	Teuchos::RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Structure-chemical interaction",true));
    bool verbose (comm->getRank() == 0);
    TimeMonitor::setStackedTimer(stackedTimer);
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
       
        ParameterListPtr_Type parameterListProblemStructure = Teuchos::getParametersFromXmlFile(xmlProblemStructureFile);
        
        ParameterListPtr_Type parameterListSolverSCI = Teuchos::getParametersFromXmlFile(xmlSolverFileSCI);

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        ParameterListPtr_Type parameterListPrecChem = Teuchos::getParametersFromXmlFile(xmlPrecFileChem);
  


 		int 		dim				= parameterListProblem->sublist("Parameter").get("Dimension",2);
        std::string		meshType    	= parameterListProblem->sublist("Parameter").get("Mesh Type","unstructured");
        int 		m				= parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string      discType        = parameterListProblem->sublist("Parameter").get("Discretization","P2");
        std::string precMethod = parameterListProblem->sublist("General").get("Preconditioner Method","Monolithic");
        int         n;

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;     
       
        ParameterListPtr_Type parameterListPrec;
        bool chemistryExplicit_ =    parameterListAll->sublist("Parameter").get("Chemistry Explicit",false);

        if(chemistryExplicit_)
            parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecCEFile);
        else
            parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
       

        parameterListAll->setParameters(*parameterListSolverSCI);
        parameterListAll->setParameters(*parameterListPrec);
                    
		parameterListAll->setParameters(*parameterListProblemStructure);
        
        ParameterListPtr_Type parameterListChemAll(new Teuchos::ParameterList(*parameterListPrecChem)) ;
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Chem") );
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );
        parameterListChemAll->setParameters(*parameterListSolverSCI);
        parameterListChemAll->setParameters(*parameterListPrecChem);
        parameterListChemAll->setParameters(*parameterListProblemStructure);


        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrec));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        parameterListStructureAll->setParameters(*parameterListPrec);
        parameterListStructureAll->setParameters(*parameterListProblem);
        parameterListStructureAll->setParameters(*parameterListProblemStructure);
		
        TimePtr_Type totalTime(TimeMonitor_Type::getNewCounter("FEDD - main - Total Time"));
        TimePtr_Type buildMesh(TimeMonitor_Type::getNewCounter("FEDD - main - Build Mesh"));

        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);

        int size = comm->getSize() - numProcsCoarseSolve;

        // #####################
        // Mesh bauen und wahlen
        // #####################
    
        if (verbose)
        {
            std::cout << "###############################################" <<std::endl;
            std::cout << "############ Starting SCI  ... ################" <<std::endl;
            std::cout << "###############################################" <<std::endl;
        }

        DomainPtr_Type domainP1chem;
        DomainPtr_Type domainP1struct;
        DomainPtr_Type domainP2chem;
        DomainPtr_Type domainP2struct;
        
        
        DomainPtr_Type domainChem;
        DomainPtr_Type domainStructure;
                
        std::string rhsType = parameterListAll->sublist("Parameter").get("RHS Type","Constant");
    
        int minNumberSubdomains=1;
        if (verbose)
        {
            std::cout << "\t  ###############################################" <<std::endl;
            std::cout << "\t  ############ sci_Test.main:: INFOs ############" <<std::endl;
        }
        
    
        TEUCHOS_TEST_FOR_EXCEPTION( size%minNumberSubdomains != 0 , std::logic_error, "Wrong number of processors for structured mesh.");
        n = (int)(std::pow( size/minNumberSubdomains, 1/3.) + 100*Teuchos::ScalarTraits<double>::eps()); // 1/H
        // Ranks beyond the n^3 subdomains would build copies of subdomain (rank % n^3)
        TEUCHOS_TEST_FOR_EXCEPTION( n*n*n*minNumberSubdomains != size, std::logic_error, "The structured cube needs n^3 ranks (1, 8, 27, ...) apart from the coarse ranks, got " << size << ".");
        std::vector<double> x(3);
        x[0]=0.0;    x[1]=0.0;	x[2]=0.0;
        domainStructure.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
        domainChem.reset(new Domain<SC,LO,GO,NO>( x, 1., 1., 1., comm));
    
        // "Square5Element" and not "Square": the 6-element subcube decomposition
        // gives tetrahedra with a negative Jacobian determinant, which the AceGen
        // SCI elements cannot assemble on (Newton diverges at the first step).
        domainStructure->buildMesh( 3,"Square5Element", dim, discType, n, m, numProcsCoarseSolve);
        domainChem->buildMesh( 3,"Square5Element", dim, discType, n, m, numProcsCoarseSolve);
		
        domainStructure->setDofs(dim);
        domainChem->setDofs(1);
     
        domainChem->setReferenceConfiguration();
        domainStructure->setReferenceConfiguration();

        // Defining diffusion tensor. Only used for explicit approach
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
 
        
        Teuchos::RCP<SmallMatrix<int>> defTS;
        if(chemistryExplicit_){ 
            defTS.reset( new SmallMatrix<int> (1) );// In case of explicit approach the SCI system is only the structure block
            // Stucture
            (*defTS)[0][0] = 1;
        }
        else {
            defTS.reset( new SmallMatrix<int> (2) );
            // Stucture
            (*defTS)[0][0] = 1;
            // Chem
            (*defTS)[1][1] = 1;
        }

        SCI<SC,LO,GO,NO> sci(domainStructure, discType,
                                domainChem, discType, diffusionTensor, reactionFunc,
                                parameterListStructureAll,
                                parameterListChemAll,
                                parameterListAll,
                                defTS);
        
        sci.info();
        
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) ); 
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryChem( new BCBuilder<SC,LO,GO,NO>( ) ); 
        
    
        // Boundary Conditions for structure
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );

        TEUCHOS_TEST_FOR_EXCEPTION( dim<3, std::logic_error, "Only 3D Test available");                               

        

        if(verbose){
            std::cout << " \t Boundary Condition type Cube " << std::endl;
            std::cout << " \t \t Point (0,0,0) held in all directions " << std::endl;
            std::cout << " \t \t Connecting x,y,z planes are held in their respective planes" << std::endl;
            std::cout << " \t \t Connecting x,y,z edges are held in x,y-direction, y,z-direction, and x,z-direction" << std::endl;
            std::cout << " \t \t Pull on z=1 and y=1 plane. " << std::endl;
        }

        bcFactory->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_X", dim); // x=0
        bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Y", dim); // y=0
        bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim); // z=0
        
        bcFactory->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet", dim);
        bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X_Y", dim); //x,y = 0
        bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Y_Z", dim); // y,z= 0
        bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_X_Z", dim); // x,z = 0
        
        bcFactoryStructure->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_X", dim);
        bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_Y", dim);
        bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
        
        bcFactoryStructure->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet", dim);
        bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X_Y", dim);
        bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_Y_Z", dim);
        bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_X_Z", dim);

            bcFactoryStructure->addBC(zeroDirichlet3D, 99, 0, domainStructure, "Dirichlet", dim);

         sci.problemStructureNonLin_->addBoundaries(bcFactoryStructure);
        
                        
        if(verbose)
            std::cout << " \t Possible rhs functions for 'Cube' boundary type are: Constant, Paper, Heart Beat" << std::endl;
        if(rhsType=="Constant")
            sci.problemStructureNonLin_->addRhsFunction( rhsYZ,0 );
        else if(rhsType=="Paper")
            sci.problemStructureNonLin_->addRhsFunction( rhsCubePaper,0 );
        else if(rhsType=="Heart Beat")
            sci.problemStructureNonLin_->addRhsFunction( rhsHeartBeatCube,0 );
        else
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "You have not selected a valid RHS Type for BC Type cube.");                               

       
        double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
        
        sci.problemStructureNonLin_->addParemeterRhs( force );
        double loadStep = parameterListAll->sublist("Parameter").get("Load Step Size",1.);
        double loadRampEnd= parameterListAll->sublist("Parameter").get("Load Ramp End",1.);
        sci.problemStructureNonLin_->addParemeterRhs( loadStep );
        sci.problemStructureNonLin_->addParemeterRhs( loadRampEnd );
        double heartBeatStart= parameterListAll->sublist("Parameter").get("Heart Beat Start",70.);
        sci.problemStructureNonLin_->addParemeterRhs( heartBeatStart );
        sci.problemStructureNonLin_->addParemeterRhs( 0. ); // degree of the load function in space: develop's surface integral reads it from the last parameter
        if(verbose){
            std::cout << "\t The following parameters were set for the RHS of your problem:" << std::endl;
            std::cout << "\t \t Force applied to interior wall (Volume Force per List): " << force << std::endl;
            std::cout << "\t \t End of Loadstepping ramp: " << loadRampEnd << std::endl;
            std::cout << "\t \t Loadstep size: " << loadStep << std::endl;
            std::cout << "\t \t Heart beat start time: " << heartBeatStart << std::endl;          
        }  
        
                   
    

        std::vector<double> parameter_vec(1, parameterListAll->sublist("Parameter").get("Inflow Start Time",0.));
        bcFactory->addBC(inflowChem, 0, 1, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
        bcFactory->addBC(inflowChem, 1, 1, domainChem, "Dirichlet", 1, parameter_vec); // inflow of Chem
        bcFactory->addBC(inflowChem, 7, 1, domainChem, "Dirichlet", 1,parameter_vec);            		
        bcFactory->addBC(inflowChem, 9, 1, domainChem, "Dirichlet", 1,parameter_vec);
        
        
        bcFactoryChem->addBC(inflowChem, 0, 0, domainChem, "Dirichlet", 1,parameter_vec); // inflow of Chem
        bcFactoryChem->addBC(inflowChem, 1, 0, domainChem, "Dirichlet", 1, parameter_vec); // inflow of Chem
        bcFactoryChem->addBC(inflowChem, 7, 0, domainChem, "Dirichlet", 1,parameter_vec);            		
        bcFactoryChem->addBC(inflowChem, 9, 0, domainChem, "Dirichlet", 1,parameter_vec);
         
       
        if (verbose)
        {
            std::cout << "\t ... setup done." <<std::endl;
            std::cout << "\t  ###############################################" <<std::endl;
            std::cout << "\t  start solving problem ..." <<std::endl;
        }
        // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        sci.problemChem_->addBoundaries(bcFactoryChem);
        
          
        // #####################
        // Zeitintegration
        // #####################
        sci.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        sci.initializeProblem();

        sci.initializeCE();
        // Matrizen assemblieren
        sci.assemble();
                    

                    
        DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

        // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
        // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
        daeTimeSolver.defineTimeStepping(*defTS);

        // Uebergebe das (nicht) lineare Problem
        daeTimeSolver.setProblem(sci);

        // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massmatrizen auf den Zeilen, welche in
        // defTS definiert worden sind.
        daeTimeSolver.setupTimeStepping();

        daeTimeSolver.advanceInTime();

        MultiVectorConstPtr_Type solutionChem ;

        if(chemistryExplicit_ )
            solutionChem = sci.getChemProblem()->getSolution()->getBlock(0);
        else
            solutionChem = sci.getSolution()->getBlock(1);

        // Restart test: a restarted run (phase 2) must reproduce the uninterrupted run
        // (phase 1): its solution at the final time is compared with the checkpoint
        // phase 1 wrote at that time.
        if(parameterListAll->sublist("Timestepping Parameter").get("Restart", false)){
        double finalTime = parameterListAll->sublist("Timestepping Parameter").get("Final time compare", 0.0);
        double tolerance = parameterListAll->sublist("Timestepping Parameter").get("Restart tolerance", 1.e-8);
        HDF5Import<SC,LO,GO,NO> importer(sci.getSolution()->getBlock(0)->getMap(),restartFile(parameterListAll, "Solutiond_s"));
        HDF5Import<SC,LO,GO,NO> importerC(solutionChem->getMap(),restartFile(parameterListAll, "Solutionc"));

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > solutionImported = importer.readVariablesHDF5(std::to_string(finalTime));
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > solutionImportedC = importerC.readVariablesHDF5(std::to_string(finalTime));

        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaResults(new ExporterParaView<SC,LO,GO,NO>());
        exParaResults->setup("Restart_Error_D", domainStructure->getMesh(), domainStructure->getFEType());

        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaResultsC(new ExporterParaView<SC,LO,GO,NO>());
        exParaResultsC->setup("Restart_Error_C", domainChem->getMesh(), domainChem->getFEType());

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionD = sci.getSolution()->getBlock(0);
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionC = solutionChem;

        // Adding solution to paraview exporter
        exParaResults->addVariable(exportSolutionD, "d", "Vector", 3, domainStructure->getMapUnique());
        exParaResults->addVariable(solutionImported, "d_import", "Vector", 3,  domainStructure->getMapUnique());

        exParaResultsC->addVariable(exportSolutionC, "c", "Scalar", 1, domainChem->getMapUnique());
        exParaResultsC->addVariable(solutionImportedC, "c_import", "Scalar", 1,  domainChem->getMapUnique());

        // Calculating the error per node
        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > errorValues = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>( sci.getSolution()->getBlock(0)->getMap() ) ); 
        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > errorValuesC = Teuchos::rcp(new MultiVector<SC,LO,GO,NO>( solutionChem->getMap() ) ); 

        //this = alpha*A + beta*B + gamma*this
        errorValues->update( 1., exportSolutionD, -1. ,solutionImported, 0.);
        errorValuesC->update( 1., exportSolutionC, -1. ,solutionImportedC, 0.);

        // Taking abs norm
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > errorValuesAbs = errorValues;
        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > errorValuesAbsC = errorValuesC;

        errorValues->abs(errorValuesAbs);
        errorValuesC->abs(errorValuesAbsC);

        Teuchos::Array<SC> norm(1),normC(1); 
        errorValues->normInf(norm);
        errorValuesC->normInf(normC);

        double res = norm[0];
        if(comm->getRank() ==0){
            std::cout << " Inf Norm of Error of Solution Displacement " << norm[0] << std::endl;
            std::cout << " Inf Norm of Error of Solution Concentration " << normC[0] << std::endl;
        }

        errorValues->norm2(norm);//const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms);
        errorValuesC->norm2(normC);//const Teuchos::ArrayView<typename Teuchos::ScalarTraits<SC>::magnitudeType> &norms);

        res = norm[0];
        if(comm->getRank() ==0){
            std::cout << " 2 Norm of Error of Solution Concentration " << norm[0] << std::endl;
            std::cout << " 2 Norm of Error of Solution Concentration " << normC[0] << std::endl;

        }
        double NormError = norm[0];
        double NormErrorC = normC[0];

        // Relative to the uninterrupted solution of each field
        Teuchos::Array<SC> reference(1), referenceC(1);
        solutionImported->norm2(reference);
        solutionImportedC->norm2(referenceC);
        double relativeError = reference[0] > 0. ? NormError/reference[0] : NormError;
        double relativeErrorC = referenceC[0] > 0. ? NormErrorC/referenceC[0] : NormErrorC;
        if(comm->getRank() ==0){
            std::cout << " Restart test: d_s at t = " << finalTime << ", relative 2-norm of the difference to the uninterrupted run: " << relativeError << " (tolerance " << tolerance << ")" << std::endl;
            std::cout << " Restart test: c at t = " << finalTime << ", relative 2-norm of the difference to the uninterrupted run: " << relativeErrorC << " (tolerance " << tolerance << ")" << std::endl;
        }
        if(!(relativeError <= tolerance) || !(relativeErrorC <= tolerance))
            testPassed = false;
        if(!(reference[0] > 0.) && !(referenceC[0] > 0.)){ // comparing identically zero solutions checks nothing
            if(comm->getRank() ==0)
                std::cout << " Restart test: the uninterrupted solution is zero at t = " << finalTime << ", the comparison checks nothing" << std::endl;
            testPassed = false;
        }

        // Adding error to paraview exporter
        exParaResults->addVariable(errorValuesAbs, "d-dImport", "Vector", 3,  domainStructure->getMapUnique());
        exParaResultsC->addVariable(errorValuesAbsC, "c-cImport", "Scalar", 1,  domainChem->getMapUnique());

        exParaResults->save(0.0);
        exParaResultsC->save(0.0);
        }

    }
    TimeMonitor_Type::report(std::cout);
        stackedTimer->stop("Structure-chemical interaction");
	StackedTimer::OutputOptions options;
	options.output_fraction = options.output_histogram = options.output_minmax = true;
	stackedTimer->report((std::cout),comm,options);
	
    return testPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}
