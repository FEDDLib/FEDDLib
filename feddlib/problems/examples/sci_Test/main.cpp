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
    res[0] = 1.;
    
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

void rhsYZ(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[0] = 0.;
    double force = parameters[1];
    double TRamp = 0.0;
       
    if(parameters[0] <= TRamp+1e-06)
        force = parameters[0] * force * 1./(TRamp);
    else
        force = parameters[1];

    if(parameters[2] == 5)
        res[1] = force;
    else
        res[1] =0.;
        
    if (parameters[2] == 4)
        res[2] = force;
    else
        res[2] = 0.;
    
    return;
}

void rhsHeartBeatCube(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[2] = 0.;
    double force = parameters[1];
    double TRamp = 1001.0;
    
	double a0    = 11.693284502463376;
	double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
	double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		         
    double Q = 0.5*a0;
    

    double t_min = parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
    double t_max = t_min + 1.0; // One heartbeat lasts 1.5 second    
    double y = M_PI * ( 2.0*( parameters[0]-t_min ) / ( t_max - t_min ) -1.0  );
    
    for(int i=0; i< 20; i++)
        Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
    
    
    // Remove initial offset due to FFT
    Q -= 0.026039341343493;
    Q = (Q - 2.85489)/(7.96908-2.85489);
    
    if(parameters[0]< TRamp){
    	Q = 0.;
    }
    
    if(parameters[0] <= 1.)
       force = force * parameters[0];
    
    if(parameters[2] == 5)
        res[1] = force+Q*0.005329;
    else
        res[1] =0.;
        
    if (parameters[2] == 4)
        res[2] = force+Q*0.005329;
    else
        res[2] = 0.;
}

void rhsHeartBeatArtery(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[2] = 0.;
    double force = parameters[1];
    double TRamp = 1001.0;
    
	double a0    = 11.693284502463376;
	double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
	double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		         
    double Q = 0.5*a0;
    

    double t_min = parameters[0] - fmod(parameters[0],1.0); //FlowConditions::t_start_unsteady;
    double t_max = t_min + 1.0; // One heartbeat lasts 1.5 second    
    double y = M_PI * ( 2.0*( parameters[0]-t_min ) / ( t_max - t_min ) -1.0  );
    
    for(int i=0; i< 20; i++)
        Q += (a[i]*std::cos((i+1.)*y) + b[i]*std::sin((i+1.)*y) ) ;
    
    
    // Remove initial offset due to FFT
    Q -= 0.026039341343493;
    Q = (Q - 2.85489)/(7.96908-2.85489);
    
    if(parameters[0]< TRamp){
    	Q = 0.;
    }
    if(parameters[0] <= 1.)
        force = force * parameters[0];
        
    if(parameters[2]==5){
        res[0] = x[0];
        res[1] = x[1];
        double r2= sqrt(res[0]*res[0]+res[1]*res[1]);
        res[0] = res[0]*(force+Q*0.005329);
        res[1] = res[1]*(force+Q*0.005329);
       
    }
    else{
        res[0] =0.;
        res[1] =0.;
    }
  
}

void rhsArteryPaper(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[2] = 0.;
    double force = parameters[1];
    double TRamp = 2001.0;
    double lambda=0.;
    
    if(parameters[0] < 1.)
        lambda = 0.875 * parameters[0];
    else if(parameters[0] <= TRamp)
    	lambda = 0.875;
    else if( parameters[0] <= 2001.5 )
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= 2001.5 && (parameters[0] - std::floor(parameters[0]))<= 0.5)
    	lambda= 0.75;
    else
        lambda = 0.875 - 0.125 * cos(4*M_PI*(parameters[0]+0.02));
    
    if(parameters[2]==5){
        res[0] = x[0];
        res[1] = x[1];
        double r2= sqrt(res[0]*res[0]+res[1]*res[1]);
        res[0] = res[0]*lambda*force;
        res[1] = res[1]*lambda*force;
       
    }
    else{
        res[0] =0.;
        res[1] =0.;
    }
            
}

void rhsCubePaper(double* x, double* res, double* parameters){
    // parameters[0] is the time, not needed here
    res[2] = 0.;
    double force = parameters[1];
    double TRamp = 2001.0;
    double lambda=0.;
    
    if(parameters[0] < 1.)
        lambda = 0.875 * parameters[0];
    else if(parameters[0] <= TRamp)
    	lambda = 0.875;
    else if( parameters[0] <= 2001.5 )
		lambda = 0.8125+0.0625*cos(2*M_PI*parameters[0]);
    else if( parameters[0] >= 2001.5 && (parameters[0] - std::floor(parameters[0]))<= 0.5)
    	lambda= 0.75;
    else
        lambda = 0.875 - 0.125 * cos(4*M_PI*(parameters[0]+0.02));
     
    if(parameters[2] == 5)
        res[1] = force*lambda;
    else
        res[1] =0.;
        
    if (parameters[2] == 4)
        res[2] = force*lambda;
    else
        res[2] = 0.;
            
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
    string xmlProblemFile = "parametersProblemSCI.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");       
    string xmlSolverFileSCI = "parametersSolverSCI.xml"; 
    myCLP.setOption("solverfileSCI",&xmlSolverFileSCI,".xml file with Inputparameters.");
    
    string xmlPrecFileStructure = "parametersPrecStructure.xml";
    myCLP.setOption("precfileStructure",&xmlPrecFileStructure,".xml file with Inputparameters.");
    string xmlPrecFileChem = "parametersPrecChem.xml";
    myCLP.setOption("precfileChem",&xmlPrecFileChem,".xml file with Inputparameters.");

    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");


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
       
        ParameterListPtr_Type parameterListSolverSCI = Teuchos::getParametersFromXmlFile(xmlSolverFileSCI);

        ParameterListPtr_Type parameterListPrecStructure = Teuchos::getParametersFromXmlFile(xmlPrecFileStructure);
        ParameterListPtr_Type parameterListPrecChem = Teuchos::getParametersFromXmlFile(xmlPrecFileChem);
  
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);


        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;     
        
        parameterListAll->setParameters(*parameterListSolverSCI);
        parameterListAll->setParameters(*parameterListPrec);

        
        ParameterListPtr_Type parameterListChemAll(new Teuchos::ParameterList(*parameterListPrecChem)) ;
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Chem") );
        sublist(parameterListChemAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter") );
        parameterListChemAll->setParameters(*parameterListPrecChem);

        
        ParameterListPtr_Type parameterListStructureAll(new Teuchos::ParameterList(*parameterListPrecStructure));
        sublist(parameterListStructureAll, "Parameter")->setParameters( parameterListProblem->sublist("Parameter Solid") );
        parameterListStructureAll->setParameters(*parameterListPrecStructure);
        parameterListStructureAll->setParameters(*parameterListProblem);

                 
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
            cout << "############ Starting SCI  ... ################" <<endl;
            cout << "###############################################" <<endl;
        }

        DomainPtr_Type domainP1chem;
        DomainPtr_Type domainP1struct;
        DomainPtr_Type domainP2chem;
        DomainPtr_Type domainP2struct;
        
        
        DomainPtr_Type domainChem;
        DomainPtr_Type domainStructure;
        
        std::string bcType = parameterListAll->sublist("Parameter").get("BC Type","Cube");
        
        std::string rhsType = parameterListAll->sublist("Parameter").get("RHS Type","Constant");
    
        domainP1chem.reset( new Domain_Type( comm, dim ) );
        domainP1struct.reset( new Domain_Type( comm, dim ) );
        domainP2chem.reset( new Domain_Type( comm, dim ) );
        domainP2struct.reset( new Domain_Type( comm, dim ) );
                                
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
        domainP1Array[0] = domainP1struct;
       // domainP1Array[1] = domainP1struct;
    
        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );                    

        pListPartitioner->set("Build Edge List",true);
        pListPartitioner->set("Build Surface List",true);
                        
        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
        
        int volumeID=10;
        if(bcType=="Artery")
        	volumeID = 15;
        partitionerP1.readAndPartition(volumeID);
                    
        if (!discType.compare("P2")){
			domainP2chem->buildP2ofP1Domain( domainP1struct );
			domainP2struct->buildP2ofP1Domain( domainP1struct );

			domainChem = domainP2chem;
			domainStructure = domainP2struct;   
		}        
		else{
			domainStructure = domainP1struct;
			domainChem = domainP1struct;
		}

       // ########################
        // Flags check
        // ########################

		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaF(new ExporterParaView<SC,LO,GO,NO>());

		Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapUnique()));
		vec_int_ptr_Type BCFlags = domainStructure->getBCFlagUnique();

		Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
		for(int i=0; i< entries.size(); i++){
			entries[i] = BCFlags->at(i);
		}

		Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

		exParaF->setup("Flags", domainStructure->getMesh(), discType);

		exParaF->addVariable(exportSolutionConst, "Flags", "Scalar", 1,domainStructure->getMapUnique(), domainStructure->getMapUniqueP2());

		exParaF->save(0.0);
		
		/*double a0    = 11.693284502463376;
		double a [20] = {1.420706949636449,-0.937457438404759,0.281479818173732,-0.224724363786734,0.080426469802665,0.032077024077824,0.039516941555861, 
		  0.032666881040235,-0.019948718147876,0.006998975442773,-0.033021060067630,-0.015708267688123,-0.029038419813160,-0.003001255512608,-0.009549531539299, 
		  0.007112349455861,0.001970095816773,0.015306208420903,0.006772571935245,0.009480436178357};
		double b [20] = {-1.325494054863285,0.192277311734674,0.115316087615845,-0.067714675760648,0.207297536049255,-0.044080204999886,0.050362628821152,-0.063456242820606,
		  -0.002046987314705,-0.042350454615554,-0.013150127522194,-0.010408847105535,0.011590255438424,0.013281630639807,0.014991955865968,0.016514327477078, 
		  0.013717154383988,0.012016806933609,-0.003415634499995,0.003188511626163};
		
		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exParaHeartBeat(new ExporterParaView<SC,LO,GO,NO>());

		
		Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolutionBeat(new MultiVector<SC,LO,GO,NO>(domainStructure->getMapUnique()));

		Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConstBeat = exportSolutionBeat;

		exParaHeartBeat->setup("HeartBeat", domainStructure->getMesh(), discType);

		exParaHeartBeat->addVariable(exportSolutionConstBeat, "Beat", "Scalar", 1,domainStructure->getMapUnique(), domainStructure->getMapUniqueP2());
		

	
		Teuchos::ArrayRCP< SC > entriesPulse  = exportSolutionBeat->getDataNonConst(0);
		double dt = 0.02;
		double lambda =0.;
		for(int i= 1; i<1500; i++){	         
			double Q = 0.5*a0;
		

			double t_min = dt * i - fmod(dt*i,1.0); //FlowConditions::t_start_unsteady;
			double t_max = t_min + 1.0; // One heartbeat lasts 1.5 second    
			double y = M_PI * ( 2.0*( dt * i-t_min ) / ( t_max - t_min ) -1.0  );
		
		
		
			for(int j=0; j< 20; j++)
				Q += (a[j]*std::cos((j+1.)*y) + b[j]*std::sin((j+1.)*y) ) ;
			
			
			// Remove initial offset due to FFT
			Q -= 0.026039341343493;
			//Q = (Q - 4.3637)/3.5427;
			entriesPulse[0] = Q ;
			
			if(dt*i < 20.)
				lambda = 0.875;
			else if( dt*i <= 20.5 )
				lambda = 0.8125+0.0625*cos(2*M_PI*dt*i);
			else if ( dt*i >= 20.5 && (dt*i - std::floor(dt*i))<= 0.5)
				lambda= 0.75;
			else
				lambda = 0.875 - 0.125 * cos(4*M_PI*(dt*i+0.02));

			entriesPulse[0] = 0.02133*lambda;
			
			
			exParaHeartBeat->save(double (dt*i));
		}*/



        if (parameterListAll->sublist("General").get("ParaView export subdomains",false) ){
            
            if (verbose)
                std::cout << "\t### Exporting subdomains ###\n";

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

        defTS.reset( new SmallMatrix<int> (2) );

        // Stucture
        (*defTS)[0][0] = 1;
        // Chem
        (*defTS)[1][1] = 1;
			

        SCI<SC,LO,GO,NO> sci(domainStructure, discType,
                                domainChem, discType, diffusionTensor, reactionFunc,
                                parameterListStructureAll,
                                parameterListChemAll,
                                parameterListAll,
                                defTS);
        
        sci.info();
        
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) ); 
            
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryChem( new BCBuilder<SC,LO,GO,NO>( ) ); 
        
    
        // Struktur-RW
        
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactoryStructure( new BCBuilder<SC,LO,GO,NO>( ) );

        if(dim == 2)
        {
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
        }
        else if(dim == 3 && bcType=="Cube")
        {

            bcFactory->addBC(zeroDirichlet, 1, 0, domainStructure, "Dirichlet_X", dim);
            bcFactory->addBC(zeroDirichlet, 2, 0, domainStructure, "Dirichlet_Y", dim);
            bcFactory->addBC(zeroDirichlet, 3, 0, domainStructure, "Dirichlet_Z", dim);
            bcFactory->addBC(zeroDirichlet, 4, 0, domainStructure, "Dirichlet_Z", dim);
            
            bcFactory->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet", dim);
            bcFactory->addBC(zeroDirichlet2D, 7, 0, domainStructure, "Dirichlet_X_Y", dim);
            bcFactory->addBC(zeroDirichlet2D, 8, 0, domainStructure, "Dirichlet_Y_Z", dim);
            bcFactory->addBC(zeroDirichlet2D, 9, 0, domainStructure, "Dirichlet_X_Z", dim);
            
            bcFactoryStructure->addBC(zeroDirichlet, 1, 0, domainStructure, "Dirichlet_X", dim);
            bcFactoryStructure->addBC(zeroDirichlet, 2, 0, domainStructure, "Dirichlet_Y", dim);
            bcFactoryStructure->addBC(zeroDirichlet, 3, 0, domainStructure, "Dirichlet_Z", dim);
           bcFactoryStructure->addBC(zeroDirichlet, 4, 0, domainStructure, "Dirichlet_Z", dim);
            
            bcFactoryStructure->addBC(zeroDirichlet3D, 0, 0, domainStructure, "Dirichlet", dim);
            bcFactoryStructure->addBC(zeroDirichlet2D, 7, 0, domainStructure, "Dirichlet_X_Y", dim);
            bcFactoryStructure->addBC(zeroDirichlet2D, 8, 0, domainStructure, "Dirichlet_Y_Z", dim);
            bcFactoryStructure->addBC(zeroDirichlet2D, 9, 0, domainStructure, "Dirichlet_X_Z", dim);

        }
        else if(dim==3 && bcType=="Artery"){
        
			bcFactory->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_Y", dim);
			bcFactory->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_X", dim);
			bcFactory->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Z", dim);

			bcFactory->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_X_Z", dim);

			bcFactory->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X", dim);
			bcFactory->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Y", dim);

			bcFactory->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactory->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_X_Z", dim);


			bcFactoryStructure->addBC(zeroDirichlet3D, 1, 0, domainStructure, "Dirichlet_Y", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 2, 0, domainStructure, "Dirichlet_X", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 3, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 4, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 13, 0, domainStructure, "Dirichlet_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 14, 0, domainStructure, "Dirichlet_Z", dim);


			bcFactoryStructure->addBC(zeroDirichlet3D, 9, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 8, 0, domainStructure, "Dirichlet_X_Z", dim);

			bcFactoryStructure->addBC(zeroDirichlet3D, 7, 0, domainStructure, "Dirichlet_X", dim);

			bcFactoryStructure->addBC(zeroDirichlet3D, 10, 0, domainStructure, "Dirichlet_Y", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 11, 0, domainStructure, "Dirichlet_Y_Z", dim);
			bcFactoryStructure->addBC(zeroDirichlet3D, 12, 0, domainStructure, "Dirichlet_X_Z", dim);
        
        }
        
        // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        if (!sci.problemStructure_.is_null())
            sci.problemStructure_->addBoundaries(bcFactoryStructure);
        else
            sci.problemStructureNonLin_->addBoundaries(bcFactoryStructure);
        
        // RHS dummy for structure
        if (dim==2) {
            if (!sci.problemStructure_.is_null())
                sci.problemStructure_->addRhsFunction( rhsX,0 );
            else
                sci.problemStructureNonLin_->addRhsFunction( rhsX,0 );
            
        }
        else if (dim==3) {
            
            if (!sci.problemStructure_.is_null()){
                if(bcType=="Cube"){
					if(rhsType=="Constant")
		    		 	sci.problemStructure_->addRhsFunction( rhsYZ,0 );
		    		if(rhsType=="Paper")
		    		 	sci.problemStructure_->addRhsFunction( rhsCubePaper,0 );
		    		if(rhsType=="Heart Beat")
		        		 	sci.problemStructure_->addRhsFunction( rhsHeartBeatCube,0 );
					
				}
				else if(bcType=="Artery"){
					if(rhsType=="Paper")
            		 	sci.problemStructure_->addRhsFunction( rhsArteryPaper,0 );
            		if(rhsType=="Heart Beat")
            		 	sci.problemStructure_->addRhsFunction( rhsHeartBeatArtery,0 );
				}
                     
                double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
                sci.problemStructure_->addParemeterRhs( force );
                double degree = 0.;
                sci.problemStructure_->addParemeterRhs( degree );

            }
            else{             
				if(bcType=="Cube"){
					if(rhsType=="Constant")
		    		 	sci.problemStructureNonLin_->addRhsFunction( rhsYZ,0 );
		    		if(rhsType=="Paper")
		        		sci.problemStructureNonLin_->addRhsFunction( rhsCubePaper,0 );
		    		if(rhsType=="Heart Beat")
		        		sci.problemStructureNonLin_->addRhsFunction( rhsHeartBeatCube,0 );
					
				}
				else if(bcType=="Artery"){
					if(rhsType=="Paper")
            		 	sci.problemStructureNonLin_->addRhsFunction( rhsArteryPaper,0 );
            		if(rhsType=="Heart Beat")
            		 	sci.problemStructureNonLin_->addRhsFunction( rhsHeartBeatArtery,0 );
				}
                double force = parameterListAll->sublist("Parameter").get("Volume force",1.);
                sci.problemStructureNonLin_->addParemeterRhs( force );
                double degree = 0.;
                sci.problemStructureNonLin_->addParemeterRhs( degree );

            }
            

        }
        if (dim==2)
        {
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Only 3D Test available");                               
                            
        }
        else if(dim==3 && bcType=="Cube")
        {

            bcFactory->addBC(inflowChem, 0, 1, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactory->addBC(inflowChem, 1, 1, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactory->addBC(inflowChem, 7, 1, domainChem, "Dirichlet", 1);            		
            //bcFactory->addBC(zeroDirichlet, 8, 1, domainChem, "Dirichlet", 1);
            bcFactory->addBC(inflowChem, 9, 1, domainChem, "Dirichlet", 1);
            /*bcFactory->addBC(zeroDirichlet, 2, 1, domainChem, "Dirichlet", 1);
            bcFactory->addBC(zeroDirichlet, 3, 1, domainChem, "Dirichlet", 1);            
            bcFactory->addBC(zeroDirichlet, 4, 1, domainChem, "Dirichlet", 1);            
            bcFactory->addBC(zeroDirichlet, 5, 1, domainChem, "Dirichlet", 1);            
           // bcFactory->addBC(zeroDirichlet, 6, 1, domainChem, "Dirichlet", 1);            
            */
            
            bcFactoryChem->addBC(inflowChem, 0, 0, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactoryChem->addBC(inflowChem, 1, 0, domainChem, "Dirichlet", 1); // inflow of Chem
            bcFactoryChem->addBC(inflowChem, 7, 0, domainChem, "Dirichlet", 1);            		
            bcFactoryChem->addBC(inflowChem, 9, 0, domainChem, "Dirichlet", 1);
           /* bcFactoryChem->addBC(zeroDirichlet, 2, 0, domainChem, "Dirichlet", 1);
            bcFactoryChem->addBC(zeroDirichlet, 3, 0, domainChem, "Dirichlet", 1);            
            bcFactoryChem->addBC(zeroDirichlet, 4, 0, domainChem, "Dirichlet", 1);            
            bcFactoryChem->addBC(zeroDirichlet, 5, 0, domainChem, "Dirichlet", 1);            
            //bcFactoryChem->addBC(zeroDirichlet, 6, 0, domainChem, "Dirichlet", 1);
            bcFactoryChem->addBC(zeroDirichlet, 8, 0, domainChem, "Dirichlet", 1);
            
            */
        }
        else if(dim==3 && bcType=="Artery"){
           bcFactory->addBC(inflowChem, 5, 1, domainChem, "Dirichlet", 1); // inflow of Chem
		   bcFactory->addBC(inflowChem, 13, 1, domainChem, "Dirichlet", 1); // inflow of Chem
		   bcFactory->addBC(inflowChem, 14, 1, domainChem, "Dirichlet", 1); // inflow of Chem
		   bcFactory->addBC(inflowChem, 7, 1, domainChem, "Dirichlet", 1); // inflow of Chem
		   bcFactory->addBC(inflowChem, 10, 1, domainChem, "Dirichlet", 1); // inflow of Chem
		   
           bcFactoryChem->addBC(inflowChem, 5, 0, domainChem, "Dirichlet", 1);
           bcFactoryChem->addBC(inflowChem, 13, 0, domainChem, "Dirichlet", 1);
           bcFactoryChem->addBC(inflowChem, 14, 0, domainChem, "Dirichlet", 1);
           bcFactoryChem->addBC(inflowChem, 7, 0, domainChem, "Dirichlet", 1);
           bcFactoryChem->addBC(inflowChem, 10, 0, domainChem, "Dirichlet", 1);
        }

        // Fuer die Teil-TimeProblems brauchen wir bei TimeProblems
        // die bcFactory; vgl. z.B. Timeproblem::updateMultistepRhs()
        sci.problemChem_->addBoundaries(bcFactoryChem);
        
          
        // #####################
        // Zeitintegration
        // #####################
        sci.addBoundaries(bcFactory); // Dem Problem RW hinzufuegen

        sci.initializeProblem();
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
    }
    TimeMonitor_Type::report(std::cout);

    return(EXIT_SUCCESS);
}
