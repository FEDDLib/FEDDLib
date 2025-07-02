#include <Tpetra_Core.hpp>

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"

#include "feddlib/problems/specific/DiffusionReaction.hpp"
#include "feddlib/problems/Solver/DAESolverInTime.hpp"

/*!
 main of unsteadyLinearDiffusion problem

 @brief unsteadyLinearDiffusion main
 @author Lea Saßmannshausen
 @version 1.0
 @copyright LS
 */

void zeroBC(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
}
void oneBC(double* x, double* res, double t, const double* parameters){
    res[0] = 1.;
}
void twoBC(double* x, double* res, double t, const double* parameters){
    res[0] = 2.;
}
void threeBC(double* x, double* res, double t, const double* parameters){
    res[0] = 3.;
}
void zeroBC2D(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
    res[1] = 0.;
}
void zeroBC3D(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
    res[1] = 0.;
    res[2] = 0.;
}

void oneFunc(double* x, double* res, double* parameters){
    res[0] = 1.;
}

void zeroFunc(double* x, double* res, double* parameters){
    res[0] = 0.;
}


void inflowChem(double* x, double* res, double t, const double* parameters)
{
    res[0] = 1.;
    
    return;
}

void reactionFunc(double* x, double* res, double* parameters){
	
    double m = 1.;	
    res[0] = m * x[0];

}

typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;

    // MPI boilerplate
    Tpetra::ScopeGuard tpetraScope (&argc, &argv); // initializes MPI
    Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;

    std::string vectorLaplace = "false";
    myCLP.setOption("vectorLaplace",&vectorLaplace,"vectorLaplace");
    std::string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    std::string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    std::string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
    double length = 4.;
    myCLP.setOption("length",&length,"length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        return EXIT_SUCCESS;
    }
    bool vL ( !vectorLaplace.compare("true") );
    
    {
        ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
        ParameterListPtr_Type parameterListPrec = Teuchos::getParametersFromXmlFile(xmlPrecFile);
        ParameterListPtr_Type parameterListSolver = Teuchos::getParametersFromXmlFile(xmlSolverFile);

        ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;
        parameterListAll->setParameters(*parameterListPrec);
        parameterListAll->setParameters(*parameterListSolver);

        // Mesh
        int dim = parameterListProblem->sublist("Parameter").get("Dimension",2);
        int level = parameterListProblem->sublist("Parameter").get("Refinement Level",2);
        int m = parameterListProblem->sublist("Parameter").get("H/h",5);
        std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization","P1");
        std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
        std::string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","");
        std::string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");

        int n;
        int size = comm->getSize();
        int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
        size -= numProcsCoarseSolve;

        int minNumberSubdomains;
        if (!meshType.compare("structured") || !meshType.compare("unstructured_struct")) {
            minNumberSubdomains = 1;
        }
        else if(!meshType.compare("structured_bfs") || !meshType.compare("unstructured_bfs")){
            minNumberSubdomains = (int) 2*length+1;
        }


        Teuchos::RCP<Domain<SC,LO,GO,NO> > domain;

        Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP1;
        Teuchos::RCP<Domain<SC,LO,GO,NO> > domainP2;
        domainP1.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
        
        MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
        domainP1Array[0] = domainP1;
        
        ParameterListPtr_Type pListPartitioner = sublist( parameterListAll, "Mesh Partitioner" );
        MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
        
        partitionerP1.readAndPartition();

	if (FEType=="P2") {
            domainP2.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
            domainP2->buildP2ofP1Domain( domainP1 );
            domain = domainP2;
        }
        else
            domain = domainP1;
        



        // ####################
        Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory(new BCBuilder<SC,LO,GO,NO>( ));

        if(dim==3){
            // Inflow through one surface. Some flags represent edge components.
            bcFactory->addBC(inflowChem, 0, 0, domain, "Dirichlet", 1); // inflow of Chem
            bcFactory->addBC(inflowChem, 1, 0, domain, "Dirichlet", 1); // inflow of Chem
            bcFactory->addBC(inflowChem, 7, 0, domain, "Dirichlet", 1);            		
            bcFactory->addBC(inflowChem, 9, 0, domain, "Dirichlet", 1);
           
		}
		else if(dim==2){
            // Inflow on one side / edge
		    bcFactory->addBC(inflowChem, 2, 0, domain, "Dirichlet", 1); // inflow of Chem
		}
	
        
        // We construct the diffusion tensor. In general it is a scale unit matrix. 
	    vec2D_dbl_Type diffusionTensor(dim,vec_dbl_Type(dim));
        double D0 = parameterListAll->sublist("Parameter").get("D0",1.);
        for(int i=0; i<dim; i++){
            diffusionTensor[i][i] =D0;
            if(i>0){
                diffusionTensor[i][i-1] = 0;
                diffusionTensor[i-1][i] = 0;
            }
            else
                diffusionTensor[i][i+1] = 0;				
        }
 
        DiffusionReaction<SC,LO,GO,NO> diffusionReaction(domain,FEType,parameterListAll, diffusionTensor,reactionFunc);
        {
       
            domain->info();
            diffusionReaction.info();
            
            {

                diffusionReaction.addBoundaries( bcFactory );

                // diffusionReaction.addRhsFunction( zeroFunc );
                // diffusionReaction.addParemeterRhs( parameterListProblem->sublist("Parameter").get("Metabolic Rate",0.0) );

                diffusionReaction.initializeProblem();
                diffusionReaction.assemble();
                                    
                diffusionReaction.setBoundaries();
                
                DAESolverInTime<SC,LO,GO,NO> daeTimeSolver(parameterListAll, comm);

                // Fuer das Strukturproblem haben wir analog, da nur d_s als Varable vorhanden:
                SmallMatrix<int> defTS(1);
                defTS[0][0] = 1;

                // Uebergebe auf welchen Bloecken die Zeitintegration durchgefuehrt werden soll
                // und Uebergabe der parameterList, wo die Parameter fuer die Zeitintegration drin stehen
                daeTimeSolver.defineTimeStepping(defTS);

                // Uebergebe das (nicht) lineare Problem
                daeTimeSolver.setProblem(diffusionReaction);

                // Setup fuer die Zeitintegration, wie z.B. Aufstellen der Massematrizen auf den Zeilen, welche in
                // defTS definiert worden sind.
                daeTimeSolver.setupTimeStepping();

                // Fuehre die komplette Zeitintegration + ggf. Newton/Fixpunkt + Loesen + Exporter durch
                daeTimeSolver.advanceInTime();



            }


       		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

            Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolution = diffusionReaction.getSolution()->getBlock(0);

            exPara->setup("solutionDiffusionReaction", domain->getMesh(), FEType);
            
            exPara->addVariable(exportSolution, "u", "Scalar", 1, domain->getMapUnique());

            exPara->save(0.0);

         }





        

       
    }
    return EXIT_SUCCESS;
}
