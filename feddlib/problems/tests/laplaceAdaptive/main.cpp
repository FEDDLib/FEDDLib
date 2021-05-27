#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"

#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FiniteElement.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructuredRefinement.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/LinearAlgebra/BlockMultiVector.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/problems/specific/Laplace.hpp"
#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/core/LinearAlgebra/BlockMultiVector.hpp"

#include "feddlib/core/Mesh/Mesh.hpp"
#include "feddlib/core/Mesh/MeshInterface.hpp"
# include "feddlib/core/Mesh/MeshFileReader.hpp"

#include <boost/function.hpp>
#include <chrono> 


/*!
 main of Laplace problem

 @brief Laplace main
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

void zeroBC(double* x, double* res, double t, const double* parameters){
    res[0] = 0.;
}

void bc2D(double* x, double* res, double t, const double* parameters){

	double r = sqrt(x[0]*x[0] + x[1]*x[1]);
    double phi;
    if(x[1] < 0.0)
		phi = 2.0*M_PI+atan2(x[1],x[0]);
    else
		phi = atan2(x[1],x[0]);
	
    res[0] =  pow(r,2/3.)*sin(2/3.*phi); 

}


// Function for the rhs
void rhs2D(double* x, double* res, double* parameters){
    
    res[0] = 0.;


    return;
}

void bc3D(double* x, double* res, double t, const double* parameters){

    res[0] =  -1/6. *(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)); 

}


// Function for the rhs
void rhs3D(double* x, double* res, double* parameters){
    
    res[0] = 1.;


    return;
}



void dummyFunc(double* x, double* res, double t, const double* parameters){

    return;
}


typedef unsigned UN;
typedef default_sc SC;
typedef default_lo LO;
typedef default_go GO;
typedef default_no NO;

using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::outArg;

using namespace FEDD;

int main(int argc, char *argv[]) {
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;
	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
	typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
	typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
	typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;
	typedef EdgeElements EdgeElements_Type;
	typedef Teuchos::RCP<EdgeElements_Type> EdgeElementsPtr_Type;
	typedef Mesh<SC,LO,GO,NO> Mesh_Type;
	typedef Teuchos::RCP<Mesh_Type > MeshPtr_Type;
	typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
	typedef Teuchos::RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
	typedef typename Mesh_Type::Elements_Type Elements_Type;
	typedef typename Mesh_Type::ElementsPtr_Type ElementsPtr_Type;
    typedef Map<LO,GO,NO> Map_Type;
    typedef typename Map_Type::MapPtr_Type MapPtr_Type;

    Teuchos::oblackholestream blackhole;
    Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    Teuchos::RCP<const Teuchos::Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

    // Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");

    std::string vectorLaplace = "false";
    myCLP.setOption("vectorLaplace",&vectorLaplace,"vectorLaplace");
    string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");
    string xmlPrecFile = "parametersPrec.xml";
    myCLP.setOption("precfile",&xmlPrecFile,".xml file with Inputparameters.");
    string xmlSolverFile = "parametersSolver.xml";
    myCLP.setOption("solverfile",&xmlSolverFile,".xml file with Inputparameters.");
    double length = 1.0;
    myCLP.setOption("length",&length,"length of domain.");

    myCLP.recogniseAllOptions(true);
    myCLP.throwExceptions(false);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = myCLP.parse(argc,argv);
    if(parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
        mpiSession.~GlobalMPISession();
        return 0;
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

		domain = domainP1;
	

		// MeshRefinement Parameters
		double tol= parameterListProblem->sublist("Mesh Refinement").get("Toleranz",0.001);
		double theta = parameterListProblem->sublist("Mesh Refinement").get("Theta",0.35);
		string strategy = parameterListProblem->sublist("Mesh Refinement").get("RefinementType","Uniform");
		int maxIter = parameterListProblem->sublist("Mesh Refinement").get("MaxIter",3);
		bool checkRestrictions = parameterListProblem->sublist("Mesh Refinement").get("Check Restrictions",true);
		string restriction = parameterListProblem->sublist("Mesh Refinement").get("Restriction Type","keepRegularity");

		vec_dbl_Type maxErrorEl(maxIter+1);
		vec_int_Type numElements(maxIter+1);
		vec_int_Type numElementsProc(maxIter+1);


	  	MeshPartitioner_Type::DomainPtrArray_Type domainP1RefinedArray(1);
		domainP1RefinedArray[0] = domainP1;
		Teuchos::RCP<Domain<SC,LO,GO,NO> > domainRefined;
				
		Teuchos::RCP<BCBuilder<SC,LO,GO,NO> > bcFactory( new BCBuilder<SC,LO,GO,NO>( ) );

		Teuchos::RCP<const MultiVector<SC,LO,GO,NO> >  exportSolution;
		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara2;
		Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

		std::chrono::duration<double> elapsed_total;
		std::chrono::duration<double> elapsed_MeshRef=(std::chrono::duration<double>) 0;
		std::chrono::duration<double> elapsed_MeshRefTmp =(std::chrono::duration<double>) 0;
		std::chrono::duration<double> elapsed_Solv=(std::chrono::duration<double>) 0;
		
		std::vector<std::chrono::duration<double>> meshTiming(maxIter);

		Teuchos::RCP<Teuchos::Time> buildMesh(Teuchos::TimeMonitor::getNewCounter("main: Refine Mesh"));

		auto startTotal = std::chrono::high_resolution_clock::now();

		int maxRank = std::get<1>(domain->getMesh()->rankRange_);
		int j=0;
		while(j<maxIter+1 ){

			if (FEType=="P2" ) {
					domainP2.reset( new Domain<SC,LO,GO,NO>( comm, dim ));
					domainP2->buildP2ofP1Domain( domainP1 );
					domain = domainP2;
			   		 }
			else 
					domain = domainP1; 

			// Initialize meshUnstrRefinement Type for domain
			domainP1->initMeshRef(domainP1);
			domain->initMeshRef(domain);

			if(dim == 2){
				bcFactory.reset( new BCBuilder<SC,LO,GO,NO>( ) );
		   		bcFactory->addBC(bc2D, 1, 0, domain, "Dirichlet", 1);
			}
			if(dim == 3){
				bcFactory.reset( new BCBuilder<SC,LO,GO,NO>( ) );
		   		bcFactory->addBC(bc3D, 1, 0, domain, "Dirichlet", 1);
			}
		   
			auto startSolv = std::chrono::high_resolution_clock::now();
			Laplace<SC,LO,GO,NO> laplace(domain,FEType,parameterListAll,vL);
				{
				laplace.addBoundaries(bcFactory);
				if(dim==2)
					laplace.addRhsFunction(rhs2D);
				if(dim==3)
					laplace.addRhsFunction(rhs3D);
		        laplace.initializeProblem();
				laplace.assemble();
				laplace.setBoundaries();
				laplace.solve();
				}

			auto finishSolv = std::chrono::high_resolution_clock::now();
			elapsed_Solv = elapsed_Solv + finishSolv- startSolv;
				
			const MultiVectorPtrConst_Type valuesSolution = laplace.getSolution()->getBlock(0);

			// Error Estimation and tagging of Elements
			vec_dbl_Type errorElement  = domainP1->errorEstimation(valuesSolution, theta, strategy);

			// Determine max error of elements on Proc
			vec_dbl_Type::iterator it;
			it = max_element(errorElement.begin(), errorElement.end());
			double maxErrorProc = errorElement.at(distance(errorElement.begin(), it)); // accumulate(errorElement.begin(), errorElement.end(),0.0);

		 	// Export Solution
			exPara.reset(new ExporterParaView<SC,LO,GO,NO>());

			exportSolution = laplace.getSolution()->getBlock(0);

			if(dim==2)
                exPara->setup("MeshRefinement2D", domain->getMesh(), FEType, parameterListAll);				
			else
                exPara->setup("MeshRefinement3D", domain->getMesh(), FEType, parameterListAll);
//
			exPara->addVariable(exportSolution, "u", "Scalar", 1, domain->getMapUnique(), domain->getMapUniqueP2());
			exPara->save(0.0);
			
	 		// Export distribution of elements among processors
			MultiVectorPtr_Type procNumTmp = Teuchos::rcp( new MultiVector_Type(domain->getMapUnique() , 1 ) );
			procNumTmp->putScalar(comm->getRank());
			MultiVectorPtrConst_Type procNum = procNumTmp;

			exPara2.reset(new ExporterParaView<SC,LO,GO,NO>());
			if(dim==2)
                exPara2->setup("ElementDistribution2D", domain->getMesh(), FEType, parameterListAll);
			else
                exPara2->setup("ElementDistribution3D", domain->getMesh(), FEType, parameterListAll);
				
			exPara2->addVariable(procNum, "Procs", "Scalar", 1, domain->getMapUnique(), domain->getMapUniqueP2());
			exPara2->save(0.0); // 
			
			numElements[j]=domain->getElementMap()->getMaxAllGlobalIndex()+1;
			numElementsProc[j] = domain->getElementsC()->numberElements();

			double errorBreak;
			// Collect max Error in Elements
			reduceAll<int, double> (*comm, REDUCE_MAX, maxErrorProc, outArg (errorBreak));

			maxErrorEl[j] = errorBreak;

			if(dim == 2){
				if(( j==maxIter || errorBreak < tol ) && j>0)
					break;
			}
			if(dim==3){
				if( j==maxIter)
					break;
			}
						
			// Refinement
			domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
			auto startRef = std::chrono::high_resolution_clock::now();
			{
				domainRefined->refineMesh(domainP1RefinedArray,j, checkRestrictions, restriction); // always use the P1 domain, P2 Domain has has lost its' edge (höhö)
			}
			auto finishRef = std::chrono::high_resolution_clock::now();
			auto timeTmp = finishRef - startRef;
			meshTiming[j] = timeTmp;

			elapsed_MeshRef = elapsed_MeshRef + finishRef - startRef;

			domainP1RefinedArray.push_back(domainRefined);

			domainP1 = domainRefined;
			domain = domainP1;
			
			j++;
			
		}	
		auto finishTotal = std::chrono::high_resolution_clock::now();

		elapsed_total = finishTotal - startTotal;

		// Collect the number of elements each processor hold after refinement
		vec_GO_Type globalProcs(0);
		for (int i=0; i<= maxRank; i++)
				globalProcs.push_back(i);

		Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

		vec_GO_Type localProc(0);
		localProc.push_back(comm->getRank());
		Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

		MapPtr_Type mapGlobalProc =
			Teuchos::rcp( new Map_Type( domain->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, comm) );

		// Global IDs of Procs
		MapPtr_Type mapProc =
			Teuchos::rcp( new Map_Type( domain->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, comm) );
		
		MultiVectorPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVector_Type( mapProc, 1 ) );

		exportLocalEntry->putScalar( (LO) numElementsProc[j] );

		MultiVectorPtr_Type elementList= Teuchos::rcp( new MultiVector_Type( mapGlobalProc, 1 ) );
		elementList->putScalar( 0 ); 
		elementList->importFromVector( exportLocalEntry, true, "Insert");

		Teuchos::ArrayRCP<const double > elementProcList = elementList->getData(0);


		if(comm->getRank() == 0){	
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Summary Mesh Refinement" << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Marking Strategy:	" << strategy << endl;
			cout << " Theta:			" << theta << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Tolerance:			" << tol << endl;
			cout << " Max number of Iterations:	" <<  maxIter << endl;
			cout << " Number of Processors:		" << maxRank+1 << endl;
			cout << " Number of Refinements:		" << j << endl;
			cout << " Checking Restrictions:		" << checkRestrictions << endl;
			cout << " Restriction Type: 		" << restriction << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Number of elements after Refinement.... " << endl;
			for(int i=1; i<=j ; i++)
				cout <<" "<< i << ":	" << numElements[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Errorestimation: max error in Elements according to error Estimator after Refinement.... " << endl;
			for (int i=1; i<=j ; i++)
				cout <<" "<< i << ":	" << maxErrorEl[i] << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << "Distribution of elements on .. " << endl;
			for(int l=0; l<maxRank+1 ; l++)
				cout <<" Processor "<< l << " carries " << elementProcList[l] << " Elements "<< endl; 
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Elapsed Time of MeshRefinement Step.. " << endl;
			for (int i=1; i<=j ; i++)
				cout <<" " << i << ":	" << meshTiming[i-1].count() << " s" << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Elapsed time of Mesh Refinement, sum of refinement steps:	" << elapsed_MeshRef.count() << " s\n";
			cout << " Elapsed time solving PDE, sum of all loops :			" << elapsed_Solv.count() << " s\n";
			cout << " Elapsed time of Refinement & Solving:				" << elapsed_total.count() << " s\n";
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			}
		
		}
	 
 return(EXIT_SUCCESS);

}






