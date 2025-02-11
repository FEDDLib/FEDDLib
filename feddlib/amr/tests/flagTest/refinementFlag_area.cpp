#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshPartitioner.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/General/ExporterParaView.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/General/BCBuilder.hpp"
#include "feddlib/amr/AdaptiveMeshRefinement.hpp"

#include <Teuchos_GlobalMPISession.hpp>
#include <Xpetra_DefaultPlatform.hpp>

/*!
 Mesh Element Flags test

 @brief  Mesh Element Flags test
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


using namespace std;
using namespace Teuchos;

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef default_go GO;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType NO;
using namespace FEDD;
int main(int argc, char *argv[]) {

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef RCP<Domain_Type > DomainPtr_Type;
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;
    typedef MeshPartitioner<SC,LO,GO,NO> MeshPartitioner_Type;
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
    
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > comm = Xpetra::DefaultPlatform::getDefaultPlatform().getComm();

 	// Command Line Parameters
    Teuchos::CommandLineProcessor myCLP;
    string ulib_str = "Tpetra";
    myCLP.setOption("ulib",&ulib_str,"Underlying lib");
	string xmlProblemFile = "parametersProblem.xml";
    myCLP.setOption("problemfile",&xmlProblemFile,".xml file with Inputparameters.");

 	ParameterListPtr_Type parameterListProblem = Teuchos::getParametersFromXmlFile(xmlProblemFile);
   
    ParameterListPtr_Type parameterListAll(new Teuchos::ParameterList(*parameterListProblem)) ;


    // Mesh
    int dim = parameterListProblem->sublist("Parameter").get("Dimension",2);
    int level = parameterListProblem->sublist("Parameter").get("Level",2);
    int m = parameterListProblem->sublist("Parameter").get("H/h",5);
    std::string FEType = parameterListProblem->sublist("Parameter").get("Discretization","P2");
    std::string meshType = parameterListProblem->sublist("Parameter").get("Mesh Type","structured");
    std::string meshName = parameterListProblem->sublist("Parameter").get("Mesh Name","");
    std::string meshDelimiter = parameterListProblem->sublist("Parameter").get("Mesh Delimiter"," ");
    int volumeID= parameterListProblem->sublist("Parameter").get("VolumeID",99);

    int n;
    int size = comm->getSize();
    int numProcsCoarseSolve = parameterListProblem->sublist("General").get("Mpi Ranks Coarse",0);
    size -= numProcsCoarseSolve;
    
    bool boolExportMesh = true;


    DomainPtr_Type domainP1;
    DomainPtr_Type domain;

    
    ParameterListPtr_Type pListPartitioner = Teuchos::rcp( new ParameterList("Mesh Partitioner") );
    pListPartitioner->set( "Mesh 1 Name", meshName );
    
    //Paralleles P1 mesh (Partitionierung mit metis)
    domainP1.reset( new Domain_Type( comm, dim ) );
    MeshPartitioner_Type::DomainPtrArray_Type domainP1Array(1);
    domainP1Array[0] = domainP1;

    MeshPartitioner<SC,LO,GO,NO> partitionerP1 ( domainP1Array, pListPartitioner, "P1", dim );
    
    partitionerP1.readAndPartition(volumeID);


    domain.reset( new Domain_Type( comm, dim ) );
                    
    // Refining mesh uniformly
    // ----------
    AdaptiveMeshRefinement<SC,LO,GO,NO> meshRefiner(parameterListProblem); // exactSolLShape
    DomainPtr_Type domainRefined;
	domainRefined.reset( new Domain<SC,LO,GO,NO>( comm, dim ) );
	vec2D_dbl_Type area(3,vec_dbl_Type(2));
	area[0][0] = 0.;
	area[0][1] = 0.5;
	area[1][0] = 0.;
	area[1][1] = 0.5;
	area[2][0] = 0.;
	area[2][1] = 0.5;


	{
		domainRefined = meshRefiner.refineArea(domainP1,area,level );
	}
	domainP1=domainRefined;
	domain = domainRefined;
	domain->buildP2ofP1Domain( domainP1 );		
			
	// ----------
                   
                    
    MultiVectorPtr_Type mvFlag = rcp(new MultiVector_Type( domain->getMapUnique() ) );
	vec_int_ptr_Type BCFlags = domain->getBCFlagUnique();
    Teuchos::ArrayRCP< SC > values = mvFlag->getDataNonConst(0);
    TEUCHOS_TEST_FOR_EXCEPTION( values.size()!= BCFlags->size(), std::runtime_error, "MultiVector and node list have different sizes.");
    for (int i=0; i<values.size(); i++)
        values[i] = BCFlags->at(i);
    
    Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());
    
    exPara->setup( "NodeFlagsP2_area", domain->getMesh(), FEType );
    MultiVectorConstPtr_Type mvFlagConst = mvFlag;
    exPara->addVariable( mvFlagConst, "Flags", "Scalar", 1, domain->getMapUnique(), domain->getMapUnique());
    
    exPara->save(0.0);
    exPara->closeExporter();

    return(EXIT_SUCCESS);
}
