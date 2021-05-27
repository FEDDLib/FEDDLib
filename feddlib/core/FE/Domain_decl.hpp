#ifndef Domain_decl_hpp
#define Domain_decl_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructuredRefinement.hpp"

/*!
 Declaration of Domain

 @brief  Domain
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class Domain {

public:
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;
    typedef std::vector<DomainPtr_Type> DomainPtrArray_Type; // Array of domains for meshrefinement
    typedef Teuchos::RCP< const Domain<SC,LO,GO,NO> > DomainConstPtr_Type;
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef Teuchos::RCP<Mesh_Type > MeshPtr_Type;

    typedef Teuchos::RCP<const Mesh_Type > MeshConstPtr_Type;
    typedef MeshStructured<SC,LO,GO,NO> MeshStr_Type;
    typedef Teuchos::RCP<MeshStr_Type> MeshStrPtr_Type;

    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef Teuchos::RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    typedef std::vector<MeshUnstrPtr_Type> MeshUnstrPtrArray_Type; // Array of meshUnstr for meshRefinement

    typedef MeshUnstructuredRefinement<SC,LO,GO,NO> MeshUnstrRef_Type;
    typedef Teuchos::RCP<MeshUnstrRef_Type> MeshUnstrRefPtr_Type;
    typedef std::vector<MeshUnstrRefPtr_Type> MeshUnstrRefPtrArray_Type; // Array of meshUnstr for meshRefinement
    
    typedef typename MeshUnstr_Type::MeshInterfacePtr_Type MeshInterfacePtr_Type;
    
    typedef typename Mesh_Type::Elements_Type Elements_Type;
    typedef typename Mesh_Type::ElementsPtr_Type ElementsPtr_Type;
            
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorPtrConst_Type;


    typedef Map<LO,GO,NO> Map_Type;
    typedef Teuchos::RCP<Map_Type> MapPtr_Type;
    typedef Teuchos::RCP<const Map_Type> MapConstPtr_Type;

    typedef Teuchos::Comm<int> Comm_Type;
    typedef Teuchos::RCP<const Comm_Type> CommConstPtr_Type;

    typedef std::vector<GO> vec_GO_Type;
    typedef std::vector<vec_GO_Type> vec2D_GO_Type;
    typedef std::vector<vec2D_GO_Type> vec3D_GO_Type;
    typedef Teuchos::RCP<vec3D_GO_Type> vec3D_GO_ptr_Type;
    
    /* ------------------------------------------------------------------------ */

    Domain();

    Domain(CommConstPtr_Type comm);
    
    Domain(CommConstPtr_Type comm, int dimension);

    Domain(vec_dbl_Type coor, double l, double h, CommConstPtr_Type comm);

    Domain(vec_dbl_Type coor, double l, double w, double h, CommConstPtr_Type comm);
    
    void initializeFEData();
    
    vec_int_ptr_Type getElementsFlag() const;

    void info();

    void buildMesh(int flags, std::string meshType, int dim, std::string FEType, int N, int M, int numProcsCoarseSolve = 0);
    
    LO getApproxEntriesPerRow() const;

    void setMeshParameterList( ParameterListPtr_Type& pl );
    
    UN getDimension() const;

    MapConstPtr_Type getMapUnique() const;

    MapConstPtr_Type getMapRepeated() const;

    MapConstPtr_Type getMapUniqueP2() const;

    MapConstPtr_Type getMapRepeatedP2() const;
    
    MapConstPtr_Type getElementMap() const;

    MapConstPtr_Type getEdgeMap() const; // edgeMap

    vec2D_dbl_ptr_Type getPointsRepeated() const;

    vec2D_dbl_ptr_Type getPointsUnique() const;

    vec_int_ptr_Type getBCFlagRepeated() const;

    vec_int_ptr_Type getBCFlagUnique() const;

    vec2D_int_ptr_Type getElements() const;

    ElementsPtr_Type getElementsC() const;
    
    MapConstPtr_Type getMapVecFieldUnique() const;

    MapConstPtr_Type getMapVecFieldRepeated() const;

    std::string getFEType() const;

    CommConstPtr_Type getComm() const;

    void readMesh(string filename, string delimiter, int dim, string FEType, int volumeID=10);

    void readMeshSize(string filename, string delimiter);
    
    void partitionMesh( bool partitionDistance = false );

    void readAndPartitionMesh( std::string filename, std::string delimiter, int dim, std::string FEType, int volumeID=10 );

    void buildP2ofP1Domain( DomainPtr_Type domainP1 );

    void refineMesh( DomainPtrArray_Type domainP1, int j, bool checkRestrictions, string restriction); // Mesh Refinement

	vec_dbl_Type errorEstimation(MultiVectorPtrConst_Type valuesSolution, double theta, string strategy);

	void initMeshRef( DomainPtr_Type domainP1 );
    
    // Baue unique node- und dof-InterfaceMap in der Interface-Nummerierung
    void buildUniqueInterfaceMaps();

    void setPartialCoupling(int flag, std::string type);
    
    // Baue die Interface-Maps in der globalen Nummerierung.
    // Dies machen wir mit Hilfe des unique partitionierten Vektors indicesGlobalMatchedUnique_.
    void buildGlobalInterfaceMaps();

    void buildInterfaceMaps();
    
    MapConstPtr_Type getInterfaceMapUnique() const;

    MapConstPtr_Type getInterfaceMapVecFieldUnique() const;
    
    MapConstPtr_Type getGlobalInterfaceMapVecFieldPartial() const{ return partialGlobalInterfaceVecFieldMap_; };

    MapConstPtr_Type getOtherGlobalInterfaceMapVecFieldPartial() const{ return otherPartialGlobalInterfaceVecFieldMap_; };

    MapConstPtr_Type getGlobalInterfaceMapUnique() const; // Brauchen wir fuer Kopplungsblock C4

    MapConstPtr_Type getGlobalInterfaceMapVecFieldUnique() const; // Brauchen wir fuer Kopplungsblock C4
    
    MapConstPtr_Type getOtherGlobalInterfaceMapVecFieldUnique() const; 

    GO getNumElementsGlobal() const;

    LO getNumElements() const;

    LO getNumPoints(std::string type="Unique") const;/*local*/

    int checkGeomentry(std::string MeshType, int dim) const;
    
    void identifyInterfaceParallelAndDistance( DomainPtr_Type domainOther, vec_int_Type interfaceID_vec );
    
    void calculateDistancesToInterface();
        
    vec_dbl_ptr_Type getDistancesToInterface() const;

    void partitionDistanceToInterface();

    void setReferenceConfiguration();

    void moveMesh(MultiVectorPtr_Type displacementUnique, MultiVectorPtr_Type displacementRepeated);

    MeshPtr_Type getMesh();

    MeshConstPtr_Type getMesh() const;

    void initializeUnstructuredMesh(int dimension, string feType, int volumeID=10);

    // Hilfsfunktionen fuer buildLocalInterfaceIDInGlobal().
    // Gibt fuer eine gegebene nodeID die entsprechende dofID und umgekehrt.
    // localDofNumber entspricht dem Rest der Division, also ob es sich um die
    // x- (=0), y- (=1) oder z-Komponente (=2) handelt.
    void toNodeID(UN dim, GO dofID, GO &nodeID, LO &localDofNumber);

    void toDofID(UN dim, GO nodeID, LO localDofNumber, GO &dofID );

    vec_long_Type getLocalInterfaceIDInGlobal() const;

    void setDummyInterfaceDomain(DomainPtr_Type domain);
    
    int findInPointsUnique(const vec_dbl_Type& x) const;

    MultiVectorPtr_Type getNodeListMV() const;
/* ----------------------------------------------------------------------------------------*/

private:

    CommConstPtr_Type 		comm_;
    MeshPtr_Type 			mesh_;
    int                     dim_;
    vec_dbl_Type            coorRec;
    double 					length;
    double		 			height;
    double 					width;
    int 					n_;
    int 					m_;
    std::string				FEType_;
    mutable MapPtr_Type mapVecFieldUnique_;
    mutable  MapPtr_Type mapVecFieldRepeated_;
    
    string_vec_ptr_Type     geometries2DVec_;
    string_vec_ptr_Type		geometries3DVec_;
    vec_dbl_ptr_Type        distancesToInterface_;

    // Unique Interface-Maps als nodes und als dofs in der Interface-Nummerierung
    MapPtr_Type             interfaceMapUnique_; // nodes
    MapPtr_Type             interfaceMapVecFieldUnique_; // dofs

    // Unique Fluid/Struktur-Interface-Maps als nodes und als dofs in der globalen Nummerierung
    MapPtr_Type             globalInterfaceMapUnique_;
    MapPtr_Type             globalInterfaceMapVecFieldUnique_;
    MapPtr_Type             partialGlobalInterfaceVecFieldMap_;
    MapPtr_Type otherGlobalInterfaceMapUnique_;
    MapPtr_Type otherGlobalInterfaceMapVecFieldUnique_;
    MapPtr_Type otherPartialGlobalInterfaceVecFieldMap_;
    // Dies ist ein (unique) partitionierter Vektor, der fuer jeden Stelle im Vektor i
    // (= lokale Interface ID; also jeder Proz. haelt nur ein Teil des Interfaces)
    // angibt, welche lokale ID dies in der globalen Nummerierung ist.
    // Beide zeigen auf denselben physischen Knoten des Interfaces!
    // TODO Fehlerhaft
    vec_long_Type           vecLocalInterfaceIDinGlobal_;

    std::string meshType_;
    int numProcsCoarseSolve_;
    int flagsOption_;

    };
}

#endif
