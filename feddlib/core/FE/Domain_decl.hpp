#ifndef Domain_decl_hpp
#define Domain_decl_hpp

#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/Mesh/MeshStructured.hpp"
#include "feddlib/core/Mesh/MeshUnstructured.hpp"

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

  /*!
    \class Domain
    \brief This class defines the finite-element domain of your finite element problem for the respective variable/component.
    

    \tparam SC The scalar type. So far, this is always double, but having it as a template parameter would allow flexibily, e.g., for using complex instead
    \tparam LO The local ordinal type. The is the index type for local indices
    \tparam GO The global ordinal type. The is the index type for global indices
    @todo This should actually be removed since the class should operate only on element level)
    \tparam NO The Kokkos Node type. This would allow for performance portability when using Kokkos. Currently, this is not used.
    
    Example: If you construct a Stokes finite element problem, you get a velocity and pressure 'P2-P1' discretization and, thus, one domain for the P2 elements and one for the P1 elements, with the respective node list etc.
    
	*/

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

	 /*!
         \brief Constructor
        
     */
    Domain();

 	/*!
         \brief Constructor
         @param[in] comm
    */
    Domain(CommConstPtr_Type comm);
    
    /*!
         \brief Constructor
         @param[in] comm
         @param[in] dimension 
    */
    Domain(CommConstPtr_Type comm, int dimension);

    /*!
         \brief Constructor for 2D structured meshes build in FEDDLib
         @param[in] coor
         @param[in] l
         @param[in] h
         @param[in] comm 
    */
    Domain(vec_dbl_Type coor, double l, double h, CommConstPtr_Type comm);

 	/*!
         \brief Constructor for 3D strucutred meshes build in FEDDLib
         @param[in] coor
         @param[in] l length
         @param[in] w width
         @param[in] h height
         @param[in] comm 
    */
    Domain(vec_dbl_Type coor, double l, double w, double h, CommConstPtr_Type comm);
    
    /*!
         \brief initializeFEData
    */
    
    void initializeFEData();
    
    /*!
         \brief Returns flags of elements as vector of type int
         \return elementFlags
    */
    vec_int_ptr_Type getElementsFlag() const;

    /*!
         \brief Information about the domain
    */
    void info() const;

 	/*!
         \brief Build structured mesh in FEDDLib 
         @param[in] flags
         @param[in] meshType
         @param[in] dim 
         @param[in] FEType discretization
         @param[in] N
         @param[in] M
         @param[in] numProcsCoarseSolve
         
    */
    
    void buildMesh(int flags, std::string meshType, int dim, std::string FEType, int N, int M, int numProcsCoarseSolve = 0);
    
     /*!
         \brief Estimate depending on FEType and dimension for the numer of entries in system matrix, as approximate number of entries is requiered for matrix initialization.
         		This serves as a upper bound.
         \return approxEntriesPerRow
    */
    LO getApproxEntriesPerRow() const;

    /*!
         \brief Get dimension
         \return dimension

    */
    UN getDimension() const;

    /*!
         \brief  Get map of all uniquely distributed nodes of the processors. 
         \return mapUnique

    */
    MapConstPtr_Type getMapUnique() const;

    /*!
         \brief Get map of all repeated (not uniquely distributed) nodes of processors
         \return mepRepeated

    */
    MapConstPtr_Type getMapRepeated() const;

    /*!
         \brief Get map of elements
         \return elementMap

    */
    MapConstPtr_Type getElementMap() const;

    /*!
         \brief Get map of edges
         \return edgeMap

    */
    MapConstPtr_Type getEdgeMap() const; // edgeMap

    /*!
         \brief Get vector of repeated points (on your processor). 2D Vector with size: numPoints x dim
         \return pointsRepeated

    */
    vec2D_dbl_ptr_Type getPointsRepeated() const;

    /*!
         \brief Get vector of unique points (on your processor). 2D Vector with size: numPoints x dim
         \return pointsUnique

    */
    vec2D_dbl_ptr_Type getPointsUnique() const;

    /*!
         \brief Get vector of the flags corrsponding to the repeated points (on your processor)
         \return BCFlagsRepeated

    */
    vec_int_ptr_Type getBCFlagRepeated() const;

    /*!
         \brief Get vector of the flags corrsponding to the unique points (on your processor)
         \return BCFlagsUnique

    */
    vec_int_ptr_Type getBCFlagUnique() const;

    /*!
         \brief Get the elements (on your processor) as vector data type. 
         \return elements

    */
    vec2D_int_ptr_Type getElements() const;

    /*!
         \brief Get the elements (on your processor) as Elements_Type 
         \return elementsC

    */
    ElementsPtr_Type getElementsC() const;
    
    /*!
         \brief Get map of all unique (not uniquely distributed) nodes of processors in a vector field point of view. Here, a node has more than one degree of freedom. For Example in with dofs=2: 1_x , 1_y , 2_x , 2_y ,... . This degree of freedom map captures this. 
         \return mepVecFieldUnique

    */
    MapConstPtr_Type getMapVecFieldUnique() const;

    /*!
         brief Get map of all repeated (not uniquely distributed) nodes of processors in a vector field point of view. Here, a node has more than one degree of freedom. For Example in with dofs=2: 1_x , 1_y , 2_x , 2_y ,... . This degree of freedom map captures this. 
         \return mepVecFieldRepeated

    */
    MapConstPtr_Type getMapVecFieldRepeated() const;

    /*!
         \brief Finite element discretization. Represented as string, i.e., P2
         \return FEType

    */
    std::string getFEType() const;

    /*!
         \brief Communicator
         \return comm

    */
    CommConstPtr_Type getComm() const;

    /*!
         \brief Reading mesh from .mesh file with name 'filename'
         @param[in] filename
         @param[in] delimiter
         @param[in] dim
         @param[in] FEType
         @param[in] volumeID element flag 

    */
    void readMesh(std::string filename, std::string delimiter, int dim, std::string FEType, int volumeID=10);

    /*!
         \brief Reading mesh size
         @param[in] filename
         @param[in] delimiter

    */
    void readMeshSize(std::string filename, std::string delimiter);
    
    /*!
         \brief Partition mesh according to number of processors. Partition with parmetis.
         @param[in] partitionDistance

    */
    void partitionMesh( bool partitionDistance = false );

    /*!
         \brief Function called when mesh is read and paritioned. In turn it calls readMesh(...) and partitionMesh(...) 
         @param[in] filename
         @param[in] delimiter
         @param[in] dim
         @param[in] FEType
         @param[in] volumeID element flag 

    */
    void readAndPartitionMesh( std::string filename, std::string delimiter, int dim, std::string FEType, int volumeID=10 );

    /*!
         \brief Building a P2 mesh from the P1 one mesh. We generally habe P1-input meshes and build the P2 mesh on top of that.
         @param[in] domainP1 the domain with the mesh we read from .mesh file

    */
    void buildP2ofP1Domain( DomainPtr_Type domainP1 );

    /*!
         \brief Initialize domain with already existing domain.
         @param[in] domainP1 the domain with the mesh we read from .mesh file

    */
    void initWithDomain(DomainPtr_Type domainsP1); 
    
    /*!
         \brief Settng mesh from domain
         @param[in] meshUnstr mesh of MeshUnstr_Type which is generally the type of meshes from .mesh files

    */
    void setMesh(MeshUnstrPtr_Type meshUnstr); 

     /*!
         \brief Initialize dummy mesh for i.e. lagrange multiplier that only represents on 'point' in that sense. i.e. for setting pressure mean value in P2-P1 stokes problem. This is necassary if a variable does not live on a mesh. 
         \brief This function might not be necassary in the long run.
         @param[in] map for this dummy mesh. Maybe the mesh only represents one point (i.e. for lagrange mp). 

    */
    void initDummyMesh(MapPtr_Type map); 
    
    /*!
         \brief  Build unique node and dof interfaceMap in interface numbering
    */
    void buildUniqueInterfaceMaps();

    /*!
         \brief  Set partial coupling
         @param[in] flag
         @param[in] type
    */
    void setPartialCoupling(int flag, std::string type);

     /*!
         \brief  Build interface maps with global numbering with indicesGlobalMatchedUnique_
    */
    void buildGlobalInterfaceMaps();

     /*!
         \brief  Build interface maps 
    */
    void buildInterfaceMaps();
    
    /*!
         \brief Get interface map unique 
         \return interfaceMapUnique
    */
    MapConstPtr_Type getInterfaceMapUnique() const;

    /*!
         \brief Get interface vector field map unique
         \return interfaceMapVecFieldUnique
    */
    MapConstPtr_Type getInterfaceMapVecFieldUnique() const;
    
    /*!
         \brief Get global interface map vector field partial
         \return partialGlobalInterfaceVecFieldMap
    */    
    MapConstPtr_Type getGlobalInterfaceMapVecFieldPartial() const{ return partialGlobalInterfaceVecFieldMap_; };

    /*!
         \brief Get other global interface map vec field partial
         \return otherPartialGlobalInterfaceVecFieldMap
    */
    MapConstPtr_Type getOtherGlobalInterfaceMapVecFieldPartial() const{ return otherPartialGlobalInterfaceVecFieldMap_; };

    /*!
         \brief Get interface map unique (for fsi coupling block c4)
         \return globalInterfaceMapUnique
    */
    MapConstPtr_Type getGlobalInterfaceMapUnique() const;

    /*!
         \brief Get interface vec field map unique (for fsi coupling block c4)
         \return globalInterfaceMapUnique
    */
    MapConstPtr_Type getGlobalInterfaceMapVecFieldUnique() const; // Brauchen wir fuer Kopplungsblock C4

 /*!
         \brief Get other interface vec field map unique (for fsi coupling block c4)
         \return globalInterfaceMapUnique
    */   
    MapConstPtr_Type getOtherGlobalInterfaceMapVecFieldUnique() const; 

    /*!
         \brief Get global number of elements
         \return numElementsGlobal
    */
    GO getNumElementsGlobal() const;

    /*!
         \brief Get local number of elements (on your processor)
         \return numElementsGlobal
    */
    LO getNumElements() const;

	/*!
         \brief Get local number of points of type 'type' (unique/repeated)
         \return numPoints
    */
    LO getNumPoints(std::string type="Unique") const;/*local*/

	/*!
         \brief Checks geometriy
         @param[in] MeshType
         @param[in] dim
         \return 
    */
    int checkGeomentry(std::string MeshType, int dim) const;
    
    /*!
         \brief Itentify interface parallal and distance
         @param[in] domainOther
         @param[in] interfaceID_vec
    */
    void identifyInterfaceParallelAndDistance( DomainPtr_Type domainOther, vec_int_Type interfaceID_vec );
    
     /*!
         \brief Calculate distance to interface
    */
    void calculateDistancesToInterface();
        
    /*!
         \brief Get distances to interface
         \return 
    */
    vec_dbl_ptr_Type getDistancesToInterface() const;

    /*!
         \brief Partition distance to interface
    */
    void partitionDistanceToInterface();

    /*!
         \brief Set reference configuration
    */
    void setReferenceConfiguration();

    /*!
         \brief Move mesh according to displacement (i.e. used in FSI)
         @param[in] displacementUnique
         @param[in] displacementRepeated
    */
    void moveMesh(MultiVectorPtr_Type displacementUnique, MultiVectorPtr_Type displacementRepeated);

    /*!
         \brief Get mesh
         \return mesh
    */
    MeshPtr_Type getMesh();

    /*!
         \brief Get mesh constant
         \return mesh
    */
    MeshConstPtr_Type getMesh() const;

    /*!
         \brief Generally the domain object holds only meshes from type 'Mesh'. If we read a mesh from file it becomes the type 'MeshUnstructured' and needs to be initialized
         @param[in] dimension
         @param[in] FEType
         @param[in] volumeID       
    */
    void initializeUnstructuredMesh(int dimension, std::string feType, int volumeID=10);

    /*!
		 \brief Hilfsfunktion fuer buildLocalInterfaceIDInGlobal().
		 Gibt fuer eine gegebene nodeID die entsprechende dofID und umgekehrt.
		 localDofNumber entspricht dem Rest der Division, also ob es sich um die
		 x- (=0), y- (=1) oder z-Komponente (=2) handelt.
    */
    void toNodeID(UN dim, GO dofID, GO &nodeID, LO &localDofNumber);

 	/*!
		 \brief Hilfsfunktion fuer buildLocalInterfaceIDInGlobal().
		 Gibt fuer eine gegebene nodeID die entsprechende dofID und umgekehrt.
		 localDofNumber entspricht dem Rest der Division, also ob es sich um die
		 x- (=0), y- (=1) oder z-Komponente (=2) handelt.
    */
    void toDofID(UN dim, GO nodeID, LO localDofNumber, GO &dofID );


    /*!
         \brief Get local interface id in global
         \return globalInterfaceID      
    */
    vec_long_Type getLocalInterfaceIDInGlobal() const;

    /*!
         \brief Set dummy interface domain
         @param[in] domain     
    */
    void setDummyInterfaceDomain(DomainPtr_Type domain);
    
    /*!
         \brief Find in points unique
         @param[in] x point to be found in unique nodes
         
    */
    int findInPointsUnique(const vec_dbl_Type& x) const;

    /*!
         \brief Get node list in multi vector point
         @param[in] nodeListMV
         
    */
    MultiVectorPtr_Type getNodeListMV() const;

    /*!
         \brief Exporting Mesh
         
    */
   void exportMesh(bool exportEdges = false, bool exportSurfaces=false, std::string exportMesh="export.mesh");

     /*!
         \brief Option of preprocessing mesh by making consistent outward normal and/or consistent element orientation, where we always have positive det of transformation to reference element
         @param correctSurfaceNormals bool for normal direction
         @param correctElementDirection bool for surface direction
    */ 
   void preProcessMesh(bool correctSurfaceNormals, bool correctElementDirection);

   /// @brief Exporting Paraview file displaying element flags of the underlying mesh
   /// @param name
   void exportElementFlags(std::string name = "default");

   /// @brief Exporting Paraview file displaying node flags of the underlying mesh
   /// @param name export suffix to identify flags
   void exportNodeFlags(std::string name = "default");

   /// @brief Exporting Paraview file displaying surface normals of the underlying mesh. As we are generally not able to plot only the surfaces, the normals are displayed in each node. This means, that at corners, the visualization is incorrect (i.e. node belongs to surfaces which are in different directions)
   /// @param name export suffix to identify flags
   void exportSurfaceNormals(std::string name = "default");

   /// @brief Exporting Paraview file displaying element volume of underlying mesh. 
   /// @param name export suffix to identify flags
   void exportElementOrientation(std::string name = "default");

   /* ----------------------------------------------------------------------------------------*/

   private:
   CommConstPtr_Type comm_; // underlying comm
   MeshPtr_Type mesh_;      // underlying mesh as base class mesh type. usually underlying mesh is either structured or unstructured
   int dim_;                // dimension
   vec_dbl_Type coorRec;
   double length;
   double height;
   double width;
   int n_;
   int m_;
   std::string FEType_; // Finite element discretization
   mutable MapPtr_Type mapVecFieldUnique_;
   mutable MapPtr_Type mapVecFieldRepeated_;

   string_vec_ptr_Type geometries2DVec_; // list with available 2D structured geometries
   string_vec_ptr_Type geometries3DVec_; // list with available 3D structured geometries
   vec_dbl_ptr_Type distancesToInterface_;

   // Unique Interface-Maps als nodes und als dofs in der Interface-Nummerierung
   MapPtr_Type interfaceMapUnique_;         // nodes
   MapPtr_Type interfaceMapVecFieldUnique_; // dofs

   // Unique Fluid/Struktur-Interface-Maps als nodes und als dofs in der globalen Nummerierung
   MapPtr_Type globalInterfaceMapUnique_;
   MapPtr_Type globalInterfaceMapVecFieldUnique_;
   MapPtr_Type partialGlobalInterfaceVecFieldMap_;
   MapPtr_Type otherGlobalInterfaceMapUnique_;
   MapPtr_Type otherGlobalInterfaceMapVecFieldUnique_;
   MapPtr_Type otherPartialGlobalInterfaceVecFieldMap_;
   // Dies ist ein (unique) partitionierter Vektor, der fuer jeden Stelle im Vektor i
   // (= lokale Interface ID; also jeder Proz. haelt nur ein Teil des Interfaces)
   // angibt, welche lokale ID dies in der globalen Nummerierung ist.
   // Beide zeigen auf denselben physischen Knoten des Interfaces!
   // TODO Fehlerhaft
   vec_long_Type vecLocalInterfaceIDinGlobal_;

   std::string meshType_;
   int numProcsCoarseSolve_;
   int flagsOption_;

    };
   
}

#endif
