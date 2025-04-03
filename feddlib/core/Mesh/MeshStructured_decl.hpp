#ifndef MeshStructured_decl_hpp
#define MeshStructured_decl_hpp
#include "feddlib/core/Utils/FEDDUtils.hpp"
#include "Mesh.hpp"

/*!
 Declaration of MeshStructured
 
 @brief  MeshStructured
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
*/

/*! 
	Extension of mesh class. 
    In MeshStructured meshes are build inside the FEDDLib from scratch and distributed among the processors - structured partition. 
    
    !!! Please note that the P1-P2 node structure differs from the 'buildP2Mesh' process in unstructured meshes !!!
    (P2 nodes are generally build on top of P1 nodes and added to the P1 element list.)

    Possible FEType are in general P1 or P2. For most cases a Q2-P1 (disc) discretization is available.
    
*/
namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class MeshStructured : public Mesh<SC,LO,GO,NO> {
	
public:
    typedef Mesh<SC,LO,GO,NO> Mesh_Type;
    typedef typename Mesh_Type::CommPtr_Type CommPtr_Type;
    typedef typename Mesh_Type::CommConstPtrConst_Type CommConstPtrConst_Type;
    typedef typename Mesh_Type::Map_Type Map_Type;
    typedef typename Mesh_Type::MapPtr_Type MapPtr_Type;
    typedef typename Mesh_Type::ElementsPtr_Type ElementsPtr_Type;
    typedef Teuchos::RCP<MeshStructured<SC,LO,GO,NO> > MeshStrPtr_Type;
    /* ###################################################################### */
    //
    MeshStructured();

    MeshStructured(CommConstPtrConst_Type& comm);

    ~MeshStructured();
    
    void setGeometry2DRectangle(std::vector<double> coordinates, double l, double h);
    
    void setGeometry3DBox(std::vector<double> coordinates, double l, double w, double h);
    

    /// @brief Building a general 2D rectangular mesh with the lenght and height as defined per 'setGeomerty2DRectangle'. Called by Domain class and different discretizations for 2D mesh are called within this functions.
    /// @param FEType Finite element discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void buildMesh2D(std::string FEType,
                    int N,
                    int M,
                    int numProcsCoarseSolve=0,
                    std::string underlyingLib="Tpetra");

    /// @brief  Building 2D TPM rectangular mesh with the lenght and height as defined per 'setGeomerty2DRectangle' - characterized by building addition line segments for boundary conditions 
    /// @param FEType Finite element discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void buildMesh2DTPM(std::string FEType,
                        int N,
                        int M,
                        int numProcsCoarseSolve=0,
                        std::string underlyingLib="Tpetra");
    
    /// @brief Building general 3D cuboid with length, width and height as defines by 'setGeometry3DBox'. Called by Domain class and different discretizations for 3D mesh are called within this functions.
    /// @param FEType Finite element discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void buildMesh3D(std::string FEType,
                     int N,
                     int M,
                     int numProcsCoarseSolve=0,
                     std::string underlyingLib="Tpetra");

    /// @brief Building general 3D cuboid with length, width and height as defines by 'setGeometry3DBox'. Called by Domain class. Only P2 discretization available. Subcubes are build with 5 elements
    /// @param FEType Finite element discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void buildMesh3D5Elements(std::string FEType, int N, int M, int numProcsCoarseSolve, std::string underlyingLib="Tpetra");

    
    /// @brief Building 2D backward facing step geometry with the lenght and height as defined per 'setGeomerty2DRectangle'. Called by Domain class and different discretizations for 2D mesh are called within this functions.
    /// @param FEType Finite element discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void buildMesh2DBFS(std::string FEType,
                     int N,
                     int M,
                     int numProcsCoarseSolve=0,
                     std::string underlyingLib="Tpetra");

    /// @brief Building general 3D backward facing step with length, width and height as defines by 'setGeometry3DBox'. Called by Domain class and different discretizations for 3D mesh are called within this functions.
    /// @param FEType Finite element discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void buildMesh3DBFS(std::string FEType,
                        int N,
                        int M,
                        int numProcsCoarseSolve=0,
                        std::string underlyingLib="Tpetra");
    
    /// @brief Building general 3D backward facing step with length, width and height as defines by 'setGeometry3DBox' with Q2-P1 discontinous discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library  
    void buildP1_Disc_Q2_3DBFS(int N,
                            int M,
                            int numProcsCoarseSolve,
                            std::string underlyingLib);

    /// @brief Building general 3D cuboid with length, width and height as defines by 'setGeometry3DBox' with Q2-P1 discontinous discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void buildP1_Disc_Q2_3DCube(int N,
                               int M,
                               int numProcsCoarseSolve,
                               std::string underlyingLib);

    /// @brief Building general 3D cuboid with length, width and height as defines by 'setGeometry3DBox' with Q1 discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void build3DQ1Cube(int N,
                       int M,
                       int numProcsCoarseSolve,
                       std::string underlyingLib );

    /// @brief Building general 3D cuboid with length, width and height as defines by 'setGeometry3DBox' with Q2 discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
   void build3DQ2Cube( int N,
                        int M,
                        int numProcsCoarseSolve,
                        std::string underlyingLib );

    /// @brief Building general 3D cuboid with length, width and height as defines by 'setGeometry3DBox' with Q2 discretization and 20?
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void build3DQ2_20Cube(int N,
                          int M,
                          int numProcsCoarseSolve,
                          std::string underlyingLib ); 


    /// @brief Building general 3D backward facing step with length, width and height as defines by 'setGeometry3DBox' with Q2 discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
    void build3DQ2BFS( int N,
                       int M,
                       int numProcsCoarseSolve,
                       std::string underlyingLib );
    
    /// @brief  Building 2D mini TPM rectangular mesh with the lenght and height as defined per 'setGeomerty2DRectangle' - characterized by building addition line segments for boundary conditions 
    /// @param FEType Finite element discretization
    /// @param N Number of subdomains
    /// @param M H/h with H subdomain diameter (length/H) and h characteristic mesh size (length/(M*N))
    /// @param numProcsCoarseSolve if we want to reserve certain processors for coarse solve
    /// @param underlyingLib underlying linear algebra library 
   void buildMesh2DMiniTPM(std::string FEType,
                            int N,
                            int M,
                            int numProcsCoarseSolve=0,                            
                            std::string underlyingLib="Tpetra" );

    /// @brief Building suface lines for TPM square mini. Empty.
    /// @param feType 
    void buildSurfaceLinesSquareMiniTPM( std::string feType );
    
    void setRankRange( int numProcsCoarseSolve );
    
    void buildElementsClass( vec2D_int_ptr_Type elements, vec_int_ptr_Type elementFlag = Teuchos::null );
    
    /// @brief Building suface lines aka edges. Empty.
    void buildSurfaceLinesSquare();
    
    GO globalID_Q2_20Cube(int r, int s , int t, int &rr, int off_x, int off_y, int off_z, int M, int N,
                          GO nmbPoints_oneDirFull, GO nmbPoints_oneDirMid);
    
    /// @brief Setting corresponding flags to structured mesh depending on underlying problem. Rectangles/Cuboids are treated as 'fluid geometrys' with inflow and outflow surface. 
    /// @param flagsOption depending on underlying geometry (Rectangle/Cuboid - BFS - TPM)
    /// @param FEType Discretization
    void setStructuredMeshFlags(int flagsOption,std::string FEType="P1");
    
    /// @brief Building element map
    void buildElementMap();

    /// @brief Building surfaces. This is useful for structural problems. Each surface gets another flag. Corresponds to flag option 3.
    /// @param flagsOption corresponds to flag option from 'setStructureMeshFlags'. Only valid option is 3 for now.
    /// @param FEType 
    void buildSurfaces(int flagsOption, std::string FEType);

    void flipSurface(vec_int_Type &surfaceElements_vec);


    /* ###################################################################### */
    
    std::vector<double> coorRec;
    double 				length; // length of geometry
    double		 		height; // height of geometry
    double 				width; // width of geometry

    /* ###################################################################### */
private:
};
}
#endif
