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
    
    virtual void dummy() {};
    
    void setGeometry2DRectangle(std::vector<double> coordinates, double l, double h);
    
    void setGeometry3DBox(std::vector<double> coordinates, double l, double w, double h);
    
    void buildMesh2D(std::string FEType,
                    int N,
                    int M,
                    int numProcsCoarseSolve=0,
                    std::string underlyingLib="Tpetra");

    void buildMesh2DTPM(std::string FEType,
                        int N,
                        int M,
                        int numProcsCoarseSolve=0,
                        std::string underlyingLib="Tpetra");
    
    void buildMesh3D(std::string FEType,
                     int N,
                     int M,
                     int numProcsCoarseSolve=0,
                     std::string underlyingLib="Tpetra");
    
    void buildMesh2DBFS(std::string FEType,
                     int N,
                     int M,
                     int numProcsCoarseSolve=0,
                     std::string underlyingLib="Tpetra");

    void buildMesh3DBFS(std::string FEType,
                        int N,
                        int M,
                        int numProcsCoarseSolve=0,
                        std::string underlyingLib="Tpetra");
    
    void buildP1_Disc_Q2_3DBFS(int N,
                            int M,
                            int numProcsCoarseSolve,
                            std::string underlyingLib);

    void buildP1_Disc_Q2_3DCube(int N,
                               int M,
                               int numProcsCoarseSolve,
                               std::string underlyingLib);

    void build3DQ1Cube(int N,
                       int M,
                       int numProcsCoarseSolve,
                       std::string underlyingLib );

    void build3DQ2Cube( int N,
                        int M,
                        int numProcsCoarseSolve,
                        std::string underlyingLib );

    void build3DQ2_20Cube(int N,
                          int M,
                          int numProcsCoarseSolve,
                          std::string underlyingLib ); 


    void build3DQ2BFS( int N,
                       int M,
                       int numProcsCoarseSolve,
                       std::string underlyingLib );
    
    void buildMesh2DMiniTPM(std::string FEType,
                            int N,
                            int M,
                            int numProcsCoarseSolve=0,                            
                            std::string underlyingLib="Tpetra" );

    void buildSurfaceLinesSquareMiniTPM( string feType );
    
    void setRankRange( int numProcsCoarseSolve );
    
    void buildElementsClass( vec2D_int_ptr_Type elements, vec_int_ptr_Type elementFlag = Teuchos::null );
    
    void buildSurfaceLinesSquare();
    
    GO globalID_Q2_20Cube(int r, int s , int t, int &rr, int off_x, int off_y, int off_z, int M, int N,
                          GO nmbPoints_oneDirFull, GO nmbPoints_oneDirMid);
    
    void setStructuredMeshFlags(int flagsOption,string FEType="P1");
    
    void buildElementMap();
    /* ###################################################################### */
    
    std::vector<double> coorRec;
    double 				length;
    double		 		height;
    double 				width;

    /* ###################################################################### */
private:    

    };

}
#endif
