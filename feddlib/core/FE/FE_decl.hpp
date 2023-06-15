#ifndef FE_DECL_hpp
#define FE_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/SmallMatrix.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "Domain.hpp"
#include "Helper.hpp"
#include "sms.hpp"
#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFENavierStokes_decl.hpp"

#include "feddlib/core/AceFemAssembly/AssembleFEFactory.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_BLAS.hpp>

#include <boost/function.hpp>

/*!
 Declaration of FE

 @brief  FE
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
class DataElement {
    public:
        DataElement();
        DataElement(int size);
        std::vector<double> getHp();
        std::vector<double> getHt();
        void setHp( double* ht);
    private:
        std::vector<double> ht_;
        std::vector<double> hp_;
};

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class FE {
  public:

    enum VarType {Std=0,Grad=1};

    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef Teuchos::RCP<Domain_Type> DomainPtr_Type;
    typedef Teuchos::RCP<const Domain_Type> DomainConstPtr_Type;
    typedef std::vector<DomainConstPtr_Type> DomainConstPtr_vec_Type;

    typedef Teuchos::RCP<Mesh<SC,LO,GO,NO> > MeshPtr_Type;
    typedef MeshUnstructured<SC,LO,GO,NO> MeshUnstr_Type;
    typedef Teuchos::RCP<MeshUnstr_Type> MeshUnstrPtr_Type;
    
    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
    typedef Teuchos::RCP<const Elements_Type> ElementsConstPtr_Type;
    
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

    typedef typename Matrix_Type::MapPtr_Type MapPtr_Type;
    typedef typename Matrix_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;

    typedef std::vector<GO> vec_GO_Type;
    typedef std::vector<vec_GO_Type> vec2D_GO_Type;
    typedef std::vector<vec2D_GO_Type> vec3D_GO_Type;
    typedef Teuchos::RCP<vec3D_GO_Type> vec3D_GO_ptr_Type;

    typedef boost::function<void(double* x, double* res, double t, const double* parameters)> BC_func_Type;
    
	typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;
    typedef Teuchos::RCP<AssembleFE_Type> AssembleFEPtr_Type;

	typedef AssembleFENavierStokes<SC,LO,GO,NO> AssembleFENavierStokes_Type;
    typedef Teuchos::RCP<AssembleFENavierStokes_Type> AssembleFENavierStokesPtr_Type;

    typedef AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO> AssembleFE_SCI_SMC_Active_Growth_Reorientation_Type;
    typedef Teuchos::RCP<AssembleFE_SCI_SMC_Active_Growth_Reorientation_Type> AssembleFE_SCI_SMC_Active_Growth_Reorientation_Ptr_Type;

    typedef std::vector<AssembleFEPtr_Type> AssembleFEPtr_vec_Type;	

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type ;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type ;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

	typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

    DomainConstPtr_vec_Type	domainVec_;
    Teuchos::RCP<ElementSpec> es_;

    /* ###################################################################### */

    FE(bool saveAssembly=false);
    
    void assemblyIdentity(MatrixPtr_Type &A);
    
    void assemblySurfaceIntegral(int dim,
                                 std::string FEType,
                                 MultiVectorPtr_Type  a,
                                 std::string fieldType,
                                 RhsFunc_Type func,
                                 std::vector<SC>& funcParameter);
                    
    void assemblySurfaceIntegralExternal(int dim,
                                    std::string FEType,
                                    MultiVectorPtr_Type f,
                                    MultiVectorPtr_Type d_rep,
                                    std::vector<SC>& funcParameter,
                                    RhsFunc_Type func,
                                    ParameterListPtr_Type params,
                                    int FEloc=0);

    void assemblyNonlinearSurfaceIntegralExternal(int dim,
                                    std::string FEType,
                                    MultiVectorPtr_Type f,
                                    MultiVectorPtr_Type d_rep,
                                    MatrixPtr_Type &Kext,
                                    std::vector<SC>& funcParameter,
                                    RhsFunc_Type func,
                                    ParameterListPtr_Type params,
                                    int FEloc = 0);
    
    void assemblySurfaceIntegralFlag(int dim,
                                    std::string FEType,
                                     MultiVectorPtr_Type  a,
                                     std::string fieldType,
                                     BC_func_Type func,
                                     std::vector<SC>& funcParameter);
    
    void assemblyAceGenTPM( MatrixPtr_Type &A00,
                            MatrixPtr_Type &A01,
                            MatrixPtr_Type &A10,
                            MatrixPtr_Type &A11,
                            MultiVectorPtr_Type &F0,
                            MultiVectorPtr_Type &F1,
                            MapPtr_Type &mapRepeated1,
                            MapPtr_Type &mapRepeated2,
                            ParameterListPtr_Type parameterList,
                            MultiVectorPtr_Type u_repeatedNewton=Teuchos::null,
                            MultiVectorPtr_Type p_repeatedNewton=Teuchos::null,
                            MultiVectorPtr_Type u_repeatedTime=Teuchos::null,
                            MultiVectorPtr_Type p_repeatedTime=Teuchos::null,
                           bool update=true,
                           bool updateHistory=true);
    
    
    void addFE(DomainConstPtr_Type domain);

    void doSetZeros(double eps = 10*Teuchos::ScalarTraits<SC>::eps());

    void assemblyEmptyMatrix(MatrixPtr_Type &A);

    void applyBTinv(vec3D_dbl_ptr_Type& dPhiIn,
                    vec3D_dbl_Type& dPhiOut,
                    SmallMatrix<SC>& Binv);

    void applyBTinv(vec3D_dbl_ptr_Type& dPhiIn,
                    vec3D_dbl_Type& dPhiOut,
                    const SmallMatrix<SC>& Binv);
    
    void assemblyLaplace(int Dimension,
                        std::string FEType,
                        int degree,
                        MatrixPtr_Type &A,
                        bool callFillComplete = true,
                         int FELocExternal = -1);

    void assemblyMass(int dim,
                      std::string FEType,
                      std::string fieldType,
                      MatrixPtr_Type &A,
                      bool callFillComplete = true);

    // Ueberladung der Assemblierung der Massematrix fuer FSI, da
    // checkFE sonst auch fuer das Strukturproblem FEloc = 1 liefert (= Fluid)
    // und somit die welche domain und Map in der Assemblierung genutzt wird.
    void assemblyMass(int dim,
                      std::string FEType,
                      std::string fieldType,
                      MatrixPtr_Type &A,
                      int FEloc,
                      bool callFillComplete = true);

    void assemblyLaplaceVecField(int dim,
                                 std::string FEType,
                                 int degree,
                                 MatrixPtr_Type &A,
                                 bool callFillComplete = true);

    void assemblyLaplaceVecFieldV2(int dim,
                                 std::string FEType,
                                 int degree,
                                 MatrixPtr_Type &A,
                                 bool callFillComplete = true);

    void assemblyNonlinearLaplace(
        int dim, std::string FEType, int degree, MultiVectorPtr_Type u_rep,
        BlockMatrixPtr_Type &A, BlockMultiVectorPtr_Type &resVec,
        ParameterListPtr_Type params, string assembleMode,
        bool callFillComplete = true, int FELocExternal = -1);

    // Assembling the reaction term of the reaction diffusion equation. Maybe add default function.
	void assemblyLinearReactionTerm(int dim,
    							std::string FEType,
                                MatrixPtr_Type &A,
                                bool callFillComplete,
                     			std::vector<SC>& funcParameter,
								RhsFunc_Type reactionFunc);	

	// Assembling the reaction term of the reaction diffusion equation. Maybe add default function.
	void assemblyReactionTerm(int dim,
    							std::string FEType,
                                MatrixPtr_Type &A,
                                MultiVectorPtr_Type u,
                                bool callFillComplete,
                     			std::vector<SC>& funcParameter,
								RhsFunc_Type reactionFunc);	

                                // Assembling the reaction term of the reaction diffusion equation. Maybe add default function.
	void assemblyDReactionTerm(int dim,
    							std::string FEType,
                                MatrixPtr_Type &A,
                                MultiVectorPtr_Type u,
                                bool callFillComplete,
                     			std::vector<SC>& funcParameter,
								RhsFunc_Type reactionFunc);

    void assemblyLinElasXDimE(int dim,
                                std::string FEType,
                                MatrixPtr_Type &A,
                                MultiVectorPtr_Type eModVec,
                                double nu,
                                bool callFillComplete=true);

    void determineEMod(std::string FEType, 
                       MultiVectorPtr_Type solution,
                       MultiVectorPtr_Type &eModVec,
                       DomainConstPtr_Type domain,
                       ParameterListPtr_Type params);
    void assemblyLaplaceDiffusion(int Dimension,
                        std::string FEType,
                        int degree,
                        MatrixPtr_Type &A,
		            	vec2D_dbl_Type diffusionTensor,
                        bool callFillComplete = true,
                        int FELocExternal = -1);

    void assemblyElasticityJacobianAndStressAceFEM(int dim,
                                                   std::string FEType,
                                                   MatrixPtr_Type &A,
                                                   MultiVectorPtr_Type &f,
                                                   MultiVectorPtr_Type u,
                                                   ParameterListPtr_Type pList,
                                                   double C,
                                                   bool callFillComplete=true);
    
    void assemblyElasticityJacobianAceFEM(int dim,
                                          std::string FEType,
                                          MatrixPtr_Type &A,
                                          MultiVectorPtr_Type u,
                                          std::string material_model,
                                          double E,
                                          double nu,
                                          double C,
                                          bool callFillComplete=true);

    void assemblyElasticityStressesAceFEM(int dim,
                                         std::string FEType,
                                         MultiVectorPtr_Type &f,
                                         MultiVectorPtr_Type u,
                                         std::string material_model,
                                         double E,
                                         double nu,
                                         double C,
                                         bool callFillComplete=true);

    void assemblyAdvectionVecField(int dim,
                                   std::string FEType,
                                   MatrixPtr_Type &A,
                                   MultiVectorPtr_Type u,
                                   bool callFillComplete);

    void assemblyAdvectionInUVecField(int dim,
                                      std::string FEType,
                                      MatrixPtr_Type &A,
                                      MultiVectorPtr_Type u,
                                      bool callFillComplete);

    void assemblyDivAndDivT( int dim,
                            std::string FEType1,
                            std::string FEType2,
                            int degree,
                            MatrixPtr_Type &Bmat,
                            MatrixPtr_Type &BTmat,
                            MapConstPtr_Type map1,
                            MapConstPtr_Type map2,
                            bool callFillComplete = true);

    void assemblyDivAndDivTFast( int dim,
                                std::string FEType1,
                                std::string FEType2,
                                int degree,
                                MatrixPtr_Type &Bmat,
                                MatrixPtr_Type &BTmat,
                                MapConstPtr_Type map1,
                                MapConstPtr_Type map2,
                                bool callFillComplete = true);

    
    /*! Bochev-Dohrmann stabilization for P1-P1 finite elements. Must be scaled with 1/nu for general Navier-Stokes problem. */
    void assemblyBDStabilization(int dim,
                                 std::string FEType,
                                 MatrixPtr_Type &A,
                                 bool callFillComplete = true);

    /*! Fuer diskret harmonische Fortsetzung mit heuristischer Skalierung mit Hilfe von DistancesToInterface_
     Aehnelt also nicht direkt dem AssemblyLaplace (mit CoeffFunc) von oben */
    void assemblyLaplaceXDim(int dim,
                            std::string FEType,
                            MatrixPtr_Type &A,
                            CoeffFuncDbl_Type func,
                            double* parameters,
                            bool callFillComplete = true);

    //! (\grad u + (\grad u)^T, \grad v ); symmetrischer Gradient, wenn func = 1.0
    void assemblyStress(int dim,
                       std::string FEType,
                       MatrixPtr_Type &A,
                       CoeffFunc_Type func,
                       int* parameters,
                       bool callFillComplete = true);

    /*! 2*\mu*(\eps(u):\eps(v)) + \lambda*tr(\eps(u))*tr(\eps(v)), wobei
     tr(\eps(u)) = div(u) */
    void assemblyLinElasXDim(int dim,
                            std::string FEType,
                            MatrixPtr_Type &A,
                            double lambda,
                            double mu,
                            bool callFillComplete = true);

    /*! Dieser Term entsteht durch schwache Formulierung der ALE-Zeitableitung
     und bleibt in nicht-conservativer Form vorhanden.
     In conservativer Form "verschwindet" der Term in die conservative Form der Konvektion.
     Das betrachten wir aber nicht.
     Der Term lautet: (\grad \cdot w) u \cdot v
     und wir bei der endgueltigen Assemblierung subtrahiert. Also -(\grad \cdot w) u \cdot v */
    void assemblyAdditionalConvection(int dim,
                                      std::string FEType,
                                      MatrixPtr_Type &A,
                                      MultiVectorPtr_Type w,
                                      bool callFillComplete = true);

    // Das ist Assemblierungsroutine fuer die Kopplungsbloecke C1, C2 und C3 und deren Transponierte.
    // Die entsprechende Skalierung wird hinterher in dem Problem gemacht. Hier wird die entsprechende
    // Einheitsmatrix (auf den Interfaceknoten) gesetzt.
    void assemblyFSICoupling(int dim,
                             std::string FEType,
                             MatrixPtr_Type &C,
                             MatrixPtr_Type &C_T,
                             int FEloc1, // 0 = Fluid, 2 = Struktur
                             int FEloc2,
                             MapConstPtr_Type map1, // DomainMap: InterfaceMapVecFieldUnique = Spalten von C_T
                             MapConstPtr_Type map2,  // RangeMap: this->getDomain(0)->getMapVecFieldUnique() = Zeilen von C_T
                             bool callFillComplete = true);

    void assemblyDummyCoupling( int dim,
                                std::string FEType,
                                MatrixPtr_Type &C,
                                int FEloc, // 0 = Fluid, 2 = Struktur
                                bool callFillComplete);

    
    // Das ist Assemblierungsroutine fuer den Kopplungsblock C4, der
    // Geometrie und Struktur koppelt. Ist nur bei geometrisch implizit von Noeten.
    void assemblyGeometryCoupling(int dim,
                                  std::string FEType,
                                  MatrixPtr_Type &C,
                                  int FEloc, // wie oben. Allerdings 4 = Geometrie (= 0)
                                  MapConstPtr_Type map1, // Fluid-Interface-Map
                                  MapConstPtr_Type map2, // DomainMap: this->getDomain(2)->getMapVecFieldUnique() = Spalten von C
                                  MapConstPtr_Type map3, // RangeMap: this->getDomain(4)->getMapVecFieldUnique() = Zeilen von C
                                  bool callFillComplete = true);

    // Shape-Derivatives fuer die Geschwindigkeit; die Hauptvariable ist \delta d und Testfunktion ist v:
    // Spannung 1: DK1 = - \mu \int [ \grad(u) \grad(\delta d) + [\grad(\delta d)]^T [\grad(u)]^T ] : \grad(v) d\Omega
    // Spannung 2: DK2 = \int [ \sigma \eta ] : \grad(v) d\Omega
    // mit \sigma = \mu ( \grad(u) + [\grad(u)]^T ) - pI und \eta = I (\nabla \cdot \delta d) - [\grad(\delta d)]^T.
    // Konvektion: DN = \rho \int (u - w)^T \eta [\grad(u)]^T v d\Omega
    // Impliziter Anteil der Konvektion: DW =  - \rho \frac{1}{k} \int (\delta d \cdot \nabla) u v d\Omega
    // aehnlich zur Newton-Matrix W aus Navier-Stokes, wobei k = dt = Zeitschrittweite
    // Zusaetzlicher Term wg. non-conservatie: DP = - \rho \int ( \grad(w) : \eta ) u \cdot v d\Omega
    // Massematrix und impliziter Anteil wg. non-conservative:
    // DM = \rho \frac{1.5}{k} \int (\nabla \cdot \delta d) u \cdot v d\Omega (Massematrix aus BDF2)
    //    - \rho \frac{1}{k} \int (\nabla \cdot \delta d) u \cdot v d\Omega
    //    = \rho \frac{0.5}{k} \int (\nabla \cdot \delta d) u \cdot v d\Omega
    // ACHTUNG: Der erste Koeffizient \frac{1.5}{k} ist durch den fuehrenden Koeffizienten
    // bei der BDF2-Zeitintegration nach Normierung gegeben (bei BDF2 also 3/2 = 1.5).
    void assemblyShapeDerivativeVelocity(int dim,
                                        std::string FEType1,
                                        std::string FEType2,
                                        MatrixPtr_Type &D,
                                        int FEloc, // 4 = Fluid (Velocity)
                                        MultiVectorPtr_Type u, // Geschwindigkeit
                                        MultiVectorPtr_Type w, // Beschleunigung Gitter
                                        MultiVectorPtr_Type p, // Druck
                                        double dt, // Zeitschrittweite
                                        double rho, // Dichte vom Fluid
                                        double nu, // Viskositaet vom Fluid
                                        bool callFillComplete = true);


    // Shape-Derivatives fuer Div-Nebenbedingung; Hauptvariable ist \delta d und Testfunktion ist q:
    // Divergenz: DB = - \int q \grad(u) : \eta d\Omega
    void assemblyShapeDerivativeDivergence(int dim,
                                           std::string FEType1,
                                           std::string FEType2,
                                           MatrixPtr_Type &DB,
                                           int FEloc1, // 1 = Fluid-Pressure
                                           int FEloc2, // 4 = Fluid-Velocity
                                           MapConstPtr_Type map1_unique, // Pressure-Map
                                           MapConstPtr_Type map2_unique, // Velocity-Map unique als VecField
                                           MultiVectorPtr_Type u, // Geschwindigkeit
                                           bool callFillComplete = true);


    void assemblyRHS(int dim,
                     std::string FEType,
                     MultiVectorPtr_Type  a,
                     std::string fieldType,
                     RhsFunc_Type func,
                     std::vector<SC>& funcParameter
                     );
    
    void assemblyRHSDegTest( int dim,
                             std::string FEType,
                             MultiVectorPtr_Type  a,
                             std::string fieldType,
                             RhsFunc_Type func,
                             std::vector<SC>& funcParameter,
                             int degree) ;

    void buildFullDPhi(vec3D_dbl_ptr_Type dPhi, Teuchos::Array<SmallMatrix<double> >& dPhiMat);

    void fillMatrixArray(SmallMatrix<double> &matIn, double* matArrayOut, std::string order, int offset=0);

    void epsilonTensor(vec_dbl_Type &basisValues, SmallMatrix<SC> &epsilonValues, int activeDof);

    void assemblyNavierStokes(int dim,
								string FETypeVelocity,
								string FETypePressure,
								int degree,
								int dofsVelocity,
								int dofsPressure,
								MultiVectorPtr_Type u_rep,
								MultiVectorPtr_Type p_rep,
								BlockMatrixPtr_Type &A,
								BlockMultiVectorPtr_Type &resVec,
								SmallMatrix_Type coeff,
								ParameterListPtr_Type params,
								bool reAssemble,
							    string assembleMode,
								bool callFillComplete = true,
								int FELocExternal=-1);

    void assemblyLaplaceAssFE(int dim,
                            string FEType,
                            int degree,
                            int dofs,
                            BlockMatrixPtr_Type &A,
                            bool callFillComplete,
                            int FELocExternal=-1);

    void assemblyAceDeformDiffu(int dim,
								string FETypeChem,
								string FETypeSolid,
								int degree,
								int dofsChem,
								int dofsSolid,
								MultiVectorPtr_Type c_rep,
								MultiVectorPtr_Type d_rep,
								BlockMatrixPtr_Type &A,
								BlockMultiVectorPtr_Type &resVec,
								ParameterListPtr_Type params,
							    string assembleMode,
								bool callFillComplete = true,
								int FELocExternal=-1);

    void assemblyAceDeformDiffuBlock(int dim,
                                string FETypeChem,
                                string FETypeSolid,
                                int degree,
                                int dofsChem,
                                int dofsSolid,
                                MultiVectorPtr_Type c_rep,
                                MultiVectorPtr_Type d_rep,
                                BlockMatrixPtr_Type &A,
                                int blockRow,
                                int blockCol,
                                BlockMultiVectorPtr_Type &resVec,
                                int block,
                                ParameterListPtr_Type params,
                                string assembleMode,
                                bool callFillComplete = true,
                                int FELocExternal=-1);

    void advanceInTimeAssemblyFEElements(double dt ,MultiVectorPtr_Type d_rep , MultiVectorPtr_Type c_rep) 
    {
        //UN FElocChem = 1; //checkFE(dim,FETypeChem); // Checks for different domains which belongs to a certain fetype
        UN FElocSolid = 0; //checkFE(dim,FETypeSolid); // Checks for different domains which belongs to a certain fetype

        //ElementsPtr_Type elementsChem= domainVec_.at(FElocChem)->getElementsC();

        ElementsPtr_Type elementsSolid = domainVec_.at(FElocSolid)->getElementsC();
        
        vec_dbl_Type solution_c;
	    vec_dbl_Type solution_d;
        for (UN T=0; T<assemblyFEElements_.size(); T++) {
		    vec_dbl_Type solution(0);

            solution_c = getSolution(elementsSolid->getElement(T).getVectorNodeList(), c_rep,1);
            solution_d = getSolution(elementsSolid->getElement(T).getVectorNodeList(), d_rep,3);
            // First Solid, then Chemistry
            solution.insert( solution.end(), solution_d.begin(), solution_d.end() );
            solution.insert( solution.end(), solution_c.begin(), solution_c.end() );
            
            assemblyFEElements_[T]->updateSolution(solution);

            assemblyFEElements_[T]->advanceInTime(dt);
        }
        
    };

	void assemblyLinearElasticity(int dim,
                                string FEType,
                                int degree,
                                int dofs,
                                MultiVectorPtr_Type d_rep,
                                BlockMatrixPtr_Type &A,
                                BlockMultiVectorPtr_Type &resVec,
                                ParameterListPtr_Type params,
                                bool reAssemble,
                                string assembleMode,
                                bool callFillComplete=true,
                                int FELocExternal=-1);

    void assemblyNonLinearElasticity(int dim,
                                    string FEType,
                                    int degree,
                                    int dofs,
                                    MultiVectorPtr_Type d_rep,
                                    BlockMatrixPtr_Type &A,
                                    BlockMultiVectorPtr_Type &resVec,
                                    ParameterListPtr_Type params,
                                    bool callFillComplete=true,
                                    int FELocExternal=-1);
                                    
    void assemblyNonLinearElasticity(int dim,
                                    string FEType,
                                    int degree,
                                    int dofs,
                                    MultiVectorPtr_Type d_rep,
                                    BlockMatrixPtr_Type &A,
                                    BlockMultiVectorPtr_Type &resVec,
                                    ParameterListPtr_Type params, 									
                                    DomainConstPtr_Type domain,
                                    MultiVectorPtr_Type eModVec,
                                    bool callFillComplete = true,
                                    int FELocExternal=-1);
                                    
    void checkMeshOrientation(int dim,string FEType);


/* ----------------------------------------------------------------------------------------*/
private:
	void addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element1,FiniteElement element2, MapConstPtr_Type mapFirstColumn,MapConstPtr_Type mapSecondColumn, tuple_disk_vec_ptr_Type problemDisk);

	void addFeBlock(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type mapFirstRow, int row, int column, tuple_disk_vec_ptr_Type problemDisk);

    void initAssembleFEElements(string elementType, tuple_disk_vec_ptr_Type problemDisk, ElementsPtr_Type elements, ParameterListPtr_Type params, vec2D_dbl_ptr_Type pointsRep, MapConstPtr_Type elementMap);

    void addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock1,FiniteElement elementBlock2, int dofs1, int dofs2 );

    void addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock, int dofs);
		
    void computeSurfaceNormal(int dim,vec2D_dbl_ptr_Type pointsRep,vec_int_Type nodeList,vec_dbl_Type &v_E, double &norm_v_E);


	AssembleFEPtr_vec_Type assemblyFEElements_;

	vec2D_dbl_Type getCoordinates(vec_LO_Type localIDs, vec2D_dbl_ptr_Type points);
	vec_dbl_Type getSolution(vec_LO_Type localIDs, MultiVectorPtr_Type u_rep, int dofsVelocity);

    //Start of AceGen code
    /*! AceGen code for 3D Neo-Hooke material model
    @param[out] v: values needed for the computaion of F, not needed after computation
    @param[in] E: E module
    @param[in] nu: Poisson ratio
    @param[in] F: deformation gradient, basis functions
    @param[out] P: stresses
    @param[out] A: strains
    */
    void nh3d(double* v, double (*E), double (*Nu), double** F , double** Pmat, double**** Amat);

    /*! AceGen code for 3D Mooney-Rivlin material model
     @param[out] v: values needed for the computaion of F, not needed after computation
     @param[in] E: E module
     @param[in] nu: Poisson ratio
     @param[in] C: material constant
     @param[in] F: deformation gradient, basis functions
     @param[out] P: stresses
     @param[out] A: strains
     */
    void mr3d(double* v, double (*E), double (*Nu), double (*C), double** F, double** Pmat, double**** Amat);

    void stvk3d(double* v,double (*lam),double (*mue),double** F,double** Pmat,double**** Amat);

    void stvk2d(double* v, double (*lam),double (*mue),double** F ,double** Pmat,double**** Amat);
    
    void SMTSetElSpecBiot(ElementSpec *es,int *idata/*not needed*/,int ic,int ng, vec_dbl_Type& paraVec);
    
    void SMTSetElSpecBiotStVK(ElementSpec *es,int *idata/*not needed*/,int ic,int ng, vec_dbl_Type& paraVec);
    
    void SMTSetElSpecBiot3D(ElementSpec *es,int *idata/*not needed*/,int ic,int ng, vec_dbl_Type& paraVec);
    
    void SKR_Biot(double* v,ElementSpec *es,ElementData *ed, NodeSpec **ns, NodeData **nd,double *rdata,int *idata,double *p,double **s);
    
    void SKR_Biot_StVK(double* v,ElementSpec *es,ElementData *ed, NodeSpec **ns, NodeData **nd,double *rdata,int *idata,double *p,double **s);
    
    void SKR_Biot3D(double* v,ElementSpec *es,ElementData *ed, NodeSpec **ns, NodeData **nd,double *rdata,int *idata,double *p,double **s);
    
    //End of AceGen code
    
    void buildTransformation(const vec_int_Type& element,
                             vec2D_dbl_ptr_Type pointsRep,
                             SmallMatrix<SC>& B,
                             std::string FEType="P");
    
    void buildTransformation(const vec_int_Type& element,
                             vec2D_dbl_ptr_Type pointsRep,
                             SmallMatrix<SC>& B,
                             vec_dbl_Type& b,
                             std::string FEType="P");

    
    void buildTransformationSurface(const vec_int_Type& element,
                                    vec2D_dbl_ptr_Type pointsRep,
                                    SmallMatrix<SC>& B,
                                    vec_dbl_Type& b,
                                    std::string FEType="P");

    void applyDiff(vec3D_dbl_Type& dPhiIn,
                   vec3D_dbl_Type& dPhiOut,
                   SmallMatrix<SC>& diffT);
    
    
    void phi(       int Dimension,
                    int intFE,
            		int i,
            		vec_dbl_Type &QuadPts,
            		double* value);


    void gradPhi(	int Dimension,
                    int intFE,
                    int i,
                    vec_dbl_Type &QuadPts,
                    vec_dbl_ptr_Type &value);
    
    /*! Most of the quadrature formulas can be found in http://code-aster.org/doc/v11/en/man_r/r3/r3.01.01.pdf 01/2021  */
    void getQuadratureValues(int Dimension,
                            int Degree,
                            vec2D_dbl_ptr_Type &QuadPts,
                            vec_dbl_ptr_Type &QuadW,
                            std::string FEType);
    
    int getPhi(	vec2D_dbl_ptr_Type &Phi,
                vec_dbl_ptr_Type &weightsPhi,
                int Dimension,
                std::string FEType,
                int Degree,
                std::string FETypeQuadPoints="");

    int getPhiGlobal(vec2D_dbl_ptr_Type &Phi,
                     vec_dbl_ptr_Type &weightsPhi,
                     int Dimension,
                     std::string FEType,
                     int Degree);

    int getDPhi(	vec3D_dbl_ptr_Type &DPhi,
                	vec_dbl_ptr_Type &weightsDPhi,
                    int Dimension,
                    std::string FEType,
                    int Degree);

    int checkFE(int Dimension,
                std::string FEType);

    /*UN determineDegree(UN dim,
                       std::string FEType1,
                       std::string FEType2,
                       VarType type1,
                       VarType type2,
                       UN extraDeg = 0);

    UN determineDegree(UN dim,
                       std::string FEType,
                       VarType type);

    UN determineDegree(UN dim,
                       std::string FEType,
                       UN degFunc);*/
    

    bool setZeros_;
    SC myeps_;
    std::vector<Teuchos::RCP<DataElement> > ed_;
    bool saveAssembly_;
};
}
#endif
