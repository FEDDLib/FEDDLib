#ifndef FE_DECL_hpp
#define FE_DECL_hpp

#include "FE_ElementAssembly.hpp"

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
class FE: public FE_ElementAssembly<SC,LO,GO,NO> {
  public:

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

    typedef AssembleFEGeneralizedNewtonian<SC,LO,GO,NO> AssembleFEGeneralizedNewtonian_Type;
    typedef Teuchos::RCP<AssembleFEGeneralizedNewtonian_Type> AssembleFEGeneralizedNewtonianPtr_Type;

    typedef AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO> AssembleFE_SCI_SMC_Active_Growth_Reorientation_Type;
    typedef Teuchos::RCP<AssembleFE_SCI_SMC_Active_Growth_Reorientation_Type> AssembleFE_SCI_SMC_Active_Growth_Reorientation_Ptr_Type;

    typedef std::vector<AssembleFEPtr_Type> AssembleFEPtr_vec_Type;	

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type ;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type ;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

	typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

    /* ###################################################################### */

    FE(bool saveAssembly=false);
    
    void assemblyIdentity(MatrixPtr_Type &A);
    
    void assemblySurfaceRobinBC(int dim, 
                                std::string FETypeP, 
                                std::string FETypeV, 
                                MultiVectorPtr_Type u, 
                                MatrixPtr_Type A, 
                                std::vector<SC> &funcParameter, 
                                RhsFunc_Type func, 
                                ParameterListPtr_Type params);

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
    
    double assemblyResistanceBoundary(int dim, 
                                std::string FEType, 
                                MultiVectorPtr_Type f, 
                                MultiVectorPtr_Type u_rep,
                                vec_dbl_Type flowRate_vec, 
                                std::vector<SC>& funcParameter, 
                                RhsFunc_Type func, 
                                ParameterListPtr_Type params, 
                                int FEloc=0);
    double assemblyAbsorbingBoundary(int dim, 
                                std::string FEType, 
                                MultiVectorPtr_Type f, 
                                MultiVectorPtr_Type u_rep, 
                                vec_dbl_Type flowRate_vec, 
                                std::vector<SC>& funcParameter, 
                                RhsFunc_Type func, 
                                double areaOutlet_init, 
                                double areaOutlet_T, 
                                ParameterListPtr_Type params, 
                                int FEloc=0);

    double assemblyAbsorbingBoundaryPaper(int dim, 
                                std::string FEType, 
                                MultiVectorPtr_Type f, 
                                MultiVectorPtr_Type u_rep, 
                                vec_dbl_Type flowRate_vec, 
                                std::vector<SC>& funcParameter, 
                                RhsFunc_Type func, 
                                double areaOutlet_init, 
                                double areaOutlet_T, 
                                ParameterListPtr_Type params, 
                                int FEloc=0);
    double assemblyAbsorbingResistanceBoundary(int dim, 
                                std::string FEType, 
                                MultiVectorPtr_Type f, 
                                MultiVectorPtr_Type u_rep, 
                                vec_dbl_Type flowRate_vec, 
                                std::vector<SC>& funcParameter, 
                                RhsFunc_Type func, 
                                double areaOutlet_init, 
                                ParameterListPtr_Type params, 
                                int FEloc=0);
    void assemblyArea(int dim,
                            double &area,
                            int inflowFlag,
                            int FEloc=0);       
                                 
    int assemblyFlowRate(int dim,
                            double &flowRateParabolic,
                            std::string FEType, 
                            int dofs,
                            int inflowFlag,
                            MultiVectorPtr_Type solution_rep,
                            int FEloc=0) ;

    void assemblyAverageVelocity(int dim, 
                                double &averageVelocity, 
                                std::string FEType, 
                                int dofs, int flag, 
                                MultiVectorPtr_Type solution_rep, 
                                int FEloc=0);

    void assemblyPressureMeanValue(int dim, std::string FEType, MultiVectorPtr_Type a);

    // void assemblyAceGenTPM( MatrixPtr_Type &A00,
    //                         MatrixPtr_Type &A01,
    //                         MatrixPtr_Type &A10,
    //                         MatrixPtr_Type &A11,
    //                         MultiVectorPtr_Type &F0,
    //                         MultiVectorPtr_Type &F1,
    //                         MapPtr_Type &mapRepeated1,
    //                         MapPtr_Type &mapRepeated2,
    //                         ParameterListPtr_Type parameterList,
    //                         MultiVectorPtr_Type u_repeatedNewton=Teuchos::null,
    //                         MultiVectorPtr_Type p_repeatedNewton=Teuchos::null,
    //                         MultiVectorPtr_Type u_repeatedTime=Teuchos::null,
    //                         MultiVectorPtr_Type p_repeatedTime=Teuchos::null,
    //                        bool update=true,
    //                        bool updateHistory=true);
    
    
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


    // // Assembling the reaction term of the reaction diffusion equation. Maybe add default function.
	// void assemblyLinearReactionTerm(int dim,
    // 							std::string FEType,
    //                             MatrixPtr_Type &A,
    //                             bool callFillComplete,
    //                  			std::vector<SC>& funcParameter,
	// 							RhsFunc_Type reactionFunc);	

	// // Assembling the reaction term of the reaction diffusion equation. Maybe add default function.
	// void assemblyReactionTerm(int dim,
    // 							std::string FEType,
    //                             MatrixPtr_Type &A,
    //                             MultiVectorPtr_Type u,
    //                             bool callFillComplete,
    //                  			std::vector<SC>& funcParameter,
	// 							RhsFunc_Type reactionFunc);	

    //                             // Assembling the reaction term of the reaction diffusion equation. Maybe add default function.
	// void assemblyDReactionTerm(int dim,
    // 							std::string FEType,
    //                             MatrixPtr_Type &A,
    //                             MultiVectorPtr_Type u,
    //                             bool callFillComplete,
    //                  			std::vector<SC>& funcParameter,
	// 							RhsFunc_Type reactionFunc);

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

    void assemblyAdvectionVecFieldScalar(int dim,
                                    std::string FEType,
                                    std::string FETypeV,
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

    Teuchos::RCP<ElementSpec> es_;

/* ----------------------------------------------------------------------------------------*/
private:

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

    /*! AceGen code for 3D Saint-Venant-Kirchhoff material model 
    @param[out] v: values needed for the computaion of F, not needed after computation
    @param[in] lam: lambda
    @param[in] mue: mu
    @param[in] F: deformation gradient, basis functions
    @param[out] P: stresses
    @param[out] A: strains
    */
    void stvk3d(double* v,double (*lam),double (*mue),double** F,double** Pmat,double**** Amat);
    
    /*! AceGen code for 2D Saint-Venant-Kirchhoff material model 
    @param[out] v: values needed for the computaion of F, not needed after computation
    @param[in] lam: lambda
    @param[in] mue: mu
    @param[in] F: deformation gradient, basis functions
    @param[out] P: stresses
    @param[out] A: strains
    */
    void stvk2d(double* v, double (*lam),double (*mue),double** F ,double** Pmat,double**** Amat);
    
    // void SMTSetElSpecBiot(ElementSpec *es,int *idata/*not needed*/,int ic,int ng, vec_dbl_Type& paraVec);
    
    // void SMTSetElSpecBiotStVK(ElementSpec *es,int *idata/*not needed*/,int ic,int ng, vec_dbl_Type& paraVec);
    
    // void SMTSetElSpecBiot3D(ElementSpec *es,int *idata/*not needed*/,int ic,int ng, vec_dbl_Type& paraVec);
    
    // void SKR_Biot(double* v,ElementSpec *es,ElementData *ed, NodeSpec **ns, NodeData **nd,double *rdata,int *idata,double *p,double **s);
    
    // void SKR_Biot_StVK(double* v,ElementSpec *es,ElementData *ed, NodeSpec **ns, NodeData **nd,double *rdata,int *idata,double *p,double **s);
    
    // void SKR_Biot3D(double* v,ElementSpec *es,ElementData *ed, NodeSpec **ns, NodeData **nd,double *rdata,int *idata,double *p,double **s);
    
    //End of AceGen code

    void applyDiff(vec3D_dbl_Type& dPhiIn,
                   vec3D_dbl_Type& dPhiOut,
                   SmallMatrix<SC>& diffT);
    
    
    std::vector<Teuchos::RCP<DataElement> > ed_;
};
}
#endif
