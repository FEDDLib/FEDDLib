#ifndef FE_ELEMENTASSEMBLY_DECL_hpp
#define FE_ELEMENTASSEMBLY_DECL_hpp


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
#include "feddlib/core/AceFemAssembly/specific/AssembleFE_SCI_SMC_Active_Growth_Reorientation.hpp"
#include "feddlib/core/AceFemAssembly/specific/AssembleFENavierStokes.hpp"

#include "feddlib/core/AceFemAssembly/AssembleFEFactory.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_BLAS.hpp>


/*!
 Declaration of FE

 @brief  FE
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {


template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class FE_ElementAssembly {
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

    FE_ElementAssembly(bool saveAssembly=false);
    
    void addFE(DomainConstPtr_Type domain);

    void doSetZeros(double eps = 10*Teuchos::ScalarTraits<SC>::eps());

    void assemblyEmptyMatrix(MatrixPtr_Type &A);

    void assemblyNonlinearLaplace(
        int dim, std::string FEType, int degree, MultiVectorPtr_Type u_rep,
        BlockMatrixPtr_Type &A, BlockMultiVectorPtr_Type &resVec,
        ParameterListPtr_Type params, std::string assembleMode,
        bool callFillComplete = true, int FELocExternal = -1);


    void assemblyNavierStokes(int dim,
								std::string FETypeVelocity,
								std::string FETypePressure,
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
							        std::string assembleMode,
								bool callFillComplete = true,
								int FELocExternal=-1);

    void assemblyLaplaceAssFE(int dim,
                            std::string FEType,
                            int degree,
                            int dofs,
                            BlockMatrixPtr_Type &A,
                            bool callFillComplete,
                            int FELocExternal=-1);

    void assemblyAceDeformDiffu(int dim,
								std::string FETypeChem,
								std::string FETypeSolid,
								int degree,
								int dofsChem,
								int dofsSolid,
								MultiVectorPtr_Type c_rep,
								MultiVectorPtr_Type d_rep,
								BlockMatrixPtr_Type &A,
								BlockMultiVectorPtr_Type &resVec,
								ParameterListPtr_Type params,
							        std::string assembleMode,
								bool callFillComplete = true,
								int FELocExternal=-1);

    void assemblyAceDeformDiffuBlock(int dim,
                                std::string FETypeChem,
                                std::string FETypeSolid,
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
                                std::string assembleMode,
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
        
    }

	void assemblyLinearElasticity(int dim,
                                std::string FEType,
                                int degree,
                                int dofs,
                                MultiVectorPtr_Type d_rep,
                                BlockMatrixPtr_Type &A,
                                BlockMultiVectorPtr_Type &resVec,
                                ParameterListPtr_Type params,
                                bool reAssemble,
                                std::string assembleMode,
                                bool callFillComplete=true,
                                int FELocExternal=-1);

    void assemblyNonLinearElasticity(int dim,
                                    std::string FEType,
                                    int degree,
                                    int dofs,
                                    MultiVectorPtr_Type d_rep,
                                    BlockMatrixPtr_Type &A,
                                    BlockMultiVectorPtr_Type &resVec,
                                    ParameterListPtr_Type params,
                                    bool callFillComplete=true,
                                    int FELocExternal=-1);
                                    
    void assemblyNonLinearElasticity(int dim,
                                    std::string FEType,
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
                                    

    /* Given a converged velocity solution this function 
       computes the averaged viscosity estimate in each cell at the center of mass 
       - CM stands for center of mass so the values at the node are averaged to obtain one value
    */
    void computeSteadyViscosityFE_CM(int dim,
	                                    std::string FETypeVelocity,
	                                    std::string FETypePressure,
										int dofsVelocity,
										int dofsPressure,
										MultiVectorPtr_Type u_rep,
										MultiVectorPtr_Type p_rep,
 										ParameterListPtr_Type params);

    // Change for all assembleFEElements the linearization type to the specified ones
    void  changeLinearizationFE(std::string linearization);

    // Write prostprocessing output fields like e.g. the viscosity based on  velocity, pressure .. solution
    // inside this BMV -> For visualization or postprocessing                                
    BlockMultiVectorPtr_Type const_output_fields;

    DomainConstPtr_vec_Type	domainVec_;


/* ----------------------------------------------------------------------------------------*/
protected:
    // Protected attributes/ functions such that derived FE classes can access them
    bool setZeros_;
    SC myeps_;
    bool saveAssembly_;

    int checkFE(int Dimension,
                std::string FEType);

    vec2D_dbl_Type getCoordinates(vec_LO_Type localIDs, vec2D_dbl_ptr_Type points);
	vec_dbl_Type getSolution(vec_LO_Type localIDs, MultiVectorPtr_Type u_rep, int dofsVelocity);


/* ----------------------------------------------------------------------------------------*/
private:
    // Functions to add element contributions to global system matrices/ vectors
	void addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element1,FiniteElement element2, MapConstPtr_Type mapFirstColumn,MapConstPtr_Type mapSecondColumn, tuple_disk_vec_ptr_Type problemDisk);

	void addFeBlock(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type mapFirstRow, int row, int column, tuple_disk_vec_ptr_Type problemDisk);

    void initAssembleFEElements(std::string elementType, tuple_disk_vec_ptr_Type problemDisk, ElementsPtr_Type elements, ParameterListPtr_Type params, vec2D_dbl_ptr_Type pointsRep, MapConstPtr_Type elementMap);

    void addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock1,FiniteElement elementBlock2, int dofs1, int dofs2 );

    void addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock, int dofs);


	AssembleFEPtr_vec_Type assemblyFEElements_;


};
}
#endif
