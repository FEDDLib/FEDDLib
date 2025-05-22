#ifndef BCBuilder_decl_hpp
#define BCBuilder_decl_hpp

#define BCBuilder_TIMER

#include "DefaultTypeDefs.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/FE.hpp"

#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/core/LinearAlgebra/BlockMultiVector.hpp"

#include <boost/function.hpp>

/*!
 Declaration of BCBuilder
 
 @brief  BCBuilder
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
      /*!
    \class BCBuilder
    \brief This class is responsible for setting the boundary conditions into the system and rhs

    \tparam SC The scalar type. So far, this is always double, but having it as a template parameter would allow flexibily, e.g., for using complex instead
    \tparam LO The local ordinal type. The is the index type for local indices
    \tparam GO The global ordinal type. The is the index type for global indices
    \tparam NO The Kokkos Node type. This would allow for performance portibility when using Kokkos. Currently, this is not used.

    When we define our test/example we add boundary conditions to the BCBuilder. These conditions match the geometry flags to boundary functions (i.e. zeroDirchlet). 

    Then, when 'setSystem()' or 'setRhs()' is called, the boundary conditiones are set according to the set boundary conditions.

    These boundary conditions are limited to dirichlet boundary conditions. 

    BC Builder contains all boundary condition information of current problem.
   
    */
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class BCBuilder {
    
public:
    
    typedef unsigned UN;
    
    typedef Teuchos::RCP<Domain<SC,LO,GO,NO> > DomainPtr_Type;

    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;
    
    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;
    
    typedef typename  Matrix_Type::Map_Type  Map_Type;
    typedef typename  Matrix_Type::MapPtr_Type  MapPtr_Type;
    typedef typename  Matrix_Type::MapConstPtr_Type MapConstPtr_Type;
    
    typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
    typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
    typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
    
    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    typedef FE<SC,LO,GO,NO> FEFac_Type;
    typedef Teuchos::RCP<FEFac_Type> FEFacPtr_Type;
    
    typedef boost::function<void(double* x, double* res, double t, const double* parameters)>    BC_func_Type;
    
    /// @brief 
    BCBuilder();

    /*! 
     @brief Adding Boundary Condition
    @param funcBC function representing bc
     @param flag flag corresponding to the funcBC
     @param block diagonal block corresponding to the bc (i.e. in block systems like stokes equation)
     @param domain finite element space correspong to block 
     @param type of bc. Dirichlet, Dirichlet_X ...
     @param dofs degrees of freedom of bc (i.e. 3 for a vector valued problem)
    */
    void addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs);

    /*! 
     @brief Adding Boundary Condition with extra parameters
    @param funcBC function representing bc
     @param flag flag corresponding to the funcBC
     @param block diagonal block corresponding to the bc (i.e. in block systems like stokes equation)
     @param domain finite element space correspong to block 
     @param type of bc. Dirichlet, Dirichlet_X ...
     @param dofs degrees of freedom of bc (i.e. 3 for a vector valued problem)
     @param parameter_vec vector containing different parameters for bc
    
    */
    void addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs, vec_dbl_Type& parameter_vec);
    
    /*! 
     @brief Adding Boundary Condition with extra parameters and vector values (i.e. when the inflow of a region is computed)
     @param funcBC function representing bc
     @param flag flag corresponding to the funcBC
     @param block diagonal block corresponding to the bc (i.e. in block systems like stokes equation)
     @param domain finite element space correspong to block 
     @param type of bc. Dirichlet, Dirichlet_X ...
     @param dofs degrees of freedom of bc (i.e. 3 for a vector valued problem)
     @param parameter_vec vector containing different parameters for bc
     @param externalSol vector with values for bc 

    */
    void addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs, vec_dbl_Type& parameter_vec, MultiVectorConstPtr_Type& externalSol);
    

    /// @brief Setting bundary condtions to problem
    /// @param blockMatrix System matrix 
    /// @param b rhs vector
    /// @param t 
    void set(const BlockMatrixPtr_Type &blockMatrix, const BlockMultiVectorPtr_Type &b, double t=0.) const;
    
//    void set(const MatrixPtr_Type &matrix, const MultiVectorPtr_Type &b, double t=0.) const;
    
    /// @brief Setting boundary conditions to (block)vector
    /// @param blockMV vector
    /// @param t 
    void setRHS(const BlockMultiVectorPtr_Type &blockMV, double t=0.) const;

    /// @brief 
    /// @param outBlockMV 
    /// @param substractBlockMV 
    /// @param t 
    void setBCMinusVector(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;
    
    /// @brief 
    /// @param outBlockMV 
    /// @param substractBlockMV 
    /// @param t 
    void setVectorMinusBC(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;
    
//    void setDirichletBCMinusVector(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;
//    
//    void setVectorMinusDirichletBC(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;

    /// @brief 
    /// @param values 
    /// @param index 
    /// @param loc 
    /// @param time 
    /// @param type 
    /// @param valuesSubstract 
    void setDirichletBoundaryFromExternal(Teuchos::ArrayRCP<SC>& values/*values will be set to this vector*/, LO index, int loc, double time, std::string type, Teuchos::ArrayRCP<SC> valuesSubstract = Teuchos::null ) const;

    /// @brief 
    /// @param matrix 
    /// @param loc 
    /// @param blockRow 
    /// @param isDiagonalBlock 
    void setDirichletBCScaled(const MatrixPtr_Type &matrix, int loc, int blockRow, bool isDiagonalBlock, double eps=1.0) const;

    /// @brief 
    /// @param matrix 
    /// @param localNode 
    /// @param dofsPerNode 
    /// @param loc 
    void setLocalRowEntry(const MatrixPtr_Type &matrix, LO localNode, UN dofsPerNode, int loc, double eps) const;
    
    /// @brief 
    /// @param blockMV 
    void setAllDirichletZero( const BlockMultiVectorPtr_Type &blockMV ) const;
//    void setRHS(const MultiVectorPtr_Type &mv, double t=0.) const;
    
    /// @brief Set boundary conditions to system
    /// @param blockMatrix 
    void setSystem(const BlockMatrixPtr_Type &blockMatrix) const;
    
    /// @brief Set boundary conditions to system
    /// @param blockMatrix 
    void setSystemScaled(const BlockMatrixPtr_Type &blockMatrix,double eps=1.0) const;
        
    /// @brief 
    /// @param matrix 
    /// @param loc 
    /// @param blockRow 
    /// @param isDiagonalBlock 
    void setDirichletBC(const MatrixPtr_Type &matrix, int loc, int blockRow, bool isDiagonalBlock) const;

    /// @brief 
    /// @param matrix 
    /// @param localNode 
    /// @param dofsPerNode 
    /// @param loc 
    void setLocalRowOne(const MatrixPtr_Type &matrix, LO localNode, UN dofsPerNode, int loc) const;
    
    /// @brief 
    /// @param matrix 
    /// @param localNode 
    /// @param dofsPerNode 
    /// @param loc 
    void setLocalRowZero(const MatrixPtr_Type &matrix, LO localNode, UN dofsPerNode, int loc) const;

    /// @brief 
    /// @param block 
    /// @return 
    bool blockHasDirichletBC(int block) const;
    
    /// @brief 
    /// @param block 
    /// @param loc 
    /// @return 
    bool blockHasDirichletBC(int block, int &loc) const;
        
    /// @brief 
    /// @param flag 
    /// @param block 
    /// @param loc 
    /// @return 
    bool findFlag(LO flag, int block, int &loc) const;
    
    /// @brief 
    /// @param block 
    /// @return 
    int dofsPerNodeAtBlock(int block);
        
//    DomainPtr_Type domainOfBlock(int block) const;
    
private:
    
    std::vector<BC_func_Type> vecBC_func_;
    vec_int_Type vecFlag_;
    vec_int_Type vecBlockID_;
    std::vector<DomainPtr_Type> vecDomain_;
    std::vector<std::string> vecBCType_;
    vec_int_Type vecDofs_;
    vec2D_dbl_Type vecBC_Parameters_;
    std::vector<MultiVectorConstPtr_Type> vecExternalSol_;
    mutable vec_dbl_ptr_Type resultPtr_;
    mutable vec_dbl_ptr_Type pointPtr_;
#ifdef BCBuilder_TIMER
    TimePtr_Type SetSystemRowTimer_;
    TimePtr_Type BlockRowHasDirichletTimer_;
    TimePtr_Type FindFlagTimer_;
#endif

    
};
}


#endif
