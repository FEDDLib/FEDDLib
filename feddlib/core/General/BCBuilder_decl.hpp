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
    
    BCBuilder();
    
    void addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs);

    void addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs, vec_dbl_Type& parameter_vec);
    
    void addBC(BC_func_Type funcBC, int flag, int block, const DomainPtr_Type &domain, std::string type, int dofs, vec_dbl_Type& parameter_vec, MultiVectorConstPtr_Type& externalSol);
    
    void set(const BlockMatrixPtr_Type &blockMatrix, const BlockMultiVectorPtr_Type &b, double t=0.) const;
    
//    void set(const MatrixPtr_Type &matrix, const MultiVectorPtr_Type &b, double t=0.) const;
    
    void setRHS(const BlockMultiVectorPtr_Type &blockMV, double t=0.) const;

    void setBCMinusVector(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;
    
    void setVectorMinusBC(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;
    
//    void setDirichletBCMinusVector(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;
//    
//    void setVectorMinusDirichletBC(const BlockMultiVectorPtr_Type &outBlockMV, const BlockMultiVectorPtr_Type &substractBlockMV, double t=0.) const;

    void setDirichletBoundaryFromExternal(Teuchos::ArrayRCP<SC>& values/*values will be set to this vector*/, LO index, int loc, double time, std::string type, Teuchos::ArrayRCP<SC> valuesSubstract = Teuchos::null ) const;

    
    void setAllDirichletZero( const BlockMultiVectorPtr_Type &blockMV ) const;
//    void setRHS(const MultiVectorPtr_Type &mv, double t=0.) const;
    
    void setSystem(const BlockMatrixPtr_Type &blockMatrix) const;
    
//    void setSystem(const MatrixPtr_Type &matrix) const;
    
    void setDirichletBC(const MatrixPtr_Type &matrix, int loc, int blockRow, bool isDiagonalBlock) const;

    void setLocalRowOne(const MatrixPtr_Type &matrix, LO localNode, UN dofsPerNode, int loc) const;
    
    void setLocalRowZero(const MatrixPtr_Type &matrix, LO localNode, UN dofsPerNode, int loc) const;

    bool blockHasDirichletBC(int block) const;
    
    bool blockHasDirichletBC(int block, int &loc) const;
        
    bool findFlag(LO flag, int block, int &loc) const;
    
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
