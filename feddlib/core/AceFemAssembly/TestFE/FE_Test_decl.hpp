#ifndef FE_TEST_DECL_hpp
#define FE_TEST_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/SmallMatrix.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/BlockMatrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/LinearAlgebra/Map.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/sms.hpp"
#include "feddlib/core/AceFemAssembly/AssembleFE.hpp"
#include "feddlib/core/AceFemAssembly/AssembleFEFactory.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_BLAS.hpp>

#include <boost/function.hpp>
#include <tuple>  
namespace FEDD {

template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class FE_Test {
  public:


    typedef Domain<SC,LO,GO,NO> Domain_Type;
    typedef Teuchos::RCP<Domain_Type> DomainPtr_Type;
    typedef Teuchos::RCP<const Domain_Type> DomainConstPtr_Type;
    typedef std::vector<DomainConstPtr_Type> DomainConstPtr_vec_Type;

    typedef Elements Elements_Type;
    typedef Teuchos::RCP<Elements_Type> ElementsPtr_Type;
    
    typedef Matrix<SC,LO,GO,NO> Matrix_Type;
    typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

	typedef SmallMatrix<SC> SmallMatrix_Type;
    typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	typedef Teuchos::RCP<const MultiVector_Type> MultiVectorConstPtr_Type;
	    
	typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;
    typedef Teuchos::RCP<AssembleFE_Type> AssembleFEPtr_Type;
    typedef std::vector<AssembleFEPtr_Type> AssembleFEPtr_vec_Type;	

    typedef typename Matrix_Type::Map_Type Map_Type;
    typedef typename Matrix_Type::MapPtr_Type MapPtr_Type;
    typedef typename Matrix_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef BlockMatrix<SC,LO,GO,NO> BlockMatrix_Type ;
    typedef Teuchos::RCP<BlockMatrix_Type> BlockMatrixPtr_Type;

    typedef BlockMultiVector<SC,LO,GO,NO> BlockMultiVector_Type ;
    typedef Teuchos::RCP<BlockMultiVector_Type> BlockMultiVectorPtr_Type;

    typedef boost::function<void(double* x, double* res, double t, const double* parameters)> BC_func_Type;

    FE_Test(bool saveAssembly=false);

    void addFE(DomainConstPtr_Type domain);

	void assemblyLaplace(int dim,
				string FEType,
				int degree,
				int dofs,
				MatrixPtr_Type &A,
				bool callFillComplete = true,
				int FELocExternal = -1 );

	void assemblyLinElas(int dim,
				string FEType,
				int degree,
				int dofs,
				MatrixPtr_Type &A,
				bool callFillComplete = true,
				int FELocExternal = -1 );
				
	void assemblyNonLinElas(int dim,
                            string FEType,
                            int degree,
                            int dofs,
                            MultiVectorPtr_Type d_rep,
                            MatrixPtr_Type &A,
                            MultiVectorPtr_Type &resVec,
                            ParameterListPtr_Type params,
                            bool reAssemble,
                            string assembleMode,
                            bool callFillComplete = true,
                            int FELocExternal = -1);

	void assemblyNavierStokes(int dim,
								string FETypeVelocity,
								string FETypePressure,
								int degree,
								int dofsVelocity,
								int dofsPressure,
								MultiVectorPtr_Type u_rep,
								BlockMatrixPtr_Type &A,
								bool callFillComplete = true,
								int FELocExternal=-1);

	void assemblyRHS(int dim,
                       string FEType,
                       MultiVectorPtr_Type a,
                       string fieldType,
                       RhsFunc_Type func,
                      vector<SC>& funcParameter
                      );

	bool saveAssembly_;
    DomainConstPtr_vec_Type	domainVec_;

	private:
		void addFeMatrix(MatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type map, int dofs);
		void addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type mapFirstColumn,MapConstPtr_Type mapSecondColumn, tuple_disk_vec_ptr_Type problemDisk);

		void addFeMv(MultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock, int dofs);
			
		void initAssembleFEElements(string elementType,tuple_disk_vec_ptr_Type problemDisk,ElementsPtr_Type elements, ParameterListPtr_Type params,vec2D_dbl_ptr_Type pointsRep);
		int checkFE(int dim, string FEType);

		AssembleFEPtr_vec_Type assemblyFEElements_;

		vec2D_dbl_Type getCoordinates(vec_LO_Type localIDs, vec2D_dbl_ptr_Type points);
		vec_dbl_Type getSolution(vec_LO_Type localIDs, MultiVectorPtr_Type u_rep, int dofsVelocity);
    
};
}
#endif










