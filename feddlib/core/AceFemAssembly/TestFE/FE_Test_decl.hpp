#ifndef FE_TEST_DECL_hpp
#define FE_TEST_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/General/SmallMatrix.hpp"
#include "feddlib/core/General/DefaultTypeDefs.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector.hpp"
#include "feddlib/core/FE/Domain.hpp"
#include "feddlib/core/FE/sms.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_BLAS.hpp>

#include <boost/function.hpp>

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


    typedef typename Matrix_Type::MapPtr_Type MapPtr_Type;
    typedef typename Matrix_Type::MapConstPtr_Type MapConstPtr_Type;

    typedef boost::function<void(double* x, double* res, double t, const double* parameters)> BC_func_Type;
    


    FE_Test(bool saveAssembly=false);

    void addFE(DomainConstPtr_Type domain);

	void assemblyLaplace(int dim,
					    string FEType,
					    int degree,
					    MatrixPtr_Type &A,
					    bool callFillComplete,
					     int FELocExternal );

	void assemblyRHS(int dim,
                       string FEType,
                       MultiVectorPtr_Type a,
                       string fieldType,
                       RhsFunc_Type func,
                      vector<SC>& funcParameter
                      );

	void addFeMatrix(MatrixPtr_Type &A, SmallMatrix<SC> elementMatrix, FiniteElement element, MapConstPtr_Type map);
	void addFeVector(VectorPtr_Type &a, vec_dbl_Type elementVector, FiniteElement element);

    bool saveAssembly_;
    DomainConstPtr_vec_Type	domainVec_;
};
}
#endif










