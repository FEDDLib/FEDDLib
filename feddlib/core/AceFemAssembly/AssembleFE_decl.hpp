#ifndef ASSEMBLEFE_DECL_hpp
#define ASSEMBLEFE_DECL_hpp


#include "feddlib/core/FEDDCore.hpp"
#include "feddlib/core/LinearAlgebra/Matrix.hpp"

namespace FEDD {
template <class SC = default_sc, class LO = default_lo, class GO = default_go, class NO = default_no>
class AssembleFE {
  public:

	typedef Matrix<SC,LO,GO,NO> Matrix_Type;
   	typedef Teuchos::RCP<Matrix_Type> MatrixPtr_Type;

	typedef MultiVector<SC,LO,GO,NO> MultiVector_Type;
	typedef Teuchos::RCP<MultiVector_Type> MultiVectorPtr_Type;
	
	typedef SmallMatrix<SC> SmallMatrix_Type;
    	typedef Teuchos::RCP<SmallMatrix_Type> SmallMatrixPtr_Type;

	typedef AssembleFE<SC,LO,GO,NO> AssembleFE_Type;

	AssembleFE(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type parameters);
		
	virtual void assemblyJacobian(SmallMatrixPtr_Type &elementMatrix) =0;
	virtual void assemblyRHS(MultiVectorPtr_Type &elementVector) =0;

	//virtual void assemblyMass(MatrixPtr_Type &A) =0;

	void checkParameter();
	void updateParams( ParameterListPtr_Type params);
	void advanceInTime(double dt);
	double getTimestep();

	void updateSolution(vec_dbl_Type solution);
	vec_dbl_Type getSolution();

	void preProcessing();
	void postProcessing(); // gespeicherte Daten rein und rausziehen. Manche Sachen nachträglich nicht mehr ändern
	
	int getDim();
	vec2D_dbl_Type getNodesRefConfig();

  protected:
	// For example if we consider an element with multiple discretizations i.e. P2-P1 Velocity-Pressure.
	// We might need more than one FEType or Degree of freedom information
	int numFEType_;
	string FEType1_; 
	string FEType2_;
	int dofs1_;
	int dofs2_;

	int dim_;

	vec2D_dbl_Type nodesRefConfig_;
	bool timeProblem_;
	int flag_;
	double timeStep_ ;
	ParameterListPtr_Type paramsMaterial_;// Including flags
	ParameterListPtr_Type params_;// Including flags
	vec_dbl_Type solution_ ; 
	

};
}
#endif

