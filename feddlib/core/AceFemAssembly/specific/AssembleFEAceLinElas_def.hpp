#ifndef ASSEMBLEFEACELINELAS_DEF_hpp
#define ASSEMBLEFEACELINELAS_DEF_hpp

#include "AssembleFEAceLinElas_decl.hpp"

namespace FEDD {


/*!

 \brief Constructor for AssembleFEAceLaplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem
@param[in} tuple Vector of tuples with Discretization information

*/
template <class SC, class LO, class GO, class NO>
AssembleFEAceLinElas<SC,LO,GO,NO>::AssembleFEAceLinElas(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{
	/// Extracting values from ParameterList params:
	mu_ = this->params_->sublist("Parameter").get("Mu",0.3571); // the last value is the dafault value, in case no parameter is set
    //lambda_ = this->params_->sublist("Parameter").get("lambda",1.);
    poissonRatio_ = this->params_->sublist("Parameter").get("Poisson Ratio",0.4e-0);

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	FEType_ = std::get<1>(this->diskTuple_->at(0)); // FEType of Disk
	dofs_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
	numNodes_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element

	dofsElement_ = dofs_*numNodes_; // "Dimension of return matrix"

}

/*!

 \brief Assembly Jacobian is simply assemblyLaplacian for Laplace Problem

@param[in] &elementMatrix

*/ 

template <class SC, class LO, class GO, class NO>
void AssembleFEAceLinElas<SC,LO,GO,NO>::assembleJacobian() {


	SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement_)); // Matrix we fill with entries.

	assemblyLinElas(elementMatrix); // Function that fills the matrix. We pass though a pointer that will be filled.

	this->jacobian_ = elementMatrix ; // We init the jacobian matrix with the matrix we just build.
}

/*!

 \brief Assembly function 

@param[in] &elementMatrix

*/
template <class SC, class LO, class GO, class NO>
void AssembleFEAceLinElas<SC,LO,GO,NO>::assemblyLinElas(SmallMatrixPtr_Type &elementMatrix) {

	/// We can access the following values we initialized/extracted in the constructor:
	// dofs_
	// FEType_
	// numNodes_
	// dofsElement_
	// mu_
	// poissonRatio_

	/// Writing entries in the element matrix for nodes 1,2..n , n=numNodes_
	/// 1_x 1_y 1_z 2_x 2_y 2_z .... n_x n_y n_z 
	
	
    for (UN i=0; i < dofsElement_; i++) {
        for (UN j=0; j < dofsElement_; j++) {
              (*elementMatrix)[i][j] = 0.; // <-- insert correct value here. For now = zero
        }
    }

}

/*!

 \brief Assembly RHS


*/
template <class SC, class LO, class GO, class NO>
void AssembleFEAceLinElas<SC,LO,GO,NO>::assembleRHS() {



}

}
#endif

