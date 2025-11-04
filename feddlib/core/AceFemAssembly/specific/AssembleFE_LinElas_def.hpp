#ifndef ASSEMBLEFE_LINELAS_DEF_hpp
#define ASSEMBLEFE_LINELAS_DEF_hpp

#include "feddlib/core/AceFemAssembly/AceInterface/NeoHookQuadraticTets.hpp"
#include <vector>
#include <iostream>

namespace FEDD {


/*!

 \brief Constructor for AssembleFE_Laplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem
@param[in} tuple Vector of tuples with Discretization information

*/
template <class SC, class LO, class GO, class NO>
AssembleFE_LinElas<SC,LO,GO,NO>::AssembleFE_LinElas(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{
	/// Extracting values from ParameterList params:
	E_ = this->params_->sublist("Parameter").get("E",3500.0); // the last value is the dafault value, in case no parameter is set
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
void AssembleFE_LinElas<SC,LO,GO,NO>::assembleJacobian() {


	SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement_)); // Matrix we fill with entries.

	assemblyLinElas(elementMatrix); // Function that fills the matrix. We pass though a pointer that will be filled.

	this->jacobian_ = elementMatrix ; // We init the jacobian matrix with the matrix we just build.
}

/*!

 \brief Assembly function 

@param[in] &elementMatrix

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_LinElas<SC,LO,GO,NO>::assemblyLinElas(SmallMatrixPtr_Type &elementMatrix) {

	/// We can access the following values we initialized/extracted in the constructor:
	// dofs_
	// FEType_
	// numNodes_
	// dofsElement_
	// mu_
	// poissonRatio_

	/// Writing entries in the element matrix for nodes 1,2..n , n=numNodes_
	/// 1_x 1_y 1_z 2_x 2_y 2_z .... n_x n_y n_z 
	
	std::vector<double> v(1002); //Working vector, size defined by AceGen-FEAP
	std::vector<double> d(2); // Material parameters
	std::vector<double> ul(30); // The solution vector(or displacement in this case)
	std::vector<double> ul0(30); // Currently unused but must be passed to match FEAP template
	std::vector<double> xl(30); // Nodal Positions in reference coordinates
	std::vector<double> s(900); // Element Stiffness Matrix [Output from skr]
	std::vector<double> p(30); // Residual vector [Output from skr]
	std::vector<double> ht(10); // History parameters currently unused
	std::vector<double> hp(10); // History parameters currently unused

	d[0] = this->E_; // TODO: Check order if there is a problem
	d[1] = this->poissonRatio_;

	for(int i=0;i<30;i++)
		ul[i] = (*this->solution_)[i]; // What is the order? I need it in the form (u1,v1,w1,u2,v2,w2,...)

	int count = 0;
	for(int i=0;i<this->numNodes_;i++)
		for(int j=0;j<this->dofs_;j++){
			xl[count] = this->getNodesRefConfig()[i][j];
			count++;}	

	for(int i=0;i<s.size();i++)
		s[i]=0.0;
	for(int i=0;i<p.size();i++)
		p[i]=0.0;

	// std::cout << "[DEBUG] SKR-Jacobian Calls after this line!" << std::endl;
	skr1(&v[0],&d[0],&ul[0],&ul0[0],&xl[0],&s[0],&p[0],&ht[0],&hp[0]); // Fortran subroutine call modifies s and p
	// std::cout << "[DEBUG] SKR-Jacobian Call successful!" << std::endl;
	// Note: FEAP/Fortran returns matrices unrolled in column major form. This must be converted for use here.

	/*std::cout << "[DEBUG] Printing FEAP output Stiffness vector" << std::endl;
	for (int i=0;i<900;i++)
		std::cout << std::setprecision(767) << s[i] << std::endl;*/

    for (UN i=0; i < this->dofsElement_; i++) {
        for (UN j=0; j < this->dofsElement_; j++) {
            (*elementMatrix)[i][j] = -s[this->dofsElement_*j+i]; // Rolling into a matrix using column major (m*j+i)
        }
    }

}

/*!

 \brief Assembly RHS


*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_LinElas<SC,LO,GO,NO>::assembleRHS() {

	// [Efficiency] Need to know which is called first: assembleRHS() or assembleJacobian(), so that multiple calls to skr() may be avoided.
	// Note skr() computes both elementMatrix_ and rhsVec_
	std::vector<double> v(1002); //Working vector, size defined by AceGen-FEAP
	std::vector<double> d(2); // Material parameters
	std::vector<double> ul(30); // The solution vector(or displacement in this case)
	std::vector<double> ul0(30); // Currently unused but must be passed to match FEAP template
	std::vector<double> xl(30); // Nodal Positions in reference coordinates
	std::vector<double> s(900); // Element Stiffness Matrix [Output from skr]
	std::vector<double> p(30); // Residual vector [Output from skr]
	std::vector<double> ht(10); // History parameters currently unused
	std::vector<double> hp(10); // History parameters currently unused

	d[0] = this->E_; // TODO: Check order if there is a problem
	d[1] = this->poissonRatio_;

	for(int i=0;i<30;i++)
		ul[i] = (*this->solution_)[i];

	int count = 0;
	for(int i=0;i<this->numNodes_;i++)
		for(int j=0;j<this->dofs_;j++){
			xl[count] = this->getNodesRefConfig()[i][j];
			count++;}

	// Initialize arrays to 0
	for(int i=0;i<s.size();i++)
		s[i]=0.0;
	for(int i=0;i<p.size();i++)
		p[i]=0.0;

	// std::cout << "[DEBUG] SKR-Rhs Calls after this line!" << std::endl;
	skr1(&v[0],&d[0],&ul[0],&ul0[0],&xl[0],&s[0],&p[0],&ht[0],&hp[0]); // Fortran subroutine call modifies s and p
	// std::cout << "[DEBUG] SKR-Rhs Call successful!" << std::endl;
	(*this->rhsVec_) = p;
}

}
#endif

