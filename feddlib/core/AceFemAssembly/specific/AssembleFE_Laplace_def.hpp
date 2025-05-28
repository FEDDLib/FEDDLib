#ifndef ASSEMBLEFE_LAPLACE_DEF_hpp
#define ASSEMBLEFE_LAPLACE_DEF_hpp

#include "AssembleFE_Laplace_decl.hpp"

namespace FEDD {


/*!

 \brief Constructor for AssembleFE_Laplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem

*/
template <class SC, class LO, class GO, class NO>
AssembleFE_Laplace<SC,LO,GO,NO>::AssembleFE_Laplace(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{

}

/*!

 \brief Assembly Jacobian is simply assemblyLaplacian for Laplace Problem

@param[in] &elementMatrix

*/ 

template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::assembleJacobian() {

	int nodesElement = this->nodesRefConfig_.size();
	int dofs = std::get<2>(this->diskTuple_->at(0));
	int dofsElement = nodesElement*dofs;
	SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement));

	assemblyLaplacian(elementMatrix);

	this->jacobian_ = elementMatrix ;
}

/*!

 \brief Assembly function for \f$ \int_T \nabla v \cdot \nabla u ~dx\f$ 

@param[in] &elementMatrix

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= std::get<3>(this->diskTuple_->at(0));//this->getNodesRefConfig().size();
	std::string FEType = std::get<1>(this->diskTuple_->at(0));
	int dofs = std::get<2>(this->diskTuple_->at(0));

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    // inner( grad(u) , grad(v) ) has twice the polyonimial degree than grad(u) or grad(v).
    UN deg = 2*Helper::determineDegree(dim,FEType,Helper::Grad);
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    buildTransformation(B);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    applyBTinv( dPhi, dPhiTrans, Binv );
    for (UN i=0; i < numNodes; i++) {
        Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
        for (UN j=0; j < numNodes; j++) {
            for (UN w=0; w<dPhiTrans.size(); w++) {
                for (UN d=0; d<dim; d++){
                    value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                }
            }
            value[j] *= absDetB;
			
			for (UN d=0; d<dofs; d++) {
              (*elementMatrix)[i*dofs +d][j*dofs+d] = value[j];
            }
        }

    }

}

/*!

 \brief Assembly function for \f$ \int_T f ~ v ~dx \f$, we need to 

@param[in] &elementVector

*/
template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::assembleRHS() {


	int dim = this->getDim();
	int numNodes= std::get<3>(this->diskTuple_->at(0));//this->getNodesRefConfig().size();
	std::string FEType = std::get<1>(this->diskTuple_->at(0));
	vec_dbl_Type elementVector(numNodes);

    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

	UN extraDegree = 1; // TODO: [JK] 2025/04 This should be user defined, based on the irregularity of f(x). Function parameter? Parameter list?
    UN deg = Helper::determineDegree(dim,FEType,Helper::Std) + extraDegree;
    Helper::getPhi(phi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    this->buildTransformation(B);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    std::vector<double> paras0(1);

    double x;

    SC value;

    //for now just const!
    std::vector<double> valueFunc(dim);
    SC* paras = &(paras0[0]);
    
    this->rhsFunc_( &x, &valueFunc[0], paras );

    for (UN i=0; i < phi->at(0).size(); i++) {
  	    value = Teuchos::ScalarTraits<SC>::zero();
		for (UN w=0; w<weights->size(); w++){
       		value += weights->at(w) * phi->at(w).at(i);
		}
        value *= absDetB *valueFunc[0];
        elementVector[i] += value;
    }

	(*this->rhsVec_) = elementVector;
}

/*!

 \brief Building Transformation

@param[in] &B

*/

template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::buildTransformation(SmallMatrix<SC>& B){

    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = 0;
    for (UN j=0; j<B.size(); j++) {
        index = j+1;
        for (UN i=0; i<B.size(); i++) {
            B[i][j] = this->nodesRefConfig_.at(index).at(i) - this->nodesRefConfig_.at(index0).at(i);
        }
    }

}

/*!

 \brief Assembly function for \f$ \int_T f ~ v ~dx \f$, we need to 

@param[in] &elementVector

*/

template <class SC, class LO, class GO, class NO>
void AssembleFE_Laplace<SC,LO,GO,NO>::applyBTinv( vec3D_dbl_ptr_Type& dPhiIn,
                                    vec3D_dbl_Type& dPhiOut,
                                    SmallMatrix<SC>& Binv){
    UN dim = Binv.size();
    for (UN w=0; w<dPhiIn->size(); w++){
        for (UN i=0; i < dPhiIn->at(w).size(); i++) {
            for (UN d1=0; d1<dim; d1++) {
                for (UN d2=0; d2<dim; d2++) {
                    dPhiOut[w][i][d1] += dPhiIn->at(w).at(i).at(d2) * Binv[d2][d1];
                }
            }
        }
    }
}

}
#endif

