#ifndef ASSEMBLEFEACELAPLACE_DEF_hpp
#define ASSEMBLEFEACELAPLACE_DEF_hpp

#include "AssembleFEAceLaplace_decl.hpp"

namespace FEDD {


/*!

 \brief Constructor for AssembleFEAceLaplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem

*/
template <class SC, class LO, class GO, class NO>
AssembleFEAceLaplace<SC,LO,GO,NO>::AssembleFEAceLaplace(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params)
{

}

/*!

 \brief Assembly Jacobian is simply assemblyLaplacian for Laplace Problem

@param[in] &elementMatrix

*/ 

template <class SC, class LO, class GO, class NO>
void AssembleFEAceLaplace<SC,LO,GO,NO>::assemblyJacobian(SmallMatrixPtr_Type &elementMatrix) {

	assemblyLaplacian(elementMatrix);
}

/*!

 \brief Assembly function for \f$ \int_T \nabla v \cdot \nabla u ~dx\f$ 

@param[in] &elementMatrix

*/
template <class SC, class LO, class GO, class NO>
void AssembleFEAceLaplace<SC,LO,GO,NO>::assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= this->getNodesRefConfig().size();

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    //UN deg = determineDegree(dim,FEType);
    //getDPhi(dPhi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    //this->buildTransformation(B);
    //detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    //applyBTinv( dPhi, dPhiTrans, Binv );
    for (UN i=0; i < numNodes; i++) {
        Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
        for (UN j=0; j < numNodes; j++) {
            for (UN w=0; w<dPhiTrans.size(); w++) {
                for (UN d=0; d<dim; d++){
                    value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                }
            }
            value[j] *= absDetB;
			(*elementMatrix)[i][j] = value[j];

		
        }

    }

}

/*!

 \brief Assembly function for \f$ \int_T f ~ v ~dx \f$, we need to 

@param[in] &elementVector

*/
template <class SC, class LO, class GO, class NO>
void AssembleFEAceLaplace<SC,LO,GO,NO>::assemblyRHS(MultiVectorPtr_Type &elementVector) {

	Teuchos::ArrayRCP< SC > elementVec = elementVector->getDataNonConst(0);

	int dim = this->getDim();
	int numNodes= this->getNodesRefConfig().size();

    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    //UN deg = determineDegree(dim,FEType,Grad);
    //getPhi(phi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    //this->buildTransformation(B);
    //detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    SC value;
    for (UN i=0; i < phi->at(0).size(); i++) {
  	    value = Teuchos::ScalarTraits<SC>::zero();
		for (UN w=0; w<weights->size(); w++){
			//func(&quadPointsTrans[w][0], &valueFunc[0] ,paras); // We need to define the function
       		value += weights->at(w) * phi->at(w).at(i);//*valueFunc[0];
		}
        value *= absDetB;
        elementVec[i] += value;
    }

}

}
#endif

