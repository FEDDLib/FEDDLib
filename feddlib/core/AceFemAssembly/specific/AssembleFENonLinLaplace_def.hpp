#ifndef ASSEMBLEFENONLINLAPLACE_DEF_hpp
#define ASSEMBLEFENONLINLAPLACE_DEF_hpp

#include "AssembleFENonLinLaplace_decl.hpp"
#include "feddlib/core/FEDDCore.hpp"
#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraitsDecl.hpp>
#include <cmath>

namespace FEDD {

void rhsFunc(double *x, double *res, double *parameters) {
    res[0] = 1; // x[0] * sin(x[1]);
}

/*!
 \brief Constructor for AssembleFENonLinLaplace

@param[in] flag Flag of element
@param[in] nodesRefConfig Nodes of element in reference configuration
@param[in] params Parameterlist for current problem
@param[in] tuple Information on the discretization
*/
template <class SC, class LO, class GO, class NO>
AssembleFENonLinLaplace<SC, LO, GO, NO>::AssembleFENonLinLaplace(
    int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,
    tuple_disk_vec_ptr_Type tuple)
    : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple) {

    this->FEType_ = std::get<1>(this->diskTuple_->at(0));
    this->dofs_ = std::get<2>(this->diskTuple_->at(0));
    // Same as this->getNodesRefConfig().size();
    this->numNodes_ = std::get<3>(this->diskTuple_->at(0));
    this->dofsElement_ = this->numNodes_ * this->dofs_;
    // Specifying rhs func here for now
    // Should be done in main when possible
    this->rhsFunc_ = rhsFunc;
}

/*!
 \brief Assemble Jacobian as per Gateaux derivative of the laplacian
@param[in] &elementMatrix
*/

template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::assembleJacobian() {

    SmallMatrixPtr_Type elementMatrix =
        Teuchos::rcp(new SmallMatrix_Type(this->dofsElement_));
    assemblyNonLinLaplacian(elementMatrix);

    this->jacobian_ = elementMatrix;
}

/*!
 \brief Assembly function for the Jacobian of the nonlinear Laplacian
 \int((1+u_0^2)\nabla\delta u+2u_0\delta u\nabla u_0)\nabla vdx

@param[in] &elementMatrix
*/
template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::assemblyNonLinLaplacian(
    SmallMatrixPtr_Type &elementMatrix) {

    int dim = this->getDim();
    UN deg = Helper::determineDegree(dim, this->FEType_, this->FEType_,
                                      Helper::Grad, Helper::Grad);

    vec3D_dbl_ptr_Type dPhi;
    vec2D_dbl_ptr_Type phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    // Get values of nodal basis at the quadrature nodes
    Helper::getDPhi(dPhi, weights, dim, this->FEType_, deg);
    Helper::getPhi(phi, weights, dim, this->FEType_, deg);

    // Built affine transformation from the reference element to the current
    // element
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    buildTransformation(B);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);
    vec3D_dbl_Type dPhiTrans(
        dPhi->size(),
        vec2D_dbl_Type(dPhi->at(0).size(), vec_dbl_Type(dim, 0.)));
    // Apple transformation to gradients of basis functions
    applyBTinv(dPhi, dPhiTrans, Binv);

    vec_dbl_Type uLoc(weights->size(), 0.);
    // At each quad node p_i gradient vector is duLoc[w] = [\partial_x phi,
    // \partial_y phi]
    vec2D_dbl_Type duLoc(weights->size(), vec_dbl_Type(dim, 0));

    // Build vector of current solution at quadrature nodes
    // Iterate over quadrature nodes
    for (int w = 0; w < phi->size(); w++) {
        // Iterate over each basis function at the current quadrature node
        for (int i = 0; i < phi->at(0).size(); i++) {
            uLoc[w] += (*this->solution_)[i] * phi->at(w).at(i);
            for (int d = 0; d < dim; d++) {
                duLoc[w][d] += (*this->solution_)[i] * dPhiTrans[w][i][d];
            }
        }
    }

    // Build local stiffness matrix
    auto value = Teuchos::ScalarTraits<SC>::zero();
    for (UN i = 0; i < this->numNodes_; i++) {
        for (UN j = 0; j < this->numNodes_; j++) {
            value = 0.;
            for (UN w = 0; w < weights->size(); w++) {
                for (UN d = 0; d < dim; d++) {
                    value += weights->at(w) *
                             ((1 + uLoc[w] * uLoc[w]) * dPhiTrans[w][i][d] +
                              2 * uLoc[w] * phi->at(w).at(j) * duLoc[w][d]) *
                             dPhiTrans[w][j][d];
                }
            }
            value *= absDetB;
            (*elementMatrix)[i][j] = value;
        }
    }
}

/*!
 \brief Assembly function for the residual of the nonlinear Laplacian
 \int fv - (1+u_0^2)\nabla u_0\nabla v dx

@param[in] &elementVector
*/
template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::assembleRHS() {

    int dim = this->getDim();
    UN deg = Helper::determineDegree(dim, this->FEType_, this->FEType_,
                                      Helper::Grad, Helper::Grad);

    vec3D_dbl_ptr_Type dPhi;
    vec2D_dbl_ptr_Type phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    Helper::getDPhi(dPhi, weights, dim, this->FEType_, deg);
    Helper::getPhi(phi, weights, dim, this->FEType_, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    buildTransformation(B);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);
    vec3D_dbl_Type dPhiTrans(
        weights->size(),
        vec2D_dbl_Type(dPhi->at(0).size(), vec_dbl_Type(dim, 0.)));
    applyBTinv(dPhi, dPhiTrans, Binv);
    vec_dbl_Type uLoc(weights->size(), 0.);
    vec_dbl_ptr_Type funcValues;
    Helper::getFuncAtQuadNodes(funcValues, this->rhsFunc_, dim, this->FEType_,
                               deg);

    // Build vector of current solution at quadrature nodes
    for (int w = 0; w < phi->size(); w++) { // quadrature nodes
        for (int i = 0; i < phi->at(0).size();
             i++) { // each basis function at the current quadrature node
            uLoc[w] += (*this->solution_)[i] * phi->at(w).at(i);
        }
    }

    this->rhsVec_.reset(new vec_dbl_Type(this->dofsElement_, 0.));
    // $fv$ term
    auto valueOne = Teuchos::ScalarTraits<SC>::zero();
    // $(1 + u_0^2)\nabla u_0 \nabla v$ term
    auto valueTwo = Teuchos::ScalarTraits<SC>::zero();
    // for storing matrix_column \cdot current_solution reduction operation
    auto reductionValue = Teuchos::ScalarTraits<SC>::zero();
    //  Build local stiffness matrx
    for (UN i = 0; i < this->numNodes_; i++) {
        reductionValue = 0.;
        // Iterate test functions (rows)
        for (UN j = 0; j < this->numNodes_; j++) {
            valueTwo = 0.;
            // Iterate quadrature nodes
            for (UN w = 0; w < dPhiTrans.size(); w++) {
                for (UN d = 0; d < dim; d++) {
                    valueTwo += weights->at(w) * (1 + uLoc[w] * uLoc[w]) *
                                dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                }
            }
            // Perform reduction of current valueTwo and solution
            reductionValue += valueTwo * (*this->solution_)[j];
        }
        reductionValue *= absDetB;
        valueOne = 0.;
        for (UN w = 0; w < dPhiTrans.size(); w++) {
            // Note that phi does not require transformation. absDetB suffices.
            valueOne += weights->at(w) * funcValues->at(w) * phi->at(w).at(i);
        }
        valueOne *= absDetB;
        (*this->rhsVec_)[i] = -(valueOne - reductionValue);
    }
}
/*!
 \brief Building Transformation
TODO This function could be outsourced to helper?
@param[in] &B The transformation matrix to be constructed
*/

template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::buildTransformation(
    SmallMatrix<SC> &B) {

    TEUCHOS_TEST_FOR_EXCEPTION((B.size() < 2 || B.size() > 3), std::logic_error,
                               "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = 0;
    for (UN j = 0; j < B.size(); j++) {
        index = j + 1;
        for (UN i = 0; i < B.size(); i++) {
            // name nodesRefConfig_ is deceptive: actually holds coords of nodes
            // of current element
            B[i][j] = this->nodesRefConfig_.at(index).at(i) -
                      this->nodesRefConfig_.at(index0).at(i);
        }
    }
}

/*!
 \brief performs matrix multiplication of two small matrices Binv * dPhiIn
TODO This function could be outsourced to helper? Does not need to be in each
assemble class
 @param[in] dPhiIn
 */
template <class SC, class LO, class GO, class NO>
void AssembleFENonLinLaplace<SC, LO, GO, NO>::applyBTinv(
    vec3D_dbl_ptr_Type &dPhiIn, vec3D_dbl_Type &dPhiOut,
    SmallMatrix<SC> &Binv) {
    UN dim = Binv.size();
    for (UN w = 0; w < dPhiIn->size(); w++) {
        for (UN i = 0; i < dPhiIn->at(w).size(); i++) {
            for (UN d1 = 0; d1 < dim; d1++) {
                for (UN d2 = 0; d2 < dim; d2++) {
                    dPhiOut[w][i][d1] +=
                        dPhiIn->at(w).at(i).at(d2) * Binv[d2][d1];
                }
            }
        }
    }
}

} // namespace FEDD
#endif
