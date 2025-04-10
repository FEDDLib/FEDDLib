#ifndef AssembleFEGeneralizedNewtonian_DEF_hpp
#define AssembleFEGeneralizedNewtonian_DEF_hpp

#include "AssembleFENavierStokes_decl.hpp"

namespace FEDD
{
    // All important things are so far defined in AssembleFENavierStokes. Please check there.
    /* Interesting paper to generalized Newtonian fluids
    @article{Poole2023,
    author = {Poole, Robert J.},
    doi = {10.1016/j.jnnfm.2023.105106},
    issn = {03770257},
    journal = {Journal of Non-Newtonian Fluid Mechanics},
    keywords = {Constitutive modelling,Flow-type,Generalised Newtonian fluids,Inelastic},
    number = {August},
    pages = {105106},
    publisher = {Elsevier B.V.},
    title = {{Inelastic and flow-type parameter models for non-Newtonian fluids}},
    url = {https://doi.org/10.1016/j.jnnfm.2023.105106},
    volume = {320},
    year = {2023}
    }

    */

    /* *******************************************************************************
    This class is for Generalized-Newtonian fluids, where we consider that the viscosity is non-constant. Because the viscosity is no longer constant, the conventional formulation with the Laplacian term cannot be considered 
    (although there is a generalized Laplacian version of the equation, see "On outflow boundary conditions in finite element simulations of non-Newtonian internal flow" 2021). 
    Instead, we use the stress-divergence formulation of the momentum equation and derive from that the element-wise entrie
    ******************************************************************************* */
    template <class SC, class LO, class GO, class NO>
    AssembleFEGeneralizedNewtonian<SC, LO, GO, NO>::AssembleFEGeneralizedNewtonian(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : AssembleFENavierStokes<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple)
    {

        ////******************* If we have an analytical formula we could also just use Paraview Postprocessing tools to compute the viscosity **********************************
        dofsElementViscosity_ = this->dofsPressure_ * this->numNodesVelocity_; // So it is a scalar quantity but as it depend on the velocity it is defined at the nodes of the velocity
        this->constOutputField_ = vec_dbl_Type(dofsElementViscosity_);         ////**********************************************************************************

        // Reading through parameterlist
        shearThinningModel = params->sublist("Material").get("ShearThinningModel", "");
        // New: We have to check which material model we use
        if (shearThinningModel == "Carreau-Yasuda")
        {
            Teuchos::RCP<CarreauYasuda<SC, LO, GO, NO>> viscosityModelSpecific(new CarreauYasuda<SC, LO, GO, NO>(params));
            viscosityModel = viscosityModelSpecific;
        }
        else if (shearThinningModel == "Power-Law")
        {
            Teuchos::RCP<PowerLaw<SC, LO, GO, NO>> viscosityModelSpecific(new PowerLaw<SC, LO, GO, NO>(params));
            viscosityModel = viscosityModelSpecific;
        }
        else if (shearThinningModel == "Dimless-Carreau")
        {
            Teuchos::RCP<Dimless_Carreau<SC, LO, GO, NO>> viscosityModelSpecific(new Dimless_Carreau<SC, LO, GO, NO>(params));
            viscosityModel = viscosityModelSpecific;
        }
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "No specific implementation for your material model request. Valid are:Carreau-Yasuda, Power-Law, Dimless-Carreau");

    }

    template <class SC, class LO, class GO, class NO>
    void AssembleFEGeneralizedNewtonian<SC, LO, GO, NO>::assembleJacobian()
    {

        // For nonlinear generalized newtonian stress tensor part
        SmallMatrixPtr_Type elementMatrixN = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
        SmallMatrixPtr_Type elementMatrixW = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));

        // For nonlinear convection
        SmallMatrixPtr_Type elementMatrixNC = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
        SmallMatrixPtr_Type elementMatrixWC = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));

        // In the first iteration step we initialize the constant matrices
        // So in the case of a newtonian fluid we would have the matrix A with the contributions of the Laplacian term
        // and the matrix B with the mixed-pressure terms. Latter one exists also in the generlized-newtonian case
        if (this->newtonStep_ == 0)
        {
            SmallMatrixPtr_Type elementMatrixA = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
            SmallMatrixPtr_Type elementMatrixB = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));

            this->constantMatrix_.reset(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
            // Construct the matrix B from FE formulation - as it is equivalent to Newtonian case we call the same function
            this->assemblyDivAndDivT(elementMatrixB); // For Matrix B
            elementMatrixB->scale(-1.);
            this->constantMatrix_->add((*elementMatrixB), (*this->constantMatrix_));
        }

        // The other element matrices are not constant so we have to update them in each step
        // As the stress tensor term, considering a generalized-newtonian constitutive equation, is nonlinear we add its contribution here

        // ANB is the FixedPoint Formulation which was named for newtonian fluids.
        // Matrix A (Laplacian term (here not occuring)), Matrix B for div-Pressure Part, Matrix N for nonlinear parts -

        this->ANB_.reset(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_)); // A + B + N
        this->ANB_->add(((*this->constantMatrix_)), ((*this->ANB_)));

        // Nonlinear advection term \rho (u \cdot \nabla) u
        // As this class is derived from NavierStokes class we can call already implemented function
        //*************** ADVECTION************************
        this->assemblyAdvection(elementMatrixNC);
        elementMatrixNC->scale(this->density_);
        this->ANB_->add((*elementMatrixNC), ((*this->ANB_)));

        // For a generalized-newtonian fluid we add additional element matrix and fill it with specific contribution
        // Remember that this term is based on the stress-divergence formulation of the momentum equation
        // \nabla \dot \tau with \tau=\eta(\gammaDot)(\nabla u + (\nabla u)^T)
        //*************** STRESS TENSOR************************
        this->assemblyStress(elementMatrixN);
        this->ANB_->add((*elementMatrixN), ((*this->ANB_)));

        this->jacobian_.reset(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
        this->jacobian_->add((*this->ANB_), (*this->jacobian_));

        // If linearization is not FixdPoint (so NOX or Newton) we add the derivative to the Jacobian matrix. Otherwise the FixedPoint formulation becomes the jacobian.
        if (this->linearization_ != "FixedPoint")
        {

            this->assemblyStressDev(elementMatrixW);     // shear stress tensor
            this->assemblyAdvectionInU(elementMatrixWC); // convection
            elementMatrixWC->scale(this->density_);

            this->jacobian_->add((*elementMatrixW), (*this->jacobian_));
            this->jacobian_->add((*elementMatrixWC), (*this->jacobian_)); // int add(SmallMatrix<T> &bMat, SmallMatrix<T> &cMat); //this+B=C elementMatrix + constantMatrix_;
        }

        //*************** BOUNDARY TERM *******************************
        /* Because we have stress-divergence form of Navier-Stokes equations in the non-newtonian case
         we have to add a extra boundary term to get the same outflow boundary condition as in the conventional formulation with
         the laplacian operator in the equations due to the fact that in the stress-divergence formulation the
         natural boundary condition is different
         We have to check whether it is an element which has edges (2D) / surfaces (3D) corresponding to an Outflow Neumann boundary
         Then we have to compute contribution
         @ToDo Add if element normal computation is integrated
        */
    }



    template <class SC, class LO, class GO, class NO>
    void AssembleFEGeneralizedNewtonian<SC,LO,GO,NO>::assembleFixedPoint() {

        SmallMatrixPtr_Type elementMatrixN = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
        SmallMatrixPtr_Type elementMatrixNC = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));

	    if(this->newtonStep_ ==0){
            SmallMatrixPtr_Type elementMatrixB = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));

            this->constantMatrix_.reset(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
            // Construct the matrix B from FE formulation - as it is equivalent to Newtonian case we call the same function
            this->assemblyDivAndDivT(elementMatrixB); // For Matrix B
            elementMatrixB->scale(-1.);
            this->constantMatrix_->add((*elementMatrixB), (*this->constantMatrix_));
        }

        this->ANB_.reset(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_)); // A + B + N
        this->ANB_->add(((*this->constantMatrix_)), ((*this->ANB_)));

        // Nonlinear advection term \rho (u \cdot \nabla) u
        // As this class is derived from NavierStokes class we can call already implemented function
        //*************** ADVECTION************************
        this->assemblyAdvection(elementMatrixNC);
        elementMatrixNC->scale(this->density_);
        this->ANB_->add((*elementMatrixNC), ((*this->ANB_)));

        // For a generalized-newtonian fluid we add additional element matrix and fill it with specific contribution
        // Remember that this term is based on the stress-divergence formulation of the momentum equation
        // \nabla \dot \tau with \tau=\eta(\gammaDot)(\nabla u + (\nabla u)^T)
        //*************** STRESS TENSOR************************
        this->assemblyStress(elementMatrixN);
        this->ANB_->add((*elementMatrixN), ((*this->ANB_)));

}

    // Extra stress term resulting from chosen non-newtonian constitutive model  - Compute element matrix entries
    template <class SC, class LO, class GO, class NO>
    void AssembleFEGeneralizedNewtonian<SC, LO, GO, NO>::assemblyStress(SmallMatrixPtr_Type &elementMatrix)
    {

        int dim = this->getDim();
        int numNodes = this->numNodesVelocity_;
        UN Grad = 1; // Needs to be fixed before 2
        string FEType = this->FETypeVelocity_;
        int dofs = this->dofsVelocity_; // For pressure it would be 1

        vec3D_dbl_ptr_Type dPhi;
        vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

        UN deg = Helper::determineDegree(dim, FEType, Grad); //  e.g. for P1 3
        Helper::getDPhi(dPhi, weights, dim, FEType, deg);    //  e.g. for deg 5 we get weight vector with 7 entries
        // Example Values: dPhi->size() = 7 so number of quadrature points, dPhi->at(0).size() = 3 number of local element points, dPhi->at(0).at(0).size() = 2 as we have dim 2 therefore we have 2 derivatives (xi/eta in natural coordinates)
        // Phi is defined on reference element

        SC detB;
        SC absDetB;
        SmallMatrix<SC> B(dim);
        SmallMatrix<SC> Binv(dim);

        this->buildTransformation(B);
        detB = B.computeInverse(Binv); // The function computeInverse returns a double value corrsponding to determinant of B   
        absDetB = std::fabs(detB);     // Absolute value of B

        // dPhiTrans are the transorfmed basisfunctions, so B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1)
        // Corresponds to \hat{grad_phi}.
        vec3D_dbl_Type dPhiTrans(dPhi->size(), vec2D_dbl_Type(dPhi->at(0).size(), vec_dbl_Type(dim, 0.)));
        Helper::applyBTinv(dPhi, dPhiTrans, Binv); // dPhiTrans corresponds now to our basisfunction in natural coordinates

        TEUCHOS_TEST_FOR_EXCEPTION(dim == 1, std::logic_error, "AssemblyStress Not implemented for dim=1");
        //***************************************************************************
        if (dim == 2)
        {
            //************************************
            // Compute shear rate gammaDot, which is a vector because it is evaluated at each gaussian quadrature point
            // for that first compute velocity gradient
            vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(weights->size(), 0.0)); // gammaDot->at(j) j=0...weights
            computeShearRate(dPhiTrans, gammaDot, dim);                        // updates gammaDot using velocity solution
            //************************************
            // Compute entries
            // Initialize some helper vectors/matrices
            double v11, v12, v21, v22, value1_j, value2_j, value1_i, value2_i, viscosity_atw;
            viscosity_atw = 0.;

            // Construct element matrices
            for (UN i = 0; i < numNodes; i++)
            {
                // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3
                for (UN j = 0; j < numNodes; j++)
                {
                    // Reset values
                    v11 = 0.0;
                    v12 = 0.0;
                    v21 = 0.0;
                    v22 = 0.0;

                    // So in general compute the components of eta*[ dPhiTrans_i : ( dPhiTrans_j + (dPhiTrans_j)^T )]
                    for (UN w = 0; w < dPhiTrans.size(); w++)
                    {

                        value1_j = dPhiTrans[w][j][0]; // so this corresponds to d\phi_j/dx
                        value2_j = dPhiTrans[w][j][1]; // so this corresponds to d\phi_j/dy

                        value1_i = dPhiTrans[w][i][0]; // so this corresponds to d\phi_i/dx
                        value2_i = dPhiTrans[w][i][1]; // so this corresponds to d\phi_i/dy

                        // viscosity function evaluated where we consider the dynamic viscosity!!
                        this->viscosityModel->evaluateMapping(this->params_, gammaDot->at(w), viscosity_atw);

                        v11 = v11 + viscosity_atw * weights->at(w) * (2.0 * value1_i * value1_j + value2_i * value2_j);
                        v12 = v12 + viscosity_atw * weights->at(w) * (value2_i * value1_j);
                        v21 = v21 + viscosity_atw * weights->at(w) * (value1_i * value2_j);
                        v22 = v22 + viscosity_atw * weights->at(w) * (2.0 * value2_i * value2_j + value1_i * value1_j);

                    } // loop end quadrature points

                    // multiply determinant from transformation
                    v11 *= absDetB;
                    v12 *= absDetB;
                    v21 *= absDetB;
                    v22 *= absDetB;

                    // Put values on the right position in element matrix - d=2 because we are in two dimensional case
                    // [v11  v12  ]
                    // [v21  v22  ]
                    (*elementMatrix)[i * dofs][j * dofs] = v11; // d=0, first dimension
                    (*elementMatrix)[i * dofs][j * dofs + 1] = v12;
                    (*elementMatrix)[i * dofs + 1][j * dofs] = v21;
                    (*elementMatrix)[i * dofs + 1][j * dofs + 1] = v22; // d=1, second dimension

                } // loop end over j node

            } // loop end over i node

        } // end if dim 2
        //***************************************************************************
        else if (dim == 3)
        {
            //************************************#

            // Compute shear rate gammaDot, which is a vector because it is evaluated at a gaussian quadrature point
            // for that compute velocity gradient
            vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(weights->size(), 0.0)); // gammaDot->at(j) j=0...weights
            computeShearRate(dPhiTrans, gammaDot, dim);                        // updates gammaDot using velcoity solution

            // Initialize some helper vectors/matrices
            double v11, v12, v13, v21, v22, v23, v31, v32, v33, value1_j, value2_j, value3_j, value1_i, value2_i, value3_i, viscosity_atw;

            viscosity_atw = 0.;

            // Construct element matrices
            for (UN i = 0; i < numNodes; i++)
            {
                // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3

                for (UN j = 0; j < numNodes; j++)
                {
                    // Reset values
                    v11 = 0.0;
                    v12 = 0.0;
                    v13 = 0.0;
                    v21 = 0.0;
                    v22 = 0.0;
                    v23 = 0.0;
                    v31 = 0.0;
                    v32 = 0.0;
                    v33 = 0.0;

                    // So in general compute the components of eta*[ dPhiTrans_i : ( dPhiTrans_j + (dPhiTrans_j)^T )]
                    for (UN w = 0; w < dPhiTrans.size(); w++)
                    {

                        value1_j = dPhiTrans.at(w).at(j).at(0); // so this corresponds to d\phi_j/dx
                        value2_j = dPhiTrans.at(w).at(j).at(1); // so this corresponds to d\phi_j/dy
                        value3_j = dPhiTrans.at(w).at(j).at(2); // so this corresponds to d\phi_j/dz

                        value1_i = dPhiTrans.at(w).at(i).at(0); // so this corresponds to d\phi_i/dx
                        value2_i = dPhiTrans.at(w).at(i).at(1); // so this corresponds to d\phi_i/dy
                        value3_i = dPhiTrans.at(w).at(i).at(2); // so this corresponds to d\phi_i/dz

                        this->viscosityModel->evaluateMapping(this->params_, gammaDot->at(w), viscosity_atw);

                        // Construct entries - we go over all quadrature points and if j is updated we set v11 etc. again to zero
                        v11 = v11 + viscosity_atw * weights->at(w) * (2.0 * value1_j * value1_i + value2_j * value2_i + value3_j * value3_i);
                        v12 = v12 + viscosity_atw * weights->at(w) * (value2_i * value1_j);
                        v13 = v13 + viscosity_atw * weights->at(w) * (value3_i * value1_j);

                        v21 = v21 + viscosity_atw * weights->at(w) * (value1_i * value2_j);
                        v22 = v22 + viscosity_atw * weights->at(w) * (value1_i * value1_j + 2.0 * value2_j * value2_i + value3_j * value3_i);
                        v23 = v23 + viscosity_atw * weights->at(w) * (value3_i * value2_j);

                        v31 = v31 + viscosity_atw * weights->at(w) * (value1_i * value3_j);
                        v32 = v32 + viscosity_atw * weights->at(w) * (value2_i * value3_j);
                        v33 = v33 + viscosity_atw * weights->at(w) * (value1_i * value1_j + value2_i * value2_j + 2.0 * value3_i * value3_j);

                    } // loop end quadrature points

                    // multiply determinant from transformation
                    v11 *= absDetB;
                    v12 *= absDetB;
                    v13 *= absDetB;
                    v21 *= absDetB;
                    v22 *= absDetB;
                    v23 *= absDetB;
                    v31 *= absDetB;
                    v32 *= absDetB;
                    v33 *= absDetB;

                    // Put values on the right position in element matrix
                    // [v11  v12  v13]
                    // [v21  v22  v23]
                    // [v31  v32  v33]
                    (*elementMatrix)[i * dofs][j * dofs] = v11; // d=0, first dimension
                    (*elementMatrix)[i * dofs][j * dofs + 1] = v12;
                    (*elementMatrix)[i * dofs][j * dofs + 2] = v13;
                    (*elementMatrix)[i * dofs + 1][j * dofs] = v21;
                    (*elementMatrix)[i * dofs + 1][j * dofs + 1] = v22; // d=1, second dimension
                    (*elementMatrix)[i * dofs + 1][j * dofs + 2] = v23; // d=1, second dimension
                    (*elementMatrix)[i * dofs + 2][j * dofs] = v31;
                    (*elementMatrix)[i * dofs + 2][j * dofs + 1] = v32; // d=2, third dimension
                    (*elementMatrix)[i * dofs + 2][j * dofs + 2] = v33; // d=2, third dimension

                } // loop end over j node
            }     // loop end over i node
        }         // end if dim==3
    }

    // Directional Derivative of shear stress term resulting from chosen nonlinear non-newtonian model  -----
    // Same structure and functions as in assemblyStress
    //  ( -2.0*deta/dgammaDot * dgammaDot/dTau * (0.5(dv^k + (dvh^k)^T): 0.5( dPhiTrans_j + (dPhiTrans_j)^T))0.5(dv^k + (dvh^k)^T): 0.5( dPhiTrans_i + (dPhiTrans_i)^T)    )
    template <class SC, class LO, class GO, class NO>
    void AssembleFEGeneralizedNewtonian<SC, LO, GO, NO>::assemblyStressDev(SmallMatrixPtr_Type &elementMatrix)
    {

        int dim = this->getDim();
        int numNodes = this->numNodesVelocity_;
        UN Grad = 2; // Needs to be fixed
        string FEType = this->FETypeVelocity_;
        int dofs = this->dofsVelocity_; // for pressure it would be 1

        vec3D_dbl_ptr_Type dPhi;
        vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

        UN deg = Helper::determineDegree(dim, FEType, Grad);
        Helper::getDPhi(dPhi, weights, dim, FEType, deg);

        SC detB;
        SC absDetB;
        SmallMatrix<SC> B(dim);
        SmallMatrix<SC> Binv(dim);

        this->buildTransformation(B);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        vec3D_dbl_Type dPhiTrans(dPhi->size(), vec2D_dbl_Type(dPhi->at(0).size(), vec_dbl_Type(dim, 0.)));
        Helper::applyBTinv(dPhi, dPhiTrans, Binv);

        TEUCHOS_TEST_FOR_EXCEPTION(dim == 1, std::logic_error, "AssemblyStress Not implemented for dim=1");

        if (dim == 2)
        {
            //************************************
            //************************************
            // Due to the extra term related to the Gaetaeux-derivative there arise prefactors which depend on the velocity gradients solutions
            // which therefore also have to be computed here therefore we compute it directly here
            vec_dbl_Type u11(dPhiTrans.size(), -1.);                           // should correspond to du/dx at each quadrature point
            vec_dbl_Type u12(dPhiTrans.size(), -1.);                           // should correspond to du/dy at each quadrature point
            vec_dbl_Type u21(dPhiTrans.size(), -1.);                           // should correspond to dv/dx at each quadrature point
            vec_dbl_Type u22(dPhiTrans.size(), -1.);                           // should correspond to dv/dy at each quadrature point
            vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(weights->size(), 0.0)); // gammaDot->at(j) j=0...weights

            vec_dbl_ptr_Type mixed_term_xy(new vec_dbl_Type(weights->size(), 0.0));
            for (UN w = 0; w < dPhiTrans.size(); w++)
            { // quads points
                // set again to zero
                u11[w] = 0.0; // du_dx
                u12[w] = 0.0; // du_dy
                u21[w] = 0.0; // dv_dx
                u22[w] = 0.0; // dv_dy
                for (UN i = 0; i < dPhiTrans[0].size(); i++)
                {                            // loop unrolling
                    LO index1 = dim * i + 0; // x
                    LO index2 = dim * i + 1; // y
                    // uLoc[d][w] += (*this->solution_)[index] * phi->at(w).at(i);
                    u11[w] += (*this->solution_)[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                    u12[w] += (*this->solution_)[index1] * dPhiTrans[w][i][1]; // because we are in 2D , 0 and 1
                    u21[w] += (*this->solution_)[index2] * dPhiTrans[w][i][0];
                    u22[w] += (*this->solution_)[index2] * dPhiTrans[w][i][1];
                }
                gammaDot->at(w) = sqrt(2.0 * u11[w] * u11[w] + 2.0 * u22[w] * u22[w] + (u12[w] + u21[w]) * (u12[w] + u21[w]));
                mixed_term_xy->at(w) = 0.5 * (u12[w] + u21[w]);
            }
            //*******************************

            // Initialize some helper vectors/matrices
            double v11, v12, v21, v22, value1_j, value2_j, value1_i, value2_i, deta_dgamma_dgamma_dtau;

            deta_dgamma_dgamma_dtau = 0.;

            // Construct element matrices
            for (UN i = 0; i < numNodes; i++)
            {
                // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3

                for (UN j = 0; j < numNodes; j++)
                {
                    // Reset values
                    v11 = 0.0;
                    v12 = 0.0;
                    v21 = 0.0;
                    v22 = 0.0;

                    // Only the part  deta/dgammaDot is different for all shear thinning models (because we make the assumption of incompressibility)? NO ALSO DIFFERENT FOR EXAMPLE FOR CASSON SO CONSIDERING YIELD STRESS
                    // but we put the two terms together because then we can multiply them together and get e.g. for carreau yasuda  : gammaDot^{a-2.0} which is for a=2.0 equals 0 and we do not have to worry about the problem what if gammaDot = 0.0
                    for (UN w = 0; w < dPhiTrans.size(); w++)
                    {

                        value1_j = dPhiTrans[w][j][0]; // so this corresponds to d\phi_j/dx
                        value2_j = dPhiTrans[w][j][1]; // so this corresponds to d\phi_j/dy

                        value1_i = dPhiTrans[w][i][0]; // so this corresponds to d\phi_i/dx
                        value2_i = dPhiTrans[w][i][1]; // so this corresponds to d\phi_i/dy

                        this->viscosityModel->evaluateDerivative(this->params_, gammaDot->at(w), deta_dgamma_dgamma_dtau);
                        /* EInfacher in unausmultiplizierter Form
                        v11 = v11 + (-2.0)*deta_dgamma_dgamma_dtau  * weights->at(w) *(u11[w]*u11[w]*value1_i*value1_j+u11[w]*mixed_terms->at(w)*(value1_i*value2_j+value2_i*value1_j)+ mixed_terms->at(w)*mixed_terms->at(w)*(value2_i*value2_j)); // xx contribution: (dv_x/dx)^2*dphi_i/dx*dphi_j/dx+dv_x/dx*f*(dphi_i/dx*dphi_j/dy+dphi_i/dy*dphi_j/dx)+f^2*dphi_i/dy*dphi_j/dy
                        v12 = v12 + (-2.0)*deta_dgamma_dgamma_dtau  * weights->at(w) *(u11[w]*mixed_terms->at(w)*value1_i*value1_j+u11[w]*u22[w]*(value1_i*value2_j)        +mixed_terms->at(w)*mixed_terms->at(w)*value2_i*value1_j+mixed_terms->at(w)*u22[w]*value2_i*value2_j ); // xy contribution:  dv_x/dx*f*dphi_i/dx*dphi_j/dx+dv_x/dx*dv_y/dy*dphi_i/dx*dphi_j/dy+f^2*dphi_i_dy*dphi_j/dx+f*dv_y/dy*dphi_i/dy*dphi_j/dy
                        v21 = v21 + (-2.0)*deta_dgamma_dgamma_dtau  * weights->at(w) *(u11[w]*mixed_terms->at(w)*value1_i*value1_j+mixed_terms->at(w)*mixed_terms->at(w)*value1_i*value2_j+u11[w]*u22[w]*value2_i*value1_j          +mixed_terms->at(w)*u22[w]*value2_i*value2_j ); // yx contribution:  dv_x/dx*f*dphi_i/dx*dphi_j/dx+dv_x/dx*dv_y/dy*dphi_i/dy*dphi_j/dx+f^2*dphi_i_dx*dphi_j/dy+f*dv_y/dy*dphi_i/dy*dphi_j/dy
                        v22 = v22 + (-2.0)*deta_dgamma_dgamma_dtau  * weights->at(w) *(u22[w]*u22[w]*value2_i*value2_j+u22[w]*mixed_terms->at(w)*(value1_i*value2_j+value2_i*value1_j)+ mixed_terms->at(w)*mixed_terms->at(w)*(value1_i*value1_j) ); // yy contribution: (dv_y/dy)^2*dphi_i/dy*dphi_j/dy+dv_y/dy*f*(dphi_i/dx*dphi_j/dy+dphi_i/dy*dphi_j/dx)+f^2*dphi_i/dx*dphi_j/dx
                        */
                        v11 = v11 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * u11[w] + value2_j * mixed_term_xy->at(w)) * (value1_i * u11[w] + value2_i * mixed_term_xy->at(w)));
                        v12 = v12 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * mixed_term_xy->at(w) + u22[w] * value2_j) * (value1_i * u11[w] + value2_i * mixed_term_xy->at(w)));

                        v21 = v21 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * u11[w] + value2_j * mixed_term_xy->at(w)) * (value1_i * mixed_term_xy->at(w) + value2_i * u22[w]));
                        v22 = v22 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * mixed_term_xy->at(w) + u22[w] * value2_j) * (value1_i * mixed_term_xy->at(w) + value2_i * u22[w]));

                    } // loop end quadrature points

                    // multiply determinant from transformation
                    v11 *= absDetB;
                    v12 *= absDetB;
                    v21 *= absDetB;
                    v22 *= absDetB;

                    // Put values on the right position in element matrix - d=2 because we are in two dimensional case
                    // [v11  v12  ]
                    // [v21  v22  ]
                    (*elementMatrix)[i * dofs][j * dofs] = v11; // d=0, first dimension
                    (*elementMatrix)[i * dofs][j * dofs + 1] = v12;
                    (*elementMatrix)[i * dofs + 1][j * dofs] = v21;
                    (*elementMatrix)[i * dofs + 1][j * dofs + 1] = v22; // d=1, second dimension

                } // loop end over j node

            } // loop end over i node

        } // end if dim 2
        //***************************************************************************
        //***************************************************************************
        else if (dim == 3)
        {
            //************************************
            // Compute shear rate gammaDot, which is a vector because it is evaluated at a gaussian quadrature point
            // for that compute velocity gradient
            vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
            vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
            vec_dbl_Type u13(dPhiTrans.size(), -1.); // should correspond to du/dz at each quadrature point

            vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
            vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point
            vec_dbl_Type u23(dPhiTrans.size(), -1.); // should correspond to dv/dz at each quadrature point

            vec_dbl_Type u31(dPhiTrans.size(), -1.); // should correspond to dw/dx at each quadrature point
            vec_dbl_Type u32(dPhiTrans.size(), -1.); // should correspond to dw/dy at each quadrature point
            vec_dbl_Type u33(dPhiTrans.size(), -1.); // should correspond to dw/dz at each quadrature point

            vec_dbl_ptr_Type gammaDot(new vec_dbl_Type(weights->size(), 0.0)); // gammaDot->at(j) j=0...weights

            vec_dbl_ptr_Type mixed_term_xy(new vec_dbl_Type(weights->size(), 0.0));
            vec_dbl_ptr_Type mixed_term_xz(new vec_dbl_Type(weights->size(), 0.0));
            vec_dbl_ptr_Type mixed_term_yz(new vec_dbl_Type(weights->size(), 0.0));

            for (UN w = 0; w < dPhiTrans.size(); w++)
            { // quads points

                u11[w] = 0.0;
                u12[w] = 0.0;
                u13[w] = 0.0;
                u21[w] = 0.0;
                u22[w] = 0.0;
                u23[w] = 0.0;
                u31[w] = 0.0;
                u32[w] = 0.0;
                u33[w] = 0.0;

                for (UN i = 0; i < dPhiTrans[0].size(); i++)
                {
                    LO index1 = dim * i + 0; // x
                    LO index2 = dim * i + 1; // y
                    LO index3 = dim * i + 2; // z
                    // uLoc[d][w] += (*this->solution_)[index] * phi->at(w).at(i);
                    u11[w] += (*this->solution_)[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                    u12[w] += (*this->solution_)[index1] * dPhiTrans[w][i][1]; // because we are in 3D , 0 and 1, 2
                    u13[w] += (*this->solution_)[index1] * dPhiTrans[w][i][2];
                    u21[w] += (*this->solution_)[index2] * dPhiTrans[w][i][0]; // v*dphi_dx
                    u22[w] += (*this->solution_)[index2] * dPhiTrans[w][i][1];
                    u23[w] += (*this->solution_)[index2] * dPhiTrans[w][i][2];
                    u31[w] += (*this->solution_)[index3] * dPhiTrans[w][i][0]; // w*dphi_dx
                    u32[w] += (*this->solution_)[index3] * dPhiTrans[w][i][1];
                    u33[w] += (*this->solution_)[index3] * dPhiTrans[w][i][2];
                }
                gammaDot->at(w) = sqrt(2.0 * u11[w] * u11[w] + 2.0 * u22[w] * u22[w] + 2.0 * u33[w] * u33[w] + (u12[w] + u21[w]) * (u12[w] + u21[w]) + (u13[w] + u31[w]) * (u13[w] + u31[w]) + (u23[w] + u32[w]) * (u23[w] + u32[w]));

                mixed_term_xy->at(w) = 0.5 * (u12[w] + u21[w]);
                mixed_term_xz->at(w) = 0.5 * (u31[w] + u13[w]);
                mixed_term_yz->at(w) = 0.5 * (u32[w] + u23[w]);
            }

            // Initialize some helper vectors/matrices
            double v11, v12, v13, v21, v22, v23, v31, v32, v33, value1_j, value2_j, value3_j, value1_i, value2_i, value3_i, deta_dgamma_dgamma_dtau;

            deta_dgamma_dgamma_dtau = 0.;

            // Construct element matrices
            for (UN i = 0; i < numNodes; i++)
            {
                // Teuchos::Array<SC> value(dPhiTrans[0].size(), 0. ); // dPhiTrans[0].size() is 3

                for (UN j = 0; j < numNodes; j++)
                {
                    // Reset values
                    v11 = 0.0;
                    v12 = 0.0;
                    v13 = 0.0;
                    v21 = 0.0;
                    v22 = 0.0;
                    v23 = 0.0;
                    v31 = 0.0;
                    v32 = 0.0;
                    v33 = 0.0;

                    // So in general compute the components of eta*[ dPhiTrans_i : ( dPhiTrans_j + (dPhiTrans_j)^T )]
                    for (UN w = 0; w < dPhiTrans.size(); w++)
                    {

                        value1_j = dPhiTrans.at(w).at(j).at(0); // so this corresponds to d\phi_j/dx
                        value2_j = dPhiTrans.at(w).at(j).at(1); // so this corresponds to d\phi_j/dy
                        value3_j = dPhiTrans.at(w).at(j).at(2); // so this corresponds to d\phi_j/dz

                        value1_i = dPhiTrans.at(w).at(i).at(0); // so this corresponds to d\phi_i/dx
                        value2_i = dPhiTrans.at(w).at(i).at(1); // so this corresponds to d\phi_i/dy
                        value3_i = dPhiTrans.at(w).at(i).at(2); // so this corresponds to d\phi_i/dz

                        this->viscosityModel->evaluateDerivative(this->params_, gammaDot->at(w), deta_dgamma_dgamma_dtau);

                        // Construct entries - we go over all quadrature points and if j is updated we set v11 etc. again to zero
                        v11 = v11 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * u11[w] + value2_j * mixed_term_xy->at(w) + value3_j * mixed_term_xz->at(w)) * (value1_i * u11[w] + value2_i * mixed_term_xy->at(w) + value3_i * mixed_term_xz->at(w)));
                        v12 = v12 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * mixed_term_xy->at(w) + u22[w] * value2_j + mixed_term_yz->at(w) * value3_j) * (value1_i * u11[w] + value2_i * mixed_term_xy->at(w) + value3_i * mixed_term_xz->at(w)));
                        v13 = v13 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value3_j * mixed_term_xz->at(w) + value2_j * mixed_term_yz->at(w) + value3_j * u33[w]) * (value1_i * u11[w] + value2_i * mixed_term_xy->at(w) + value3_i * mixed_term_xz->at(w)));

                        v21 = v21 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * u11[w] + value2_j * mixed_term_xy->at(w) + value3_j * mixed_term_xz->at(w)) * (value1_i * mixed_term_xy->at(w) + value2_i * u22[w] + value3_i * mixed_term_yz->at(w)));
                        v22 = v22 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * mixed_term_xy->at(w) + u22[w] * value2_j + mixed_term_yz->at(w) * value3_j) * (value1_i * mixed_term_xy->at(w) + value2_i * u22[w] + value3_i * mixed_term_yz->at(w)));
                        v23 = v23 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value3_j * mixed_term_xz->at(w) + value2_j * mixed_term_yz->at(w) + value3_j * u33[w]) * (value1_i * mixed_term_xy->at(w) + value2_i * u22[w] + value3_i * mixed_term_yz->at(w)));

                        v31 = v31 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * u11[w] + value2_j * mixed_term_xy->at(w) + value3_j * mixed_term_xz->at(w)) * (value1_i * mixed_term_xz->at(w) + mixed_term_yz->at(w) * value2_i + u33[w] * value3_i));
                        v32 = v32 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value1_j * mixed_term_xy->at(w) + u22[w] * value2_j + mixed_term_yz->at(w) * value3_j) * (value1_i * mixed_term_xz->at(w) + mixed_term_yz->at(w) * value2_i + u33[w] * value3_i));
                        v33 = v33 + (-2.0) * deta_dgamma_dgamma_dtau * weights->at(w) * ((value3_j * mixed_term_xz->at(w) + value2_j * mixed_term_yz->at(w) + value3_j * u33[w]) * (value1_i * mixed_term_xz->at(w) + mixed_term_yz->at(w) * value2_i + u33[w] * value3_i));

                    } // loop end quadrature points

                    // multiply determinant from transformation
                    v11 *= absDetB;
                    v12 *= absDetB;
                    v13 *= absDetB;
                    v21 *= absDetB;
                    v22 *= absDetB;
                    v23 *= absDetB;
                    v31 *= absDetB;
                    v32 *= absDetB;
                    v33 *= absDetB;

                    // Put values on the right position in element matrix
                    // [v11  v12  v13]
                    // [v21  v22  v23]
                    // [v31  v32  v33]
                    (*elementMatrix)[i * dofs][j * dofs] = v11; // d=0, first dimension
                    (*elementMatrix)[i * dofs][j * dofs + 1] = v12;
                    (*elementMatrix)[i * dofs][j * dofs + 2] = v13;
                    (*elementMatrix)[i * dofs + 1][j * dofs] = v21;
                    (*elementMatrix)[i * dofs + 1][j * dofs + 1] = v22; // d=1, second dimension
                    (*elementMatrix)[i * dofs + 1][j * dofs + 2] = v23; // d=1, second dimension
                    (*elementMatrix)[i * dofs + 2][j * dofs] = v31;
                    (*elementMatrix)[i * dofs + 2][j * dofs + 1] = v32; // d=2, third dimension
                    (*elementMatrix)[i * dofs + 2][j * dofs + 2] = v33; // d=2, third dimension

                } // loop end over j node

            } // loop end over i node

        } // end if dim = 3
    }

    
    //ToDo: Using new functionalty integrate again assembly of Neumann boundary term

    // "Fixpunkt"- Matrix without jacobian for calculating Ax
    // Here update please to unlinearized System Matrix accordingly.
    template <class SC, class LO, class GO, class NO>
    void AssembleFEGeneralizedNewtonian<SC, LO, GO, NO>::assembleRHS()
    {

        SmallMatrixPtr_Type elementMatrixN = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
        SmallMatrixPtr_Type elementMatrixNC = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));

        this->ANB_.reset(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_)); // A + B + N
        this->ANB_->add((*this->constantMatrix_), (*this->ANB_));

        // Nonlinear shear stress tensor *******************************
        this->assemblyStress(elementMatrixN);
        this->ANB_->add((*elementMatrixN), (*this->ANB_));

        // Nonlinear convection term *******************************
        this->assemblyAdvection(elementMatrixNC);
        elementMatrixNC->scale(this->density_);
        this->ANB_->add((*elementMatrixNC), (*this->ANB_));

        // @ToDo If underlying element is an outflow boundary element - will be added when func nonlinear boundar term *******************************
        /*if (this->surfaceElement == true)
        {
            SmallMatrixPtr_Type elementMatrixNB = Teuchos::rcp(new SmallMatrix_Type(this->dofsElementVelocity_ + this->numNodesPressure_));
            this->assemblyNeumannBoundaryTerm(elementMatrixNB);
            this->ANB_->add((*elementMatrixNB), ((*this->ANB_)));
        }
        */

        this->rhsVec_.reset(new vec_dbl_Type(this->dofsElement_, 0.));
        // Multiplying ANB_ * solution // ANB Matrix without nonlinear part.
        int s = 0, t = 0;
        for (int i = 0; i < this->ANB_->size(); i++)
        {
            if (i >= this->dofsElementVelocity_)
                s = 1;
            for (int j = 0; j < this->ANB_->size(); j++)
            {
                if (j >= this->dofsElementVelocity_)
                    t = 1;
                (*this->rhsVec_)[i] += (*this->ANB_)[i][j] * (*this->solution_)[j] * this->coeff_[s][t];
            }
            t = 0;
        }
    }

    /*!
    Building Transformation
    */
    /* In 2D B=[ nodesRefConfig(1).at(0)-nodesRefConfig(0).at(0)    nodesRefConfig(2).at(0)-nodesRefConfig(0).at(0)    ]
               [ nodesRefConfig(1).at(1)-nodesRefConfig(0).at(1)    nodesRefConfig(2).at(1)-nodesRefConfig(0).at(1)   ]
    /*!
        - Triangle numbering

                        2
                        * *
                        *   *
                        5     4
                        *       *
                        *         *
                        0 * * 3 * * 1

        In 2D B=[ 1(x)-0(x)    2(x)-0(x)  ]
                [ 1(y)-0(y)    2(y)-0(y)  ]
    */


    // Compute Shear Rate on quadrature points depending on gradient of velocity solution at nodes
    template <class SC, class LO, class GO, class NO>
    void AssembleFEGeneralizedNewtonian<SC, LO, GO, NO>::computeShearRate(vec3D_dbl_Type dPhiTrans,
                                                                          vec_dbl_ptr_Type &gammaDot, int dim)
    {

        //****************** TWO DIMENSIONAL *********************************
        if (dim == 2)
        {

            vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
            vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
            vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
            vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point

            for (UN w = 0; w < dPhiTrans.size(); w++)
            { // quads points
                // set again to zero
                u11[w] = 0.0;
                u12[w] = 0.0;
                u21[w] = 0.0;
                u22[w] = 0.0;
                for (UN i = 0; i < dPhiTrans[0].size(); i++)
                {                            // loop unrolling
                    LO index1 = dim * i + 0; // x
                    LO index2 = dim * i + 1; // y
                    // uLoc[d][w] += (*this->solution_)[index] * phi->at(w).at(i);
                    u11[w] += (*this->solution_)[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                    u12[w] += (*this->solution_)[index1] * dPhiTrans[w][i][1]; // because we are in 2D , 0 and 1
                    u21[w] += (*this->solution_)[index2] * dPhiTrans[w][i][0];
                    u22[w] += (*this->solution_)[index2] * dPhiTrans[w][i][1];
                }
                gammaDot->at(w) = sqrt(2.0 * u11[w] * u11[w] + 2.0 * u22[w] * u22[w] + (u12[w] + u21[w]) * (u12[w] + u21[w])); 
            }
        } // end if dim == 2
        //****************** THREE DIMENSIONAL *********************************
        else if (dim == 3)
        {

            vec_dbl_Type u11(dPhiTrans.size(), -1.); // should correspond to du/dx at each quadrature point
            vec_dbl_Type u12(dPhiTrans.size(), -1.); // should correspond to du/dy at each quadrature point
            vec_dbl_Type u13(dPhiTrans.size(), -1.); // should correspond to du/dz at each quadrature point

            vec_dbl_Type u21(dPhiTrans.size(), -1.); // should correspond to dv/dx at each quadrature point
            vec_dbl_Type u22(dPhiTrans.size(), -1.); // should correspond to dv/dy at each quadrature point
            vec_dbl_Type u23(dPhiTrans.size(), -1.); // should correspond to dv/dz at each quadrature point

            vec_dbl_Type u31(dPhiTrans.size(), -1.); // should correspond to dw/dx at each quadrature point
            vec_dbl_Type u32(dPhiTrans.size(), -1.); // should correspond to dw/dy at each quadrature point
            vec_dbl_Type u33(dPhiTrans.size(), -1.); // should correspond to dw/dz at each quadrature point

            for (UN w = 0; w < dPhiTrans.size(); w++)
            { // quads points
                // set again to zero
                u11[w] = 0.0;
                u12[w] = 0.0;
                u13[w] = 0.0;
                u21[w] = 0.0;
                u22[w] = 0.0;
                u23[w] = 0.0;
                u31[w] = 0.0;
                u32[w] = 0.0;
                u33[w] = 0.0;

                for (UN i = 0; i < dPhiTrans[0].size(); i++)
                {
                    LO index1 = dim * i + 0; // x
                    LO index2 = dim * i + 1; // y
                    LO index3 = dim * i + 2; // z
                    // uLoc[d][w] += (*this->solution_)[index] * phi->at(w).at(i);
                    u11[w] += (*this->solution_)[index1] * dPhiTrans[w][i][0]; // u*dphi_dx
                    u12[w] += (*this->solution_)[index1] * dPhiTrans[w][i][1]; // because we are in 3D , 0 and 1, 2
                    u13[w] += (*this->solution_)[index1] * dPhiTrans[w][i][2];
                    u21[w] += (*this->solution_)[index2] * dPhiTrans[w][i][0]; // v*dphi_dx
                    u22[w] += (*this->solution_)[index2] * dPhiTrans[w][i][1];
                    u23[w] += (*this->solution_)[index2] * dPhiTrans[w][i][2];
                    u31[w] += (*this->solution_)[index3] * dPhiTrans[w][i][0]; // w*dphi_dx
                    u32[w] += (*this->solution_)[index3] * dPhiTrans[w][i][1];
                    u33[w] += (*this->solution_)[index3] * dPhiTrans[w][i][2];
                }
                gammaDot->at(w) = sqrt(2.0 * u11[w] * u11[w] + 2.0 * u22[w] * u22[w] + 2.0 * u33[w] * u33[w] + (u12[w] + u21[w]) * (u12[w] + u21[w]) + (u13[w] + u31[w]) * (u13[w] + u31[w]) + (u23[w] + u32[w]) * (u23[w] + u32[w]));
            }
        } // end if dim == 3
    }

    /* Based on the current solution (velocity, pressure etc.) we want to be able to compute postprocessing fields
    like here the viscosity inside an element.
    */
    template <class SC, class LO, class GO, class NO>
    void AssembleFEGeneralizedNewtonian<SC, LO, GO, NO>::computeLocalconstOutputField()
    {
        int dim = this->getDim();
        string FEType = this->FETypeVelocity_;

        SC detB;
        SmallMatrix<SC> B(dim);
        SmallMatrix<SC> Binv(dim);

        this->buildTransformation(B);
        detB = B.computeInverse(Binv);

        vec3D_dbl_ptr_Type dPhiAtCM;

        // Compute viscosity at center of mass using nodal values and shape function **********************************************************************************
        TEUCHOS_TEST_FOR_EXCEPTION(dim == 1, std::logic_error, "computeLocalconstOutputField Not implemented for dim=1");

        Helper::getDPhiAtCM(dPhiAtCM, dim, FEType); // These are the original coordinates of the reference element
        vec3D_dbl_Type dPhiTransAtCM(dPhiAtCM->size(), vec2D_dbl_Type(dPhiAtCM->at(0).size(), vec_dbl_Type(dim, 0.)));
        Helper::applyBTinv(dPhiAtCM, dPhiTransAtCM, Binv); // We need transformation because of velocity gradient in shear rate equation

        vec_dbl_ptr_Type gammaDoti(new vec_dbl_Type(dPhiAtCM->size(), 0.0)); // Only one value because size is one
        computeShearRate(dPhiTransAtCM, gammaDoti, dim);                     // updates gammaDot using velocity solution
        this->viscosityModel->evaluateMapping(this->params_, gammaDoti->at(0), this->constOutputField_.at(0));
    }

}
#endif
