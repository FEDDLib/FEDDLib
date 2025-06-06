#ifndef ASSEMBLEFENAVIERSTOKES_DEF_hpp
#define ASSEMBLEFENAVIERSTOKES_DEF_hpp

#include "AssembleFENavierStokes_decl.hpp"

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AssembleFENavierStokes<SC,LO,GO,NO>::AssembleFENavierStokes(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{
	int locVelocity=0;
	int locPressure=0;		
	if(std::get<0>(this->diskTuple_->at(0))=="Velocity"){
		locVelocity=0;
		locPressure=1;
	}
	else if(std::get<0>(this->diskTuple_->at(1))=="Velocity"){
		locVelocity=1;
		locPressure=0;
	}
	else
    	TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "No discretisation Information for Velocity in Navier Stokes Element." );
		

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	FETypeVelocity_ = std::get<1>(this->diskTuple_->at(locVelocity));
	FETypePressure_ =std::get<1>(this->diskTuple_->at(locPressure));

	dofsVelocity_ = std::get<2>(this->diskTuple_->at(locVelocity));
	dofsPressure_ =std::get<2>(this->diskTuple_->at(locPressure));

	numNodesVelocity_ = std::get<3>(this->diskTuple_->at(locVelocity));
	numNodesPressure_=std::get<3>(this->diskTuple_->at(locPressure));

	dofsElementVelocity_ = dofsVelocity_*numNodesVelocity_;
	dofsElementPressure_  = dofsPressure_*numNodesPressure_;	

	//this->solution_ = vec_dbl_Type(dofsElementVelocity_); //dofsElementPressure_+
	this->solutionVelocity_ = vec_dbl_Type(dofsElementVelocity_);
	this->solutionPressure_ = vec_dbl_Type(dofsElementPressure_);

 	viscosity_ = this->params_->sublist("Parameter").get("Viscosity",1.);
    density_ = this->params_->sublist("Parameter").get("Density",1.);

	dofsElement_ = dofsElementVelocity_+ dofsElementPressure_;

	SmallMatrix_Type coeff(2);
	coeff[0][0]=1.; coeff[0][1] = 1.; coeff[1][0] = 1.; coeff[1][1] = 1.; // we keep it constant like this for now. For BDF time disc. okay.
	coeff_ = coeff;

    linearization_ = this->params_->sublist("General").get("Linearization","FixedPoint"); // Information to assemble Jacobian accordingly
}

template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::setCoeff(SmallMatrix_Type coeff) {
	// We only substitute the coefficients if the matrix has the same 
	// size. In some non timedepenent cases the coeff matrix can be empty. 
	// We prevent that case.
	if(coeff.size() == 2)
		coeff_ = coeff;		

}

template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assembleJacobian() {

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));
	SmallMatrixPtr_Type elementMatrixW =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

	if(this->newtonStep_ ==0){

		SmallMatrixPtr_Type elementMatrixA =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));
		SmallMatrixPtr_Type elementMatrixB =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

		constantMatrix_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

		assemblyLaplacian(elementMatrixA);

		elementMatrixA->scale(viscosity_);
		elementMatrixA->scale(density_);

		constantMatrix_->add( (*elementMatrixA),(*constantMatrix_));

		assemblyDivAndDivT(elementMatrixB); // For Matrix B

		elementMatrixB->scale(-1.);

		constantMatrix_->add( (*elementMatrixB),(*constantMatrix_));
    }

	ANB_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_)); // A + B + N
	ANB_->add( (*constantMatrix_),(*ANB_));

	assemblyAdvection(elementMatrixN);
	elementMatrixN->scale(density_);
	ANB_->add( (*elementMatrixN),(*ANB_));
    if(linearization_ != "FixedPoint"){
	    assemblyAdvectionInU(elementMatrixW);
	    elementMatrixW->scale(density_);
    }

	//elementMatrix->add((*constantMatrix_),(*elementMatrix));
	this->jacobian_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

	this->jacobian_->add((*ANB_),(*this->jacobian_));
    // If the linearization is Newtons Method we need to add W-Matrix
    if(linearization_ != "FixedPoint"){
    	this->jacobian_->add((*elementMatrixW),(*this->jacobian_));  // int add(SmallMatrix<T> &bMat, SmallMatrix<T> &cMat); //this+B=C elementMatrix + constantMatrix_;
    }
}

template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assembleFixedPoint() {

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

	if(this->newtonStep_ ==0){
		SmallMatrixPtr_Type elementMatrixA =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));
		SmallMatrixPtr_Type elementMatrixB =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

		constantMatrix_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

		assemblyLaplacian(elementMatrixA);

		elementMatrixA->scale(viscosity_);
		elementMatrixA->scale(density_);

		constantMatrix_->add( (*elementMatrixA),(*constantMatrix_));

		assemblyDivAndDivT(elementMatrixB); // For Matrix B

		elementMatrixB->scale(-1.);

		constantMatrix_->add( (*elementMatrixB),(*constantMatrix_));
    }

	ANB_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_)); // A + B + N
	ANB_->add( (*constantMatrix_),(*ANB_));

	assemblyAdvection(elementMatrixN);
	elementMatrixN->scale(density_);
	ANB_->add( (*elementMatrixN),(*ANB_));

}


template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= numNodesVelocity_;
	std::string FEType = FETypeVelocity_;
	int dofs = dofsVelocity_;

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    // inner( grad(u) , grad(v) ) has twice the polyonimial degree than grad(u) or grad(v).
    UN deg = 2*Helper::determineDegree(dim,FEType,Helper::Grad);
    // cout << " Degree Laplace " << deg << " Grad " << Grad << " FeType " << FEType << endl;
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
  
    buildTransformation(B);

    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );
  	
    for (UN i=0; i < numNodes; i++) {
        Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
        for (UN j=0; j < numNodes; j++) {
            for (UN w=0; w<dPhiTrans.size(); w++) {
                for (UN d=0; d<dim; d++){
                    value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                }
            }
            value[j] *= absDetB;
			 /*if (std::fabs(value[j]) < std::pow(10,-14)) {
		            value[j] = 0.;
		        }*/
			for (UN d=0; d<dofs; d++) {
              (*elementMatrix)[i*dofs +d][j*dofs+d] = value[j];
            }
        }

    }
}

// Assemble RHS with updated solution coming from Fixed Point Iter or der Newton.
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assembleRHS(){

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_));

	ANB_.reset(new SmallMatrix_Type( dofsElementVelocity_+numNodesPressure_)); // A + B + N
	ANB_->add( (*constantMatrix_),(*ANB_));

	assemblyAdvection(elementMatrixN);
	elementMatrixN->scale(density_);
	ANB_->add( (*elementMatrixN),(*ANB_));

	this->rhsVec_.reset( new vec_dbl_Type ( dofsElement_,0.) );
	// Multiplying ANB_ * solution // ANB Matrix without nonlinear part.
	int s=0,t=0;
	for(int i=0 ; i< ANB_->size();i++){
		if (i >= dofsElementVelocity_)
			s=1;
		for(int j=0; j < ANB_->size(); j++){
			if(j >= dofsElementVelocity_)
				t=1;
			(*this->rhsVec_)[i] += (*ANB_)[i][j]*(*this->solution_)[j]*coeff_[s][t];
		}
		t=0;
	}
}


template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assemblyAdvection(SmallMatrixPtr_Type &elementMatrix){

	int dim = this->getDim();
	int numNodes= numNodesVelocity_;
	std::string FEType = FETypeVelocity_;
	int dofs = dofsVelocity_;

	vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type  phi;
	vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    // inner( inner( uh , nabla )phi_j , phi_i )
    // For uh: degPhi, for nabla phi_j: degGradPhi, for phi_i: degPhi.
	UN degPhi = Helper::determineDegree(dim,FEType,Helper::Std);
	UN degGradPhi = Helper::determineDegree(dim,FEType,Helper::Grad);
    UN deg = degPhi + degGradPhi + degPhi;

	Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    Helper::getPhi(phi, weights, dim, FEType, deg);

	SC detB;
	SC absDetB;
	SmallMatrix<SC> B(dim);
	SmallMatrix<SC> Binv(dim);

    vec2D_dbl_Type uLoc( dim, vec_dbl_Type( weights->size() , -1. ) );

    buildTransformation(B);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );

    for (int w=0; w<phi->size(); w++){ //quads points
        for (int d=0; d<dim; d++) {
            uLoc[d][w] = 0.;
            for (int i=0; i < phi->at(0).size(); i++) {
                LO index = dim * i + d;
                uLoc[d][w] += (*this->solution_)[index] * phi->at(w).at(i);
            }
        }

    }

    for (UN i=0; i < phi->at(0).size(); i++) {
        Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
        Teuchos::Array<GO> indices( dPhiTrans[0].size(), 0 );
        for (UN j=0; j < value.size(); j++) {
            for (UN w=0; w<dPhiTrans.size(); w++) {
                for (UN d=0; d<dim; d++){
                    value[j] += weights->at(w) * uLoc[d][w] * (*phi)[w][i] * dPhiTrans[w][j][d];
				}
            }
            value[j] *= absDetB;

            /*if (setZeros_ && std::fabs(value[j]) < myeps_) {
                value[j] = 0.;
            }*/

     
		}
 		for (UN d=0; d<dim; d++) {
    	    for (UN j=0; j < indices.size(); j++)
    	   		 (*elementMatrix)[i*dofs +d][j*dofs+d] = value[j];
			
        }
        
    }

}


template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assemblyAdvectionInU(SmallMatrixPtr_Type &elementMatrix){

	int dim = this->getDim();
	int numNodes= numNodesVelocity_;
	std::string FEType = FETypeVelocity_;
	int dofs = dofsVelocity_;


	vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type  phi;
	vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    // inner( inner( phi_j , nabla )uh , phi_i )
    // For phi_j: degPhi, for nabla uh: degGradPhi, for phi_i: degPhi.
	UN degPhi = Helper::determineDegree(dim,FEType,Helper::Std);
	UN degGradPhi = Helper::determineDegree(dim,FEType,Helper::Grad);
    UN deg = degPhi + degGradPhi + degPhi;

	Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    Helper::getPhi(phi, weights, dim, FEType, deg);

	SC detB;
	SC absDetB;
	SmallMatrix<SC> B(dim);
	SmallMatrix<SC> Binv(dim);

    vec2D_dbl_Type uLoc( dim, vec_dbl_Type( weights->size() , -1. ) );

    buildTransformation(B);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );
    //UN FEloc = checkFE(dim,FEType);


    std::vector<SmallMatrix<SC> > duLoc( weights->size(), SmallMatrix<SC>(dim) ); //for all quad points p_i each matrix is [u_x * grad Phi(p_i), u_y * grad Phi(p_i), u_z * grad Phi(p_i) (if 3D) ], duLoc[w] = [[phixx;phixy],[phiyx;phiyy]] (2D)

    for (int w=0; w<dPhiTrans.size(); w++){ //quads points
        for (int d1=0; d1<dim; d1++) {
            for (int i=0; i < dPhiTrans[0].size(); i++) {
                LO index = dim *i+ d1;
                for (int d2=0; d2<dim; d2++)
                    duLoc[w][d2][d1] += (*this->solution_)[index] * dPhiTrans[w][i][d2];
            }
        }
    }

    for (UN i=0; i < phi->at(0).size(); i++) {
        for (UN d1=0; d1<dim; d1++) {
            Teuchos::Array<SC> value( dim*phi->at(0).size(), 0. ); //These are value (W_ix,W_iy,W_iz)
            for (UN j=0; j < phi->at(0).size(); j++) {
                for (UN d2=0; d2<dim; d2++){
                    for (UN w=0; w<phi->size(); w++) {
                        value[ dim * j + d2 ] += weights->at(w) * duLoc[w][d2][d1] * (*phi)[w][i] * (*phi)[w][j];
                    }
                    value[ dim * j + d2 ] *= absDetB;
                }
            }
            for (UN j=0; j < phi->at(0).size(); j++){
                for (UN d=0; d<dofs; d++) {
         	    	(*elementMatrix)[i*dofs +d1][j*dofs+d] = value[j*dofs+d];
				}
            }
          
        }
    }
}
 


template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assemblyDivAndDivT(SmallMatrixPtr_Type &elementMatrix) {

    vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
	int dim = this->dim_;

    // psi_k * div(phi_j)
    // for pressure function psi_k: degPhi_pres, for velocity term div(phi_j): degGradPhi_vel
	UN degGradPhi_vel = Helper::determineDegree(dim,FETypeVelocity_,Helper::Grad);
    UN degPhi_pres = Helper::determineDegree(dim,FETypePressure_,Helper::Std);
    UN deg = degGradPhi_vel + degPhi_pres;

    Helper::getDPhi(dPhi, weights, dim, FETypeVelocity_, deg);

    //if (FETypePressure_=="P1-disc-global")
      //  Helper::getPhiGlobal(phi, weights, dim, FETypePressure_, deg);
    if (FETypePressure_=="P1-disc" && FETypeVelocity_=="Q2" )
        Helper::getPhi(phi, weights, dim, FETypePressure_, deg, FETypeVelocity_);
    else
        Helper::getPhi(phi, weights, dim, FETypePressure_, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);


	buildTransformation(B);
    detB = B.computeInverse(Binv);
    absDetB = std::fabs(detB);

    vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
    Helper::applyBTinv( dPhi, dPhiTrans, Binv );

    Teuchos::Array<GO> rowIndex( 1, 0 );
    Teuchos::Array<SC> value(1, 0.);


    for (UN i=0; i < phi->at(0).size(); i++) {
        if (FETypePressure_=="P0")
            rowIndex[0] = GO ( 0 );
        else
            rowIndex[0] = GO ( i );

        for (UN j=0; j < dPhiTrans[0].size(); j++) {
            for (UN d=0; d<dim; d++){
                value[0] = 0.;
                for (UN w=0; w<dPhiTrans.size(); w++)
                    value[0] += weights->at(w) * phi->at(w)[i] * dPhiTrans[w][j][d];
                value[0] *= absDetB;


				(*elementMatrix)[rowIndex[0]+dofsVelocity_*numNodesVelocity_][dofsVelocity_ * j + d] +=value[0];
				(*elementMatrix)[dofsVelocity_ * j + d][dofsVelocity_*numNodesVelocity_+rowIndex[0]] +=value[0];
            }
        }
    }
	//elementMatrix->print();
    // We compute value twice, maybe we should change this
    /*for (UN i=0; i < dPhiTrans[0].size(); i++) {

        Teuchos::Array<Teuchos::Array<SC> >valueVec( dim, Teuchos::Array<SC>( phi->at(0).size(), 0. ) );
        Teuchos::Array<GO> indices( phi->at(0).size(), 0 );
        for (UN j=0; j < valueVec[0].size(); j++) {
            for (UN w=0; w<dPhiTrans.size(); w++) {
                for (UN d=0; d<dim; d++)
                    valueVec[d][j] += weights->at(w) * phi->at(w)[j] * dPhiTrans[w][i][d];
            }
            for (UN d=0; d<dim; d++){
                valueVec[d][j] *= absDetB;
                if (setZeros_ && std::fabs(valueVec[d][j]) < myeps_) {
                    valueVec[d][j] = 0.;
                }
            }
        }

        for (UN j=0; j < indices.size(); j++){
            if (FEType2=="P0")
                indices[j] = GO ( mapping2->getGlobalElement( T ) );
            else
                indices[j] = GO ( mapping2->getGlobalElement( elements2->getElement(T).getNode(j) ) );
        }
        for (UN d=0; d<dim; d++) {
            GO row = GO ( dim * mapping1->getGlobalElement( elements1->getElement(T).getNode(i) ) + d );
            BTmat->insertGlobalValues( row, indices(), valueVec[d]() );
        }

    }*/


}


/*!

 \brief Building Transformation

@param[in] &B

*/

template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::buildTransformation(SmallMatrix<SC>& B){

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

}
#endif

