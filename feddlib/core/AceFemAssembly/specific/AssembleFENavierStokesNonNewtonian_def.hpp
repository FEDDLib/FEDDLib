#ifndef ASSEMBLEFENAVIERSTOKESNONNEWTONIAN_DEF_hpp
#define ASSEMBLEFENAVIERSTOKESNONNEWTONIAN_DEF_hpp

#include "AssembleFENavierStokes_decl.hpp"

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::AssembleFENavierStokesNonNewtonian(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFENavierStokes<SC,LO,GO,NO>(flag, nodesRefConfig, params,tuple)
{
	// All important things are so far defined in AssembleFENavierStokes. Please check there.
}


template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assembleJacobian() {

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
	SmallMatrixPtr_Type elementMatrixW =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

	if(this->newtonStep_ ==0){
		SmallMatrixPtr_Type elementMatrixA =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));
		SmallMatrixPtr_Type elementMatrixB =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

		this->constantMatrix_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

		assemblyLaplacian(elementMatrixA); // Here inser new routine

		elementMatrixA->scale(this->viscosity_);
		elementMatrixA->scale(this->density_);

		this->constantMatrix_->add( (*elementMatrixA),(*this->constantMatrix_));

		this->assemblyDivAndDivT(elementMatrixB); // For Matrix B

		elementMatrixB->scale(-1.);

		this->constantMatrix_->add( (*elementMatrixB),(*this->constantMatrix_));
    }

    // ANB is the FixedPoint Formulation. Matrix A + N for advection part and B for div-Pressure Part.
	this->ANB_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_)); // A + B + N
	this->ANB_->add( ((*this->constantMatrix_)),((*this->ANB_)));

	this->assemblyAdvection(elementMatrixN);
	elementMatrixN->scale(this->density_);
	this->ANB_->add( (*elementMatrixN),((*this->ANB_)));

    // If linearization is not FixdPoint (so NOX or Newton) we add the derivative to the Jacobian matrix. Otherwise the FixedPoint formulation becomes the jacobian.
    if(this->linearization_ != "FixedPoint"){
	    this->assemblyAdvectionInU(elementMatrixW);
	    elementMatrixW->scale(this->density_);
    }

	//elementMatrix->add((*constantMatrix_),(*elementMatrix));
	this->jacobian_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

	this->jacobian_->add((*this->ANB_),(*this->jacobian_));
    // If the linearization is Newtons Method we need to add W-Matrix
    if(this->linearization_ != "FixedPoint"){
    	this->jacobian_->add((*elementMatrixW),(*this->jacobian_));  // int add(SmallMatrix<T> &bMat, SmallMatrix<T> &cMat); //this+B=C elementMatrix + constantMatrix_;
    }
}

// Laplacian in  ' \Delta u ' sense. Here apply new assembly :D  
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assemblyLaplacian(SmallMatrixPtr_Type &elementMatrix) {

	int dim = this->getDim();
	int numNodes= this->numNodesVelocity_;
	int Grad =2; // Needs to be fixed	
	string FEType = this->FETypeVelocity_;
	int dofs = this->dofsVelocity_;

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,Grad);
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
			 /*if (std::fabs(value[j]) < pow(10,-14)) {
		            value[j] = 0.;
		        }*/
			for (UN d=0; d<dofs; d++) {
              (*elementMatrix)[i*dofs +d][j*dofs+d] = value[j];
            }
        }

    }
}

// Here update please to unlinearized System Matrix accordingly.
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::assembleRHS(){

	SmallMatrixPtr_Type elementMatrixN =Teuchos::rcp( new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_));

	this->ANB_.reset(new SmallMatrix_Type( this->dofsElementVelocity_+this->numNodesPressure_)); // A + B + N
	this->ANB_->add( (*this->constantMatrix_),(*this->ANB_));

	this->assemblyAdvection(elementMatrixN);
	elementMatrixN->scale(this->density_);
	this->ANB_->add( (*elementMatrixN),(*this->ANB_));

	this->rhsVec_.reset( new vec_dbl_Type ( dofsElement_,0.) );
	// Multiplying ANB_ * solution // System Matrix times solution
	int s=0,t=0;
	for(int i=0 ; i< this->ANB_->size();i++){
		if (i >= this->dofsElementVelocity_)
			s=1;
		for(int j=0; j < this->ANB_->size(); j++){
			if(j >= this->dofsElementVelocity_)
				t=1;
			(*this->rhsVec_)[i] += (*this->ANB_)[i][j]*(*this->solution_)[j]*this->coeff_[s][t];
			//cout <<"Solution["<<j <<"]" << this->solution_[i] << endl;
		}
		t=0;
		//cout <<"RHS["<<i <<"]" << this->rhsVec_[i] << endl;
	}
}

/*
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assemblyAdvection(SmallMatrixPtr_Type &elementMatrix){

	int dim = this->getDim();
	int numNodes= numNodesVelocity_;
	int Grad =2; // Needs to be fixed	
	string FEType = FETypeVelocity_;
	int dofs = dofsVelocity_;


	vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type  phi;
	vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

	UN deg = 5; //Helper::determineDegree(dim,FEType,Grad); // Not complete
	//UN extraDeg = determineDegree( dim, FEType, Std); //Elementwise assembly of grad u
    //UN deg = determineDegree( dim, FEType, FEType, Grad, Std, extraDeg);

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
    applyBTinv( dPhi, dPhiTrans, Binv );

    for (int w=0; w<phi->size(); w++){ //quads points
        for (int d=0; d<dim; d++) {
            uLoc[d][w] = 0.;
            for (int i=0; i < phi->at(0).size(); i++) {
                LO index = dim * i + d;
                uLoc[d][w] += this->solution_[index] * phi->at(w).at(i);
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
         

     
		}
 		for (UN d=0; d<dim; d++) {
    	    for (UN j=0; j < indices.size(); j++)
    	   		 (*elementMatrix)[i*dofs +d][j*dofs+d] = value[j];
			
        }
        
    }

}*/

/*
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assemblyAdvectionInU(SmallMatrixPtr_Type &elementMatrix){

	int dim = this->getDim();
	int numNodes= numNodesVelocity_;
	int Grad =2; // Needs to be fixed	
	string FEType = FETypeVelocity_;
	int dofs = dofsVelocity_;


	vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type  phi;
	vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

	UN deg = 5 ;// Helper::determineDegree(dim,FEType,Grad); // Not complete
	//UN extraDeg = determineDegree( dim, FEType, Std); //Elementwise assembly of grad u
    //UN deg = determineDegree( dim, FEType, FEType, Grad, Std, extraDeg);

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
    applyBTinv( dPhi, dPhiTrans, Binv );
    //UN FEloc = checkFE(dim,FEType);


    std::vector<SmallMatrix<SC> > duLoc( weights->size(), SmallMatrix<SC>(dim) ); //for all quad points p_i each matrix is [u_x * grad Phi(p_i), u_y * grad Phi(p_i), u_z * grad Phi(p_i) (if 3D) ], duLoc[w] = [[phixx;phixy],[phiyx;phiyy]] (2D)

    for (int w=0; w<dPhiTrans.size(); w++){ //quads points
        for (int d1=0; d1<dim; d1++) {
            for (int i=0; i < dPhiTrans[0].size(); i++) {
                LO index = dim *i+ d1;
                for (int d2=0; d2<dim; d2++)
                    duLoc[w][d2][d1] += this->solution_[index] * dPhiTrans[w][i][d2];
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
}*/
 

/*
template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokes<SC,LO,GO,NO>::assemblyDivAndDivT(SmallMatrixPtr_Type &elementMatrix) {

    vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
	int dim = this->dim_;

    UN deg =2; // Helper::determineDegree2( dim, FETypeVelocity_, FETypePressure_, Grad, Std);

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
    applyBTinv( dPhi, dPhiTrans, Binv );

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

    }


}*/


/*!

 \brief Building Transformation

@param[in] &B

*/

template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::buildTransformation(SmallMatrix<SC>& B){

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

template <class SC, class LO, class GO, class NO>
void AssembleFENavierStokesNonNewtonian<SC,LO,GO,NO>::applyBTinv( vec3D_dbl_ptr_Type& dPhiIn,
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

