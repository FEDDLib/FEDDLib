#ifndef AssembleFE_SCI_SMC_MLCK_DEF_hpp
#define AssembleFE_SCI_SMC_MLCK_DEF_hpp

#include "AssembleFE_SCI_SMC_MLCK_decl.hpp"

#ifdef FEDD_HAVE_ACEGENINTERFACE
#include <aceinterface.hpp>
#endif

#include <vector>
//#include <iostream>

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AssembleFE_SCI_SMC_MLCK<SC,LO,GO,NO>::AssembleFE_SCI_SMC_MLCK(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params, tuple)
{
		/*
		fA -Fibre angle_1  							30e0, 
		$[Lambda]$C50 -LambdaC50_2				 	0.12e1
		$[Gamma]$3 -Gamma3_3 						0.9e0,
		$[Lambda]$BarCDotMax -LambdaBarCDotMax_4 	0.3387e-1
		$[Lambda]$BarCDotMin -LambdaBarCDotMin_5 	-0.3387e-1
		$[Gamma]$2 -Gamma2_6 						50e0,
		$[Gamma]$1 -Gamma1_7 						0.50247e0, 
		$[Eta]$1 -Eta1_8 							0.18745e0, 
		Ca50 -Ca50_9 								0.4e0
		k2 -K2_10 									0.2e0
		k5 -K5_11 									0.2e0,
		k3 -K3_12 									0.134e0
		k4 -K4_13 									0.166e-2
		k7 -K7_14 									0.66e-4
		$[Kappa]$C -KappaC_15 						0.14636600000000002e3
		$[Beta]$1 -Beta1_16							0.10097e-2
		$[Mu]$a -MuA_17 							0.9291e1
		$[Alpha]$ -Alpha_18 						0.2668e2
		$[Epsilon]$1 -Epsilon1_19					0.15173775e3
		$[Epsilon]$2 -Epsilon2_20 					0.27566199999999996e1
		c1 -C1_21 									0.1152507e2
		$[Alpha]$1 -Alpha1_22 						0.127631e1
		$[Alpha]$2 -Alpha2_23 						0.308798e1
		p1 -P1_24 									0.3e0
		p3 -P3_25 									0.2e0
		c50 -C50_26 								0.5e0
		d0 -D0_27 									0.6e-4
		m -M_28 								    0e0
		startTime -StartTime_29  					1000e0
		$[Rho]$0 -Density_30 						1e0
													
		*/

	fA_= this->params_->sublist("Parameter Solid").get("FA",30.e0); // ??
	lambdaC50_ = this->params_->sublist("Parameter Solid").get("LambdaC50",0.12e1); // ??
	gamma3_= this->params_->sublist("Parameter Solid").get("Gamma3",0.9e0);
	lambdaBarCDotMax_= this->params_->sublist("Parameter Solid").get("LambdaBarCDotMax",0.3387e-1); // ??
	lambdaBarCDotMin_= this->params_->sublist("Parameter Solid").get("LambdaBarCDotMin",-0.3387e-1); // ?? 
	gamma2_ = this->params_->sublist("Parameter Solid").get("Gamma2",50.0e0); // ??
	gamma1_ = this->params_->sublist("Parameter Solid").get("Gamma1",0.50247e0); 
	eta1_ = this->params_->sublist("Parameter Solid").get("Eta1",0.18745e0); // ??
	ca50_ = this->params_->sublist("Parameter Solid").get("Ca50",0.4e0); // ??
	k2_ = this->params_->sublist("Parameter Solid").get("K2",0.2e0); 
	k5_ = this->params_->sublist("Parameter Solid").get("K5",0.2e0);
	k3_ = this->params_->sublist("Parameter Solid").get("K3",0.134e0); // ??
	k4_ = this->params_->sublist("Parameter Solid").get("K4",0.166e-2); // ??
	k7_= this->params_->sublist("Parameter Solid").get("K7",0.66e-4); // ?? 
	kappaC_ = this->params_->sublist("Parameter Solid").get("KappaC",146.36600000000002e0);
	beta1_ = this->params_->sublist("Parameter Solid").get("Beta1",0.10097e-2); // ??
	muA_ = this->params_->sublist("Parameter Solid").get("MuA",0.9291e1); 
	alpha_ = this->params_->sublist("Parameter Solid").get("Alpha",0.2668e2); 
	epsilon1_ = this->params_->sublist("Parameter Solid").get("Epsilon1", 0.15173775e3);
	epsilon2_ = this->params_->sublist("Parameter Solid").get("Epsilon2",0.27566199999999996e1); // ??
	c1_ = this->params_->sublist("Parameter Solid").get("C1",11.52507e0);
	alpha1_ = this->params_->sublist("Parameter Solid").get("Alpha1",1.27631e0);
	alpha2_ = this->params_->sublist("Parameter Solid").get("Alpha2",0.308798e1); // ?? 
	p1_ = this->params_->sublist("Parameter Solid").get("P1",0.3e0);
	p3_ = this->params_->sublist("Parameter Solid").get("P3",0.2e0);
	c50_ = this->params_->sublist("Parameter Solid").get("C50",0.5e0);
	d0_ = this->params_->sublist("Parameter Diffusion").get("D0",6.e-05);
	m_ = this->params_->sublist("Parameter Solid").get("m",0.e0);
	startTime_ = this->params_->sublist("Parameter Solid").get("ActiveStartTime",1001.e0); // At Starttime 1000 the diffused drug influences the material model. -> Active response at T=starttime
	rho_ = this->params_->sublist("Parameter Solid").get("Rho",1.e0);

	// iCode_ = this->params_->sublist("Parameter Solid").get("Intergration Code",18);
	iCode_=18; //Only works for 18 currently!!

    FEType_ = std::get<1>(this->diskTuple_->at(0)); // FEType of Disk
	dofsSolid_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
	dofsChem_ = std::get<2>(this->diskTuple_->at(1)); // Degrees of freedom per node

	numNodesSolid_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element
	numNodesChem_ = std::get<3>(this->diskTuple_->at(1)); // Number of nodes of element

	dofsElement_ = dofsSolid_*numNodesSolid_ + dofsChem_*numNodesChem_; // "Dimension of return matrix"

	// Einlesen durch Parameterdatei irgendwann cool
	history_ ={1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
	//.resize(48);
	
	historyUpdated_.resize(48,0.);

	solutionC_n_.resize(10,0.);
	solutionC_n1_.resize(10,0.);

	this->solution_.reset( new vec_dbl_Type ( dofsElement_,0.) );
	
	/*timeParametersVec_.resize(0, vec_dbl_Type(2));
    numSegments_ = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

 	for(int i=1; i <= numSegments_; i++){

        double startTime = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
        double dtTmp = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
        
        vec_dbl_Type segment = {startTime,dtTmp};
        timeParametersVec_.push_back(segment);
    }*/

}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_MLCK<SC,LO,GO,NO>::assembleJacobian() {

    SmallMatrixPtr_Type elementMatrix = Teuchos::rcp( new SmallMatrix_Type(dofsElement_,0.));

    assemble_SCI_SMC_MLCK(elementMatrix);

    this->jacobian_ = elementMatrix;
	
}
template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_MLCK<SC,LO,GO,NO>::advanceInTime( double dt){

	// If we have a time segment setting we switch to the demanded time increment
	/*for(int i=0; i<numSegments_ ; i++){
		if(this->timeStep_ +1.0e-12 > timeParametersVec_[i][0])
			this->timeIncrement_=timeParametersVec_[i][1];
	}*/
	this->timeStep_ = this->timeStep_ + this->timeIncrement_;
	
	this->timeIncrement_ = dt;

	for(int i=0; i< 48; i++){
	//	if(this->timeStep_  > startTime_ +dt )
			history_[i] = historyUpdated_[i];
		
	}

	for(int i=0; i< 10 ; i++)
		solutionC_n_[i]=(*this->solution_)[i+30]; // this is the LAST solution of newton iterations

	
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_MLCK<SC,LO,GO,NO>::assembleRHS(){

	this->rhsVec_.reset( new vec_dbl_Type ( dofsElement_,0.) );

#ifdef FEDD_HAVE_ACEGENINTERFACE


	double deltaT=this->getTimeIncrement();

	double positions[30];
	int count = 0;
	//cout << "Positions " << endl;
	for(int i=0;i<10;i++){
		for(int j=0;j<3;j++){
			positions[count] = this->getNodesRefConfig()[i][j];
			count++;
	//		cout << " i " << count-1 << " " <<  positions[count-1] << endl;
		}
	}
	//cout << endl;

	double displacements[30];
	for(int i = 0; i < 30; i++)
	{
		displacements[i]=(*this->solution_)[i];		
	}

	double history[48];// = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
    for(int i = 0; i < history_.size(); i++)
		history[i] = history_[i];
   
   
    double concentrations[10];
    for(int i = 0; i < 10; i++)
	{
		concentrations[i]= (*this->solution_)[i+30];	
		solutionC_n1_[i] = (*this->solution_)[i+30];		// in each newtonstep solution for n+1 is updated.
	}	
	
    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   	
    double domainData[] = {fA_, lambdaC50_, gamma3_, lambdaBarCDotMax_, lambdaBarCDotMin_, 
							gamma2_ , gamma1_, eta1_, ca50_, k2_, k5_, k3_, k4_, k7_ , kappaC_ , 
							beta1_ , muA_ , alpha_, epsilon1_,epsilon2_,c1_,alpha1_,alpha2_,
							p1_,p3_,c50_,d0_,m_,startTime_,rho_};

    double rates[10];

	for(int i=0; i<10 ; i++){
		rates[i] =(solutionC_n1_[i]-solutionC_n_[i])/deltaT;//(solutionC_n1_[i])/deltaT; //-solutionC_n_[i](solutionC_n1_[i]-solutionC_n_[i])/deltaT;//
	}

    
    double time = this->getTimeStep()+deltaT;
    double subIterationTolerance = 1.e-7;
 
    // immer speicher und wenn es konvergiert, dann zur history machen
    double historyUpdated[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
 	
	AceGenInterface::DeformationDiffusionSmoothMuscleMLCKOnlyTetrahedra3D10 elem(positions, displacements, concentrations, accelerations, rates, domainData, history, subIterationTolerance, deltaT, time, iCode_);
	elem.compute();

	double *residuumRint = elem.getResiduumVectorRint();

	double *residuumRDyn = elem.getResiduumVectorRdyn();

	// getResiduumVectorRint(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &domainData[0], &history[0], subIterationTolerance, deltaT, time, iCode_, &historyUpdated[0], residuumRint);
    
	// getResiduumVectorRdyn(&positions[0], &displacements[0], &concentrations[0], &accelerations[0],&rates[0], &domainData[0], &history[0], subIterationTolerance, deltaT, time, iCode_, &historyUpdated[0], residuumRDyn);
	for(int i=0; i< 30 ; i++){
		(*this->rhsVec_)[i] = residuumRint[i]; //+residuumRDyn[i];
	}
	double *residuumRc = elem.getResiduumVectorRc();
	// getResiduumVectorRc(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &domainData[0], &history[0], subIterationTolerance, deltaT, time, iCode_, &historyUpdated[0], residuumRc);


	for(int i=0; i< 10 ; i++){		
		(*this->rhsVec_)[i+30] = residuumRc[i];
	}

	// free(residuumRc);
	// free(residuumRint);
	// free(residuumRDyn);
	

	
#endif

}

/**
 * @brief Example for stiffness matrix Kuc
 * @param nodalPositionReference [in] The nodal positions in the reference coordinates in order (x1,y1,z1,x2,y2,z2..)
 * @param displacements [in] The nodal displacements
 * @param concentrations [in] The nodal concentrations
 * @param accelerations [in] The nodal accelerations
 * @param rates [in] The nodal concentrations rates (dC/dt)
 * @param domainData [in] Vector of domain data (Use function getDomainDataNames() to get the names and order)
 * @param history  [in] Vector of history variables [Order: LambdaBarC1, LambdaBarC2, nA1, nA2, nB1, nB2, nC1, nC2, nD1, nD2, LambdaA1, LambdaA2] (The length must be equal to number of history variables per gauss point * number of gauss points)
 * @param subIterationTolerance [in] Tolernace for the local Newton iterations (Gauss Point values). Recommended 1e-7.
 * @param timeIncrement [in] Time increment
 * @param time [in] Current time
 * @param integrationScheme [in] Integration code (18 for 4 point, 19 for 5 point, 40 for 14 point, 43 for 8 point)
 * @param historyUpdated [out] Updated values of the history variables after the Newton-Raphson iteration
 * @param stiffnessMatrixKuu [out] The stiffness matrix Kuu
 */


template <class SC,class LO, class GO, class NO>
void AssembleFE_SCI_SMC_MLCK<SC,LO,GO,NO>::assemble_SCI_SMC_MLCK(SmallMatrixPtr_Type &elementMatrix){
	// double deltat=this->getTimeIncrement();
	// std::vector<double> deltat(1);
	// deltat[0]=this->getTimeIncrement();
	
#ifdef FEDD_HAVE_ACEGENINTERFACE

	double deltaT=this->getTimeIncrement();
 	
	double positions[30];
	int count = 0;
	for(int i=0;i<10;i++){
		for(int j=0;j<3;j++){
			positions[count] =  this->getNodesRefConfig()[i][j];
			count++;
		}
	}
	double displacements[30];
	for(int i = 0; i < 30; i++)
	{
		displacements[i]= (*this->solution_)[i];	
		
	}
    double history[48];// = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
    for(int i = 0; i < history_.size(); i++)
		history[i] = history_[i];   //if (integrationCode==19)


    double concentrations[10];
    for(int i = 0; i < 10; i++)
	{
		concentrations[i]=(*this->solution_)[i+30];	
		solutionC_n1_[i]=(*this->solution_)[i+30];		// in each newtonstep solution for n+1 is updated.
	}

    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   
    double domainData[] = {fA_, lambdaC50_, gamma3_, lambdaBarCDotMax_, lambdaBarCDotMin_, 
							gamma2_ , gamma1_, eta1_, ca50_, k2_, k5_, k3_, k4_, k7_ , kappaC_ , 
							beta1_ , muA_ , alpha_, epsilon1_,epsilon2_,c1_,alpha1_,alpha2_,
							p1_,p3_,c50_,d0_,m_,startTime_,rho_};

	double rates[10];
	for(int i=0; i<10 ; i++){
		rates[i] =(solutionC_n1_[i]-solutionC_n_[i]) / deltaT;//
	}
	// ##########################
    
    double time = this->getTimeStep()+deltaT;
    double subIterationTolerance = 1.e-7;
    	
    // immer speicher und wenn es konvergiert, dann zur history machen
    // double historyUpdated[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
	
	AceGenInterface::DeformationDiffusionSmoothMuscleMLCKOnlyTetrahedra3D10 elem(positions, displacements, concentrations, accelerations, rates, domainData, history, subIterationTolerance, deltaT, time, iCode_);
	elem.compute();

	double *historyUpdated = elem.getHistoryUpdated();

	double** stiffnessMatrixKuu = elem.getStiffnessMatrixKuu();
	double** stiffnessMatrixKuc = elem.getStiffnessMatrixKuc();
	double** stiffnessMatrixKcu = elem.getStiffnessMatrixKcu();
	double** stiffnessMatrixKcc = elem.getStiffnessMatrixKcc();
	double** massMatrixMc = elem.getMassMatrixMc();

    
	for(int i=0; i< 48; i++){
		historyUpdated_[i] = historyUpdated[i];
	}
	
	for(int i=0; i< 30; i++){
		for(int j=0; j<30; j++){
			//if(std::fabs(stiffnessMatrixKuu[i][j]) > 1e7)
			//	cout << " !!! Sus entry Kuu [" << i << "][" << j << "] " << stiffnessMatrixKuu[i][j] << endl; 
			
			(*elementMatrix)[i][j]=stiffnessMatrixKuu[i][j];
		}
	}
	for(int i=0; i< 30; i++){
		for(int j=0; j<10; j++){
			//if(std::fabs(stiffnessMatrixKuc[i][j]) > 1e7)
			//	cout << " !!! Sus entry Kuc [" << i << "][" << j << "] " << stiffnessMatrixKuc[i][j] << endl; 
			
			(*elementMatrix)[i][j+30]=stiffnessMatrixKuc[i][j];
		}
	}
	for(int i=0; i< 10; i++){
		for(int j=0; j<30; j++){
			//if(std::fabs(stiffnessMatrixKcu[i][j]) > 1e7)
			//	cout << " !!! Sus entry Kcu [" << i << "][" << j << "] " << stiffnessMatrixKcu[i][j] << endl; 
			
			(*elementMatrix)[i+30][j]=stiffnessMatrixKcu[i][j];
		}
	}
	for(int i=0; i< 10; i++){
		for(int j=0; j<10; j++){
			//if(std::fabs(massMatrixMc[i][j]) > 1e5 || std::fabs(stiffnessMatrixKcc[i][j]) > 1e5 )
			//	cout << " !!! Sus entry Mass [" << i << "][" << j << "] " << massMatrixMc[i][j] << " or stiff Kcc " << stiffnessMatrixKcc[i][j] << endl; 
			 
			(*elementMatrix)[i+30][j+30] =stiffnessMatrixKcc[i][j] +(1./deltaT)*massMatrixMc[i][j]; //
		}
	}

    


#endif

	
}


} // namespace FEDD
#endif // AssembleFE_SCI_SMC_MLCK_DEF_hpp
