#ifndef AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp
#define AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp

#include "AssembleFE_SCI_SMC_Active_Growth_Reorientation_decl.hpp"

#ifdef FEDD_HAVE_ACEGENINTERFACE
#include "aceinterface.hpp"
#endif

#include <vector>
//#include <iostream>

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO>::AssembleFE_SCI_SMC_Active_Growth_Reorientation(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
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
		$[Gamma]$6 -Gamma6_24 									0.15e1
		$[Lambda]$P50 -LambdaP50_25								1.e0
		kDotMin -KDotMin_26										-0.118863e-2
		$[Zeta]$1 -Zeta1_27										100.e0
		kDotMax -KDotMax_28										0.51028e-3
		$[Gamma]$4 -Gamma4_29									200.e0	
		$[Lambda]$BarDotPMin -LambdaBarDotPMin_30 				-0.17689e-3
		$[Lambda]$BarDotPMax -LambdaBarDotPMax_31 				0.13957e-3
		$[Gamma]$5 -Gamma5_32									50.e0
		$[Zeta]$2 -Zeta2_33										1000.e0
		$[CapitalDelta]$$[Lambda]$BarPMin -DeltaLambdaBarPMin_34	-0.1e-4
		p1 -P1_35												0.3e0
		p3 -P3_36												0.2e0
		c50 -C50_37												0.5e0
		d0 -D0_38												0.6e-4					
		m -M_39													0e0
		activeStartTime -ActiveStartTime_40						1000.e0
		k$[Eta]$Plus -kEtaPlus_41								0.6e0
		m$[Eta]$Plus -mEtaPlus_42								5.e0
		growthStartTime -growthStartTime_43						1.e0
		reorientationStartTime -reorientationStartTime_44 		1.e0
		growthEndTime -growthEndTime_45							100.e0
		reorientationEndTime -reorientationEndTime_46			100.e0
		k$[Theta]$Plus -KThetaPlus_47							1.e0
		k$[Theta]$Minus -KThetaMinus_48							1.e0	
		m$[Theta]$Plus -MThetaPlus_49							3.e0
		m$[Theta]$Minus -MThetaMinus_50							3.e0
		$[Theta]$Plus1 -ThetaPlus1_51							0.1000882e1
		$[Theta]$Plus2 -ThetaPlus2_52							0.1234826e1					
		$[Theta]$Plus3 -ThetaPlus3_53							0.11414189999999999e1						
		$[Theta]$Minus1 -ThetaMinus1_54							0.98e0
		$[Theta]$Minus2 -ThetaMinus2_55							0.98e0
		$[Theta]$Minus3 -ThetaMinus3_56							0.98e0
		$[Rho]$ -Density_57										1.e0
										
	*/

	
	// -------------------- Parameter ---------------------
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
	kappaC_ = this->params_->sublist("Parameter Solid").get("KappaC",146.36600000000002e-0); 
	beta1_ = this->params_->sublist("Parameter Solid").get("Beta1",0.10097e-2); // ??
	muA_ = this->params_->sublist("Parameter Solid").get("MuA",0.9291e1); 
	alpha_ = this->params_->sublist("Parameter Solid").get("Alpha",0.2668e2); 
	epsilon1_ = this->params_->sublist("Parameter Solid").get("Epsilon1",151.73775e0); 
	epsilon2_ = this->params_->sublist("Parameter Solid").get("Epsilon2",0.27566199999999996e1); // ??
	c1_ = this->params_->sublist("Parameter Solid").get("C1",11.52507e0);
	alpha1_ = this->params_->sublist("Parameter Solid").get("Alpha1",1.27631e0);
	alpha2_ = this->params_->sublist("Parameter Solid").get("Alpha2",0.308798e1); // ?? 
	gamma6_ = this->params_->sublist("Parameter Solid").get("Gamma6",0.15e1);
	lambdaP50_ = this->params_->sublist("Parameter Solid").get("LambdaP50",1.0e0);
	kDotMin_ = this->params_->sublist("Parameter Solid").get("KDotMin",-0.118863e-2);
	zeta1_ = this->params_->sublist("Parameter Solid").get("Zeta1",100.e0);
	kDotMax_ = this->params_->sublist("Parameter Solid").get("KDotMax",0.51028e-3);
	gamma4_ = this->params_->sublist("Parameter Solid").get("Gamma4",200.e0);
	lambdaBarDotPMin_ = this->params_->sublist("Parameter Solid").get("LambdaBarDotMin",-0.17689e-3);
	lambdaBarDotPMax_ = this->params_->sublist("Parameter Solid").get("LambdaBarDotMax",0.13957e-3);
	gamma5_ = this->params_->sublist("Parameter Solid").get("Gamma5", 50.e0 );
	zeta2_ = this->params_->sublist("Parameter Solid").get("Zeta2", 1000.e0 );
	DeltaLambdaBarPMin_ =this->params_->sublist("Parameter Solid").get("DeltaLambdaBarPMin", -0.1e-4);
	p1_ = this->params_->sublist("Parameter Solid").get("P1",0.3e0);
	p3_ = this->params_->sublist("Parameter Solid").get("P3",0.2e0);
	c50_ = this->params_->sublist("Parameter Solid").get("C50",0.5e0);
	d0_ = this->params_->sublist("Parameter Diffusion").get("D0",6.e-05);
	m_ = this->params_->sublist("Parameter Solid").get("m",0.e0);
	activeStartTime_ = this->params_->sublist("Parameter Solid").get("ActiveStartTime",1000.e0); // At Starttime 1000 the diffused drug influences the material model. -> Active response at T=starttime	
	kEtaPlus_ = this->params_->sublist("Parameter Solid").get("KEtaPlus",0.6e0);
	mEtaPlus_ = this->params_->sublist("Parameter Solid").get("MEtaPlus",5.0e0);
	growthStartTime_ = this->params_->sublist("Parameter Solid").get("GrowthStartTime",1.e0);
	reorientationStartTime_ = this->params_->sublist("Parameter Solid").get("ReorientationStartTime",1.e0);
	growthEndTime_ = this->params_->sublist("Parameter Solid").get("GrowthEndTime",100.e0);
	reorientationEndTime_ = this->params_->sublist("Parameter Solid").get("ReorientationEndTime",100.e0);
	kThetaPlus_ = this->params_->sublist("Parameter Solid").get("KThetaPlus",1.e0);
	kThetaMinus_ = this->params_->sublist("Parameter Solid").get("KThetaMinus",1.e0);
	mThetaPlus_ = this->params_->sublist("Parameter Solid").get("MThetaPlus",3.e0);
	mThetaMinus_ = this->params_->sublist("Parameter Solid").get("MThetaMinus",3.e0);
	thetaPlus1_ = this->params_->sublist("Parameter Solid").get("ThetaPlus1",0.1000882e1);
	thetaPlus2_ = this->params_->sublist("Parameter Solid").get("ThetaPlus2",0.1234826e1);
	thetaPlus3_ = this->params_->sublist("Parameter Solid").get("ThetaPlus3",0.11414189999999999e1);
	thetaMinus1_ = this->params_->sublist("Parameter Solid").get("ThetaMinus1",0.98e0);
	thetaMinus2_ = this->params_->sublist("Parameter Solid").get("ThetaMinus2",0.98e0);
	thetaMinus3_ = this->params_->sublist("Parameter Solid").get("ThetaMinus3",0.98e0);
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
	this->history_ ={1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                     1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                     1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                     1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	
	// 34 * 4 
	this->historyUpdated_.resize(136,0.);

	solutionC_n_.resize(10,0.);
	solutionC_n1_.resize(10,0.);

	this->solution_.reset( new vec_dbl_Type ( dofsElement_,0.) );

	timeParametersVec_.resize(0, vec_dbl_Type(2));
    numSegments_ = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").get("Number of Segments",0);

 	for(int i=1; i <= numSegments_; i++){

        double startTime = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("Start Time",0.);
        double dtTmp = this->params_->sublist("Timestepping Parameter").sublist("Timestepping Intervalls").sublist(std::to_string(i)).get("dt",0.1);
        
        vec_dbl_Type segment = {startTime,dtTmp};
        timeParametersVec_.push_back(segment);
    }


}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO>::assembleJacobian() {

    SmallMatrixPtr_Type elementMatrix = Teuchos::rcp( new SmallMatrix_Type(this->dofsElement_,0.));

    assemble_SCI_SMC_Active_Growth_Reorientation(elementMatrix);

    this->jacobian_ = elementMatrix;
	
}
template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO>::advanceInTime( double dt){

	//cout << " advanced in time for this element with dt " << dt << endl;
	this->timeIncrement_ = dt;

	// If we have a time segment setting we switch to the demanded time increment
	for(int i=0; i<numSegments_ ; i++){
		if(this->timeStep_ >= timeParametersVec_[i][0])
			this->timeIncrement_=timeParametersVec_[i][1];
	}

	//cout << " Changed to timeincrement " << this->timeIncrement_<< endl;
	this->timeStep_ = this->timeStep_ + this->timeIncrement_;
	
	for(int i=0; i< 136; i++){
		if(this->timeStep_  > this->activeStartTime_ +dt )
			this->history_[i] = this->historyUpdated_[i];
		
	}

	for(int i=0; i< 10 ; i++)
		this->solutionC_n_[i]=(*this->solution_)[i+30]; // this is the LAST solution of newton iterations

	
}

template <class SC, class LO, class GO, class NO>
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO>::assembleRHS(){

	this->rhsVec_.reset( new vec_dbl_Type ( this->dofsElement_,0.) );

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

	double history[136];// = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
    for(int i = 0; i < this->history_.size(); i++)
		history[i] = this->history_[i];
   
   
    double concentrations[10];
    for(int i = 0; i < 10; i++)
	{
		concentrations[i]= (*this->solution_)[i+30];	
		this->solutionC_n1_[i] = (*this->solution_)[i+30];		// in each newtonstep solution for n+1 is updated.
	}	
	
    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   	

    double domainData[] = {this->fA_, this->lambdaC50_, this->gamma3_, this->lambdaBarCDotMax_, this->lambdaBarCDotMin_, 
							this->gamma2_ , this->gamma1_, this->eta1_, this->ca50_, this->k2_, this->k5_, this->k3_, this->k4_, this->k7_ , this->kappaC_ , 
							this->beta1_ , this->muA_ , this->alpha_, this->epsilon1_,this->epsilon2_,this->c1_,this->alpha1_,this->alpha2_,
							gamma6_ ,lambdaP50_ , kDotMin_ , zeta1_ , kDotMax_ ,gamma4_ ,lambdaBarDotPMin_ ,lambdaBarDotPMax_ , gamma5_, zeta2_, DeltaLambdaBarPMin_ ,
							this->p1_,this->p3_,this->c50_,this->d0_,this->m_,this->activeStartTime_,
		 					kEtaPlus_ , mEtaPlus_ , growthStartTime_ ,reorientationStartTime_ , growthEndTime_ , reorientationEndTime_ ,kThetaPlus_ ,
							kThetaMinus_ , mThetaPlus_ , mThetaMinus_ , thetaPlus1_ , thetaPlus2_ ,thetaPlus3_ , thetaMinus1_ ,thetaMinus2_ ,thetaMinus3_ ,	this->rho_,0.};

    double rates[10];

	for(int i=0; i<10 ; i++){
		rates[i] =(this->solutionC_n1_[i]-this->solutionC_n_[i])/deltaT;//(solutionC_n1_[i])/deltaT; //-solutionC_n_[i](solutionC_n1_[i]-solutionC_n_[i])/deltaT;//
	}

    
    double time = this->getTimeStep();
    double subIterationTolerance = 1.e-7;
 
    // immer speicher und wenn es konvergiert, dann zur history machen
    double historyUpdated[] = {1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    		   1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                    		   1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                     		   1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 1., 1., 1.512656, 1.512656, 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	
	AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 elem(positions, displacements, concentrations, accelerations, rates, domainData, history, subIterationTolerance, deltaT, time, this->iCode_);
	elem.compute();

	double *residuumRint = elem.getResiduumVectorRint();

	double *residuumRDyn = elem.getResiduumVectorRdyn();

	for(int i=0; i< 30 ; i++){
		(*this->rhsVec_)[i] = residuumRint[i]; //+residuumRDyn[i];
	}
	double *residuumRc = elem.getResiduumVectorRc();



	for(int i=0; i< 10 ; i++){		
		(*this->rhsVec_)[i+30] = residuumRc[i];
	}


	
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
void AssembleFE_SCI_SMC_Active_Growth_Reorientation<SC,LO,GO,NO>::assemble_SCI_SMC_Active_Growth_Reorientation(SmallMatrixPtr_Type &elementMatrix){
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
    double history[136];// = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
    for(int i = 0; i < this->history_.size(); i++)
		history[i] = this->history_[i];   //if (integrationCode==19)


    double concentrations[10];
    for(int i = 0; i < 10; i++)
	{
		concentrations[i]=(*this->solution_)[i+30];	
		solutionC_n1_[i]=(*this->solution_)[i+30];		// in each newtonstep solution for n+1 is updated.
	}

    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   
   
    double domainData[] = {this->fA_, this->lambdaC50_, this->gamma3_, this->lambdaBarCDotMax_, this->lambdaBarCDotMin_, 
							this->gamma2_ , this->gamma1_, this->eta1_, this->ca50_, this->k2_, this->k5_, this->k3_, this->k4_, this->k7_ , this->kappaC_ , 
							this->beta1_ , this->muA_ , this->alpha_, this->epsilon1_,this->epsilon2_,this->c1_,this->alpha1_,this->alpha2_,
							gamma6_ ,lambdaP50_ , kDotMin_ , zeta1_ , kDotMax_ ,gamma4_ ,lambdaBarDotPMin_ ,lambdaBarDotPMax_ , gamma5_, zeta2_, DeltaLambdaBarPMin_ ,
							this->p1_,this->p3_,this->c50_,this->d0_,this->m_,this->activeStartTime_,
		 					kEtaPlus_ , mEtaPlus_ , growthStartTime_ ,reorientationStartTime_ , growthEndTime_ , reorientationEndTime_ ,kThetaPlus_ ,
							kThetaMinus_ , mThetaPlus_ , mThetaMinus_ , thetaPlus1_ , thetaPlus2_ ,thetaPlus3_ , thetaMinus1_ ,thetaMinus2_ ,thetaMinus3_ ,	this->rho_,0.};


	double rates[10];
	for(int i=0; i<10 ; i++){
		rates[i] =(this->solutionC_n1_[i]-this->solutionC_n_[i]) / deltaT;//
	}
	// ##########################

    // history  [in] Vector of history variables [Order: LambdaBarC1, LambdaBarC2, nA1, nA2, nB1, nB2, nC1, nC2, nD1, nD2, LambdaA1, LambdaA2] (The length must be equal to number of history variables per gauss point * number of gauss points)
    
    double time = this->getTimeStep();
    double subIterationTolerance = 1.e-7;
    	
    // immer speicher und wenn es konvergiert, dann zur history machen
    // double historyUpdated[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
	
	AceGenInterface::DeformationDiffusionSmoothMuscleActiveGrowthReorientationTetrahedra3D10 elem(positions, displacements, concentrations, accelerations, rates, domainData, history, subIterationTolerance, deltaT, time, this->iCode_);
	elem.compute();

	double *historyUpdated = elem.getHistoryUpdated();

	double** stiffnessMatrixKuu = elem.getStiffnessMatrixKuu();
	double** stiffnessMatrixKuc = elem.getStiffnessMatrixKuc();
	double** stiffnessMatrixKcu = elem.getStiffnessMatrixKcu();
	double** stiffnessMatrixKcc = elem.getStiffnessMatrixKcc();
	double** massMatrixMc = elem.getMassMatrixMc();

	for(int i=0; i< 136; i++){
		this->historyUpdated_[i] = historyUpdated[i];
	}


	for(int i=0; i< 30; i++){
		for(int j=0; j<30; j++){
			//if(fabs(stiffnessMatrixKuu[i][j]) > 1e7)
			//	cout << " !!! Sus entry Kuu [" << i << "][" << j << "] " << stiffnessMatrixKuu[i][j] << endl; 
			
			(*elementMatrix)[i][j]=stiffnessMatrixKuu[i][j];
		}
	}
	for(int i=0; i< 30; i++){
		for(int j=0; j<10; j++){
			//if(fabs(stiffnessMatrixKuc[i][j]) > 1e7)
			//	cout << " !!! Sus entry Kuc [" << i << "][" << j << "] " << stiffnessMatrixKuc[i][j] << endl; 
			
			(*elementMatrix)[i][j+30]=stiffnessMatrixKuc[i][j];
		}
	}
	for(int i=0; i< 10; i++){
		for(int j=0; j<30; j++){
			//if(fabs(stiffnessMatrixKcu[i][j]) > 1e7)
			//	cout << " !!! Sus entry Kcu [" << i << "][" << j << "] " << stiffnessMatrixKcu[i][j] << endl; 
			
			(*elementMatrix)[i+30][j]=stiffnessMatrixKcu[i][j];
		}
	}
	for(int i=0; i< 10; i++){
		for(int j=0; j<10; j++){
			//if(fabs(massMatrixMc[i][j]) > 1e5 || fabs(stiffnessMatrixKcc[i][j]) > 1e5 )
			//	cout << " !!! Sus entry Mass [" << i << "][" << j << "] " << massMatrixMc[i][j] << " or stiff Kcc " << stiffnessMatrixKcc[i][j] << endl; 
			 
			(*elementMatrix)[i+30][j+30] =stiffnessMatrixKcc[i][j] +(1./deltaT)*massMatrixMc[i][j]; //
		}
	}
    
    


#endif

	
}


} // namespace FEDD
#endif // AssembleFE_SCI_SMC_Active_Growth_Reorientation_DEF_hpp
