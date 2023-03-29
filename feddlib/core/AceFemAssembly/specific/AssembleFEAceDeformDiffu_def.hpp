#ifndef ASSEMBLEFEACEDEFORMDIFFU_DEF_hpp
#define ASSEMBLEFEACEDEFORMDIFFU_DEF_hpp

#include "AssembleFEAceDeformDiffu_decl.hpp"
#ifdef FEDD_HAVE_ACEGENINTERFACE
#include "aceinterface.h"
#endif

#include <vector>
// #include <iostream>

namespace FEDD
{

	template <class SC, class LO, class GO, class NO>
	AssembleFEAceDeformDiffu<SC, LO, GO, NO>::AssembleFEAceDeformDiffu(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params, tuple_disk_vec_ptr_Type tuple) : AssembleFE<SC, LO, GO, NO>(flag, nodesRefConfig, params, tuple)
	{
		// Extracting values from ParameterList
		E0_ = this->params_->sublist("Parameter Solid").get("E", 379.95e-6);
		E1_ = this->params_->sublist("Parameter Solid").get("E1", 300.0e-6);
		poissonRatio_ = this->params_->sublist("Parameter Solid").get("Poisson Ratio", 0.49e-0);
		c1_ = this->params_->sublist("Parameter Solid").get("c1", 0.25e-0);
		D0_ = this->params_->sublist("Parameter Diffusion").get("D0", 6.0e-5);
		m_ = this->params_->sublist("Parameter Diffusion").get("m", 0.0);
		dofOrdering_ = this->params_->sublist("Parameter").get("Ordering", 2);

		FEType_ = std::get<1>(this->diskTuple_->at(0));	   // FEType of Disk
		dofsSolid_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
		dofsChem_ = std::get<2>(this->diskTuple_->at(1));  // Degrees of freedom per node

		numNodesSolid_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element
		numNodesChem_ = std::get<3>(this->diskTuple_->at(1));  // Number of nodes of element

		dofsElement_ = dofsSolid_ * numNodesSolid_ + dofsChem_ * numNodesChem_; // "Dimension of return matrix"

		solution_n_.resize(60, 0.);
		solution_n1_.resize(60, 0.);
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFEAceDeformDiffu<SC, LO, GO, NO>::assembleJacobian()
	{

		SmallMatrixPtr_Type elementMatrix = Teuchos::rcp(new SmallMatrix_Type(dofsElement_, 0.));

		assembleDeformationDiffusionNeoHook(elementMatrix);
		// elementMatrix->print();
		this->jacobian_ = elementMatrix;
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFEAceDeformDiffu<SC, LO, GO, NO>::advanceInTime(double dt)
	{

		// cout << " advanced in time for this element " << endl;
		if (this->timeStep_ < 1.)
			this->timeIncrement_ = 0.05;
		if (this->timeStep_ >= 1.)
			this->timeIncrement_ = 1.0;
		if (this->timeStep_ >= 10.)
			this->timeIncrement_ = 10.0;
		if (this->timeStep_ >= 1000.)
			this->timeIncrement_ = 0.01;

		//	this->timeIncrement_ = dt;
		this->timeStep_ = this->timeStep_ + this->timeIncrement_;

		for (int i = 0; i < 40; i++)
		{
			if (i < 30)
				solution_n_[i] = (*this->solution_)[i];
			else
				solution_n_[30 + 3 * (i - 30)] = (*this->solution_)[i];
		}
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFEAceDeformDiffu<SC, LO, GO, NO>::assembleRHS()
	{

		this->rhsVec_.reset(new vec_dbl_Type(dofsElement_, 0.));

		std::vector<double> positions(30);
		int count = 0;
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				positions[count] = this->getNodesRefConfig()[i][j];
				count++;
			}
		}

		std::vector<double> displacements(30);
		for (int i = 0; i < 30; i++)
			displacements[i] = (*this->solution_)[i];

		std::vector<double> concentrations(10);
		for (int i = 0; i < 10; i++)
			concentrations[i] = (*this->solution_)[i + 30];

		std::vector<double> concentrationsLastConverged(10);
		for (int i = 0; i < 10; i++)
			concentrationsLastConverged[i] = (solution_n_)[i + 30];

		std::vector<double> domainData(6);
		domainData[0] = this->E0_;
		domainData[1] = this->E1_;
		domainData[2] = this->poissonRatio_;
		domainData[3] = this->c1_;
		domainData[4] = this->D0_;
		domainData[5] = this->m_;

		double timeIncrement = this->getTimeIncrement();

		int integrationCode = 19;

		std::vector<double> residuum(40);

		std::vector<double> stiffnessMatrixFlat(40*40);
		double **stiffnessMatrix = (double**)calloc(40, sizeof(double *));
		for (int i = 0; i < 40; i++)
			stiffnessMatrix[i] = &stiffnessMatrixFlat[i * 40];

		getResiduumAndStiffnessMatrixNeoHooke(&positions[0], &displacements[0], &concentrations[0], &concentrationsLastConverged[0], &domainData[0], timeIncrement, integrationCode, &residuum[0], stiffnessMatrix);

		for (int i = 0; i < 40; i++)
			(*this->rhsVec_)[i] = residuum[i];

		free(stiffnessMatrix);
		// std::vector<double> v(2238);
		// std::vector<double> d(6);
		// std::vector<double> ul(60);
		// std::vector<double> ul0(60);
		// std::vector<double> xl(60);
		// std::vector<double> s(3600);
		// std::vector<double> p(60);
		// std::vector<double> ht(10);
		// std::vector<double> hp(10);

		// std::vector<double> deltat(1);

		// deltat[0]=this->getTimeIncrement();

		// d[0] = this->E0_;
		// d[1] = this->E1_;
		// d[2] = this->poissonRatio_;
		// d[3] = this->c1_;
		// d[4] = this->D0_;
		// d[5] = this->m_;

		// if(dofOrdering_ == 1)
		// {
		// 	for(int i = 0; i < 40; i++)
		// 		if((i+1)%4 == 0)
		// 			solution_n1_[27 + 3*(i+1)/4] = (*this->solution_)[i];
		// 		else
		// 			solution_n1_[i - ((i+1) - (i+1)%4)/4] = (*this->solution_)[i];
		// }
		// else if(dofOrdering_ == 2)
		// {
		// 	for(int i = 0; i < 40; i++)
		// 	{
		// 		if(i<30)
		// 			solution_n1_[i]=(*this->solution_)[i];
		// 		else
		// 			solution_n1_[30 + 3*(i-30)] = (*this->solution_)[i];
		// 	}
		// }
		// else
		// 	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown DOF ordering sequence. Known identifiers: 1 and 2. Check parameters file!");

		// int count = 0;
		// for(int i=0;i<10;i++)
		// 	for(int j=0;j<3;j++){
		// 		xl[count] = this->getNodesRefConfig()[i][j];
		// 		count++;}

		// for(int i=0;i<s.size();i++)
		// 	s[i]=0.0;
		// for(int i=0;i<p.size();i++)
		// 	p[i]=0.0;

		// std::cout << "From RHS call: " << std::endl;

		// std::cout << "Nodes ref config: " << std::endl;
		// for(int i=0;i<30;i++)
		// 	std::cout << xl[i] << " ";
		// std::cout << std::endl;

		// std::cout << "solution_n1_ " << std::endl;
		// for(int i=0;i<30;i++)
		// 	std::cout << solution_n1_[i] << " ";
		// std::cout << std::endl;

		// std::cout << "solution_n_ " << std::endl;
		// for(int i=0;i<30;i++)
		// 	std::cout << solution_n_[i] << " ";
		// std::cout << std::endl;

		// skr_DDNH(&v[0],&d[0],&solution_n1_[0],&solution_n_[0],&xl[0],&s[0],&p[0],&ht[0],&hp[0],&deltat[0]);

		// std::cout << "p " << std::endl;
		// for(int i=0;i<40;i++)
		// 	std::cout << p[i] << " ";
		// std::cout << std::endl;

		// if(dofOrdering_ == 1)
		// {
		// 	for(int i=0;i<40;i++)
		// 	{
		// 		if((i+1)%4==0)
		// 			(*this->rhsVec_)[i] = p[29+(i+1)/4];
		// 		else
		// 			(*this->rhsVec_)[i] = p[i - ((i+1) - (i+1)%4)/4];
		// 	}
		// }
		// else if(dofOrdering_ == 2)
		// {
		// 	for(int i = 0; i < 40; i++){
		// 		if(fabs(p[i])>1e-14)
		// 			(*this->rhsVec_)[i] = p[i];
		// 		else
		// 			(*this->rhsVec_)[i] = 0.0;

		// 	}
		// }
		// else
		// 	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown DOF ordering sequence. Known identifiers: 1 and 2. Check parameters file!");
	}

	template <class SC, class LO, class GO, class NO>
	void AssembleFEAceDeformDiffu<SC, LO, GO, NO>::assembleDeformationDiffusionNeoHook(SmallMatrixPtr_Type &elementMatrix)
	{

		std::vector<double> positions(30);
		int count = 0;
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				positions[count] = this->getNodesRefConfig()[i][j];
				count++;
			}
		}

		std::vector<double> displacements(30);
		for (int i = 0; i < 30; i++)
			displacements[i] = (*this->solution_)[i];

		std::vector<double> concentrations(10);
		for (int i = 0; i < 10; i++)
			concentrations[i] = (*this->solution_)[i + 30];

		std::vector<double> concentrationsLastConverged(10);
		for (int i = 0; i < 10; i++)
			concentrationsLastConverged[i] = (solution_n_)[i + 30];

		std::vector<double> domainData(6);
		domainData[0] = this->E0_;
		domainData[1] = this->E1_;
		domainData[2] = this->poissonRatio_;
		domainData[3] = this->c1_;
		domainData[4] = this->D0_;
		domainData[5] = this->m_;

		double timeIncrement = this->getTimeIncrement();

		int integrationCode = 19;

		std::vector<double> residuum(40);

		std::vector<double> stiffnessMatrixFlat(40*40);
		double **stiffnessMatrix = (double**)calloc(40, sizeof(double *));
		for (int i = 0; i < 40; i++)
			stiffnessMatrix[i] = &stiffnessMatrixFlat[i * 40];

		std::cout << "From stiffness matrix call: " << std::endl;

		std::cout << "Nodes ref config: " << std::endl;
		for (int i = 0; i < 30; i++)
			std::cout << positions[i] << ",";
		std::cout << std::endl;

		std::cout << "displacements: " << std::endl;
		for (int i = 0; i < 30; i++)
			std::cout << displacements[i] << ",";
		std::cout << std::endl;

		std::cout << "concentrations: " << std::endl;
		for (int i = 0; i < 10; i++)
			std::cout << concentrations[i] << ",";
		std::cout << std::endl;

		std::cout << "concentrationsLastConverged: " << std::endl;
		for (int i = 0; i < 10; i++)
			std::cout << concentrationsLastConverged[i] << ",";
		std::cout << std::endl;

		std::cout << "domainData: " << std::endl;
		for (int i = 0; i < 6; i++)
			std::cout << domainData[i] << ",";
		std::cout << std::endl;

		std::cout << "timeIncrement: " << timeIncrement << std::endl;

		std::cout << "integrationCode: " << integrationCode << std::endl;

		getResiduumAndStiffnessMatrixNeoHooke(&positions[0], &displacements[0], &concentrations[0], &concentrationsLastConverged[0], &domainData[0], timeIncrement, integrationCode, &residuum[0], stiffnessMatrix);

		std::cout << "residuum: " << std::endl;

		for (int i = 0; i < 40; i++)
			std::cout << residuum[i] << ",";
		std::cout << std::endl;
		

		for (UN i = 0; i < this->dofsElement_; i++)
		{
			for (UN j = 0; j < this->dofsElement_; j++)
			{
				(*elementMatrix)[i][j] = stiffnessMatrix[i][j];
			}
		}

		free(stiffnessMatrix);
		// double deltat=this->getTimeIncrement();
		//  std::vector<double> deltat(1);

		// deltat[0]=this->getTimeIncrement();

		// d[0] = this->E0_;
		// d[1] = this->E1_;
		// d[2] = this->poissonRatio_;
		// d[3] = this->c1_;
		// d[4] = this->D0_;
		// d[5] = this->m_;

		// if(dofOrdering_ == 1)
		// {
		// 	for(int i = 0; i < 40; i++)
		// 		if((i+1)%4 == 0)
		// 			solution_n1_[27 + 3*(i+1)/4] = (*this->solution_)[i];
		// 		else
		// 			solution_n1_[i - ((i+1) - (i+1)%4)/4] = (*this->solution_)[i];
		// }
		// else if(dofOrdering_ == 2)
		// {
		// 	for(int i = 0; i < 40; i++)
		// 	{
		// 		if(i<30)
		// 			solution_n1_[i]=(*this->solution_)[i];
		// 		else
		// 			solution_n1_[30 + 3*(i-30)] = (*this->solution_)[i];
		// 	}
		// }
		// else
		// 	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown DOF ordering sequence. Known identifiers: 1 and 2. Check parameters file!");

		// int count = 0;
		// for(int i=0;i<10;i++)
		// 	for(int j=0;j<3;j++){
		// 		xl[count] = this->getNodesRefConfig()[i][j];
		// 		count++;}

		// std::cout << "Nodes ref config: " << std::endl;
		// for(int i=0;i<30;i++)
		// 	std::cout << xl[i] << " ";
		// std::cout << std::endl;

		// std::cout << "solution_n1_ " << std::endl;
		// for(int i=0;i<30;i++)
		// 	std::cout << solution_n1_[i] << " ";
		// std::cout << std::endl;

		// std::cout << "solution_n_ " << std::endl;
		// for(int i=0;i<30;i++)
		// 	std::cout << solution_n_[i] << " ";
		// std::cout << std::endl;

		// for(int i=0;i<s.size();i++)
		// 	s[i]=0.0;
		// for(int i=0;i<p.size();i++)
		// 	p[i]=0.0;

		// skr_DDNH(&v[0],&d[0],&solution_n1_[0],&solution_n1_[0],&xl[0],&s[0],&p[0],&ht[0],&hp[0],&deltat[0]);

		// if(dofOrdering_ == 2)
		// {

		// 	std::vector<double> s_temp(1600);
		// 	int counter=0;
		// 	for(int i=0;i<40;i++)
		// 	{
		// 		for(int j=0;j<40;j++)
		// 		{
		// 			s_temp[counter] = s[60*i+j];
		// 			counter++;
		// 		}
		// 	}

		// 	for (UN i=0; i < this->dofsElement_; i++) {
		// 		for (UN j=0; j < this->dofsElement_; j++) {
		// 			if(fabs(s_temp[40*j+i]) > 1e-14)
		// 				(*elementMatrix)[i][j] = s_temp[40*j+i]; // Rolling into a matrix using column major (m*j+i)
		// 		}
		// 	}

		// }
		// else if(dofOrdering_ == 1) // Caution: Untested!
		// {
		// 	for(int i = 0; i < 40; i++)
		// 	{
		// 		if((i+1)%4==0)
		// 		{
		// 			int l = 29 + (i+1)/4;
		// 			for(int j = 0; j < 40; j++)
		// 			{
		// 				if((j+1)%4==0)
		// 				{
		// 					int m = 29 + (j+1)/4;
		// 					(*elementMatrix)[i][j] = s[40*m+l];
		// 				}
		// 				else
		// 				{
		// 					int m = j - ((j+1)-(j+1)%4)/4;
		// 					(*elementMatrix)[i][j] = s[40*m+l];
		// 				}
		// 			}
		// 		}
		// 		else
		// 		{
		// 			int l = i - ((i+1)-(i+1)%4)/4;
		// 			for(int j = 0; j < 40; j++)
		// 			{
		// 				if((j+1)%4==0)
		// 				{
		// 					int m = 29 + (j+1)/4;
		// 					(*elementMatrix)[i][j] = s[40*m+l];
		// 				}
		// 				else
		// 				{
		// 					int m = j - ((j+1)-(j+1)%4)/4;
		// 					(*elementMatrix)[i][j] = s[40*m+l];
		// 				}
		// 			}
		// 		}
		// 	}
		//}
		// else
		//	TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown DOF ordering sequence. Known identifiers: 1 and 2. Check parameters file!");

		// Need to modify the above based on dof ordering flag selected
	}

} // namespace FEDD
#endif // ASSEMBLEFEACEDEFORMDIFFU_DEF_hpp
