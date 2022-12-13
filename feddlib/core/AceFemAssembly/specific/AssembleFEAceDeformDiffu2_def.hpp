#ifndef AssembleFEAceDeformDiffu2_DEF_hpp
#define AssembleFEAceDeformDiffu2_DEF_hpp

#include "AssembleFEAceDeformDiffu2_decl.hpp"

#ifdef FEDD_HAVE_ACEGENINTERFACE
#include "aceinterface.h"
#include "ace2.h"
#endif

#include <vector>
//#include <iostream>

namespace FEDD {

template <class SC, class LO, class GO, class NO>
AssembleFEAceDeformDiffu2<SC,LO,GO,NO>::AssembleFEAceDeformDiffu2(int flag, vec2D_dbl_Type nodesRefConfig, ParameterListPtr_Type params,tuple_disk_vec_ptr_Type tuple):
AssembleFE<SC,LO,GO,NO>(flag, nodesRefConfig, params, tuple)
{
    // Extracting values from ParameterList
    E0_ = this->params_->sublist("Parameter Solid").get("E",379.95e-6);
    E1_ = this->params_->sublist("Parameter Solid").get("E1",300.0e-6);
    poissonRatio_ = this->params_->sublist("Parameter Solid").get("Poisson Ratio",0.49e-0);
    c1_ = this->params_->sublist("Parameter Solid").get("c1",0.25e-0);
    D0_ = this->params_->sublist("Parameter Diffusion").get("D0",6.0e-5);
    m_ = this->params_->sublist("Parameter Diffusion").get("m",0.0);
    dofOrdering_ = this->params_->sublist("Parameter").get("Ordering", 2);

    FEType_ = std::get<1>(this->diskTuple_->at(0)); // FEType of Disk
	dofsSolid_ = std::get<2>(this->diskTuple_->at(0)); // Degrees of freedom per node
	dofsChem_ = std::get<2>(this->diskTuple_->at(1)); // Degrees of freedom per node

	numNodesSolid_ = std::get<3>(this->diskTuple_->at(0)); // Number of nodes of element
	numNodesChem_ = std::get<3>(this->diskTuple_->at(1)); // Number of nodes of element

	dofsElement_ = dofsSolid_*numNodesSolid_ + dofsChem_*numNodesChem_; // "Dimension of return matrix"

}

template <class SC, class LO, class GO, class NO>
void AssembleFEAceDeformDiffu2<SC,LO,GO,NO>::assembleJacobian() {

    SmallMatrixPtr_Type elementMatrix = Teuchos::rcp( new SmallMatrix_Type(dofsElement_));

    assembleDeformationDiffusionNeoHook(elementMatrix);

    this->jacobian_ = elementMatrix;
}

template <class SC, class LO, class GO, class NO>
void AssembleFEAceDeformDiffu2<SC,LO,GO,NO>::assembleRHS(){

	this->rhsVec_ = vec_dbl_Type(dofsElement_,0);

	std::vector<double> v(2238);
	std::vector<double> d(6);
	std::vector<double> ul(60); 
	std::vector<double> ul0(60);
	std::vector<double> xl(60);
	std::vector<double> s(3600);
	std::vector<double> p(60);
	std::vector<double> ht(10);
	std::vector<double> hp(10);
	
	std::vector<double> deltat(1);

	deltat[0]=this->getTimeIncrement();

	d[0] = this->E0_;
	d[1] = this->E1_;
	d[2] = this->poissonRatio_;
	d[3] = this->c1_;
	d[4] = this->D0_;
	d[5] = this->m_;

	if(dofOrdering_ == 1)
	{
		for(int i = 0; i < 40; i++)
			if((i+1)%4 == 0)
				ul[27 + 3*(i+1)/4] = this->getSolution()[i];
			else
				ul[i - ((i+1) - (i+1)%4)/4] = this->getSolution()[i];
	}
	else if(dofOrdering_ == 2)
	{
		for(int i = 0; i < 40; i++)
		{
			if(i<30)
				ul[i]=this->getSolution()[i];
			else
				ul[30 + 3*(i-30)] = this->getSolution()[i];
		}
	}
	else
		TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown DOF ordering sequence. Known identifiers: 1 and 2. Check parameters file!");

	int count = 0;
	for(int i=0;i<10;i++)
		for(int j=0;j<3;j++){
			xl[count] = this->getNodesRefConfig()[i][j];
			count++;}

	for(int i=0;i<s.size();i++)
		s[i]=0.0;
	for(int i=0;i<p.size();i++)
		p[i]=0.0;

	//skr_DDNH(&v[0],&d[0],&ul[0],&ul0[0],&xl[0],&s[0],&p[0],&ht[0],&hp[0],&deltat[0]);

}

/**
 * @brief This function returns the stiffness matrix Kuc
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
 * param historyUpdated [out] Updated values of the history variables after the Newton-Raphson iteration
 * @param stiffnessMatrixKuc [out] The stiffness matrix Kuc
 */


template <class SC,class LO, class GO, class NO>
void AssembleFEAceDeformDiffu2<SC,LO,GO,NO>::assembleDeformationDiffusionNeoHook(SmallMatrixPtr_Type &elementMatrix){

	//double deltat=this->getTimeIncrement();
	std::vector<double> deltat(1);
	deltat[0]=this->getTimeIncrement();
	
	d[0] = this->E0_;
	d[1] = this->E1_;
	d[2] = this->poissonRatio_;
	d[3] = this->c1_;
	d[4] = this->D0_;
	d[5] = this->m_;

#ifdef FEDD_HAVE_ACEGENINTERFACE
 	double **stiffnessMatrixKuu = (double **)malloc(30 * sizeof(double *));
    for (int i = 0; i < 30; i++)
        stiffnessMatrixKuu[i] = (double *)malloc(30 * sizeof(double));
        
    double **stiffnessMatrixKuc = (double **)malloc(30 * sizeof(double *));
    for (int i = 0; i < 30; i++)
        stiffnessMatrixKuc[i] = (double *)malloc(10 * sizeof(double));
        
    double **stiffnessMatrixKcu = (double **)malloc(10 * sizeof(double *));
    for (int i = 0; i < 10; i++)
        stiffnessMatrixKcu[i] = (double *)malloc(30 * sizeof(double));
        
    double **stiffnessMatrixKcc = (double **)malloc(10 * sizeof(double *));
    for (int i = 0; i < 10; i++)
        stiffnessMatrixKcc[i] = (double *)malloc(10 * sizeof(double));


	double positions[30];
	int count = 0;
	for(int i=0;i<10;i++)
		for(int j=0;j<3;j++){
			positions[count] = this->getNodesRefConfig()[i][j];
			count++;
	}
	double displacements[30];
	for(int i = 0; i < 30; i++)
	{
		displacements[i]=this->getSolution()[i];			
	}
	
    // find number of entries in positions
    int n1 = sizeof(positions) / sizeof(positions[0]);
    // find number of entries in displacements
    int n2 = sizeof(displacements) / sizeof(displacements[0]);
    // printf("n1=%d\n", n1);
    // printf("n2=%d\n", n2);
    char **a = getDataNames();
    
   
    double concentrations[10];
    for(int i = 30; i < 40; i++)
	{
		concentrations[i]=this->getSolution()[i];			
	}
	
    double accelerations[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
   
    // Print all entries in a
     int i = 0;
     while (a[i] != NULL)
     {
         printf("%s\n", a[i]);
         i++;
     }

   
    double domainData[] = {30e0, 0.12e1, 0.9e0, 0.3387e-1, -0.3387e-1, 50e0,
                     0.50247e0, 0.18745e0, 0.4e0, 0.2e0, 0.2e0, 0.134e0,
                     0.166e-2, 0.66e-4, 0.14636600000000002e3, 0.10097e-2, 0.9291e1, 0.2668e2,
                     0.15173775e3, 0.27566199999999996e1, 0.1152507e2, 0.127631e1, 0.308798e1, 0.3e0,
                     0.2e0, 0.5e0, 0.6e-4, 0e0, 1000e0, 1e0,
                     0e0}; // I have to check every time
                     
    double rates[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    // history  [in] Vector of history variables [Order: LambdaBarC1, LambdaBarC2, nA1, nA2, nB1, nB2, nC1, nC2, nD1, nD2, LambdaA1, LambdaA2] (The length must be equal to number of history variables per gauss point * number of gauss points)
    
    
    double deltaT = this->getTimeIncrement();
    double time = this->getTimeStep();
    double subIterationTolerance = 1e-7;
    // integrationScheme [in] Integration code (18 for 4 point, 19 for 5 point, 40 for 14 point, 43 for 8 point)
    int integrationCode = 18;
    
    double history[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
    //if (integrationCode==19)
    //   history = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 60 values, 12 variables, 5 gausspoints
  
    	
    double historyUpdated[] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0}; // 48 values, 12 variables, 4 gausspoints
   
   
    // find number of entries in concentrations
    int n3 = sizeof(concentrations) / sizeof(concentrations[0]);
    // find number of entries in rates
    int n4 = sizeof(rates) / sizeof(rates[0]);
    // find number of entries in accelerations
    int n5 = sizeof(accelerations) / sizeof(accelerations[0]);


    getStiffnessMatrixKuu(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &domainData[0], &history[0], subIterationTolerance, deltaT, time, integrationCode, &historyUpdated[0], stiffnessMatrixKuu);
	getStiffnessMatrixKuc(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &domainData[0], &history[0], subIterationTolerance, deltaT, time, integrationCode,  stiffnessMatrixKuc);
	getStiffnessMatrixKcu(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &domainData[0], &history[0], subIterationTolerance, deltaT, time, integrationCode,  stiffnessMatrixKcu);
	getStiffnessMatrixKcc(&positions[0], &displacements[0], &concentrations[0], &accelerations[0], &rates[0], &domainData[0], &history[0], subIterationTolerance, deltaT, time, integrationCode, stiffnessMatrixKcc);



	for(int i=0; i< 30; i++){
		for(int j=0; j<30; j++){
			(*elementMatrix)[i][j]=stiffnessMatrixKuu[i][j];
		}
	}
	for(int i=0; i< 30; i++){
		for(int j=0; j<10; j++){
			(*elementMatrix)[i][j+30]=stiffnessMatrixKuc[i][j];
		}
	}
	for(int i=0; i< 10; i++){
		for(int j=0; j<30; j++){
			(*elementMatrix)[i+30][j]=stiffnessMatrixKcu[i][j];
		}
	}
	for(int i=0; i< 10; i++){
		for(int j=0; j<10; j++){
			(*elementMatrix)[i+30][j+30]=stiffnessMatrixKcc[i][j];
		}
	}
	
	elementMatrix->print();
#endif

	
}

// Need to modify the above based on dof ordering flag selected


} // namespace FEDD
#endif // AssembleFEAceDeformDiffu2_DEF_hpp
