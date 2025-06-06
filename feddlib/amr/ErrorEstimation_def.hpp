#ifndef ErrorEstimation_def_hpp
#define ErrorEstimation_def_hpp

#ifndef MESH_TIMER_START
#define MESH_TIMER_START(A,S) Teuchos::RCP<Teuchos::TimeMonitor> A = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("Mesh Refinement") + std::string(S))));
#endif

#ifndef MESH_TIMER_STOP
#define MESH_TIMER_STOP(A) A.reset();
#endif

#include "ErrorEstimation_decl.hpp"
#include "feddlib/core/LinearAlgebra/MultiVector_def.hpp"
#include <chrono> 
/*!
 Definition of ErrorEstimation
 
  \brief  ErrorEstimation
 @author Lea Saßsmannshausen

*/

using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::outArg;
using std::cout;
using std::endl;
using std::setprecision;

namespace FEDD {


/*!
\brief ErrorEstimation is initiated with the dimension and problemType, as it is necessary to fit errorEstimation to the problem at hand.

@param[in] dim Dimension of problem.
@param[in] problemType Type of problem. Choose between Laplace, Stokes and NavierStokes

*/
template <class SC, class LO, class GO, class NO>
ErrorEstimation<SC,LO,GO,NO>:: ErrorEstimation(int dim, std::string problemType, bool writeMeshQuality)
{
	this->dim_ = dim;
	this->problemType_ = problemType;	
	this->writeMeshQuality_ = writeMeshQuality;
}

template <class SC, class LO, class GO, class NO>
ErrorEstimation<SC,LO,GO,NO>::~ErrorEstimation(){
}

/*!
\brief Identifying the problem with respect to the degrees of freedom and whether we calculate pressure. By telling how many block the MultiVector has, we can tell whether we calculate pressure or not. Depending on the numbers of entries within the solution vector, we can tell how many degreesOfFreedom (dofs) we have.

@param[in] valuesSolution Blocks that contains solution for velocity and as the case may be pressure.

*/
template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::identifyProblem(BlockMultiVectorConstPtr_Type valuesSolution){

	// dofs can be determined by checking the solutions' blocks an comparing the length of the vectors to the number of unique nodes. 
	// If the length is the same: dofs =1 , if it's twice as long: dofs =2 .. and so on

	// P_12 generally represents the velocity domain:
	if(inputMesh_->getMapUnique()->getNodeNumElements() == valuesSolution->getBlock(0)->getDataNonConst(0).size())
		dofs_ = 1;
	else if(2*(inputMesh_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(0)->getDataNonConst(0).size())
		dofs_ = 2;
	else if(3*(inputMesh_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(0)->getDataNonConst(0).size())
		dofs_ = 3;
	
	// If ValuesSolution has more than one block its typically because the have also a pressure solution
	if(valuesSolution->size() > 1){
		if(inputMeshP1_->getMapUnique()->getNodeNumElements() == valuesSolution->getBlock(1)->getDataNonConst(0).size())
			dofsP_ = 1;
		else if(2*(inputMeshP1_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(1)->getDataNonConst(0).size())
			dofsP_ = 2;
		else if(3*(inputMeshP1_->getMapUnique()->getNodeNumElements()) == valuesSolution->getBlock(1)->getDataNonConst(0).size())
			dofsP_ = 3;
		
		calculatePressure_ = true;
	}
	else
		dofsP_=0;

}


/*!
\brief We split the solution from the BlockMultiVector valuesSolution into one or more seperate blocks, where the blocks represent the different dimensions

@param[in] valuesSolution Block that contains solution for velocity and as the case may be pressure.

*/

template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::makeRepeatedSolution(BlockMultiVectorConstPtr_Type valuesSolution){
	// We extraxt so solution into dofs-many vectors an use them to make the repeated solution an put it back together in a block multivector later
	BlockMultiVectorPtr_Type valuesSolutionVel = Teuchos::rcp( new BlockMultiVector_Type(dofs_ ) );

	Teuchos::ArrayRCP< SC > valuesVel  = valuesSolution->getBlock(0)->getDataNonConst(0);

	MultiVectorPtr_Type mvValues =  Teuchos::rcp( new MultiVector_Type(inputMesh_->getMapUnique(), 1 ) );
	Teuchos::ArrayRCP< SC > mvValuesA  = mvValues->getDataNonConst(0);	

	for(int d=0; d< dofs_ ; d++){
		MultiVectorPtr_Type mvValuesRep =  Teuchos::rcp( new MultiVector_Type(inputMesh_->getMapRepeated(), 1 ) );	
		for(int i=0; i< mvValuesA.size(); i++){
			mvValuesA[i] = valuesVel[i*dofs_+d];

		}
		mvValuesRep->importFromVector(mvValues,false,"Insert");
		valuesSolutionVel->addBlock(mvValuesRep,d);		
		
	}
	if(calculatePressure_ == true){

		BlockMultiVectorPtr_Type valuesSolutionPre = Teuchos::rcp( new BlockMultiVector_Type(dofsP_ ) );

		Teuchos::ArrayRCP< SC > valuesPre  = valuesSolution->getBlock(1)->getDataNonConst(0);

		MultiVectorPtr_Type mvValues =  Teuchos::rcp( new MultiVector_Type(inputMeshP1_->getMapUnique(), 1 ) );
		Teuchos::ArrayRCP< SC > mvValuesA  = mvValues->getDataNonConst(0);	

		for(int d=0; d< dofsP_ ; d++){
			MultiVectorPtr_Type mvValuesRep =  Teuchos::rcp( new MultiVector_Type(inputMeshP1_->getMapRepeated(), 1 ) );	
			for(int i=0; i< mvValuesA.size(); i++){
				mvValuesA[i] = valuesPre[i*dofsP_+d];

			}
			mvValuesRep->importFromVector(mvValues,false,"Insert");
			valuesSolutionPre->addBlock(mvValuesRep,d);		
			
		}
		
		valuesSolutionRepPre_ = valuesSolutionPre;

	}

	valuesSolutionRepVel_ = valuesSolutionVel;

}
/*!
\brief Main Function for a posteriori error estimation

\brief depending on the problem the the error estimation is calculated accordingly.

@param[in] inputMeshP1 The P1 Mesh that is used for later refinement.
@param[in] inputMeshP12 The possible P2 Mesh, if one of the solutions is of P2 Discretisation, otherwise both meshes are P1.
@param[in] solution Solution of the PDE in BlockMultiVector Format (Block 0: Velocity, Block 1: Pressure).
@param[in] rhs The right hand side function of the pde.
@param[in] FETypeV Finite element type as the maximum FEType for the Velocity, pressure is assumed to be P1 always.
*/ 
template <class SC, class LO, class GO, class NO>
typename ErrorEstimation<SC,LO,GO,NO>::MultiVectorPtr_Type ErrorEstimation<SC,LO,GO,NO>::estimateError(MeshUnstrPtr_Type inputMeshP12, MeshUnstrPtr_Type inputMeshP1, BlockMultiVectorConstPtr_Type valuesSolution, RhsFunc_Type rhsFunc, std::string FETypeV){
	
	// Setting the InputMeshes
	inputMesh_ = inputMeshP12;
	inputMeshP1_ = inputMeshP1;

	// Identifying the Problem with respect to the dofs_ and splitting the Solution
	this->identifyProblem(valuesSolution);
	this->makeRepeatedSolution(valuesSolution);


	// Setting FETypes
	FEType1_ = "P1"; // Pressure FEType
	FEType2_ = FETypeV; // Velocity or u FEType

	if(FEType2_ != "P1" && FEType2_ != "P2"){ 
   		TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Error Estimation only works for Triangular P1 or P2 Elements");
	}


	ElementsPtr_Type elements = inputMesh_->getElementsC();
   	EdgeElementsPtr_Type edgeElements = inputMesh_->getEdgeElements();
	MapConstPtr_Type elementMap = inputMesh_->getElementMap();
	vec2D_dbl_ptr_Type points = inputMesh_->getPointsRepeated();
	int myRank = inputMesh_->comm_->getRank();
	surfaceElements_ = inputMesh_->getSurfaceTriangleElements();

	this->dim_ = inputMesh_->dim_;

    
	// 3D Residual Based Error Estimation work when using Surfaces
	// - New Surface Set
	// - 
	MESH_TIMER_START(errorEstimation," Error Estimation ");
 
	edgeElements->matchEdgesToElements(elementMap);

	vec2D_GO_Type elementsOfEdgesGlobal = edgeElements->getElementsOfEdgeGlobal();
	vec2D_LO_Type elementsOfEdgesLocal = edgeElements->getElementsOfEdgeLocal();

	// Vector size of elements for the resdual based error estimate
	MultiVectorPtr_Type errorElementMv = Teuchos::rcp(new MultiVector_Type( elementMap) ); 
	Teuchos::ArrayRCP<SC> errorElement = errorElementMv->getDataNonConst(0);

	double maxErrorEl;
	double maxErrorElLoc=0.0;


	if(this->dim_ == 2){ 

		// Edge Numbers of Element
		vec_int_Type edgeNumbers(3);

		
		// The Jump is over edges or surfaces is calculated beforehand, as it involves communication
		// In the function itselfe is decided which kind of jump is calculated.
		vec_dbl_Type u_Jump = calculateJump();

		
		// Calculating diameter of elements	
		vec_dbl_Type areaTriangles(elements->numberElements());
		vec_dbl_Type rho_T(elements->numberElements());
		vec_dbl_Type C_T(elements->numberElements());
		vec_dbl_Type h_T =  calcDiamTriangles(elements,points, areaTriangles, rho_T, C_T);

		// The divU Part and residual of the Element are calculated elementwise, as the are independet of other processors
		double divUElement=0;
		double resElement=0;
	
		for(int k=0; k< elements->numberElements();k++){
			edgeNumbers = edgeElements->getEdgesOfElement(k); // edges of Element k
			resElement = determineResElement(elements->getElement(k), rhsFunc); // calculating the residual of element k
			if(problemType_ == "Stokes" || problemType_ == "NavierStokes"){ // If the Problem is a (Navier)Stokes Problem we calculate divU (maybe start a hierarchy
				divUElement = determineDivU(elements->getElement(k));
			}

			errorElement[k] = std::sqrt(1./2*(u_Jump[edgeNumbers[0]]+u_Jump[edgeNumbers[1]]+u_Jump[edgeNumbers[2]])+divUElement+ h_T[k]*h_T[k]*resElement ); //divUElement ; 

			if(maxErrorElLoc < errorElement[k] )
					maxErrorElLoc = errorElement[k];
		}

		reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxErrorElLoc, outArg (maxErrorElLoc));

		if(	writeMeshQuality_ ){
			// We asses the Mesh Quality 
			double maxh_T, minh_T;
			double maxC_T, minC_T;
			double maxrho_T, minrho_T;
			double maxArea_T, minArea_T;

			auto it = max_element(h_T.begin(), h_T.end()); // 
			maxh_T =  h_T[distance(h_T.begin(), it)];
			it = min_element(h_T.begin(), h_T.end()); // 
			minh_T =  h_T[distance(h_T.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxh_T, outArg (maxh_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minh_T, outArg (minh_T));

			it = max_element(rho_T.begin(), rho_T.end()); // 
			maxrho_T =  rho_T[distance(rho_T.begin(), it)];
			it = min_element(rho_T.begin(), rho_T.end()); // 
			minrho_T =  rho_T[distance(rho_T.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxrho_T, outArg (maxrho_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minrho_T, outArg (minrho_T));

			it = max_element(areaTriangles.begin(), areaTriangles.end()); // 
			maxArea_T = areaTriangles[distance(areaTriangles.begin(), it)];
			it = min_element(areaTriangles.begin(), areaTriangles.end()); // 
			minArea_T =  areaTriangles[distance(areaTriangles.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxArea_T, outArg (maxArea_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minArea_T, outArg (minArea_T));


			it = max_element(C_T.begin(), C_T.end()); // 
			maxC_T =  C_T[distance(C_T.begin(), it)];
			it = min_element(C_T.begin(), C_T.end()); // 
			minC_T =  C_T[distance(C_T.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxC_T, outArg (maxC_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minC_T, outArg (minC_T));

			if(inputMesh_->getComm()->getRank() == 0){
				cout << "\tMesh Quality Assessment 2D of current mesh" << endl;
				cout << "\t__________________________________________________________________________________________________________ " << endl;
				cout << " " << endl;
				cout << " \tCircumdiameter h_T: \t\t" <<"max. = " << setprecision(5)  << maxh_T << "\tmin. = " << setprecision(5)  << minh_T  << endl;
				cout << " \tIncircumdiameter rho_T:\t\t" <<"max. = " << setprecision(5)  <<maxrho_T << "\tmin. = " <<setprecision(5)  << minrho_T  << endl;
				cout << " \tArea of Triangles:\t \t" <<"max. = " << setprecision(5)  << maxArea_T << "\tmin. = " << setprecision(5)  << minArea_T  << endl;
				cout << " \tShape parameter:\t \t" <<"max. = " << setprecision(5)  << maxC_T << "\tmin. = " << setprecision(5)  <<minC_T << endl;
				cout << " \tThe maximal Error of Elements is "  << maxErrorElLoc << endl;
				cout << "\t__________________________________________________________________________________________________________ " << endl;
			}
		}


				
	}

	else if(this->dim_ == 3){

		this->updateElementsOfSurfaceLocalAndGlobal(edgeElements, this->surfaceElements_);
		this->buildTriangleMap();
		this->surfaceElements_->matchSurfacesToElements(elementMap);
		SurfaceElementsPtr_Type surfaceElements = this->surfaceElements_;
	 
		// Jump per Element over edges
		vec_dbl_Type errorSurfacesInterior(4); 
		// Edge Numbers of Element
		vec_int_Type surfaceNumbers(4);
		// tmp values u_1,u_2 of gradient u of two elements corresponing with surface 
		vec_dbl_Type u1(3),u2(3);

		// gradient of u elementwise 
		vec_dbl_Type u_Jump = calculateJump();

		// h_E,min defined as in "A posteriori error estimation for anisotropic tetrahedral and triangular finite element meshes"
		// This is the minimal per element, later we also need the adjacent elements' minimal height in order to determine h_E
		vec_dbl_Type volTetraeder = determineVolTet(elements, points);

	    // Calculating diameter of elements	
		vec_dbl_Type areaTriangles(surfaceElements->numberElements());
		vec_dbl_Type rho_Tri(surfaceElements->numberElements());
		vec_dbl_Type C_Tri(surfaceElements->numberElements());
		vec_dbl_Type h_Tri =  calcDiamTriangles3D(surfaceElements,points, areaTriangles, rho_Tri, C_Tri);

		vec_dbl_Type h_T = calcDiamTetraeder(elements,points, volTetraeder);
		vec_dbl_Type rho_T = calcRhoTetraeder(elements,surfaceElements, volTetraeder, areaTriangles);

		// necessary entities
		vec_dbl_Type p1(3),p2(3); // normal Vector of Surface
		vec_dbl_Type v_E(3); // normal Vector of Surface
		double norm_v_E;

		int elTmp1,elTmp2;

	
		// Adjustment for 3D Implememtation	
		// In case we have more than one proc we need to exchange information via the interface. 
		// We design a multivector consisting of u's x and y values, to import and export it easier to other procs only for interface elements

		// The divU Part and residual of the Element are calculated elementwise, as the are independet of other processors
		double divUElement=0;
		double resElement=0;
	
		// Then we determine the jump over the edges, if the element we need for this is not on our proc, we import the solution u
		for(int k=0; k< elements->numberElements();k++){
			surfaceNumbers = surfaceElements->getSurfacesOfElement(k); // surfaces of Element k

			resElement = determineResElement(elements->getElement(k), rhsFunc); // calculating the residual of element k
			if(problemType_ == "Stokes" || problemType_ == "NavierStokes") // If the Problem is a (Navier)Stokes Problem we calculate divU (maybe start a hierarchy
				divUElement = determineDivU(elements->getElement(k));
			errorElement[k] = std::sqrt(1./2*(u_Jump[surfaceNumbers[0]]+u_Jump[surfaceNumbers[1]]+u_Jump[surfaceNumbers[2]]+u_Jump[surfaceNumbers[3]])+divUElement +h_T[k]*h_T[k]*resElement );
			if(maxErrorElLoc < errorElement[k] )
					maxErrorElLoc = errorElement[k];			
		}	

		reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxErrorElLoc, outArg (maxErrorElLoc));
		
		if(	writeMeshQuality_ ){
			double maxh_T, minh_T, maxh_Tri, minh_Tri;
			double maxC_T, minC_T, maxC_Tri, minC_Tri;
			double maxrho_T, minrho_T, maxrho_Tri, minrho_Tri;
			double maxArea_T, minArea_T;
			double maxVol_T, minVol_T;

			// h_Triangles
			auto it = max_element(h_Tri.begin(), h_Tri.end()); // 
			maxh_Tri =  h_Tri[distance(h_Tri.begin(), it)];
			it = min_element(h_Tri.begin(), h_Tri.end()); // 

			minh_Tri =  h_Tri[distance(h_Tri.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxh_Tri, outArg (maxh_Tri));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minh_Tri, outArg (minh_Tri));

			// h_Tetraeder
			it = max_element(h_T.begin(), h_T.end()); // 
			maxh_T =  h_T[distance(h_T.begin(), it)];
			it = min_element(h_T.begin(), h_T.end()); // 

			minh_T =  h_T[distance(h_T.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxh_T, outArg (maxh_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minh_T, outArg (minh_T));

			// rho_Tri
			it = max_element(rho_Tri.begin(), rho_Tri.end()); // 
			maxrho_Tri =  rho_Tri[distance(rho_Tri.begin(), it)];
			it = min_element(rho_Tri.begin(), rho_Tri.end()); // 
			minrho_Tri =  rho_Tri[distance(rho_Tri.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxrho_Tri, outArg (maxrho_Tri));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minrho_Tri, outArg (minrho_Tri));

			// rho_Tetraeder
			it = max_element(rho_T.begin(), rho_T.end()); // 
			maxrho_T =  rho_T[distance(rho_T.begin(), it)];
			it = min_element(rho_T.begin(), rho_T.end()); // 
			minrho_T =  rho_T[distance(rho_T.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxrho_T, outArg (maxrho_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minrho_T, outArg (minrho_T));

			// Area Triangles
			it = max_element(areaTriangles.begin(), areaTriangles.end()); // 
			maxArea_T = areaTriangles[distance(areaTriangles.begin(), it)];
			it = min_element(areaTriangles.begin(), areaTriangles.end()); // 
			minArea_T =  areaTriangles[distance(areaTriangles.begin(), it)];
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxArea_T, outArg (maxArea_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minArea_T, outArg (minArea_T));

			// Volume Tetraeder
			it = max_element(volTetraeder.begin(), volTetraeder.end()); // 
			maxVol_T = volTetraeder[distance(volTetraeder.begin(), it)];

			it = min_element(volTetraeder.begin(), volTetraeder.end()); // 
			minVol_T =  volTetraeder[distance(volTetraeder.begin(), it)];

			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxVol_T, outArg (maxVol_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minVol_T, outArg (minVol_T));

			// C_Tri
			it = max_element(C_Tri.begin(), C_Tri.end()); // 
			maxC_Tri =  C_Tri[distance(C_Tri.begin(), it)];

			it = min_element(C_Tri.begin(), C_Tri.end()); // 
			minC_Tri =  C_Tri[distance(C_Tri.begin(), it)];

			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxC_Tri, outArg (maxC_Tri));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minC_Tri, outArg (minC_Tri));

			// C_T
			vec_dbl_Type C_T(elements->numberElements());
			for(int i=0; i< h_T.size(); i++){
				C_T[i] = h_T[i] / rho_T[i];
			}
			it = max_element(C_T.begin(), C_T.end()); // 
			maxC_T =  C_T[distance(C_T.begin(), it)];

			it = min_element(C_T.begin(), C_T.end()); // 
			minC_T =  C_T[distance(C_T.begin(), it)];

			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxC_T, outArg (maxC_T));
			reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, minC_T, outArg (minC_T));


			if(inputMesh_->getComm()->getRank() == 0){
				cout << " \t-- Mesh Quality Assessment 3D of current mesh level -- \t" << endl;
				cout << "\t__________________________________________________________________________________________________________ " << endl;
				cout << " " << endl;
				cout << " \tCircumdiameter h_T:\t\t\t" <<"max. = "  << setprecision(5)  << maxh_T << " min. = " << setprecision(5)  <<minh_T  << endl;
				cout << " \tIncircumdiameter rho_T:\t\t\t" <<"max. = " <<setprecision(5)  << maxrho_T << " min. = " << setprecision(5)  <<minrho_T  << endl;
				cout << " \tCircumdiameter h_Tri:\t\t\t" <<"max. = " <<setprecision(5)  << maxh_Tri << " min. = " << setprecision(5)  <<minh_Tri  << endl;
				cout << " \tIncircumdiameter rho_Tri:\t\t" <<"max. = " << setprecision(5)  <<maxrho_Tri << " min. = " << setprecision(5)  <<minrho_Tri  << endl;
				cout << " \tArea of Triangles: \t\t\t" <<"max. = " << setprecision(5)  <<maxArea_T << " min. = " << setprecision(5)  <<minArea_T  << endl;
				cout << " \tVolume of Tetraeder: \t\t\t" <<"max. = " << setprecision(5)  <<maxVol_T << " min. = " << setprecision(5)  <<minVol_T  << endl;
				cout << " \tShape parameter Tetraeder: \t\t" <<"max. = " << setprecision(5)  << maxC_T << " min. = " << setprecision(5)  <<minC_T << endl;
				cout << " \tShape parameter Triangles: \t\t" <<"max. = " << setprecision(5)  << maxC_Tri << " min. = " << setprecision(5)  <<minC_Tri << endl;
				cout << " \tThe maximal Error of Elements is \t"  << maxErrorElLoc << endl;
				cout << "\t__________________________________________________________________________________________________________ " << endl;
			}
		}


	}

	errorEstimation_ = errorElementMv;
	
	MESH_TIMER_STOP(errorEstimation);

	return errorElementMv;

}


/*!
\brief Tags only a certain Area for refinement and is independent of any error estimation.

@param[in] inputMeshP1 The P1 Mesh that is used for later refinement.
@param[in] area  Area that is suppose to be refined. If is a vector defining the area as follows: row1:[x_0,x_1] x-limits, row2: [y_0,y_1] y-limits, row3: [z_0,z_1] z-limits .
*/
template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::tagArea( MeshUnstrPtr_Type inputMeshP1,vec2D_dbl_Type area){

	inputMesh_ = inputMeshP1;

	ElementsPtr_Type elements = inputMesh_->getElementsC();
   	EdgeElementsPtr_Type edgeElements = inputMesh_->getEdgeElements();
	MapConstPtr_Type elementMap = inputMesh_->getElementMap();
	int myRank = inputMesh_->comm_->getRank();

	int numberEl = elements->numberElements();
	int numberPoints =dim_+1;

	vec2D_dbl_ptr_Type vecPoints= inputMesh_->getPointsRepeated();

	int taggedElements=0;

	if(inputMesh_->getComm()->getRank() == 0){
		cout << "\t__________________________________________________________________________________________________________ " << endl;
		cout << " " << endl;
		cout << " \tThe area you requested for Refinement is :" << endl ;
			cout << "\tx in [" << area[0][0] << ", " << area[0][1] << "] " << endl;
		if(this->dim_>1)
			cout << "\ty in [" << area[1][0] << ", " << area[1][1] << "] " << endl;
		if(this->dim_>2)
			cout << "\tz in [" << area[2][0] << ", " << area[2][1] << "] " << endl;
		cout << "\t__________________________________________________________________________________________________________ " << endl;
	}
	vec_int_Type edgeNum(6); 
	edgeElements->matchEdgesToElements(elementMap);

	vec_dbl_Type point(this->dim_);
	int edgeNumber = 3*(this->dim_-1);

	LO p1ID, p2ID; 
	vec_dbl_Type P1,P2; 

	for(int i=0; i<numberEl; i++){
		edgeNum = edgeElements->getEdgesOfElement(i);
		// Checking whether Nodes of Element are part of the to be tagged area
		for(int k=0; k<numberPoints; k++){
				if(vecPoints->at(elements->getElement(i).getNode(k)).at(0) >= area[0][0] && vecPoints->at(elements->getElement(i).getNode(k)).at(0) <= area[0][1]){
					if(vecPoints->at(elements->getElement(i).getNode(k)).at(1) >= area[1][0] && vecPoints->at(elements->getElement(i).getNode(k)).at(1) <= area[1][1]){
						if(this->dim_>2){
							if(vecPoints->at(elements->getElement(i).getNode(k)).at(2) >= area[2][0] && vecPoints->at(elements->getElement(i).getNode(k)).at(2) <= area[2][1]){
									elements->getElement(i).tagForRefinement();
									k=numberPoints;
									taggedElements++;
							}
						}
						else if(this->dim_ == 2){
							elements->getElement(i).tagForRefinement();
							k=numberPoints;
							taggedElements++;		
						}
					
				}
			}
		}
		for(int k=0; k<edgeNumber ; k++){
			LO p1ID =edgeElements->getElement(edgeNum[k]).getNode(0);
			LO p2ID =edgeElements->getElement(edgeNum[k]).getNode(1);
			P1 = vecPoints->at(p1ID);
			P2 = vecPoints->at(p2ID);

			for (int d=0; d<this->dim_; d++){
				point[d]= ( (P1)[d] + (P2)[d] ) / 2.;
			} 
	  
			if(point.at(0) >= area[0][0] && point.at(0) <= area[0][1] && !elements->getElement(i).isTaggedForRefinement()){
					if(point.at(1) >= area[1][0] && point.at(1) <= area[1][1]){
						if(this->dim_>2){
							if(point.at(2) >= area[2][0] && point.at(2) <= area[2][1]){
									elements->getElement(i).tagForRefinement();
									taggedElements++;
							}
						}
						else if(this->dim_ == 2){
							elements->getElement(i).tagForRefinement();
							taggedElements++;		
						}
					
				}
			}

		}

	}
	reduceAll<int, int> (*inputMesh_->getComm(), REDUCE_MAX, taggedElements, outArg (taggedElements));

	if(inputMesh_->getComm()->getRank()==0){
		cout << "\t__________________________________________________________________________________________________________ " << endl;
		cout << " " << endl;
		cout << "\t With the 'tagArea tool' " << taggedElements << " Elements were tagged for Refinement " << endl;
		cout << "\t__________________________________________________________________________________________________________ " << endl;
	}

}

/*!
\brief Tags only a certain Area for refinement and is independent of any error estimation.

@param[in] inputMeshP1 The P1 Mesh that is used for later refinement.
@param[in] area  Area that is suppose to be refined. If is a vector defining the area as follows: row1:[x_0,x_1] x-limits, row2: [y_0,y_1] y-limits, row3: [z_0,z_1] z-limits .
*/

template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::tagAll( MeshUnstrPtr_Type inputMeshP1){

	inputMesh_ = inputMeshP1;

	ElementsPtr_Type elements = inputMesh_->getElementsC();
   	EdgeElementsPtr_Type edgeElements = inputMesh_->getEdgeElements();
	MapConstPtr_Type elementMap = inputMesh_->getElementMap();
	int myRank = inputMesh_->comm_->getRank();

	int numberEl = elements->numberElements();
	int numberPoints =dim_+1;

	vec2D_dbl_ptr_Type vecPoints= inputMesh_->getPointsRepeated();

	int taggedElements=0;

	for(int i=0; i<numberEl; i++){
		elements->getElement(i).tagForRefinement();
		taggedElements++;						
	}
		
	reduceAll<int, int> (*inputMesh_->getComm(), REDUCE_MAX, taggedElements, outArg (taggedElements));

	if(inputMesh_->getComm()->getRank()==0){
		cout << "\t__________________________________________________________________________________________________________ " << endl;
		cout << " " << endl;
		cout << "\tWith tag all:' " << taggedElements << " Elements were tagged for Refinement " << endl;
		cout << "\t__________________________________________________________________________________________________________ " << endl;
	}

}

// Tagging all elements adjacent to the flag.
template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::tagFlag( MeshUnstrPtr_Type inputMeshP1,int flag){

	inputMesh_ = inputMeshP1;

	ElementsPtr_Type elements = inputMesh_->getElementsC();
   	EdgeElementsPtr_Type edgeElements = inputMesh_->getEdgeElements();
	MapConstPtr_Type elementMap = inputMesh_->getElementMap();
	int myRank = inputMesh_->comm_->getRank();

	int numberEl = elements->numberElements();
	int numberPoints =dim_+1;

	vec2D_dbl_ptr_Type vecPoints= inputMesh_->getPointsRepeated();

	int taggedElements=0;

	for(int i=0; i<numberEl; i++){
		FiniteElement fe = elements->getElement( i );
		if(fe.getFlag() == flag){
			elements->getElement(i).tagForRefinement();
			taggedElements++;
		}
		else{
			ElementsPtr_Type subEl = fe.getSubElements(); // might be null
			for (int surface=0; surface<fe.numSubElements(); surface++) {
				FiniteElement feSub = subEl->getElement( surface  );
				if(feSub.getFlag() == flag){
					elements->getElement(i).tagForRefinement();
					taggedElements++;
				}
			}
		}
	}
		
	reduceAll<int, int> (*inputMesh_->getComm(), REDUCE_MAX, taggedElements, outArg (taggedElements));

	if(inputMesh_->getComm()->getRank()==0){
		cout << "\t__________________________________________________________________________________________________________ " << endl;
		cout << " " << endl;
		cout << "\tWith tag all:' " << taggedElements << " Elements were tagged for Refinement " << endl;
		cout << "\t__________________________________________________________________________________________________________ " << endl;
	}

}

/*!
\brief Part of the error estimator that calculates the jump part of the estimation. What kind of jump is calculated depends on the problemType we have at hand.
*/
template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::calculateJump(){

	int surfaceNumbers = dim_+1 ; // Number of (dim-1) - dimensional surfaces of Element (triangle has 3 edges)	

	// Necessary mesh objects
	ElementsPtr_Type elements = inputMesh_->getElementsC();
    EdgeElementsPtr_Type edgeElements = inputMesh_->getEdgeElements();
	MapConstPtr_Type elementMap = inputMesh_->getElementMap();
	int myRank = inputMesh_->comm_->getRank();
	SurfaceElementsPtr_Type surfaceTriangleElements = this->surfaceElements_;
	vec2D_dbl_ptr_Type points = inputMesh_->pointsRep_;

	int numberSurfaceElements=0;
	if(dim_==2){
		numberSurfaceElements = edgeElements->numberElements();
	}
	else if(dim_==3){
		numberSurfaceElements = surfaceTriangleElements->numberElements();
		
	}

	vec_dbl_Type surfaceElementsError(numberSurfaceElements);

	// We calculate the Gradient for all edges or surfaces in direction of N. This is done for the gradient of the velocity u or the pressure p.
	vec3D_dbl_Type u_El = calcNPhi("Gradient",  dofs_, FEType2_);

	// We generally do not need to calculate this part of the jump as the pressure is steady over the edges resulting in a zero pressure Jump.
	vec3D_dbl_Type p_El;
	//if(calculatePressure_ == true)
		//p_El = calcNPhi("None",  dofsP_, FEType1_);

	double quadWeightConst =1.;
	
	if(this->dim_ == 2){

		// necessary entities
		// calculating diameter of elements	
		vec_dbl_Type p1_2(2); // normal Vector of Surface

		double h_E ; // something with edge
		vec_dbl_Type v_E(2); // normal Vector of edges
		double norm_v_E;

		for(int k=0; k< edgeElements->numberElements();k++){
			// Normalenvektor bestimmen:
			p1_2[0] = points->at(edgeElements->getElement(k).getNode(0)).at(0) - points->at(edgeElements->getElement(k).getNode(1)).at(0);
			p1_2[1] = points->at(edgeElements->getElement(k).getNode(0)).at(1) - points->at(edgeElements->getElement(k).getNode(1)).at(1);
	
			double jumpQuad=0;
			for(int j=0; j<dofs_ ; j++){
				for(int i=0; i<u_El[k].size();i++){
					if(	calculatePressure_ == false)
						jumpQuad += std::pow(u_El[k][i][j],2);
					if(calculatePressure_ == true){
						jumpQuad += std::pow( u_El[k][i][j] ,2) ; //u_El[k][i][j] - p_El[k][i][0],2)
					}
				}
			}

			h_E =  std::sqrt(std::pow(p1_2[0],2)+std::pow(p1_2[1],2));

			if(this->FEType2_ =="P2")
				quadWeightConst = h_E /6.;
			else 	
				quadWeightConst = h_E;

			surfaceElementsError[k] =h_E *quadWeightConst*jumpQuad;


		}
	}

	if(this->dim_==3){
		// Edge Numbers of Element
		vec_int_Type surfaceElementsIds(surfaceNumbers);

		vec_dbl_Type areaTriangles(numberSurfaceElements);
		vec_dbl_Type tmpVec1(numberSurfaceElements);
		vec_dbl_Type tmpVec2(numberSurfaceElements);
		// necessary entities
		vec_dbl_Type p1(3),p2(3); // normal Vector of Surface
		vec_dbl_Type v_E(3); // normal Vector of Surface
		double norm_v_E;

		vec_dbl_Type h_Tri =  calcDiamTriangles3D(surfaceTriangleElements,points, areaTriangles, tmpVec1, tmpVec2);

		for(int k=0; k< surfaceTriangleElements->numberElements();k++){
			FiniteElement surfaceTmp = surfaceTriangleElements->getElement(k);
			double jumpQuad=0.;
			for(int j=0; j<dofs_ ; j++){
				for(int i=0; i<u_El[k].size();i++){
					if(	calculatePressure_ == false)
						jumpQuad += std::pow(u_El[k][i][j],2);
					if(calculatePressure_ == true)
						jumpQuad += std::pow(u_El[k][i][j],2); // -p_El[k][i][0]
				}
			}
			if(this->FEType2_ =="P2")
				quadWeightConst = areaTriangles[k];
			else 	
				quadWeightConst = areaTriangles[k];


			surfaceElementsError[k] = h_Tri[k]*jumpQuad*quadWeightConst; 
		}
			
	}
	return surfaceElementsError;

}

/*!
\brief Function that calculates the jump part for nabla u or p 

@param[in] phiDerivative phiDerivative is either 'Gradient' or 'None' and what kind of jump is calculated depends on the problemType we have at hand. If phiDerivative is 'Gradient' the nabla u jump part is caluculated and if its 'None' then the pressure jump.
@param[in] dofsSol Degree of freedom of the caluclated jump part. The pressure dof is always 1 whereas velocity dof can vary depending on the problem.
@param[in] FEType Finite element type of the calculated jump part. 
 
*/
template <class SC, class LO, class GO, class NO>
vec3D_dbl_Type ErrorEstimation<SC,LO,GO,NO>::calcNPhi(std::string phiDerivative, int dofsSol, std::string FEType){

	int surfaceNumbers = dim_+1 ; // Number of (dim-1) - dimensional surfaces of Element (triangle has 3 edges)		

	// necessary mesh objects
	vec2D_dbl_ptr_Type points = inputMesh_->pointsRep_;
	ElementsPtr_Type elements = inputMesh_->getElementsC();	
    EdgeElementsPtr_Type edgeElements = inputMesh_->getEdgeElements();
	MapConstPtr_Type elementMap = inputMesh_->getElementMap();
	SurfaceElementsPtr_Type surfaceTriangleElements = this->surfaceElements_;
	vec_int_Type surfaceElementsIds(surfaceNumbers);
	MapConstPtr_Type surfaceMap;

	// number of surface elements depending on dimension
	int numberSurfaceElements=0;
	if(dim_==2){
		numberSurfaceElements = edgeElements->numberElements();
		surfaceMap= inputMesh_->edgeMap_;
	}
	else if(dim_==3){
		numberSurfaceElements = surfaceTriangleElements->numberElements();
		surfaceMap= this->surfaceTriangleMap_;
	}

	// the int-version of the fe discretisation
	int intFE = 1;
	int numNodes= dim_+1;
	int quadPSize = 1;

	if(FEType2_ == "P2"){
		quadPSize=3; // Number of Surface Quadrature Points
	}
	if(FEType == "P2"){
		numNodes=6;
		intFE =2;
		if(dim_ ==3){
			numNodes=10;
		}
	}

	// We determine u for each Quad Point, so we need to determine how many Points we have
	vec3D_dbl_Type vec_El(numberSurfaceElements,vec2D_dbl_Type(quadPSize,vec_dbl_Type(dofsSol))); // return value
	vec3D_dbl_Type vec_El1(numberSurfaceElements,vec2D_dbl_Type(quadPSize,vec_dbl_Type(dofsSol))); // as we calculate a jump over a surface we generally have two solutions for each side
	vec3D_dbl_Type vec_El2(numberSurfaceElements,vec2D_dbl_Type(quadPSize,vec_dbl_Type(dofsSol)));

	vec3D_dbl_Type vec_El_Exp(numberSurfaceElements,vec2D_dbl_Type(quadPSize,vec_dbl_Type(dofsSol))); // the values of our elements as procs we communicate

	// vectors for late map building
	vec_GO_Type elementImportIDs(0); 
	vec_GO_Type elementExportIDs(0);
	vec_LO_Type surfaceIDsLocal(0);

	// gradient of u elementwise
	vec_LO_Type kn1(numNodes);

    SC detB1;
    SC absDetB1;
    SmallMatrix<SC> B1(dim_);
    SmallMatrix<SC> Binv1(dim_);  
	int index0,index;

	edgeElements->matchEdgesToElements(elementMap);

	// elementsOfSurfaceGlobal and -Local for determining the communication
	vec2D_GO_Type elementsOfSurfaceGlobal; 
	vec2D_LO_Type elementsOfSurfaceLocal;

	if(dim_==2){
		elementsOfSurfaceGlobal = edgeElements->getElementsOfEdgeGlobal();	
 		elementsOfSurfaceLocal = edgeElements->getElementsOfEdgeLocal();
	}
	else if(dim_ ==3){
		elementsOfSurfaceGlobal = surfaceTriangleElements->getElementsOfSurfaceGlobal();
		elementsOfSurfaceLocal	= surfaceTriangleElements->getElementsOfSurfaceLocal();
	}

	// the normal vector and its norm
	vec_dbl_Type v_E(dim_);
	double norm_v_E=1.;

	// We loop through all surfaces in order to calculate nabla u or p on each surface depending on which element we look at.
	// The jump is calculated via two surfaces. Consequently we have two values per edge/surface which we either have or need to import/export.
	for(int k=0; k< numberSurfaceElements;k++){
		// We only need to calculate the jump for interior egdes/surfaces, which are characetrized by the fact that they are connected to two elements
		bool interiorElement= false;
		if(dim_ == 2){
			if((edgeElements->getElement(k).getFlag() == 10 || edgeElements->getElement(k).getFlag() == 0))//&& elementsOfSurfaceLocal.at(k).size()>1)
				interiorElement=true;
		}  
		else if(dim_ == 3){
			if((surfaceTriangleElements->getElement(k).getFlag() == 10 || surfaceTriangleElements->getElement(k).getFlag() == 0))//&& elementsOfSurfaceLocal.at(k).size()>1)
				interiorElement=true;
		}  
		if(interiorElement == true && elementsOfSurfaceLocal.at(k).size() > 1){  
		// Per edge we have depending on discretization quadrature points and weights
			vec_dbl_Type quadWeights(quadPSize);
			vec2D_dbl_Type quadPoints; 
		
			if(dim_ == 2){ 
				quadPoints = getQuadValues(dim_, FEType2_ , "Surface", quadWeights, edgeElements->getElement(k));
				v_E.at(0) = points->at(edgeElements->getElement(k).getNode(0)).at(1) - points->at(edgeElements->getElement(k).getNode(1)).at(1);
				v_E.at(1) = -(points->at(edgeElements->getElement(k).getNode(0)).at(0) - points->at(edgeElements->getElement(k).getNode(1)).at(0));
				norm_v_E = std::sqrt(std::pow(v_E[0],2)+std::pow(v_E[1],2));	
			}
			else if(dim_ == 3){
				quadPoints = getQuadValues(dim_, FEType2_ , "Surface", quadWeights, surfaceTriangleElements->getElement(k));
				// Normalenvektor bestimmen:
				vec_dbl_Type p1(3),p2(3);
				p1[0] = points->at(surfaceTriangleElements->getElement(k).getNode(0)).at(0) - points->at(surfaceTriangleElements->getElement(k).getNode(1)).at(0);
				p1[1] = points->at(surfaceTriangleElements->getElement(k).getNode(0)).at(1) - points->at(surfaceTriangleElements->getElement(k).getNode(1)).at(1);
				p1[2] = points->at(surfaceTriangleElements->getElement(k).getNode(0)).at(2) - points->at(surfaceTriangleElements->getElement(k).getNode(1)).at(2);

				p2[0] = points->at(surfaceTriangleElements->getElement(k).getNode(0)).at(0) - points->at(surfaceTriangleElements->getElement(k).getNode(2)).at(0);
				p2[1] = points->at(surfaceTriangleElements->getElement(k).getNode(0)).at(1) - points->at(surfaceTriangleElements->getElement(k).getNode(2)).at(1);
				p2[2] = points->at(surfaceTriangleElements->getElement(k).getNode(0)).at(2) - points->at(surfaceTriangleElements->getElement(k).getNode(2)).at(2);

				v_E[0] = p1[1]*p2[2] - p1[2]*p2[1];
				v_E[1] = p1[2]*p2[0] - p1[0]*p2[2];
				v_E[2] = p1[0]*p2[1] - p1[1]*p2[0];

				norm_v_E = std::sqrt(std::pow(v_E[0],2)+std::pow(v_E[1],2)+std::pow(v_E[2],2));	  
			}
	
			vec_LO_Type elementsIDs(0);
			// Case that both necessary element are on the same Proc
			if(elementsOfSurfaceLocal.at(k).at(0) != -1){
				elementsIDs.push_back(elementsOfSurfaceLocal.at(k).at(0));
			}
			if(elementsOfSurfaceLocal.at(k).at(1) != -1){
				elementsIDs.push_back(elementsOfSurfaceLocal.at(k).at(1));
			}
			for(int i=0; i<elementsIDs.size(); i++){

				// We extract the nodes of the elements the surface is connected to

			kn1= elements->getElement(elementsIDs[i]).getVectorNodeListNonConst();
			
			// Transformation Matrices
			// We need this inverse Matrix to also transform the quadrature points of our surface back to the reference element
			index0 = kn1[0];
			for (int s=0; s<dim_; s++) {
				index = kn1[s+1];
				for (int t=0; t<dim_; t++) {
					B1[t][s] = points->at(index).at(t) -points->at(index0).at(t);
				}
			}

			detB1 = B1.computeInverse(Binv1);
			detB1 = std::fabs(detB1);

			vec2D_dbl_Type quadPointsT1(quadPSize,vec_dbl_Type(dim_));
			for(int l=0; l< quadPSize; l++){
				for(int p=0; p< dim_ ; p++){
					for(int q=0; q< dim_; q++){
						quadPointsT1[l][p] += Binv1[p][q]* (quadPoints[l][q] - points->at(elements->getElement(elementsIDs[i]).getNode(0)).at(q))  ; 
					}
				}
									
			}
				// We make the distinction between a gradient jump calculation or a simple jump calculation 
				for(int l=0; l< quadPSize; l++){
					if(phiDerivative == "Gradient"){
						vec2D_dbl_Type deriPhi1 = gradPhi(dim_, intFE, quadPointsT1[l]);
						vec2D_dbl_Type deriPhiT1(numNodes,vec_dbl_Type(dim_));
						for(int q=0; q<dim_; q++){
							for(int p=0;p<numNodes; p++){
								for(int s=0; s< dim_ ; s++)
									deriPhiT1[p][q] += (deriPhi1[p][s]*Binv1[s][q]);
							}
						}
						vec2D_dbl_Type u1_Tmp(dofsSol, vec_dbl_Type(dim_));
						Teuchos::ArrayRCP< SC > valuesSolRep;	
						for( int d =0; d < dofsSol ; d++){
							valuesSolRep = valuesSolutionRepVel_->getBlock(d)->getDataNonConst(0);
							for(int s=0; s<dim_; s++){
								for(int t=0; t< numNodes ; t++){
					 				u1_Tmp[d][s] += valuesSolRep[kn1[t]]*deriPhiT1[t][s] ;
								}
							}

						}
						if(i==0){
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El1[k][l][d] += (u1_Tmp[d][s])*(v_E[s]/norm_v_E);
								}
								vec_El1[k][l][d] = std::sqrt(quadWeights[l])* vec_El1[k][l][d];
							}
						}
						else {
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El2[k][l][d] += (u1_Tmp[d][s])*(v_E[s]/norm_v_E);
								}
								vec_El2[k][l][d] = std::sqrt(quadWeights[l])* vec_El2[k][l][d];
							}
						}	


					}

					if(phiDerivative == "None"){
						vec_dbl_Type phiV = phi(dim_, intFE, quadPointsT1[l]);
						vec_dbl_Type vec_Tmp(dofsSol);
						Teuchos::ArrayRCP< SC > valuesSolRep;	
						for( int d =0; d < dofsSol ; d++){
							valuesSolRep = valuesSolutionRepPre_->getBlock(d)->getDataNonConst(0);
							for(int t=0; t< phiV.size() ; t++){
				 				vec_Tmp[d] += valuesSolRep[kn1[t]]*phiV[t] ;
							}
						}
						if(i==0){
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El1[k][l][d] += (vec_Tmp[d])*(v_E[s]/norm_v_E);
								}
								vec_El1[k][l][d] = std::sqrt(quadWeights[l])* vec_El1[k][l][d];
							}
						}
						else {
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El2[k][l][d] += (vec_Tmp[d])*(v_E[s]/norm_v_E);
								}
								vec_El2[k][l][d] = std::sqrt(quadWeights[l])* vec_El2[k][l][d];
							}
						}	

					}
	
				}
			}
			// If one of out local elementsOfSurface is -1 it means, that one element connected to the surface is on another processor an we later need to exchange information
			if(elementsOfSurfaceLocal.at(k).at(0) == -1  || elementsOfSurfaceLocal.at(k).at(1) == -1){

				elementImportIDs.push_back(surfaceMap->getGlobalElement(k));
		
				for(int l=0; l< vec_El1[k].size() ;l++){
					for( int d =0; d < dofsSol ; d++)
						vec_El_Exp[k][l][d] = vec_El1[k][l][d];
				}	
			}	
		}
	}
	inputMesh_->getComm()->barrier();

	
	// We construct map from the previously extracted elementImportIDs
	sort(elementImportIDs.begin(),elementImportIDs.end());

	vec_GO_Type::iterator ip = unique( elementImportIDs.begin(), elementImportIDs.end());

	elementImportIDs.resize(distance(elementImportIDs.begin(), ip)); 
	Teuchos::ArrayView<GO> globalElementArrayImp = Teuchos::arrayViewFromVector( elementImportIDs);

	// Map which represents the surface ids that need to import element information
	MapPtr_Type mapElementImport =
		Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), globalElementArrayImp, 0, inputMesh_->getComm()) );


	int maxRank = std::get<1>(inputMesh_->rankRange_);
	MapPtr_Type mapElementsUnique = mapElementImport;


	if(maxRank>0 && mapElementsUnique->getMaxAllGlobalIndex() > 0)
		mapElementsUnique = mapElementImport->buildUniqueMap( inputMesh_->rankRange_ );

	// In case of a vector values solution we need to import/export dofs-many entries and all of those entries for each of the quadpoints

	MultiVectorPtr_Type valuesU_x_expU = Teuchos::rcp( new MultiVector_Type( mapElementsUnique, 1 ) );	
	MultiVectorPtr_Type valuesU_y_expU = Teuchos::rcp( new MultiVector_Type( mapElementsUnique, 1 ) );
	MultiVectorPtr_Type valuesU_z_expU = Teuchos::rcp( new MultiVector_Type( mapElementsUnique, 1 ) );

	MultiVectorPtr_Type valuesU_x_impU = Teuchos::rcp( new MultiVector_Type( mapElementsUnique, 1 ) );	
	MultiVectorPtr_Type valuesU_y_impU = Teuchos::rcp( new MultiVector_Type( mapElementsUnique, 1 ) );
	MultiVectorPtr_Type valuesU_z_impU = Teuchos::rcp( new MultiVector_Type( mapElementsUnique, 1 ) );


	MultiVectorPtr_Type valuesU_x_imp = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );	
	MultiVectorPtr_Type valuesU_y_imp = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );
	MultiVectorPtr_Type valuesU_z_imp = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );

	MultiVectorPtr_Type valuesU_x_impF = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );	
	MultiVectorPtr_Type valuesU_y_impF = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );
	MultiVectorPtr_Type valuesU_z_impF = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );

	// Array Pointers which will contain the to be exchanges information
	Teuchos::ArrayRCP< SC > entriesU_x_imp  = valuesU_x_imp->getDataNonConst(0);
	Teuchos::ArrayRCP< SC > entriesU_y_imp  = valuesU_y_imp->getDataNonConst(0);
	Teuchos::ArrayRCP< SC > entriesU_z_imp  = valuesU_z_imp->getDataNonConst(0);

	Teuchos::ArrayRCP< SC > entriesU_x_impF  = valuesU_x_impF->getDataNonConst(0);
	Teuchos::ArrayRCP< SC > entriesU_y_impF  = valuesU_y_impF->getDataNonConst(0);
	Teuchos::ArrayRCP< SC > entriesU_z_impF  = valuesU_z_impF->getDataNonConst(0);

	// We exchange values per quad point
	for(int i=0; i<quadPSize; i++){
		for(int j=0; j<elementImportIDs.size(); j++){
			entriesU_x_imp[j] = vec_El_Exp[surfaceMap->getLocalElement(elementImportIDs[j])][i][0] ;
			if(dofsSol == 2)
				entriesU_y_imp[j] = vec_El_Exp[surfaceMap->getLocalElement(elementImportIDs[j])][i][1] ;
			else if(dofsSol == 3)
				entriesU_z_imp[j] = vec_El_Exp[surfaceMap->getLocalElement(elementImportIDs[j])][i][2] ;
		}

		valuesU_x_impU->importFromVector(valuesU_x_imp, false, "Insert");
		valuesU_y_impU->importFromVector(valuesU_y_imp, false, "Insert");
		valuesU_z_impU->importFromVector(valuesU_z_imp, false, "Insert");

		valuesU_x_expU->exportFromVector(valuesU_x_imp, false, "Insert");
		valuesU_y_expU->exportFromVector(valuesU_y_imp, false, "Insert");
		valuesU_z_expU->exportFromVector(valuesU_z_imp, false, "Insert");

		valuesU_x_impF->importFromVector(valuesU_x_impU, false, "Insert");
		valuesU_y_impF->importFromVector(valuesU_y_impU, false, "Insert");
		valuesU_z_impF->importFromVector(valuesU_z_impU, false, "Insert");

		valuesU_x_impF->exportFromVector(valuesU_x_expU, false, "Insert");
		valuesU_y_impF->exportFromVector(valuesU_y_expU, false, "Insert");
		valuesU_z_impF->exportFromVector(valuesU_z_expU, false, "Insert");

		for(int j=0; j<elementImportIDs.size(); j++){
			vec_El2[surfaceMap->getLocalElement(elementImportIDs[j])][i][0] =entriesU_x_impF[j];
			if(dofsSol == 2)
				vec_El2[surfaceMap->getLocalElement(elementImportIDs[j])][i][1] =entriesU_y_impF[j];
			else if(dofsSol ==3)
				vec_El2[surfaceMap->getLocalElement(elementImportIDs[j])][i][2] = entriesU_z_impF[j];
		}	
	}
	for(int i=0; i< numberSurfaceElements ; i++){
		for(int j=0; j<quadPSize; j++){
			for(int k=0; k< dofsSol ; k++){
				vec_El[i][j][k] = std::fabs(vec_El1[i][j][k]) - std::fabs(vec_El2[i][j][k]); 
			
			}
		}
	}


	return vec_El;

}

/*!
\brief Function that marks the elements for refinement.

@param[in] errorElementMv MultiVector that contains the estimated error for each element.
@param[in] theta Parameter determining for marking strategies.
@param[in] markingStrategy Strategy with which the elements are marked. Implemented Strategies 'Doerfler' or 'Maximum'. 
@param[in] meshP1 P1 mesh which is used for later refinement and has to be the one beeing marked.

 \brief !! it is essential that the meshP1 mesh inserted here is the mesh that will be used for mesh refinement, as it contains the elementwise-information determining refinement. !!
*/

template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::markElements(MultiVectorPtr_Type errorElementMv, double theta, std::string strategy,  MeshUnstrPtr_Type meshP1){


	Teuchos::ArrayRCP<SC> errorEstimation = errorElementMv->getDataNonConst(0);

	ElementsPtr_Type elements = meshP1->getElementsC();	

	this->markingStrategy_ = strategy;

	theta_ = theta;
	// As we decide which element to tag based on the maximum error in the elements globally, we need to communicated this maxErrorEl
	auto it = std::max_element(errorEstimation.begin(), errorEstimation.end()); // c++11
	double maxErrorEl= errorEstimation[std::distance(errorEstimation.begin(), it)];

	// Maximum-strategy for marking the elements as proposed by Verfürth in "A posterio error estimation"
    reduceAll<int, double> (*meshP1->getComm(), REDUCE_MAX, maxErrorEl, outArg (maxErrorEl));
	int flagCount=0;
	if( markingStrategy_ == "Maximum"){
		for(int k=0; k< elements->numberElements() ; k++){
			if( errorEstimation[k] > theta_ * maxErrorEl){
				elements->getElement(k).tagForRefinement();
				flagCount ++;
				}
			}
	}
	// Equilibirum/Doerfler-strategy  for marking the elements as proposed by Verfürth in "A posterio error estimation"
	else if(markingStrategy_ == "Doerfler"){
		double thetaSumTmp=0.;
		double thetaSum=0.;
		double muMax=0.;
		double sigmaT=0.;
		vec_bool_Type flagged(elements->numberElements());
		for(int k=0; k< elements->numberElements(); k++){
			thetaSumTmp = thetaSumTmp + std::pow(errorEstimation[k],2);
			flagged[k]=false;
		}
    	reduceAll<int, double> (*meshP1->getComm(), REDUCE_SUM, thetaSumTmp, outArg (thetaSum));
		while(sigmaT < theta_*thetaSum){
			muMax=0.;
			for(int k=0; k< elements->numberElements(); k++){
				if(muMax < errorEstimation[k] && flagged[k]==false ){
					muMax = errorEstimation[k];
				}
			}
    		reduceAll<int, double> (*meshP1->getComm(), REDUCE_MAX, muMax, outArg (muMax));
			for(int k=0; k< elements->numberElements(); k++){
				if(muMax == errorEstimation[k] && flagged[k]==false ){
					flagged[k]=true;
					sigmaT = sigmaT + std::pow(errorEstimation[k],2);
					}
			}
    	reduceAll<int, double> (*meshP1->getComm(), REDUCE_MAX, sigmaT, outArg (sigmaT));
		}
			
	   	for(int k=0; k< elements ->numberElements() ; k++){
			if( flagged[k]==true){
				elements->getElement(k).tagForRefinement();
				flagCount++;
				}
		}  

	}	
		
 	// If no strategy is chosen or we choose uniform refinement ourselves, all elements are marked for refinement
	else{ 
		 for(int k=0; k< elements->numberElements() ; k++){
				elements->getElement(k).tagForRefinement();
				flagCount++;
			}
	}
    reduceAll<int, int> (*meshP1->getComm(), REDUCE_SUM, flagCount , outArg (flagCount));

	if(meshP1->getComm()->getRank() == 0){
	cout << " " << endl;
	cout << " \tA-posteriori Error Estimation for " << problemType_ << "-problem tagged " << flagCount << " Elements for adaptive Refinement with " <<  markingStrategy_ << "-Strategy and Theta= " << theta_ << endl;
	cout << " " << endl;
	}
  
}

/*!
\brief Function that that determines || div(u) ||_T for a Element T.

@param[in] element FiniteElement element where ||div(u)||_T is calculated on.

@param[out] divElement Divergence of u_h on element

*/
template <class SC, class LO, class GO, class NO>
double ErrorEstimation<SC,LO,GO,NO>::determineDivU(FiniteElement element){

// Quad Points and weights of third order
	double divElement =0.;

	int dim = this->dim_;
    SC* paras ; //= &(funcParameter[0]);

	int t=1;
	if(dim == 2)
		t =3;
	else if(dim == 3)
		t=5;

	vec_dbl_Type QuadW(t);
	vec2D_dbl_Type QuadPts = getQuadValues(dim, FEType2_, "Element", QuadW, element);


	vec2D_dbl_ptr_Type points = inputMesh_->getPointsRepeated();

	// Transformation Matrices
	// We determine deltaU for the Elements. If FEType=P1 deltaU equals 0. If FEType = P2 we get a constant solution:
	double deltaU =0.;
	vec_LO_Type nodeList = element.getVectorNodeListNonConst();

	SC detB1;
    SmallMatrix<SC> B1(dim);
    SmallMatrix<SC> Binv1(dim);  
 
	// We need this inverse Matrix to also transform the quad Points back to the reference element
	int index0 = nodeList[0];
	 for (int s=0; s<dim; s++) {
		int index = nodeList[s+1];
		for (int t=0; t<dim; t++) {
			B1[t][s] = points->at(index).at(t) - points->at(index0).at(t);
		}
	}


	
	detB1 = B1.computeInverse(Binv1);
	detB1 = std::fabs(detB1);

	int intFE = 1;
	if(this->FEType2_ == "P2")
		intFE =2;


	double divElementTmp=0.;
	for (UN w=0; w< QuadW.size(); w++){
		vec2D_dbl_Type gradPhiV =gradPhi(dim, intFE, QuadPts[w]);
		vec2D_dbl_Type gradPhiT(gradPhiV.size(),vec_dbl_Type(dim));
		for(int q=0; q<dim; q++){
			for(int p=0;p<nodeList.size(); p++){
				for(int s=0; s< dim ; s++)
					gradPhiT[p][q] += (gradPhiV[p][s]*Binv1[s][q]);			
			}
		}
		vec_dbl_Type uTmp(gradPhiV.size());
		for(int j=0; j< dofs_ ; j++){
			Teuchos::ArrayRCP<SC> valuesSolutionRep = valuesSolutionRepVel_->getBlock(j)->getDataNonConst(0);
			for(int i=0; i< nodeList.size(); i++){
				uTmp[i] += valuesSolutionRep[nodeList[i]]*gradPhiT[i][j];
			}
		}
		divElementTmp=0.;
		for(int i=0; i< nodeList.size(); i++){
			divElementTmp += uTmp[i];
		}
		divElementTmp = std::pow(divElementTmp,2);
		divElementTmp *= QuadW[w];

		divElement += divElementTmp;
	}
   	divElement =  divElement *detB1;



	return divElement;


}

/*!
\brief Function that that determines ||\Delta u_h + f ||_(L2(T)), || \Delta u_h + f - \nabla p_h ||_T or || \Delta u_h + f - \nabla p_h - (u_h \cdot \nabla)u_h||_T for an Element T

@param[in] element FiniteElement element where ||div(u)||_T is calculated on.
@param[in] rhsFunc The right hand side function of the pde.

@param[out] resElement Residual of element according to problemType 
*/

template <class SC, class LO, class GO, class NO>
double ErrorEstimation<SC,LO,GO,NO>::determineResElement(FiniteElement element, RhsFunc_Type rhsFunc){

// Quad Points and weights of third order

	vec_dbl_Type resElement(dofs_);

	int dim = this->dim_;
    SC* paras ;

	int t=1;
	if(dim == 2)
		t =3;
	else if(dim == 3)
		t=5;

	vec_dbl_Type QuadW(t);
	vec2D_dbl_Type QuadPts = getQuadValues(dim, FEType2_, "Element", QuadW, element);

	int intFE = 1;
	if(this->FEType2_ == "P2")
		intFE =2;

	vec2D_dbl_ptr_Type points = inputMesh_->getPointsRepeated();

	// Transformation Matrices
	// We determine deltaU for the Elements. If FEType=P1 deltaU equals 0. If FEType = P2 we get a constant solution:
	vec_LO_Type nodeList = 	element.getVectorNodeListNonConst();

	SC detB1;
    SmallMatrix<SC> B1(dim);
    SmallMatrix<SC> Binv1(dim);  
 
	// We need this inverse Matrix to also transform the quad Points of on edge back to the reference element
	int index0 = nodeList[0];
	 for (int s=0; s<dim; s++) {
		int index = nodeList[s+1];
		for (int t=0; t<dim; t++) {
			B1[t][s] = points->at(index).at(t) - points->at(index0).at(t);
		}
	}


	detB1 = B1.computeInverse(Binv1);
	detB1 = std::fabs(detB1);

	vec_dbl_Type valueFunc(dofs_);

	vec2D_dbl_Type quadPointsTrans(QuadW.size(),vec_dbl_Type(dim));

	for(int i=0; i< QuadW.size(); i++){
		 for(int j=0; j< dim; j++){
			for(int k=0; k< dim; k++){
		 		quadPointsTrans[i][j] += B1[j][k]* QuadPts[i].at(k) ; 
			} 
			quadPointsTrans[i][j] += points->at(element.getNode(0)).at(j);
		 }
	}
	
	vec_dbl_Type nablaP(dim);
	if(problemType_ == "Stokes" || problemType_ == "NavierStokes"){
		Teuchos::ArrayRCP<SC> valuesSolutionRepPre = valuesSolutionRepPre_->getBlock(0)->getDataNonConst(0);
		vec2D_dbl_Type deriPhi = gradPhi(dim, 1, QuadPts[0]);
		vec2D_dbl_Type deriPhiT(dim+1,vec_dbl_Type(dim));
		for(int q=0; q<dim; q++){
			for(int p=0;p<dim+1; p++){
				for(int s=0; s< dim ; s++)
					deriPhiT[p][q] += (deriPhi[p][s]*Binv1[s][q]);
				
			}
		}
		for(int s=0; s<dim; s++){
			for(int t=0; t< dim+1 ; t++){
				nablaP[s] += valuesSolutionRepPre[nodeList[t]]*deriPhiT[t][s] ;
			}
		}
	}
	// ####################################################################
	vec2D_dbl_Type nablaU(dim,vec_dbl_Type(QuadW.size()));
	if(problemType_ == "NavierStokes"){
		for (UN w=0; w< QuadW.size(); w++){
			vec2D_dbl_Type deriPhi = gradPhi(dim, intFE, QuadPts[w]);
			vec2D_dbl_Type deriPhiT(nodeList.size(),vec_dbl_Type(dim));
			for(int q=0; q<dim; q++){
				for(int p=0;p<nodeList.size(); p++){
					for(int s=0; s< dim ; s++)
						deriPhiT[p][q] += (deriPhi[p][s]*Binv1[s][q]);
					
				}
			}
			Teuchos::ArrayRCP< SC > valuesSolRep;	
			vec2D_dbl_Type vec_Tmp(nodeList.size(),vec_dbl_Type(dofs_));
			for( int d =0; d < dofs_ ; d++){
				valuesSolRep = valuesSolutionRepVel_->getBlock(d)->getDataNonConst(0);
				for(int t=0; t< nodeList.size() ; t++){
	 				vec_Tmp[t][d] = valuesSolRep[nodeList[t]];
				}	
			}

			vec2D_dbl_Type u1_Tmp(dofs_, vec_dbl_Type(dim_));

			for( int d =0; d < dofs_ ; d++){
				valuesSolRep = valuesSolutionRepVel_->getBlock(d)->getDataNonConst(0);
				for(int s=0; s<dim; s++){
					for(int t=0; t< nodeList.size() ; t++){
		 				u1_Tmp[d][s] += valuesSolRep[nodeList[t]]*deriPhiT[t][s] ;
					}
				}
			}
			vec_dbl_Type phiV = phi(dim, intFE, QuadPts[w]);
			for( int d =0; d < dofs_ ; d++){
				for(int t=0; t< nodeList.size() ; t++){
					for(int s=0; s<dim; s++){
						nablaU[d][w] += vec_Tmp[t][s]*phiV[t]*u1_Tmp[s][d];
					}
				}
			}

		}
	}

	vec_dbl_Type deltaU(dofs_);
	if(this->FEType2_ == "P2"){
		vec_dbl_Type deltaPhi(nodeList.size());
		if(this->dim_ == 2){			
			deltaPhi={8, 4, 4, -8, 0, -8 };
		}
		else if(this->dim_ == 3)
			deltaPhi={12, 4, 4, 4, -8, 0, -8, -8, 0, 0};

		for(int j=0 ; j< dofs_; j++){
			Teuchos::ArrayRCP<SC> valuesSolutionRep = valuesSolutionRepVel_->getBlock(j)->getDataNonConst(0);
			for(int i=0; i< nodeList.size() ; i++){
				deltaU[j] += deltaPhi[i]*valuesSolutionRep[nodeList[i]];	
			}
		}
	}
	for (UN w=0; w< QuadW.size(); w++){
			rhsFunc(&quadPointsTrans[w][0], &valueFunc[0] ,paras);
			for(int j=0 ; j< dofs_; j++){
		   		resElement[j] += QuadW[w] * std::pow(valueFunc[j] + deltaU[j] - nablaU[j][w] - nablaP[j] ,2); 
		}
	}
	double resElementValue =0.;
	for(int j=0; j< dofs_ ; j++)
		resElementValue += resElement[j] *detB1;

	return resElementValue;


}
/*!

\brief Returns neccesary quadrature Values. Is distinguishes between needing Element or Surface information.

@param[in] dim Dimension for which the quadrature points are needed.
@param[in] FEType Finite element type for which the quadrature points are needed.
@param[in] Type Type of quadrature points are need. Either 'Element' if you integrate over an element or 'Surface' if you need to integrate over a surface (i.e. for calculating the jump)
@param[in] QuadW Vector to be filled with the quadrature weights accordingly
@param[in] FiniteElement surface for which you need the quadrature points in case if 'Surface' type, as it is needed for figuring out the quadrature points

@param[out] QuadPts Quadrature points
@param[out] QuadW Quadrature weights
\brief Keep in mind that elementwise quadPoint are defined on reference element whereas surface quadPoints are defined on the input surface, which is typically not the reference Element. 

*/
template <class SC, class LO, class GO, class NO>
vec2D_dbl_Type ErrorEstimation<SC,LO,GO,NO>::getQuadValues(int dim, std::string FEType, std::string Type, vec_dbl_Type &QuadW, FiniteElement surface){

	vec2D_dbl_Type QuadPts(QuadW.size(), vec_dbl_Type(dim));
	vec2D_dbl_ptr_Type points = inputMesh_->pointsRep_;
	if(Type == "Element"){
		if(this->dim_ == 2){
	 		
		    double a = 1/6.;
		    QuadPts[0][0] 	= 0.5;
		    QuadPts[0][1]    = 0.5;

		    QuadPts[1][0] 	= 0.;
		    QuadPts[1][1] 	= 0.5;

		    QuadPts[2][0] 	= 0.5;
		    QuadPts[2][1] 	= 0.;

		    QuadW[0]		= a;
		    QuadW[1]        = a;
		    QuadW[2]        = a;
		}		
		else if(this->dim_ == 3){
	 		

			double a = .25;
			double b = 1./6.;
			double c = .5;
			QuadPts[0][0] = a;
			QuadPts[0][1] = a;
			QuadPts[0][2]= a;
			
			QuadPts[1][0] = b;
			QuadPts[1][1] = b;
			QuadPts[1][2] = b;
			
			QuadPts[2][0] = b;
			QuadPts[2][1] = b;
			QuadPts[2][2]=  c;
			
			QuadPts[3][0] = b;
			QuadPts[3][1] = c;
			QuadPts[3][2] = b;
			
			QuadPts[4][0] = c;
			QuadPts[4][1] = b;
			QuadPts[4][2] = b;
			
			QuadW[0] = -2./15.;
			QuadW[1] = 3./40.;
			QuadW[2] = 3./40.;
			QuadW[3] = 3./40.;
			QuadW[4] = 3./40.;
			}
	}
	else if (Type =="Surface"){
		if(dim==2){
			double x0 = points->at(surface.getNode(0)).at(0);
			double y0 = points->at(surface.getNode(0)).at(1);
			double x1 = points->at(surface.getNode(1)).at(0);
			double y1 = points->at(surface.getNode(1)).at(1);
			

			if(FEType == "P1"){
				
				QuadPts[0][0] =  (x0+x1)/2.;
				QuadPts[0][1] =  (y0+y1)/2.;

				QuadW[0] = 1.;
			}
			else if(FEType == "P2"){

				QuadPts[0][0] =  x0;
				QuadPts[0][1] =  y0;
				QuadPts[1][0] =  (x0+x1)/2.;
				QuadPts[1][1] =  (y0+y1)/2.;
				QuadPts[2][0] =   x1;
				QuadPts[2][1] =  y1;

				QuadW[0] = 1.;
				QuadW[1] = 4.;
				QuadW[2] = 1.;
			}
			
		}	
		else if(dim==3){
			// Here we choose as quadpoints the midpoints of the triangle sides
			double x0 = points->at(surface.getNode(0)).at(0);
			double y0 = points->at(surface.getNode(0)).at(1);
			double z0 = points->at(surface.getNode(0)).at(2);
			double x1 = points->at(surface.getNode(1)).at(0);
			double y1 = points->at(surface.getNode(1)).at(1);
			double z1 = points->at(surface.getNode(1)).at(2);
			double x2 = points->at(surface.getNode(2)).at(0);
			double y2 = points->at(surface.getNode(2)).at(1);
			double z2 = points->at(surface.getNode(2)).at(2);

			if(FEType == "P1"){
				// As nabla phi is a constant function, quad points don't really matter in that case 
				QuadPts[0][0] =   1/3.;
				QuadPts[0][1] =   1/3.;
				QuadPts[0][2] =   1/3.;

				QuadW[0] = 1.;
			}
			else if(FEType == "P2"){
				QuadPts[0][0] =  (x0+x1)/2.;
				QuadPts[0][1] =  (y0+y1)/2.;
				QuadPts[0][2] =  (z0+z1)/2.;
				QuadPts[1][0] =  (x0+x2)/2.;
				QuadPts[1][1] =  (y0+y2)/2.;
				QuadPts[1][2] =  (z0+z2)/2.;
				QuadPts[2][0] =  (x1+x2)/2.;
				QuadPts[2][1] =  (y1+y2)/2.;
				QuadPts[2][2] =  (z1+z2)/2.;

				QuadW[0] = 1/3.;
				QuadW[1] = 1/3.;
				QuadW[2] = 1/3.;
			}

		}
	}

	return QuadPts;	

}
/*!
 \brief function, that determines volume of tetrahedra.

@param[in] elements Elements.
@param[in] edgeElements Edges.
@param[in] points Points.

@param[out] volumeTetrahedra Volume of tetrahedras
*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::determineVolTet(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points){
	
	vec_dbl_Type volumeTetrahedra(elements->numberElements());

	vec_dbl_Type p1(3),p2(3), p3(3), v_K(3);

	vec2D_dbl_Type p(4,vec_dbl_Type(3));

	for(int k=0; k< elements->numberElements() ; k++){

		// Calculating edges of Tetraeder
		vec_LO_Type nodeList = 	elements->getElement(k).getVectorNodeListNonConst();
		
		for(int i=0; i<4; i++){
			p[i] = points->at(nodeList[i]);
		}

		p1[0] = p[0][0] - p[1][0];
		p1[1] = p[0][1] - p[1][1];
		p1[2] =	p[0][2] - p[1][2];

		p2[0] = p[0][0] - p[2][0];
		p2[1] = p[0][1] - p[2][1];
		p2[2] =	p[0][2] - p[2][2];

		p3[0] = p[0][0] - p[3][0];
		p3[1] = p[0][1] - p[3][1];
		p3[2] =	p[0][2] - p[3][2];

		
		v_K[0] = p1[1]*p2[2] - p1[2]*p2[1];
		v_K[1] = p1[2]*p2[0] - p1[0]*p2[2];
		v_K[2] = p1[0]*p2[1] - p1[1]*p2[0];

		
		volumeTetrahedra[k] = std::fabs(p3[0] * v_K[0] + p3[1] * v_K[1] +p3[2] * v_K[2]) / 6. ;
	}

		
		
	return volumeTetrahedra;
}

/*!
 \brief Calculating the circumdiameter of tetraeder.

@param[in] elements Elements.
@param[in] points Points.
@param[in] volTet Volume of tetrahedra.

@param[out] diamElements Uncircumdiameter of tetrahedra.

*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::calcDiamTetraeder(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points, vec_dbl_Type volTet){
	
	vec_dbl_Type diamElements(elements->numberElements());

	double a,b,c,A,B,C;

	vec2D_dbl_Type p(4,vec_dbl_Type(this->dim_));
	for(int k=0; k< elements->numberElements() ; k++){	
		
		// Calculating edges of Tetraeder
		vec_LO_Type nodeList = 	elements->getElement(k).getVectorNodeListNonConst();
		
		for(int i=0; i<4; i++){
			p[i] = points->at(nodeList[i]);
		}

		a = std::sqrt(std::pow(p[0][0] - p[1][0],2)+ std::pow(p[0][1] - p[1][1],2) +std::pow(p[0][2] - p[1][2],2));
		b = std::sqrt(std::pow(p[0][0] - p[2][0],2)+ std::pow(p[0][1] - p[2][1],2) +std::pow(p[0][2] - p[2][2],2));
		c = std::sqrt(std::pow(p[0][0] - p[3][0],2)+ std::pow(p[0][1] - p[3][1],2) +std::pow(p[0][2] - p[3][2],2));

		A = std::sqrt(std::pow(p[3][0] - p[2][0],2)+ std::pow(p[3][1] - p[2][1],2) +std::pow(p[3][2] - p[2][2],2));
		B = std::sqrt(std::pow(p[1][0] - p[3][0],2)+ std::pow(p[1][1] - p[3][1],2) +std::pow(p[1][2] - p[3][2],2));
		C = std::sqrt(std::pow(p[1][0] - p[2][0],2)+ std::pow(p[1][1] - p[2][1],2) +std::pow(p[1][2] - p[2][2],2));		


		diamElements[k] = std::sqrt(((a*A+b*B+c*C)*(-a*A+b*B+c*C)*(a*A-b*B+c*C)*(a*A+b*B-c*C))) / (12*volTet[k]) ;	

	}


	return diamElements;
}

/*!
 \brief Calcutlating the incircumdiameter of tetrahedra.

@param[in] elements Elements.
@param[in] surfaceTriangleElements TriangleElements.
@param[in] volTet Volume of tetrahedra.
@param[in] areaTriangles Area of faces of tetrahedra.

@param[out] rhoElements Incircumdiameter of tetrahedra.

*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::calcRhoTetraeder(ElementsPtr_Type elements,SurfaceElementsPtr_Type surfaceTriangleElements, vec_dbl_Type volTet, vec_dbl_Type areaTriangles){
	vec_dbl_Type rhoElements(elements->numberElements());
	
	
	vec_LO_Type surfaceOfEl(4);
	
	for(int k=0; k< elements->numberElements() ; k++){	
		// Calculating edges of Tetraeder
		surfaceOfEl = surfaceTriangleElements->getSurfacesOfElement(k);

		rhoElements[k] = (6*volTet[k]) / (areaTriangles[surfaceOfEl[0]] + areaTriangles[surfaceOfEl[1]] +areaTriangles[surfaceOfEl[2]] +areaTriangles[surfaceOfEl[3]]);

	}


	return rhoElements;
}
/*!
 \brief Calculating the diameter of triangles.

@param[in] elements Elements.
@param[in] points Points.

@param[out] diamElements Diameter of triangles.

*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::calcDiamTriangles3D(SurfaceElementsPtr_Type surfaceTriangleElements,vec2D_dbl_ptr_Type points,vec_dbl_Type& areaTriangles, vec_dbl_Type& rho_T, vec_dbl_Type& C_T){
	
	vec_dbl_Type diamElements(surfaceTriangleElements->numberElements());
	
	FiniteElement elementTmp;

	vec_dbl_Type vecTmp(3),vecTmp1(3),vecTmp2(3);
	for(int k=0; k< surfaceTriangleElements->numberElements() ; k++){
		double lengthA, lengthB, lengthC,s1;

		elementTmp = surfaceTriangleElements->getElement(k);

		vecTmp[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(1)).at(0);
		vecTmp[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(1)).at(1);
		vecTmp[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(1)).at(2);
		lengthA = std::sqrt(std::pow(vecTmp[0],2)+std::pow(vecTmp[1],2)+std::pow(vecTmp[2],2));
	
		vecTmp1[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp1[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		vecTmp1[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(2)).at(2);
		lengthB = std::sqrt(std::pow(vecTmp1[0],2)+std::pow(vecTmp1[1],2)+std::pow(vecTmp1[2],2));

		vecTmp2[0] = points->at(elementTmp.getNode(1)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp2[1] = points->at(elementTmp.getNode(1)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		vecTmp2[2] = points->at(elementTmp.getNode(1)).at(2) - points->at(elementTmp.getNode(2)).at(2);
		lengthC = std::sqrt(std::pow(vecTmp2[0],2)+std::pow(vecTmp2[1],2)+std::pow(vecTmp2[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		areaTriangles[k] = std::sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));	

		diamElements[k] = 2*(lengthA *lengthB *lengthC)/(4*areaTriangles[k]);

		rho_T[k] = 4*(areaTriangles[k]) / (lengthA+lengthB+lengthC);

		C_T[k] = diamElements[k] / rho_T[k];		
	}
	return diamElements;
}


/*!
 \brief Calculating the area of the triangle elements of tetrahedra

@param[in] elements Elements.
@param[in] edgeElements Edges.
@param[in] surfaceElements Triangle Elements.
@param[in] points Points.

@param[out] areaTriangles Area of triangles.
*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::determineAreaTriangles(ElementsPtr_Type elements,EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceElements, vec2D_dbl_ptr_Type points){
	
	vec_dbl_Type areaTriangles(surfaceElements->numberElements());
	vec_dbl_Type vecTmp(3),vecTmp1(3),vecTmp2(3);
	 
	FiniteElement surfaceTmp;
	for(int k=0; k< surfaceElements->numberElements() ; k++){

		double lengthA, lengthB, lengthC,s1;

		surfaceTmp= surfaceElements->getElement(k);
		vecTmp[0] = points->at(surfaceTmp.getNode(0)).at(0) - points->at(surfaceTmp.getNode(1)).at(0);
		vecTmp[1] = points->at(surfaceTmp.getNode(0)).at(1) - points->at(surfaceTmp.getNode(1)).at(1);
		vecTmp[2] = points->at(surfaceTmp.getNode(0)).at(2) - points->at(surfaceTmp.getNode(1)).at(2);
		lengthA = std::sqrt(std::pow(vecTmp[0],2)+std::pow(vecTmp[1],2)+std::pow(vecTmp[2],2));
	

		vecTmp1[0] = points->at(surfaceTmp.getNode(0)).at(0) - points->at(surfaceTmp.getNode(2)).at(0);
		vecTmp1[1] = points->at(surfaceTmp.getNode(0)).at(1) - points->at(surfaceTmp.getNode(2)).at(1);
		vecTmp1[2] = points->at(surfaceTmp.getNode(0)).at(2) - points->at(surfaceTmp.getNode(2)).at(2);
		lengthB = std::sqrt(std::pow(vecTmp1[0],2)+std::pow(vecTmp1[1],2)+std::pow(vecTmp1[2],2));

		vecTmp2[0] = points->at(surfaceTmp.getNode(1)).at(0) - points->at(surfaceTmp.getNode(2)).at(0);
		vecTmp2[1] = points->at(surfaceTmp.getNode(1)).at(1) - points->at(surfaceTmp.getNode(2)).at(1);
		vecTmp2[2] = points->at(surfaceTmp.getNode(1)).at(2) - points->at(surfaceTmp.getNode(2)).at(2);
		lengthC = std::sqrt(std::pow(vecTmp2[0],2)+std::pow(vecTmp2[1],2)+std::pow(vecTmp2[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		areaTriangles[k] = std::sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));	
	}
	return areaTriangles;
}
/*!
 \brief Calculating the diameter of triangles. 

@param[in] elements Elements
@param[in] points Points

@param[out] diamElements Diameter of triangles (Uncirclediameter).
@param[out] rho_T Incirclediameter.
@param[out] C_T Shape parameter.

*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>:: calcDiamTriangles(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points, vec_dbl_Type& areaTriangles, vec_dbl_Type& rho_T, vec_dbl_Type& C_T){
	
	vec_dbl_Type diamElements(elements->numberElements());
	
	 
	FiniteElement elementTmp;

	vec_dbl_Type vecTmp(2),vecTmp1(2),vecTmp2(2);
	for(int k=0; k< elements->numberElements() ; k++){


		double lengthA, lengthB, lengthC,s1;
		elementTmp = elements->getElement(k);
		vecTmp[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(1)).at(0);
		vecTmp[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(1)).at(1);
		lengthA = std::sqrt(std::pow(vecTmp[0],2)+std::pow(vecTmp[1],2));
	

		vecTmp1[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp1[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		lengthB = std::sqrt(std::pow(vecTmp1[0],2)+std::pow(vecTmp1[1],2));

		vecTmp2[0] = points->at(elementTmp.getNode(1)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp2[1] = points->at(elementTmp.getNode(1)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		lengthC = std::sqrt(std::pow(vecTmp2[0],2)+std::pow(vecTmp2[1],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		double area = std::sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));	

		areaTriangles[k] = area;

		diamElements[k] = 2*(lengthA *lengthB *lengthC)/(4*area);
		rho_T[k] = 4*(area) / (lengthA+lengthB+lengthC);

		C_T[k] = diamElements[k] / rho_T[k];
		
	}
	return diamElements;
}





/*!
 \brief DetermineCoarseningError is the essential part of the mesh coarsening process.
 \brief Instead of calulating an error for mesh level k, we redestribute it to lower mesh levels and refining those. We execute this function with an estimated error from level k. With this calculated error, we mark the elements according to that error and refine afterwards.
 If we decide to coarsen a certain mesh level, we take that level, look at the k-m level and refine that to the point where we are at the same level we wanted to perform the coarsening on.


@param[in] mesh_k Current mesh of level k.
@param[in] mesh_k_m Mesh of refinement level k-m.
@param[in] errorElementMv_k The error estimation of mesh level k.
@param[in] distribution Either 'forwards' or 'backwards'. We determine the error estimate in level k-m with redistributing backwards. if we are in level k-m we calculate the k-m+1 mesh level error estimation via redistributing the k-m error forward. 
@param[in] markingStrategy The strategy with which element are marked.
@param[in] theta Marking threshold.

*/
template <class SC, class LO, class GO, class NO>
typename ErrorEstimation<SC,LO,GO,NO>::MultiVectorPtr_Type ErrorEstimation<SC,LO,GO,NO>::determineCoarseningError(MeshUnstrPtr_Type mesh_k, MeshUnstrPtr_Type mesh_k_m, MultiVectorPtr_Type errorElementMv_k,  std::string distribution, std::string markingStrategy, double theta){

	theta_ =theta;
	markingStrategy_ = markingStrategy;

	// We determine which meshes we need to focus on.
	// Mesh of level k ist mesh_k and the mesh at the beginning of 'iteration'
	ElementsPtr_Type elements_k = mesh_k->getElementsC();
	// Mesh of level k-m is then mesh_k_m, an the mesh that helps us determine the new coarsening error
	ElementsPtr_Type elements_k_m = mesh_k_m->getElementsC();
	
	MultiVectorPtr_Type errorElementMv_k_m = Teuchos::rcp(new MultiVector_Type( mesh_k_m->getElementMap()) ); 
	Teuchos::ArrayRCP<SC> errorElement_k_m = errorElementMv_k_m->getDataNonConst(0);
	
	Teuchos::ArrayRCP<SC> errorElement_k = errorElementMv_k->getDataNonConst(0);

	// We determine the error of mesh level k-m with the error of mesh level k
	if(distribution == "backwards"){
		for(int i=0; i< errorElement_k.size(); i++){
			errorElement_k_m[elements_k->getElement(i).getPredecessorElement()] += errorElement_k[i];
		}
	}
	if(distribution == "forwards"){
		for(int i=0; i< errorElement_k_m.size(); i++){
			errorElement_k_m[i] = errorElement_k[elements_k_m->getElement(i).getPredecessorElement()];
		}
	}
	return errorElementMv_k_m;

	// As a result we coarsened the mesh at level k= iteration with the factor 'm'


}

/*!
 \brief UpdateElementsOfSurfaceLocalAndGlobal is performed here instead of in meshRefinement, as the information is only needed in case of error estimation. Equivalent function as updateElementsOfEdgeGlobal but with the important destinction, that it runs without communication and only relies on local information.

@param[in] edgeElements Edes.
@param[in] surfaceTriangleElements Triangle elements.

*/

template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::updateElementsOfSurfaceLocalAndGlobal(EdgeElementsPtr_Type edgeElements, SurfaceElementsPtr_Type surfaceTriangleElements){

	// The vector contains the Information which elements contain the edge
	vec2D_GO_Type elementsOfEdgeGlobal = edgeElements->getElementsOfEdgeGlobal();
	vec2D_LO_Type elementsOfEdgeLocal = edgeElements->getElementsOfEdgeLocal();

	vec2D_GO_Type elementsOfSurfaceGlobal = surfaceTriangleElements->getElementsOfSurfaceGlobal();

	GO nEl;
	vec_LO_Type edgeTmp1, edgeTmp2, edgeTmp3;
	LO id1, id2, id3;

	vec2D_LO_Type edgeList(0,vec_LO_Type(2));	

	for(int i=0; i<edgeElements->numberElements(); i++){
		edgeList.push_back({edgeElements->getElement(i).getNode(0),edgeElements->getElement(i).getNode(1)});
	}

	for(int i=0; i < surfaceTriangleElements->numberElements() ; i++){
		FiniteElement surfaceTmp = surfaceTriangleElements->getElement(i);
		if(surfaceTmp.isInterfaceElement() && (surfaceTmp.getFlag() == 0 || surfaceTmp.getFlag() == 10)){

			vec_int_Type tmpElements(0);
			//this->surfaceTriangleElements_->setElementsOfSurfaceLocalEntry(i,-1);

			edgeTmp1={surfaceTmp.getNode(0),surfaceTmp.getNode(1)};
			edgeTmp2={surfaceTmp.getNode(1),surfaceTmp.getNode(2)};
			//edgeTmp3={surfaceTmp.getNode(1),surfaceTmp.getNode(2)};

			sort(edgeTmp1.begin(),edgeTmp1.end());
			sort(edgeTmp2.begin(),edgeTmp2.end());
			//sort(edgeTmp3.begin(),edgeTmp3.end());
			
			auto it1 = find( edgeList.begin(), edgeList.end() ,edgeTmp1 );
            id1 = distance( edgeList.begin() , it1 );

			auto it2 = find(edgeList.begin(), edgeList.end() ,edgeTmp2 );
            id2 = distance(edgeList.begin() , it2 );

			for(int j=0 ; j< elementsOfEdgeGlobal[id1].size() ; j++)
				tmpElements.push_back( elementsOfEdgeGlobal[id1][j]);

			for(int j=0 ; j< elementsOfEdgeGlobal[id2].size() ; j++)
				tmpElements.push_back( elementsOfEdgeGlobal[id2][j]);


			sort(tmpElements.begin(),tmpElements.end());
			bool found =false;

			//cout << "Surface Element da " << elementsOfSurfaceGlobal[i][0] << " tmp elements in question "  << tmpElements[0] << " " ;
			for(int j=0; j< tmpElements.size()-1; j++){
				if((tmpElements[j] == tmpElements[j+1] )&& (tmpElements[j] != elementsOfSurfaceGlobal[i][0]) && (found==false)) { 
					nEl = tmpElements[j];
  					surfaceTriangleElements->setElementsOfSurfaceGlobalEntry(i,nEl);
					found = true;
				}
			}
			if(found == false){
				cout << " No Element Found for edges " << id1 << " " << id2 << " on Proc " << inputMesh_->getComm()->getRank() << " und element schon da=" << elementsOfSurfaceGlobal[i][0] << endl;
				cout << " Tmp1= ";
				for(int j=0; j< tmpElements.size(); j++){
					cout << " " << tmpElements[j];
				}
				cout<< " Flag Surface " << surfaceTmp.getFlag()   << endl;
			}			
		}
	}

	vec2D_LO_Type elementsOfSurfaceLocal = surfaceTriangleElements->getElementsOfSurfaceLocal();
	elementsOfSurfaceGlobal = surfaceTriangleElements->getElementsOfSurfaceGlobal();
}

/*!
 \brief Build Surface Map.
 Contrary to building the edge map, building the surface map is somewhat simpler as elementsOfSurfaceGlobal and elementsOfSurfaceLocal already exist.
 Via elementsOfSurface global each surface can be uniquely determined by the two elements it connects.

The surfacemap is only used for error estimation. Thus it is only build here and not in the refinementFactory.

*/

template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::buildTriangleMap(){

		int maxRank = std::get<1>(inputMesh_->rankRange_);
		const int myRank = inputMesh_->getComm()->getRank();

		vec_GO_Type globalProcs(0);
		for (int i=0; i<= maxRank; i++)
			globalProcs.push_back(i);

		Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

		vec_GO_Type localProc(0);
		localProc.push_back(inputMesh_->getComm()->getRank());
		Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

		MapPtr_Type mapGlobalProc =
			Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, inputMesh_->getComm()) );

		MapPtr_Type mapProc =
			Teuchos::rcp( new Map_Type(  Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, inputMesh_->getComm()) );



		vec2D_int_Type interfaceSurfacesLocalId(1,vec_int_Type(1));


		MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );

		// (A) First we determine a Map only for the interface Nodes
		// This will reduce the size of the Matrix we build later significantly if only look at the interface edges
		int numSurfaces= surfaceElements_->numberElements();
		vec2D_GO_Type inzidenzIndices(0,vec_GO_Type(2)); // Vector that stores global IDs of each edge (in Repeated Sense)
		vec_LO_Type localSurfaceIndex(0); // stores the local ID of surfaces in question 
		vec_GO_Type id(2);
		int surfacesUnique=0;

		vec2D_dbl_ptr_Type points = inputMesh_->getPointsRepeated();

		vec2D_GO_Type elementsOfSurfaceGlobal = surfaceElements_->getElementsOfSurfaceGlobal();

		vec_GO_Type elementRep(0);

		int interfaceNum=0;
		for(int i=0; i<numSurfaces; i++ ){
			if(surfaceElements_->getElement(i).isInterfaceElement()){

				id[0] = elementsOfSurfaceGlobal[i][0]; 
				id[1] = elementsOfSurfaceGlobal[i][1];
			 	


				sort(id.begin(),id.end());
				inzidenzIndices.push_back(id);

				localSurfaceIndex.push_back(i);
				interfaceNum++;

			}
	
			else{
				surfacesUnique++;
			}
			
			for(int j=0; j < elementsOfSurfaceGlobal[i].size() ; j++)
				elementRep.push_back(elementsOfSurfaceGlobal[i][j]);

		 }

		sort(elementRep.begin(),elementRep.end());
		vec_GO_Type::iterator ip = unique( elementRep.begin(), elementRep.end());
		elementRep.resize(distance(elementRep.begin(), ip)); 

		Teuchos::ArrayView<GO> elementRepArray = Teuchos::arrayViewFromVector( elementRep);

		MapPtr_Type elementMapRep =
			Teuchos::rcp( new Map_Type( Teuchos::OrdinalTraits<GO>::invalid(), elementRepArray, 0, inputMesh_->getComm()) );

	

		// This Matrix is row based, where the row is based on mapInterfaceNodesUnqiue
		// We then add a '1' Entry when two global Node IDs form an edge
		MatrixPtr_Type inzidenzMatrix = Teuchos::rcp( new Matrix_Type(inputMesh_->getElementMap(), 40 ) );
		Teuchos::Array<GO> index(1);
		Teuchos::Array<GO> col(1);
		Teuchos::Array<SC> value(1, Teuchos::ScalarTraits<SC>::one() );

		for(int i=0; i<inzidenzIndices.size(); i++ ){
			index[0] = inzidenzIndices[i][0];
			col[0] = inzidenzIndices[i][1];
			inzidenzMatrix->insertGlobalValues(index[0], col(), value());
		 }
   		inzidenzMatrix->fillComplete();

	
		// ---------------------------------------------------
		// 2. Set unique edges IDs ---------------------------
		// Setting the IDs of Edges that are uniquely on one
		// Processor
		// ---------------------------------------------------
		exportLocalEntry->putScalar( (LO) surfacesUnique );

		MultiVectorLOPtr_Type newSurfacesUniqueGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newSurfacesUniqueGlobal->putScalar( (LO) 0 ); 
		newSurfacesUniqueGlobal->importFromVector( exportLocalEntry, false, "Insert");

		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > newSurfacesList = newSurfacesUniqueGlobal->getData(0);

		GO procOffsetSurface=0;
		for(int i=0; i< myRank; i++)
			procOffsetSurface= procOffsetSurface + newSurfacesList[i];

		// global IDs for map
		vec_GO_Type vecGlobalIDsSurfaces(numSurfaces); 
	
		// Step 1: adding unique global edge IDs
		int count=0;
		for(int i=0; i< numSurfaces; i++){
			if(!surfaceElements_->getElement(i).isInterfaceElement()){
				vecGlobalIDsSurfaces.at(i) = procOffsetSurface+count;
				count++;
			}
		}	
		
		// Now we add the repeated ids, by first turning interfaceEdgesTag into a map
		// Offset for interface IDS:
		GO offsetInterface=0;
		for(int i=0; i< maxRank+1; i++)
			 offsetInterface=  offsetInterface + newSurfacesList[i];
		
		// (D) Now we count the row entries on each processor an set global IDs

		Teuchos::ArrayView<const LO> indices;
		Teuchos::ArrayView<const SC> values;
		vec2D_GO_Type inzidenzIndicesUnique(0,vec_GO_Type(2)); // Vector that stores only both global IDs if the first is part of my unique Interface Nodes
		MapConstPtr_Type colMap = inzidenzMatrix->getMap("col");
		MapConstPtr_Type rowMap = inzidenzMatrix->getMap("row");
		int numRows = rowMap->getNodeNumElements();
		int uniqueSurfaces =0;
		for(int i=0; i<numRows; i++ ){
			inzidenzMatrix->getLocalRowView(i, indices,values); 
			uniqueSurfaces = uniqueSurfaces+indices.size();
			vec_GO_Type surfaceTmp(2);
			for(int j=0; j<indices.size(); j++){
				surfaceTmp[0] = rowMap->getGlobalElement(i);
				surfaceTmp[1] = colMap->getGlobalElement(indices[j]);
				inzidenzIndicesUnique.push_back(surfaceTmp);
			}
		}
	
		exportLocalEntry->putScalar( uniqueSurfaces);
		MultiVectorLOPtr_Type newSurfaceInterfaceGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
		newSurfaceInterfaceGlobal->putScalar( (LO) 0 ); 
		newSurfaceInterfaceGlobal->importFromVector( exportLocalEntry,true, "Insert");

		// offset EdgesUnique for proc and globally
		Teuchos::ArrayRCP< const LO > numUniqueInterface = newSurfaceInterfaceGlobal->getData(0);

		procOffsetSurface=0;
		for(int i=0; i< myRank; i++)
			procOffsetSurface= procOffsetSurface + numUniqueInterface[i];

		int numInterfaceSurface=0;
		
		vec_GO_Type uniqueInterfaceIDsList_(inzidenzIndicesUnique.size());
		for(int i=0; i< uniqueInterfaceIDsList_.size(); i++)
			uniqueInterfaceIDsList_[i] = procOffsetSurface + i;

		MatrixPtr_Type indMatrix = Teuchos::rcp( new Matrix_Type(inputMesh_->getElementMap(), 40 ) );

		for(int i=0; i<inzidenzIndicesUnique.size(); i++ ){
			index[0] = inzidenzIndicesUnique[i][0];
			col[0] = inzidenzIndicesUnique[i][1];
			Teuchos::Array<SC> value2(1,uniqueInterfaceIDsList_[i]);
			indMatrix->insertGlobalValues(index[0], col(), value2());
		}
   		indMatrix->fillComplete();

		MatrixPtr_Type importMatrix = Teuchos::rcp( new Matrix_Type(elementMapRep, 40 ) );
   		
		importMatrix->importFromVector(indMatrix,false,"Insert");
		importMatrix->fillComplete(); 

		// Determine global indices
		GO surfaceID=0;
		colMap = importMatrix->getMap("col");
		rowMap = importMatrix->getMap("row");
		LO valueID=0;
		bool found = false;
		GO entry =0;
		for(int i=0; i<inzidenzIndices.size(); i++ ){
			
			importMatrix->getLocalRowView(rowMap->getLocalElement(inzidenzIndices[i][0]), indices,values); // Indices and values connected to node i / row i in Matrix
			// Entries in 'indices' represent the local entry in 'colmap
			// with 'getGlobalElement' we know the global Node ID that belongs to the first Node that form an edge
			// vector in with entries only for edges belonging to node i;
			vec2D_GO_Type indicesTmp(indices.size(),vec_GO_Type(2));
			vec_GO_Type indTmp(2);
			for(int j=0; j<indices.size(); j++){
				indTmp[0] = colMap->getGlobalElement(indices[j]);
				indTmp[1] = values[j];
				indicesTmp.push_back(indTmp);	// vector with the indices and values belonging to node i
			}

			found = false;
			for(int k=0; k<indicesTmp.size();k++){
				if(inzidenzIndices[i][1] == indicesTmp[k][0]){
					entry =k;
					k = indicesTmp.size();
					surfaceID = indicesTmp[entry][1];
					vecGlobalIDsSurfaces.at(localSurfaceIndex[i]) = offsetInterface + surfaceID;
					found =true;
				}
			}
			if(found == false)
				cout << " Asking for row " << rowMap->getLocalElement(inzidenzIndices[i][0]) << " for Edge [" << inzidenzIndices[i][0] << ",  " << inzidenzIndices[i][1] << "], on Proc " << myRank << " but no Value found " <<endl;
		 }


		Teuchos::RCP<std::vector<GO>> surfacesGlobMapping = Teuchos::rcp( new std::vector<GO>( vecGlobalIDsSurfaces ) );
		Teuchos::ArrayView<GO> surfacesGlobMappingArray = Teuchos::arrayViewFromVector( *surfacesGlobMapping);

		this->surfaceTriangleMap_.reset(new Map<LO,GO,NO>(Teuchos::OrdinalTraits<GO>::invalid(), surfacesGlobMappingArray, 0, inputMesh_->getComm()) );
}

/*!
 \brief Calcutlating the gradient of phi depending on quad points p.

@param[in] dim Dimension.
@param[in] intFE Integer value of discretisation P1 (1) or P2(2).
@param[in] p Quadpoints or other points for which phi is suppose to be evaluated.

@param[out] gradPhi Gradient of phi with evaluated values p.
*/

template <class SC, class LO, class GO, class NO>
vec2D_dbl_Type ErrorEstimation<SC,LO,GO,NO>::gradPhi(int dim,
                int intFE,
                vec_dbl_Type &p){

	int numNodes = dim+1;
	if(intFE == 2){
		numNodes=6;
		if(dim==3)
			numNodes=10;
	}
		
	vec2D_dbl_Type value(numNodes,vec_dbl_Type(dim));

    if (dim==2) {
        switch (intFE) {
            case 1://P1
                value[0][0]= -1.;
                value[0][1]= -1.;

                value[1][0]= 1.;
                value[1][1]= 0.;
               
                value[2][0]= 0.;
                value[2][1]= 1.;
                
                break;
            case 2://P2
                value[0][0]= 1. - 4.*(1 - p[0] - p[1]);
                value[0][1]= 1. - 4.*(1 - p[0] - p[1]);
          
                value[1][0]= 4.*p[0] - 1;
                value[1][1]= 0.;
              
                value[2][0]= 0.;
                value[2][1]= 4.*p[1] - 1;
              
                value[3][0]= 4 * (1. - 2*p[0] - p[1]);
                value[3][1]= -4 * p[0];
          
                value[4][0]= 4.*p[1];
                value[4][1]= 4.*p[0];
               
                value[5][0]= - 4.*p[1];
                value[5][1]= 4 * (1. - p[0] - 2*p[1]);
		
				break;
        }
    }
    else if(dim==3) {
        switch (intFE) {
            case 1://P1
               
                value[0][0]= -1.;
                value[0][1]= -1.;
                value[0][2]= -1.;
               
                value[1][0]= 1.;
                value[1][1]= 0.;
                value[1][2]= 0.;
               
                value[2][0]= 0.;
                value[2][1]= 1.;
                value[2][2]= 0.;
               
                value[3][0]= 0.;
                value[3][1]= 0.;
                value[3][2]= 1.;
              
      			break;

            case 2://P2
              
		        value[0][0]= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
		        value[0][1]= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
		        value[0][2]= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];

		        value[1][0]= 4.*p[0] - 1;
		        value[1][1]= 0.;
		        value[1][2]= 0.;

		        value[2][0]= 0.;
		        value[2][1]= 4.*p[1] - 1;
		        value[2][2]= 0.;

		        value[3][0]= 0.;
		        value[3][1]= 0.;
		        value[3][2]= 4.*p[2] - 1;

		        value[4][0]= 4. - 8.*p[0] - 4.*p[1] - 4.*p[2];
		        value[4][1]= - 4.*p[0];
		        value[4][2]= - 4.*p[0];

		        value[5][0]= 4.*p[1];
		        value[5][1]= 4.*p[0];
		        value[5][2]= 0.;
		        
		        value[6][0]= - 4.*p[1];
		        value[6][1]= 4. - 4.*p[0] - 8.*p[1] - 4.*p[2];
		        value[6][2]= - 4.*p[1];
		       
		        value[7][0]= - 4.*p[2];
		        value[7][1]= - 4.*p[2];
		        value[7][2]= 4. - 4.*p[0] - 4.*p[1] - 8.*p[2];

		        value[8][0]= 4.*p[2];
		        value[8][1]= 0.;
		        value[8][2]= 4.*p[0];
		       
		        value[9][0]= 0.;
		        value[9][1]= 4.*p[2];
		        value[9][2]= 4.*p[1];
			
			
			break;
			}
		}
	return value;

}

/*!
 \brief Calcutlating phi depending on quad points p.

@param[in] dim.
@param[in] intFE integer value of discretisation P1 (1) or P2(2).
@param[in] p Quadpoints or other points for which phi is to be evaluated.

@param[out] phi Phi with evaluated values p.
*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::phi(int dim,
                int intFE,
                vec_dbl_Type &p){

	int numNodes = dim+1;
	if(intFE == 2){
		numNodes=6;
		if(dim==3)
			numNodes=10;
	}
		
	vec_dbl_Type value(numNodes);

    if (dim==2) {
        switch (intFE) {
            case 1://P1
                value[0]= p[0];
                value[1]= p[1];
                value[2]= 1-p[0]-p[1];
                
                break;
          
            case 2://P2
               
                value[0]= -(1. - p[0]-p[1]) * (1 - 2.*(1-p[0] - p[1]));             
                value[1] = -p[0] *  (1 - 2*p[0]);               
                value[2] = -p[1] *  (1 - 2*p[1]);               
                value[3] = 4*p[0] * (1 - p[0]-p[1]);               
                value[4] = 4*p[0]*p[1];                
                value[5] = 4*p[1] * (1 - p[0]-p[1]);                       
                break;
        }
    }

    else if(dim==3){
        switch (intFE) {
            case 1://P1
                
                value[0] = (1. - p[0]-p[1]-p[2]);                       
                value[1] = p[0];                      
                value[2] = p[1];                       
                value[3] = p[2];
                
                break;
            case 2: //P2               
                value[0] = (1. - p[0]-p[1]-p[2]) * (1 - 2*p[0] - 2*p[1] - 2*p[2]);               
               	value[1] = p[0] * (2*p[0] - 1);               
                value[2] = p[1] * (2*p[1] - 1);                
                value[3] = p[2] * (2*p[2] - 1);               
                value[4] = 4*p[0] * (1 - p[0]-p[1]-p[2]);                
                value[5] = 4*p[0]*p[1];                
                value[6] = 4*p[1] * (1 - p[0]-p[1]-p[2]);                
                value[7] = 4*p[2] * (1 - p[0]-p[1]-p[2]);                
                value[8] = 4*p[0]*p[2];               
                value[9] = 4*p[1]*p[2];
                    break;
                }
		}
	return value;

}


}
#endif

