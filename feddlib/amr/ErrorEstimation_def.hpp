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
 
 @brief  ErrorEstimation
 @version 1.0

 */

using namespace std;
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::outArg;

namespace FEDD {

template <class SC, class LO, class GO, class NO>
ErrorEstimation<SC,LO,GO,NO>:: ErrorEstimation(int dim, string problemType)
{
	this->dim_ = dim;
	this->problemType_ = problemType;	

}

template <class SC, class LO, class GO, class NO>
ErrorEstimation<SC,LO,GO,NO>::~ErrorEstimation(){
}

///
/// @brief We need to identify the problem with respect to the DegreesOfFreedom(dofs) of u or velocity and if necessary p or pressure
///
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

}
///
/// @brief We split the solution from the BlockMultiVector valuesSolution into one or two seperate blocks, where the blocks represen the different dimensions
///

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
\brief Main Function for A-posteriori Error Estimation

\brief depending on the problem the the error estimation is calculated accordingly

@param[in] inputMeshP1 the P1 Mesh that is used for later refinement
@param[in] inputMeshP12 the possible P2 Mesh, if one of the solutions is of P2 Discretisation, otherwise both meshes are P1
@param[in] solution of the PDE in BlockMultiVector Format (Block 0: Velocity, Block 1: Pressure)
@param[in] rhs Function 
@param[in] FETypeV as the maximum FEType for the Velocity, pressure is assumed to be P1 

*/ 
template <class SC, class LO, class GO, class NO>
typename ErrorEstimation<SC,LO,GO,NO>::MultiVectorPtr_Type ErrorEstimation<SC,LO,GO,NO>::estimateError(MeshUnstrPtr_Type inputMeshP12, MeshUnstrPtr_Type inputMeshP1, BlockMultiVectorConstPtr_Type valuesSolution, RhsFunc_Type rhsFunc, string FETypeV){
	
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
		vec_dbl_Type h_T = calcDiamElements(elements,points);

		// The divU Part and residual of the Element are calculated elementwise, as the are independet of other processors
		double divUElement=0;
		double resElement=0;
	
		for(int k=0; k< elements->numberElements();k++){
			edgeNumbers = edgeElements->getEdgesOfElement(k); // edges of Element k
			resElement = determineResElement(elements->getElement(k), rhsFunc); // calculating the residual of element k
			if(problemType_ == "Stokes") // If the Problem is a Stokes Problem we calculate divU (maybe start a hierarchy
				divUElement = determineDivU(elements->getElement(k));
			errorElement[k] = sqrt(1./2*(u_Jump[edgeNumbers[0]]+u_Jump[edgeNumbers[1]]+u_Jump[edgeNumbers[2]])+divUElement +h_T[k]*h_T[k]*resElement );

			if(maxErrorElLoc < errorElement[k] )
					maxErrorElLoc = errorElement[k];
			//cout << " Error Element [k] " << errorElement[k] << " with resElement " << resElement << " divU " << divUElement << endl;
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
		vec_dbl_Type areaTriangles = determineAreaTriangles(elements,edgeElements, surfaceElements,  points);
		vec_dbl_Type volTetraeder(elements->numberElements());
		vec_dbl_Type h_T_min = determineH_T_min(elements,edgeElements, points, volTetraeder);

		vec_dbl_Type h_E_min = vec_dbl_Type(surfaceElements->numberElements());
		vec_dbl_Type h_E(surfaceElements->numberElements());
		//MultiVectorPtr_Type h_E_minMv = Teuchos::rcp( new MultiVector_Type( elementMap, 1 ) );	
		//Teuchos::ArrayRCP< SC > valuesSolutionRep  = valuesSolutionRepeated->getDataNonConst(0);

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
			if(problemType_ == "Stokes") // If the Problem is a Stokes Problem we calculate divU (maybe start a hierarchy
				divUElement = determineDivU(elements->getElement(k));
		
			errorElement[k] = h_T_min[k]*sqrt((u_Jump[surfaceNumbers[0]]+u_Jump[surfaceNumbers[1]]+u_Jump[surfaceNumbers[2]]+u_Jump[surfaceNumbers[3]]+resElement+divUElement)); 

			if(maxErrorElLoc < errorElement[k] )
					maxErrorElLoc = errorElement[k];

			cout << " Error Element " << k << " " << errorElement[k] << " Jumps: " << u_Jump[surfaceNumbers[0]] << " " << u_Jump[surfaceNumbers[1]] << " " << u_Jump[surfaceNumbers[2]] << " " << u_Jump[surfaceNumbers[3]]  << " h_T " << h_T_min[k] << "divUElement " << divUElement << " res " << resElement << endl;
			
		}	
		
		double maxArea, minArea;
		double maxhTmin, minhTmin, maxhEmin, minhEmin, maxhE, minhE;

		auto it = max_element(areaTriangles.begin(), areaTriangles.end()); // c++11
		maxArea = areaTriangles[distance(areaTriangles.begin(), it)];

		it = min_element(areaTriangles.begin(), areaTriangles.end()); // c++11
		minArea = areaTriangles[distance(areaTriangles.begin(), it)];
		if(inputMesh_->getComm()->getRank() == 0){
			cout << "__________________________________________________________________________________________________________ " << endl;
			cout << " " << endl;
			cout << " Mesh Quality Assesment 3D " << endl;
			cout << " Area Triangles- Max: " << maxArea << " Min " << minArea  << endl;
			cout << " The maximal Error of Elements is "  << maxErrorElLoc << endl;
			cout << "__________________________________________________________________________________________________________ " << endl;
		}
	}
	errorEstimation_ = errorElementMv;
	//this->markElements(errorElement);

	//cout << " done " << endl;
	
	MESH_TIMER_STOP(errorEstimation);

	return errorElementMv;

}


/*!
@brief Tags only a certain Area for refinement and is independent of any error estimation
@param[in] inputMeshP1 the P1 Mesh that is used for later refinement
@param[in] the area, that is suppose to be refined. If is a vector defining the area as follows: row1:[x_0,x_1] x-limits, row2: [y_0,y_1] y-limits, row3: [z_0,z_1] z-limits 
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
		cout << "__________________________________________________________________________________________________________ " << endl;
		cout << " " << endl;
		cout << " The area you requested for Refinement is :" << endl ;
			cout << " x in [" << area[0][0] << ", " << area[0][1] << "] " << endl;
		if(this->dim_>1)
			cout << " y in [" << area[1][0] << ", " << area[1][1] << "] " << endl;
		if(this->dim_>2)
			cout << " z in [" << area[2][0] << ", " << area[2][1] << "] " << endl;
		cout << "__________________________________________________________________________________________________________ " << endl;
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
		cout << "__________________________________________________________________________________________________________ " << endl;
		cout << " " << endl;
		cout << " With the 'tagArea tool' " << taggedElements << " Elements were tagged for Refinement " << endl;
		cout << "__________________________________________________________________________________________________________ " << endl;
	}

}

/*!
@brief Part of the error estimator that calculates the jump part of the estimation

@param[in] none, as all necessary parameters for the calculation are already part of the Error estimation class.

@brief What kind of jump is calculated depends on the problemType we have at hand.
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

	// For each Edge or Triangle we can determine deriPhi with the necessary dim-1 Quad Points.
	//vec3D_dbl_Type deriPhi(numberSurfaceElements,vec2D_dbl_Type(numRowsDPhi, vec_dbl_Type(dim))) // calcDPih();

	vec3D_dbl_Type u_El = calcNPhi("Gradient",  dofs_, FEType2_);
	
	vec3D_dbl_Type p_El;
	if(calculatePressure_ == true)
		p_El = calcNPhi("None",  dofsP_, FEType1_);

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
						jumpQuad += pow(u_El[k][i][j],2);
					if(calculatePressure_ == true)
						jumpQuad += pow(u_El[k][i][j] - p_El[k][i][0],2);
				}
			}

			h_E =  sqrt(pow(p1_2[0],2)+pow(p1_2[1],2));

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

		// h_E,min defined as in "A posteriori error estimation for anisotropic tetrahedral and triangular finite element meshes"
		// This is the minimal per element, later we also need the adjacent elements' minimal height in order to determine h_E
		vec_dbl_Type areaTriangles = determineAreaTriangles(elements,edgeElements, surfaceTriangleElements,  points);
		vec_dbl_Type volTetraeder(elements->numberElements());
		vec_dbl_Type h_T_min = determineH_T_min(elements,edgeElements, points, volTetraeder);

		vec_dbl_Type h_E_min(surfaceTriangleElements->numberElements());
		vec_dbl_Type h_E(surfaceTriangleElements->numberElements());

		int elTmp1,elTmp2;

		vec2D_GO_Type elementsOfSurfaceGlobal = surfaceTriangleElements->getElementsOfSurfaceGlobal();
		vec2D_LO_Type elementsOfSurfaceLocal  = surfaceTriangleElements->getElementsOfSurfaceLocal();

		vec_GO_Type elementImportMap(0);
		
		int sizetmp=0;
		for(int i=0; i<surfaceTriangleElements->numberElements(); i++){
			sizetmp = elementsOfSurfaceLocal.at(i).size();
			if(sizetmp > 1){
				if(elementsOfSurfaceLocal.at(i).at(0) == -1){
					elementImportMap.push_back(elementsOfSurfaceGlobal.at(i).at(0));
				}
				if(elementsOfSurfaceLocal.at(i).at(1) == -1){
					elementImportMap.push_back(elementsOfSurfaceGlobal.at(i).at(1));
				}
			}
		}
		sort(elementImportMap.begin(),elementImportMap.end());
		vec_GO_Type::iterator ip = unique( elementImportMap.begin(), elementImportMap.end());
		elementImportMap.resize(distance(elementImportMap.begin(), ip)); 
		Teuchos::ArrayView<GO> globalElementArrayImp = Teuchos::arrayViewFromVector( elementImportMap);

		// global Ids of Elements' Nodes
		MapPtr_Type mapElementImport =
			Teuchos::rcp( new Map_Type( elementMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalElementArrayImp, 0, inputMesh_->getComm()) );

		MultiVectorPtr_Type h_T_min_MV = Teuchos::rcp( new MultiVector_Type(elementMap, 1 ) );	
		MultiVectorPtr_Type volTet_MV = Teuchos::rcp( new MultiVector_Type( elementMap, 1 ) );
		Teuchos::ArrayRCP< SC > entriesh_T_min  = h_T_min_MV->getDataNonConst(0);
		Teuchos::ArrayRCP< SC > entriesVolTet  = volTet_MV->getDataNonConst(0);

		for(int l=0; l<elements->numberElements(); l++){
			entriesh_T_min[l] = h_T_min[l];
			entriesVolTet[l] = volTetraeder[l];
		}

		MultiVectorPtr_Type h_T_min_MVImport = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );	
		MultiVectorPtr_Type volTet_MVImport = Teuchos::rcp( new MultiVector_Type( mapElementImport, 1 ) );
	
		h_T_min_MVImport->putScalar(0);
		volTet_MVImport->putScalar(0);
		
		h_T_min_MVImport->importFromVector(h_T_min_MV,true,"Insert");
		volTet_MVImport->importFromVector(volTet_MV,true,"Insert");

	
		Teuchos::ArrayRCP< SC > h_T_min_Entries  = h_T_min_MVImport->getDataNonConst(0);
		Teuchos::ArrayRCP< SC > volTet_Entries  = volTet_MVImport->getDataNonConst(0);

		// Then we determine the jump over the edges, if the element we need for this is not on our proc, we import the solution u
		for(int k=0; k< elements->numberElements();k++){
			surfaceElementsIds = surfaceTriangleElements->getSurfacesOfElement(k); // surfaces of Element k
			for(int i=0;i<4;i++){
				h_E_min[surfaceElementsIds[i]] = h_T_min[k];
				h_E[surfaceElementsIds[i]] = 3. * (volTetraeder[k] ) /areaTriangles[surfaceElementsIds[i]];
				if(elementsOfSurfaceGlobal.at(surfaceElementsIds[i]).size() > 1){  // not a boundary edge
					// Case that both necessary element are on the same Proc
					if(elementsOfSurfaceLocal.at(surfaceElementsIds[i]).at(0) != -1   &&  elementsOfSurfaceLocal.at(surfaceElementsIds[i]).at(1) != -1){
						elTmp1 = elementsOfSurfaceLocal.at(surfaceElementsIds[i]).at(0);
						elTmp2 = elementsOfSurfaceLocal.at(surfaceElementsIds[i]).at(1);

						h_E_min[surfaceElementsIds[i]] = 1./2* (h_T_min[elTmp1] + h_T_min[elTmp2]);
						h_E[surfaceElementsIds[i]] = 3./2. * (volTetraeder[elTmp1] + volTetraeder[elTmp2] ) /areaTriangles[surfaceElementsIds[i]];

					}
					// if one of the necessary elements is not on our proc, we need to import the corresponding value of u
					else{
						vec_GO_Type elementInq(2); // the element in question is the one missing on our proc at entry 0
						elTmp1 = elementsOfSurfaceGlobal.at(surfaceElementsIds[i]).at(0);
						elTmp2 = elementsOfSurfaceGlobal.at(surfaceElementsIds[i]).at(1);

						if(elementMap->getLocalElement(elTmp1) !=  -1){
								elementInq[0]=elTmp2;
								elementInq[1]=elTmp1;
			
							}
						else {
								elementInq[0]=elTmp1;
								elementInq[1]=elTmp2;
							}

						h_E_min[surfaceElementsIds[i]] = 1./2* (h_T_min[elementMap->getLocalElement(elementInq[1])] + h_T_min_Entries[mapElementImport->getLocalElement(elementInq[0])]);
						h_E[surfaceElementsIds[i]] = 3./2. * (volTetraeder[elementMap->getLocalElement(elementInq[1])] + volTet_Entries[mapElementImport->getLocalElement(elementInq[0])] ) /areaTriangles[surfaceElementsIds[i]];

					}
				}
			}
		}
		//MultiVectorPtr_Type h_E_minMv = Teuchos::rcp( new MultiVector_Type( elementMap, 1 ) );	
		//Teuchos::ArrayRCP< SC > valuesSolutionRep  = valuesSolutionRepeated->getDataNonConst(0);

		// necessary entities
		vec_dbl_Type p1(3),p2(3); // normal Vector of Surface
		vec_dbl_Type v_E(3); // normal Vector of Surface
		double norm_v_E;



		for(int k=0; k< surfaceTriangleElements->numberElements();k++){
			double jumpQuad=0.;
			for(int j=0; j<dofs_ ; j++){
				for(int i=0; i<u_El[k].size();i++){
					if(	calculatePressure_ == false)
						jumpQuad += pow(u_El[k][i][j],2);
					if(calculatePressure_ == true)
						jumpQuad += pow(u_El[k][i][j] - p_El[k][i][0],2);
				}
			}
			if(this->FEType2_ =="P2")
				quadWeightConst = areaTriangles[k];
			else 	
				quadWeightConst = areaTriangles[k];


			surfaceElementsError[k] = (1./h_E[k])*jumpQuad*quadWeightConst;

			cout << " Jump Quad " << " k " << k << " " << jumpQuad << " h_E " << h_E[k] << endl;
		}
						
	}
	return surfaceElementsError;

}

/*!
@brief Function that calculates the jump part for nabla u or p 

@param[in] phiDerivative is either 'Gradient' or 'None' and what kind of jump is calculated depends on the problemType we have at hand. If phiDerivative is 'Gradient' the nabla u jump part is caluculated and if its 'None' then the p-
@param[in] dofsSol is the degree of freedom of the caluclated jump part. p's dof is typically 1 whereas u's dof can vary depending on the problem
@param[in] FEType of the calculated jump part. 
 
*/
template <class SC, class LO, class GO, class NO>
vec3D_dbl_Type ErrorEstimation<SC,LO,GO,NO>::calcNPhi(string phiDerivative, int dofsSol, string FEType){

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
		if(elementsOfSurfaceGlobal.at(k).size() > 1){  
		// Per edge we have depending on discretization quadrature points and weights
			vec_dbl_Type quadWeights(quadPSize);
			vec2D_dbl_Type quadPoints; 
		
			if(dim_ == 2){ 
				quadPoints = getQuadValues(dim_, FEType2_ , "Surface", quadWeights, edgeElements->getElement(k));
				v_E.at(0) = points->at(edgeElements->getElement(k).getNode(0)).at(1) - points->at(edgeElements->getElement(k).getNode(1)).at(1);
				v_E.at(1) = -(points->at(edgeElements->getElement(k).getNode(0)).at(0) - points->at(edgeElements->getElement(k).getNode(1)).at(0));
				norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2));	
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

				norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2)+pow(v_E[2],2));	  
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
				for (int l=0; l< numNodes; l++){
					kn1[l]= elements->getElement(elementsIDs[i]).getNode(l);
				}
				
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
				// We make the destinction between a gradient jump calculation or a simple jump calculation 

				for(int l=0; l< quadPSize; l++){
					if(phiDerivative == "Gradient"){
						vec2D_dbl_Type deriPhi1 = gradPhi(dim_, intFE, quadPointsT1[l]);
						vec2D_dbl_Type deriPhiT1(numNodes,vec_dbl_Type(dofsSol));
						for(int q=0; q<dim_; q++){
							for(int p=0;p<numNodes; p++){
								for(int s=0; s< dim_ ; s++)
									deriPhiT1[p][q] += (deriPhi1[p][s]*Binv1[s][q]);
							}
						}
						vec2D_dbl_Type u1_Tmp(dim_, vec_dbl_Type(dofsSol));
						Teuchos::ArrayRCP< SC > valuesSolRep;	
						for( int d =0; d < dofsSol ; d++){
							valuesSolRep = valuesSolutionRepVel_->getBlock(d)->getDataNonConst(0);
							for(int s=0; s<dim_; s++){
								for(int t=0; t< numNodes ; t++){
					 				u1_Tmp[s][d] += valuesSolRep[kn1[t]]*deriPhiT1[t][s] ;
								}
							}
						}
						if(i==0){
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El1[k][l][d] += (u1_Tmp[s][d])*(v_E[s]/norm_v_E);
								}
								vec_El1[k][l][d] = sqrt(quadWeights[l])* vec_El1[k][l][d];

							}
						}
						else {
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El2[k][l][d] += (u1_Tmp[s][d])*(v_E[s]/norm_v_E);
								}
								vec_El2[k][l][d] = sqrt(quadWeights[l])* vec_El2[k][l][d];
							}
						}			

					}

					if(phiDerivative == "None"){
						vec_dbl_Type phiV = phi(dim_, intFE, quadPointsT1[l]);
						vec_dbl_Type vec_Tmp(dofsSol);
						Teuchos::ArrayRCP< SC > valuesSolRep;	
						for( int d =0; d < dofsSol ; d++){
							valuesSolRep = valuesSolutionRepPre_->getBlock(d)->getDataNonConst(0);
							for(int t=0; t< numNodes ; t++){
				 				vec_Tmp[d] += valuesSolRep[kn1[t]]*phiV[t] ;
							}
						}
						if(i==0){
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El1[k][l][d] += (vec_Tmp[d])*(v_E[s]/norm_v_E);
								}
								vec_El1[k][l][d] = sqrt(quadWeights[l])* vec_El1[k][l][d];

							}
						}
						else {
							for( int d =0; d < dofsSol ; d++){
								for(int s=0; s<dim_; s++){
									vec_El2[k][l][d] += (vec_Tmp[d])*(v_E[s]/norm_v_E);
								}
								vec_El2[k][l][d] = sqrt(quadWeights[l])* vec_El2[k][l][d];
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
	
	// We construct map from the previously extracted elementImportIDs
	sort(elementImportIDs.begin(),elementImportIDs.end());

	vec_GO_Type::iterator ip = unique( elementImportIDs.begin(), elementImportIDs.end());
	elementImportIDs.resize(distance(elementImportIDs.begin(), ip)); 
	Teuchos::ArrayView<GO> globalElementArrayImp = Teuchos::arrayViewFromVector( elementImportIDs);


	// Map which represents the surface ids that need to import element information
	MapPtr_Type mapElementImport =
		Teuchos::rcp( new Map_Type( elementMap->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalElementArrayImp, 0, inputMesh_->getComm()) );

	int maxRank = std::get<1>(inputMesh_->rankRange_);
	MapPtr_Type mapElementsUnique = mapElementImport;

	if(maxRank>0)
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
			entriesU_y_imp[j] = vec_El_Exp[surfaceMap->getLocalElement(elementImportIDs[j])][i][1] ;
			if(dofsSol == 3)
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
			vec_El2[surfaceMap->getLocalElement(elementImportIDs[j])][i][1] =entriesU_y_impF[j];
			if(dofsSol ==3)
				vec_El2[surfaceMap->getLocalElement(elementImportIDs[j])][i][2] = entriesU_z_impF[j];
		}	
	}
	
	for(int i=0; i< numberSurfaceElements ; i++){
		for(int j=0; j<quadPSize; j++){
			for(int k=0; k< dofsSol ; k++)
				vec_El[i][j][k] = fabs(vec_El1[i][j][k]) - fabs(vec_El2[i][j][k]); 

		}
	}


	return vec_El;

}

/*!
@brief Function that marks the elements for refinement 

@param[in] errorElementMv is the MultiVector that contains the estimated error for each element
@param[in] theta is a parameter determining a certain error bound for marking
@param[in] markingStrategy is the strategy with which the elements are marked. Implemented Strategies 'Doerfler' or 'Maximum'  
@param[in] meshP1 is the P1 mesh which is used for later refinement and has to be the one beeing marked

@brief !! it is essential that the meshP1 mesh inserted here is the mesh that will be used for mesh refinement, as it contains the elementwise-information determining refinement. !!
*/

template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::markElements(MultiVectorPtr_Type errorElementMv, double theta, string strategy,  MeshUnstrPtr_Type meshP1){


	Teuchos::ArrayRCP<SC> errorEstimation = errorElementMv->getDataNonConst(0);

	ElementsPtr_Type elements = meshP1->getElementsC();	

	this->markingStrategy_ = strategy;

	theta_ = theta;
	// As we decide which element to tag based on the maximum error in the elements globally, we need to communicated this maxErrorEl
	auto it = max_element(errorEstimation.begin(), errorEstimation.end()); // c++11
	double maxErrorEl= errorEstimation[distance(errorEstimation.begin(), it)];

	// Maximum-strategy for marking the elements as proposed by Verfürth in "A posterio error estimation"
    reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, maxErrorEl, outArg (maxErrorEl));
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
			thetaSumTmp = thetaSumTmp + pow(errorEstimation[k],2);
			flagged[k]=false;
		}
    	reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_SUM, thetaSumTmp, outArg (thetaSum));
		while(sigmaT < theta_*thetaSum){
			muMax=0.;
			for(int k=0; k< elements->numberElements(); k++){
				if(muMax < errorEstimation[k] && flagged[k]==false ){
					muMax = errorEstimation[k];
				}
			}
    		reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, muMax, outArg (muMax));
			for(int k=0; k< elements->numberElements(); k++){
				if(muMax == errorEstimation[k] && flagged[k]==false ){
					flagged[k]=true;
					sigmaT = sigmaT + pow(errorEstimation[k],2);
					}
			}
    	reduceAll<int, double> (*inputMesh_->getComm(), REDUCE_MAX, sigmaT, outArg (sigmaT));
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
    reduceAll<int, int> (*inputMesh_->getComm(), REDUCE_SUM, flagCount , outArg (flagCount));

	if(inputMesh_->getComm()->getRank() == 0){
	cout << "__________________________________________________________________________________________________________ " << endl;
	cout << " " << endl;
	cout << " The A-posteriori Error Estimation tagged " << flagCount << " Elements for adaptive Refinement with " << markingStrategy_ << "-Strategy and Theta= " << theta_ << endl;
	cout << "__________________________________________________________________________________________________________ " << endl;
	}
  
}

/*!
@brief Function that that determines || div(u) ||_T for a Element T

@param[in] FiniteElement element where ||div(u)||_T is calculated on

*/
template <class SC, class LO, class GO, class NO>
double ErrorEstimation<SC,LO,GO,NO>::determineDivU(FiniteElement element){
	//cout << "Calculating residual " << endl;
// Quad Points and weights of third order

	double resElement =0.;

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
	//cout << " DetB1 " << detB1 << endl;
	vec_dbl_Type valueFunc(dofs_);

	vec2D_dbl_Type quadPointsTrans(QuadW.size(),vec_dbl_Type(dim));

	for(int i=0; i< QuadW.size(); i++){
		 for(int j=0; j< dim; j++){
			for(int k=0; k< dim; k++){
		 		quadPointsTrans[i][j] += B1[j][k]* QuadPts[i].at(k) ; 
			} 
			quadPointsTrans[i][j] += points->at(element.getNode(0)).at(j)	;
		 }
	}

	
	
	int intFE = 1;
	if(this->FEType2_ == "P2")
		intFE =2;

	vec_dbl_Type divPhiV;
	for (UN w=0; w< QuadW.size(); w++){
		divPhiV = divPhi(dim, intFE, quadPointsTrans[w]);
		vec_dbl_Type uTmp(dofs_);
		for(int j=0; j< dofs_ ; j++){
			Teuchos::ArrayRCP<SC> valuesSolutionRep = valuesSolutionRepVel_->getBlock(j)->getDataNonConst(0);
			for(int i=0; i< divPhiV.size(); i++){
				uTmp[j] += valuesSolutionRep[nodeList[i]]*divPhiV[i];
			}
	   		resElement += pow((Binv1[0][0]+Binv1[1][1])*uTmp[j] ,2);
		}
		resElement *= QuadW[w];

	}
   	resElement =  resElement *detB1;

	//cout<< " ... done " << endl;

	return resElement;


}

/*!
@brief Function that that determines ||\Delta u_h + f ||_(L2(T)) or || \Delta u_h + f - \nabla p ||_T for a Element T

@param[in] FiniteElement element where ||div(u)||_T is calculated on
@param[in] RhsFunc_Type rhsFunc which is the function used for the rhs of the pde

*/

template <class SC, class LO, class GO, class NO>
double ErrorEstimation<SC,LO,GO,NO>::determineResElement(FiniteElement element, RhsFunc_Type rhsFunc){
	//cout << "Calculating residual " << endl;
// Quad Points and weights of third order

	double resElement =0.;

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
	//cout << " DetB1 " << detB1 << endl;
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
	if(calculatePressure_ == true){
		Teuchos::ArrayRCP<SC> valuesSolutionRepPre = valuesSolutionRepPre_->getBlock(0)->getDataNonConst(0);
		vec2D_dbl_Type deriPhi = gradPhi(dim, 1, quadPointsTrans[0]);
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
		   		resElement += QuadW[w] * pow(valueFunc[j] + deltaU[j] + nablaP[j] ,2);
		}
	}
   	resElement =  resElement *detB1;

	//cout<< " ... done " << endl;

	return resElement;


}
/*!

@brief Returns neccesary quadrature Values. Is distinguishes between needing Element or Surface information

@param[in] dim for which the quadrature points are needed
@param[in] FEType for which the quadrature points are needed
@param[in] Type of quadrature points are need. Either 'Element' if you integrate over an element or 'Surface' if you need to integrate over a surface (i.e. for calculating the jump)
@param[in] QuadW Vector to be filled with the quadrature weights accordingly
@param[in] FiniteElement surface for which you need the quadrature points in case if 'Surface' type, as it is needed for figuring out the quadrature points

*/
template <class SC, class LO, class GO, class NO>
vec2D_dbl_Type ErrorEstimation<SC,LO,GO,NO>::getQuadValues(int dim, string FEType, string Type, vec_dbl_Type &QuadW, FiniteElement surface){

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
				QuadPts[2][0] = 	x1;
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
				QuadPts[2][0] = 	(x1+x2)/2.;
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
@brief function, that determines h_T as the shortest vector inside a tetraeder as propose in...

@param[in] elements
@param[in] edgeelements
@param[in] points
@param[in] volTetraeder is calulated along the way and also usefull at another part of 3D jump calculation

*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::determineH_T_min(ElementsPtr_Type elements,EdgeElementsPtr_Type edgeElements, vec2D_dbl_ptr_Type points, vec_dbl_Type& volTetraeder){
	
	vec_dbl_Type h_T_min(elements->numberElements());
	vec_int_Type edgeNumbers(4);

	// We define P0P1 as the length the longest edge
	int P0,P1;
	int P0P1;
	// The triangle containing P0P1 with the largest area is P0P1P2
	int P2, P2Tmp1, P2Tmp2;

	double area1,area2,area;
	// P0P2 is the shortest edge of P0P1P2, the remaining node is P3
	int P3;
	// We define three vectors in our tetrahedra
	// rho1 = P0P1;
	double rho1 =0;
	// rho2 is the orthogonal Vector on P0P1 connecting it with P2
	double rho2,rho2Tmp1,rho2Tmp2;
	// rho3 is the vector to P3 that is perpendicular to P0P1P2 
	double rho3; 

	vec_dbl_Type vecTmp(3),vecTmp0(3),vecTmp1(3),vecTmp2(3),vecTmp11(3), vecTmp12(3), vecTmp21(3), vecTmp22(3),v_E(3);
	double norm_v_E;
	 
	FiniteElement edgeTmp;
	for(int k=0; k< elements->numberElements() ; k++){
		bool found=false;
		edgeNumbers = edgeElements->getEdgesOfElement(k); // edges of Element k
		rho1 =0;
		vec_dbl_Type length(6);
		P0= -1, P1 =-1, P2 =-1, P3 =-1;
		P2Tmp1=-1, P2Tmp2=-1;
		// First determine longest Edge
		for(int i=0; i<6; i++){
			edgeTmp= edgeElements->getElement(edgeNumbers[i]);
			vecTmp[0] = points->at(edgeTmp.getNode(0)).at(0) - points->at(edgeTmp.getNode(1)).at(0);
			vecTmp[1] = points->at(edgeTmp.getNode(0)).at(1) - points->at(edgeTmp.getNode(1)).at(1);
			vecTmp[2] = points->at(edgeTmp.getNode(0)).at(2) - points->at(edgeTmp.getNode(1)).at(2);
			length[i] = sqrt(pow(vecTmp[0],2)+pow(vecTmp[1],2)+pow(vecTmp[2],2));
			if(length[i] > rho1){
				rho1 = length[i];
				P0P1 = i;
				P0 = edgeTmp.getNode(0);
				P1 = edgeTmp.getNode(1);
				vecTmp0=vecTmp;

			}
		}
		// Calculating the areas
		for(int i=0; i<6; i++){
			if(i != P0P1){
				edgeTmp= edgeElements->getElement(edgeNumbers[i]);
				if(!found ){
					if(edgeTmp.getNode(0) == P0 ) {				
						P2Tmp1 = edgeTmp.getNode(1);
						found=true;
					}
					else if(edgeTmp.getNode(1) == P0){
						P2Tmp1 = edgeTmp.getNode(0);
						found=true;

					}
					else if(edgeTmp.getNode(0) == P1){
						P2Tmp1 = edgeTmp.getNode(1);
						found=true;
					}
					else if(edgeTmp.getNode(1) == P1){
						P2Tmp1 = edgeTmp.getNode(0);
						found=true;

					}
				}
				else if(found){
					if(edgeTmp.getNode(0) == P0 && edgeTmp.getNode(1) != P2Tmp1  ) {				
						P2Tmp2 = edgeTmp.getNode(1);
					}
					else if(edgeTmp.getNode(1) == P0 && edgeTmp.getNode(0) != P2Tmp1){
						P2Tmp2 = edgeTmp.getNode(0);
					}
					else if(edgeTmp.getNode(0) == P1 && edgeTmp.getNode(1) != P2Tmp1){
						P2Tmp2 = edgeTmp.getNode(1);
					}
					else if(edgeTmp.getNode(1) == P1 && edgeTmp.getNode(0) != P2Tmp1){
						P2Tmp2 = edgeTmp.getNode(0);
					}

				}
				
			}
		}
		double lengthA, lengthB, lengthC,s1,s2,s;
		lengthA = rho1;
		edgeTmp = edgeElements->getElement(edgeNumbers[P0P1]);

		vecTmp11[0] = points->at(P0).at(0) - points->at(P2Tmp1).at(0);
		vecTmp11[1] = points->at(P0).at(1) - points->at(P2Tmp1).at(1);
		vecTmp11[2] = points->at(P0).at(2) - points->at(P2Tmp1).at(2);
		lengthB = sqrt(pow(vecTmp11[0],2)+pow(vecTmp11[1],2)+pow(vecTmp11[2],2));

		vecTmp12[0] = points->at(P1).at(0) - points->at(P2Tmp1).at(0);
		vecTmp12[1] = points->at(P1).at(1) - points->at(P2Tmp1).at(1);
		vecTmp12[2] = points->at(P1).at(2) - points->at(P2Tmp1).at(2);
		lengthC = sqrt(pow(vecTmp12[0],2)+pow(vecTmp12[1],2)+pow(vecTmp12[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		area1 = sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));
		rho2Tmp1 = (2./lengthA)*area1;

		vecTmp21[0] = points->at(P0).at(0) - points->at(P2Tmp2).at(0);
		vecTmp21[1] = points->at(P0).at(1) - points->at(P2Tmp2).at(1);
		vecTmp21[2] = points->at(P0).at(2) - points->at(P2Tmp2).at(2);
		lengthB = sqrt(pow(vecTmp21[0],2)+pow(vecTmp21[1],2)+pow(vecTmp21[2],2));

		vecTmp22[0] = points->at(P1).at(0) - points->at(P2Tmp2).at(0);
		vecTmp22[1] = points->at(P1).at(1) - points->at(P2Tmp2).at(1);
		vecTmp22[2] = points->at(P1).at(2) - points->at(P2Tmp2).at(2);
		lengthC = sqrt(pow(vecTmp22[0],2)+pow(vecTmp22[1],2)+pow(vecTmp22[2],2));

		s2 = (lengthA+lengthB+lengthC)/2.;
		area2 = sqrt(s2*(s2-lengthA)*(s2-lengthB)*(s2-lengthC));
		rho2Tmp2 = (2./lengthA)*area2;

		if(area1 >= area2){
			P2 = P2Tmp1;
			P3 = P2Tmp2;
			rho2 = rho2Tmp1;
			vecTmp1=vecTmp11;
			vecTmp2=vecTmp12;
			area = area1;

		}

		else if(area1 < area2){
			P2 = P2Tmp2;
			P3 = P2Tmp1;
			rho2 = rho2Tmp2;
			vecTmp1=vecTmp21;
			vecTmp2=vecTmp22;
			area=area2;

		}

		// Now we have to determine the distance between P3 and P0P1P2:	
		v_E[0] = vecTmp1[1]*vecTmp0[2] - vecTmp1[2]*vecTmp0[1];
		v_E[1] = vecTmp1[2]*vecTmp0[0] - vecTmp1[0]*vecTmp0[2];
		v_E[2] = vecTmp1[0]*vecTmp0[1] - vecTmp1[1]*vecTmp0[0];

		norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2)+pow(v_E[2],2));	 
  
		rho3 = fabs((v_E[0] *(points->at(P3)[0]-points->at(P0)[0] ) + v_E[1] *(points->at(P3)[1] - points->at(P0)[1]) +v_E[2] *(points->at(P3)[2] -points->at(P0)[2]) ))/norm_v_E ;		

		volTetraeder[k] = 1./3. * (area * rho3);

		h_T_min[k] = min(rho1,rho2);
		h_T_min[k] = min(h_T_min[k],rho3);
		//cout << "h_e_min["<< k << "] " << h_E_min[k] << " aus rho1=" << rho1 << " rho2=" << rho2 << " rho3="<< rho3 << endl;
	
	}
	return h_T_min;
}
/*!
@brief Calculating the diameter of elements. This is necessary for 2D A-posteriori error estimation

@param[in] elements
@param[in] points

*/

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::calcDiamElements(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points){
	
	vec_dbl_Type diamElements(elements->numberElements());
	
	vec_dbl_Type areaTriangles(elements->numberElements());
	vec_dbl_Type P1(2),P2(2),P3(2);
	 
	FiniteElement elementTmp;
	vec_dbl_Type circPoints(2);
	for(int k=0; k< elements->numberElements() ; k++){
		elementTmp = elements->getElement(k);
		P1[0] = points->at(elementTmp.getNode(0)).at(0);
		P1[1] = points->at(elementTmp.getNode(0)).at(1);
		P2[0] = points->at(elementTmp.getNode(1)).at(0);
		P2[1] = points->at(elementTmp.getNode(1)).at(1);
		P3[0] = points->at(elementTmp.getNode(2)).at(0);
		P3[1] = points->at(elementTmp.getNode(2)).at(1);
		
		double d = 2*(P1[0]*(P2[1]-P3[1])+P2[0]*(P3[1]-P1[1])+P3[0]*(P1[1]-P2[1]));

		circPoints[0] = ((P1[0]*P1[0]+P1[1]*P1[1])*(P2[1]-P3[1])+(P2[0]*P2[0]+P2[1]*P2[1])*(P3[1]-P1[1])+(P3[0]*P3[0]+P3[1]*P3[1])*(P1[1]-P2[1]))/d;

		circPoints[1] = ((P1[0]*P1[0]+P1[1]*P1[1])*(P3[0]-P2[0])+(P2[0]*P2[0]+P2[1]*P2[1])*(P1[0]-P3[0])+(P3[0]*P3[0]+P3[1]*P3[1])*(P2[0]-P1[0]))/d;

		diamElements[k] = 2*sqrt(pow(P1[0]-circPoints[0],2)+pow(P1[1]-circPoints[1],2));

		
	}
	return diamElements;
}

/*!
@brief Calculating the area of the triangle elements of tetrahedra

@param[in] elements
@param[in] edgeElements
@param[in] surfaceElements
@param[in] points


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
		lengthA = sqrt(pow(vecTmp[0],2)+pow(vecTmp[1],2)+pow(vecTmp[2],2));
	

		vecTmp1[0] = points->at(surfaceTmp.getNode(0)).at(0) - points->at(surfaceTmp.getNode(2)).at(0);
		vecTmp1[1] = points->at(surfaceTmp.getNode(0)).at(1) - points->at(surfaceTmp.getNode(2)).at(1);
		vecTmp1[2] = points->at(surfaceTmp.getNode(0)).at(2) - points->at(surfaceTmp.getNode(2)).at(2);
		lengthB = sqrt(pow(vecTmp1[0],2)+pow(vecTmp1[1],2)+pow(vecTmp1[2],2));

		vecTmp2[0] = points->at(surfaceTmp.getNode(1)).at(0) - points->at(surfaceTmp.getNode(2)).at(0);
		vecTmp2[1] = points->at(surfaceTmp.getNode(1)).at(1) - points->at(surfaceTmp.getNode(2)).at(1);
		vecTmp2[2] = points->at(surfaceTmp.getNode(1)).at(2) - points->at(surfaceTmp.getNode(2)).at(2);
		lengthC = sqrt(pow(vecTmp2[0],2)+pow(vecTmp2[1],2)+pow(vecTmp2[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		areaTriangles[k] = sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));	
	}
	return areaTriangles;
}



/*!
@brief determineCoarseningError is the essential part of the mesh coarsening process.
@brief instead of calulating a error of mesh level k, we redestribute it to lower mesh levels and defining those.// We execute this function with an estimated error from the above 'estimateCoarseningError' function. With this error, we mark the elements according to that error and refine afterwards
 If we decide to coarsen a certain mesh level, we take that level, look at the k-m level and refine that to the point where we are at the same level we wanted to perform the coarsening on


@param[in] mesh_k the current mesh of level k
@param[in] mesh_k_m the mesh of refinement level k-m
@param[in] errorElementMv_k as the error estimation of mesh level k
@param[in] distribution is either 'forwards' or 'backwards'. We determine the error estimate in level k-m with redistributing backwards. if we are in level k-m we calculate the k-m+1 mesh level error estimation via redistributing the k-m error forward. 
@param[in] markingStrategy the strategy with which element are marked
@param[in] theta as the a marking threshold


*/
template <class SC, class LO, class GO, class NO>
typename ErrorEstimation<SC,LO,GO,NO>::MultiVectorPtr_Type ErrorEstimation<SC,LO,GO,NO>::determineCoarseningError(MeshUnstrPtr_Type mesh_k, MeshUnstrPtr_Type mesh_k_m, MultiVectorPtr_Type errorElementMv_k,  string distribution, string markingStrategy, double theta){

	MultiVectorPtr_Type errorElementMv = Teuchos::rcp(new MultiVector_Type( mesh_k->getElementMap()) ); 
	Teuchos::ArrayRCP<SC> errorEstimation = errorElementMv->getDataNonConst(0);

	// We determine which meshes we need to focus on.
	// Mesh of level k ist mesh_k and the mesh at the beginning of 'iteration'
	ElementsPtr_Type elements_k = Teuchos::rcp( new Elements(*(mesh_k->getElementsC())));

	theta_ =theta;
	markingStrategy_ = markingStrategy;

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

	this->markElements(errorElementMv_k_m,theta, markingStrategy, mesh_k_m);


	cout << " Finished Coarsening Error Estimation " << endl;

	return errorElementMv;

	// As a result we coarsened the mesh at level k= iteration with the factor 'm'


}

/*!
@brief updateElementsOfSurfaceLocalAndGlobal is performed here instead of in meshRefinement, as the information is only needed in case of error estimation

@param[in] edgeElements
@param[in] surfaceTriangleElements

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
		if(surfaceTmp.isInterfaceElement()){
			//cout << " Surface " << i << ": " <<  surfaceTmp.getNode(0) << " " << surfaceTmp.getNode(1) << " " << surfaceTmp.getNode(2) << endl;
			vec_int_Type tmpElements(0);
			//this->surfaceTriangleElements_->setElementsOfSurfaceLocalEntry(i,-1);

			edgeTmp1={surfaceTmp.getNode(0),surfaceTmp.getNode(1)};
			edgeTmp2={surfaceTmp.getNode(0),surfaceTmp.getNode(2)};
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
				//cout << tmpElements[j+1] << " " ;
				if((tmpElements[j] == tmpElements[j+1] )&& (tmpElements[j] != elementsOfSurfaceGlobal[i][0]) && (found==false)) { 
					//cout << " Update: tmp1=" << tmpElements[j] << " tmp2=" << tmpElements[j+1] << " und element schon da=" << elementsOfSurfaceGlobal[i][0] << endl;
					nEl = tmpElements[j];
  					surfaceTriangleElements->setElementsOfSurfaceGlobalEntry(i,nEl);
					found = true;
				}
			}
			//cout << endl;
			if(found == false)
				cout << " No Element Found for edges " << id1 << " " << id2 << " on Proc " << inputMesh_->getComm()->getRank() << endl;
		}
	}

	vec2D_LO_Type elementsOfSurfaceLocal = surfaceTriangleElements->getElementsOfSurfaceLocal();
	elementsOfSurfaceGlobal = surfaceTriangleElements->getElementsOfSurfaceGlobal();
}


// Build Surface Map
// Contrary to building the edge map, building the surface map is somewhat simpler as elementsOfSurfaceGlobal and elementsOfSurfaceLocal already exist.
// Via elementsOfSurface global each surface can be uniquely determined by the two elements it connects.
template <class SC, class LO, class GO, class NO>
void ErrorEstimation<SC,LO,GO,NO>::buildTriangleMap(){

		cout << " ---- Building Triangle Surface Map ----- " << endl;

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
			Teuchos::rcp( new Map_Type( inputMesh_->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, inputMesh_->getComm()) );

		MapPtr_Type mapProc =
			Teuchos::rcp( new Map_Type( inputMesh_->getEdgeMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, inputMesh_->getComm()) );



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
			Teuchos::rcp( new Map_Type( inputMesh_->getElementMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), elementRepArray, 0, inputMesh_->getComm()) );

	

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
   		inzidenzMatrix->fillComplete(); //mapInterfaceNodesUnique,mapInterfaceNodesUnique);

	
		// ---------------------------------------------------
		// 2 .Set unique edges IDs ---------------------------
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
   		indMatrix->fillComplete(); //mapUniqueInterfaceNodes,mapUniqueInterfaceNodes);
		//sindMatrix->print();
		//indMatrix->writeMM();

		MatrixPtr_Type importMatrix = Teuchos::rcp( new Matrix_Type(elementMapRep, 40 ) );
   		
		importMatrix->importFromVector(indMatrix,false,"Insert");
		importMatrix->fillComplete(); // mapInterfaceNodesRep, mapInterfaceNodesUnique); //mapScaledInterfaceNodesGlobalID,mapScaledInterfaceNodesGlobalID);
		
		//importMatrix->print(Teuchos::VERB_EXTREME);

		//importMatrix->writeMM();

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
			//sort(indicesTmp.begin(),indicesTmp.end());
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


		Teuchos::RCP<std::vector<GO>> surfacesGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDsSurfaces ) );
		Teuchos::ArrayView<GO> surfacesGlobMappingArray = Teuchos::arrayViewFromVector( *surfacesGlobMapping);

		this->surfaceTriangleMap_.reset(new Map<LO,GO,NO>(inputMesh_->getElementMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), surfacesGlobMappingArray, 0, inputMesh_->getComm()) );
		//this->surfaceTriangleMap_->print();

		cout << "--- done" << endl;
}
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
                value[0]=  p[0];
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

template <class SC, class LO, class GO, class NO>
vec_dbl_Type ErrorEstimation<SC,LO,GO,NO>::divPhi(int dim,
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
                value[0]=  1;
                value[1]= 1;

                value[2]= -2;
                
                break;
          
            case 2://P2
               
                value[0]= (1. - 4.*(1 - p[0] - p[1]))*2;             
                value[1] = 4.*p[0] - 1;               
                value[2] = 4.*p[1] - 1;               
                value[3] = 4 * (1. - 2*p[0] - p[1])-4 * p[0];               
                value[4] = 4*p[0]+4*p[1]; //done                
                value[5] = -4*p[1]+4-4*p[0]-8*p[1];   // done                    
                break;

              
                
 
          
             
        }
    }
	else if(dim ==3){
    	switch (intFE) {
            case 1://P1
               
                value[0]= -3.;            
                value[1]= 1.;               
                value[2]= 1.;
                value[3]= 1.;
              
      			break;

            case 2://P2
              
		        value[0]= -9. + 12.*p[0] + 12.*p[1] + 12.*p[2];
		        value[1]= 4.*p[0] - 1;
		        value[2]= 4.*p[1] - 1;
		        value[3] = 4.*p[2] - 1;
		        value[4]= 4. - 16.*p[0] - 4.*p[1] - 4.*p[2];
		        value[5]= 4.*p[1]+ 4.*p[0];
		        value[6]= 4. - 4.*p[0] - 16.*p[1] - 4.*p[2];
		        value[7]= 4. - 4.*p[0] - 4.*p[1] - 16.*p[2];
		        value[8]= 4.*p[2]+ 4.*p[0];		       
		        value[9]= 4.*p[2] + 4.*p[1];
			
				break;
			}

		}
			
	return value;

}




}
#endif

