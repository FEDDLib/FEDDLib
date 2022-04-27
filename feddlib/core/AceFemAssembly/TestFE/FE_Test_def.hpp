#ifndef FE_TEST_DEF_hpp
#define FE_TEST_DEF_hpp
#include "FE_Test_decl.hpp"
namespace FEDD {
/*!
 \brief Constructor
*/
template <class SC, class LO, class GO, class NO> 
FE_Test<SC,LO,GO,NO>::FE_Test(bool saveAssembly):
domainVec_(0),
assemblyFEElements_(0),
saveAssembly_(saveAssembly)
{
}
/*!
 \brief Adding underlying domains.
@param[in] domain Finite element function space
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::addFE(DomainConstPtr_Type domain){
    
    if (saveAssembly_){
        DomainPtr_Type domainNC = Teuchos::rcp_const_cast<Domain_Type>( domain );
        domainNC->initializeFEData();
    }
    domainVec_.push_back(domain);
    
}
/*!
 \brief Assembly of constant stiffness matix for laplacian operator \f$ \Delta \f$
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::assemblyLaplace(int dim,
                                        string FEType,
                                        int degree,
                                        int dofs,
                                        MatrixPtr_Type &A,
                                        bool callFillComplete,
                                        int FELocExternal){
    ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemLaplace.xml");
    
    UN FEloc = checkFE(dim,FEType);
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec2D_dbl_Type nodes;
    int numNodes=dim+1;
    if(FEType == "P2"){
        numNodes= 6;
        if(dim==3)
            numNodes=10;
    }
    tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
    tuple_ssii_Type vel ("Laplace",FEType,dofs,numNodes); 
    problemDisk->push_back(vel);
    if(assemblyFEElements_.size()== 0)
        initAssembleFEElements("Laplace",problemDisk,elements, params,pointsRep);
    else if(assemblyFEElements_.size() != elements->numberElements())
         TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );
    for (UN T=0; T<elements->numberElements(); T++) {
        assemblyFEElements_[T]->assembleJacobian();
        SmallMatrixPtr_Type elementMatrix = assemblyFEElements_[T]->getJacobian(); 
        addFeMatrix(A,elementMatrix, elements->getElement(T), map,dofs);
        
    }
    if (callFillComplete)
        A->fillComplete();
}
/*!
 \brief Assembly of constant stiffness matix for linear elasticity 
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::assemblyLinElas(int dim,
                                        string FEType,
                                        int degree,
                                        int dofs,
                                        MatrixPtr_Type &A,
                                        bool callFillComplete,
                                        int FELocExternal){
    ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemLinElas.xml");
    
    UN FEloc = checkFE(dim,FEType);
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec2D_dbl_Type nodes;
    int numNodes=dim+1;
    if(FEType == "P2"){
        numNodes= 6;
        if(dim==3)
            numNodes=10;
    }
    tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
    tuple_ssii_Type vel ("LinElas",FEType,dofs,numNodes); 
    problemDisk->push_back(vel);
    if(assemblyFEElements_.size()== 0)
        initAssembleFEElements("LinearElasticity",problemDisk,elements, params,pointsRep);
    else if(assemblyFEElements_.size() != elements->numberElements())
         TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );
    vec_dbl_Type solution(dofs*numNodes,0);
    
    for (UN T=0; T<elements->numberElements(); T++) {
        assemblyFEElements_[T]->updateSolution(solution);
        assemblyFEElements_[T]->assembleJacobian();
        SmallMatrixPtr_Type elementMatrix = assemblyFEElements_[T]->getJacobian(); 
        addFeMatrix(A,elementMatrix, elements->getElement(T), map,dofs);
        
        // elementMatrix->print();
    }
    if (callFillComplete)
        A->fillComplete();
}
/*!
 \brief Assembly of constant stiffness matix for nonlinear elasticity 
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::assemblyNonLinElas(int dim,
                                    string FEType,
                                    int degree,
                                    int dofs,
                                    MultiVectorPtr_Type d_rep,
                                    MatrixPtr_Type &A,
                                    MultiVectorPtr_Type &resVec,
                                    ParameterListPtr_Type params,
                                    bool reAssemble,
                                    string assembleMode,
                                    bool callFillComplete,
                                    int FELocExternal){
                                    
    ElementsPtr_Type elements = domainVec_.at(0)->getElementsC();

	int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(0)->getPointsRepeated();

	MapConstPtr_Type map = domainVec_.at(0)->getMapRepeated();

	vec_dbl_Type solution(0);
	vec_dbl_Type solution_d;

	vec_dbl_Type rhsVec;

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numNodes=6;
	if(dim==3){
		numNodes=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type displacement ("Displacement",FEType,dofs,numNodes);
	problemDisk->push_back(displacement);

	if(assemblyFEElements_.size()== 0)
	 	initAssembleFEElements("NonLinearElasticity",problemDisk,elements, params,pointsRep);
	else if(assemblyFEElements_.size() != elements->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );

	MultiVectorPtr_Type resVecRep = Teuchos::rcp( new MultiVector_Type( domainVec_.at(0)->getMapVecFieldRepeated(), 1 ) );

 	SmallMatrixPtr_Type elementMatrix;
	for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);

		solution_d = getSolution(elements->getElement(T).getVectorNodeList(), d_rep,dofs);

		solution.insert( solution.end(), solution_d.begin(), solution_d.end() );

		assemblyFEElements_[T]->updateSolution(solution);

        if(assembleMode == "Jacobian"){
            assemblyFEElements_[T]->assembleJacobian();
            elementMatrix = assemblyFEElements_[T]->getJacobian();              
            assemblyFEElements_[T]->advanceNewtonStep();
            addFeMatrix(A, elementMatrix, elements->getElement(T), map, dofs);
        }
        if(assembleMode == "Rhs"){
		    assemblyFEElements_[T]->assembleRHS();
		    rhsVec = assemblyFEElements_[T]->getRHS(); 
			addFeMv(resVecRep, rhsVec, elements->getElement(T),  dofs);
		}


	}
	if (callFillComplete && assembleMode == "Jacobian")
	    A->fillComplete( domainVec_.at(0)->getMapVecFieldUnique(),domainVec_.at(0)->getMapVecFieldUnique());
	
    if(assembleMode == "Rhs"){
        resVec->importFromVector(resVecRep, true,"Add");

    }
        
}
/*!
 \brief Assembly of Jacobian for NavierStokes 
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::assemblyNavierStokes(int dim,
                                        string FETypeVelocity,
                                        string FETypePressure,
                                        int degree,
                                        int dofsVelocity,
                                        int dofsPressure,
                                        MultiVectorPtr_Type u_rep,
                                        BlockMatrixPtr_Type &A,
                                        bool callFillComplete,
                                        int FELocExternal){
    ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemNavierStokes.xml");
    
    UN FElocVel = checkFE(dim,FETypeVelocity); // Checks for different domains which belongs to a certain fetype
    UN FElocPres = checkFE(dim,FETypePressure); // Checks for different domains which belongs to a certain fetype
    ElementsPtr_Type elements = domainVec_.at(FElocVel)->getElementsC();
    int dofsElement = elements->getElement(0).getVectorNodeList().size();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FElocVel)->getPointsRepeated();
    MapConstPtr_Type mapVel = domainVec_.at(FElocVel)->getMapRepeated();
    MapConstPtr_Type mapPres = domainVec_.at(FElocPres)->getMapRepeated();
    vec_dbl_Type solution;
    /// Tupel construction follows follwing pattern:
    /// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
    int dofs;
    int numVelo=6;
    if(dim==3){
        numVelo=10;
    }
    tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
    tuple_ssii_Type vel ("Velocity","P2",dofsVelocity,numVelo);
    tuple_ssii_Type pres ("Pressure","P1",dofsPressure,dim+1);
    problemDisk->push_back(vel);
    problemDisk->push_back(pres);
    if(assemblyFEElements_.size()== 0)
        initAssembleFEElements("NavierStokes",problemDisk,elements, params,pointsRep);
    else if(assemblyFEElements_.size() != elements->numberElements())
         TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );
    //SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement));
    for (UN T=0; T<assemblyFEElements_.size(); T++) {
        
        solution = getSolution(elements->getElement(T).getVectorNodeList(), u_rep,dofsVelocity);
        assemblyFEElements_[T]->updateSolution(solution);
        assemblyFEElements_[T]->assembleJacobian();
        assemblyFEElements_[T]->advanceNewtonStep();
        SmallMatrixPtr_Type elementMatrix = assemblyFEElements_[T]->getJacobian(); 
        addFeBlockMatrix(A, elementMatrix, elements->getElement(T), mapVel, mapPres, problemDisk);
    }
    if (callFillComplete){
        A->getBlock(0,0)->fillComplete();
        A->getBlock(1,0)->fillComplete(domainVec_.at(FElocVel)->getMapVecFieldUnique(),domainVec_.at(FElocPres)->getMapUnique());
        A->getBlock(0,1)->fillComplete(domainVec_.at(FElocPres)->getMapUnique(),domainVec_.at(FElocVel)->getMapVecFieldUnique());
        A->getBlock(1,1)->fillComplete();
    }
}
/*!
 \brief Inserting element stiffness matrices into global stiffness matrix
@todo column indices pre determine
@param[in] &A Global Block Matrix
@param[in] elementMatrix Stiffness matrix of one element
@param[in] element Corresponding finite element
@param[in] map Map that corresponds to repeated nodes of first block
@param[in] map Map that corresponds to repeated nodes of second block
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type mapFirstRow,MapConstPtr_Type mapSecondRow, tuple_disk_vec_ptr_Type problemDisk){
        
        int numDisk = problemDisk->size();
        int dofs1 = std::get<2>(problemDisk->at(0));
        int dofs2 = std::get<2>(problemDisk->at(1));
        int numNodes1 = std::get<3>(problemDisk->at(0));
        int numNodes2=std::get<3>(problemDisk->at(1));
        int dofsBlock1 = dofs1*numNodes1;
        int dofsBlock2  = dofs2*numNodes2;
        Teuchos::Array<SC> value1( numNodes1, 0. );
        Teuchos::Array<GO> columnIndices1( numNodes1, 0 );
        for (UN i=0; i < numNodes1 ; i++) {
            for(int di=0; di<dofs1; di++){
                GO row =GO (dofs1* mapFirstRow->getGlobalElement( element.getNode(i) )+di);
                for(int d=0; d<dofs1; d++){
                    for (UN j=0; j < columnIndices1.size(); j++){
                        columnIndices1[j] = GO ( dofs1 * mapFirstRow->getGlobalElement( element.getNode(j) ) + d );
                        value1[j] = (*elementMatrix)[dofs1*i+di][dofs1*j+d];    
                    }
                    A->getBlock(0,0)->insertGlobalValues( row, columnIndices1(), value1() ); // Automatically adds entries if a value already exists 
                }          
            }
        }
        Teuchos::Array<SC> value2( 1, 0. );
        Teuchos::Array<GO> columnIndex( 1, 0. );
        Teuchos::Array<GO> rowIndex( 1, 0. );
        int offset= numNodes1*dofs1;
        //Teuchos::ArrayView<const LO> indices;
        //Teuchos::ArrayView<const SC> values;
        for (UN j=0; j < numNodes2; j++){
            rowIndex[0] = GO ( mapSecondRow->getGlobalElement( element.getNode(j) ) );
			for (UN i=0; i < numNodes1 ; i++) {
				for(int d=0; d<dofs1; d++){				
					value2[0] = (*elementMatrix)[i*dofs1+d][offset+j];			    				    		
					columnIndex[0] =GO (dofs1* mapFirstRow->getGlobalElement( element.getNode(i) )+d);

			  		A->getBlock(1,0)->insertGlobalValues( rowIndex[0], columnIndex(), value2() ); // Automatically adds entries if a value already exists   
			  		A->getBlock(0,1)->insertGlobalValues( columnIndex[0], rowIndex(), value2() ); // Automatically adds entries if a value already exists        
				}
			}      
		}
}
/*!
 \brief Initialization of vector consisting of the assembleFE Elements. Follows structure of 'normal' elements, i.e. elementMap also applicable
@param[in] elementType i.e. Laplace, Navier Stokes..
@param[in] problemDisk Tuple of specific problem Information
@param[in] elements
@param[in] params parameter list
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::initAssembleFEElements(string elementType,tuple_disk_vec_ptr_Type problemDisk,ElementsPtr_Type elements, ParameterListPtr_Type params,vec2D_dbl_ptr_Type pointsRep){
    
    vec2D_dbl_Type nodes;
    for (UN T=0; T<elements->numberElements(); T++) {
        
        nodes = getCoordinates(elements->getElement(T).getVectorNodeList(), pointsRep);
        AssembleFEFactory<SC,LO,GO,NO> assembleFEFactory;
        AssembleFEPtr_Type assemblyFE = assembleFEFactory.build(elementType,elements->getElement(T).getFlag(),nodes, params,problemDisk);
        assemblyFEElements_.push_back(assemblyFE);
    }
}
/*!
 \brief Assembly of Rhs.
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] &a Resulting vector
@param[in] fieldType Vectorvalued or scalar rhs
@param[in] func Rhs function
@param[in] funcParameter Parameter for rhs function 
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::assemblyRHS(int dim,
                       string FEType,
                       MultiVectorPtr_Type a,
                       string fieldType,
                       RhsFunc_Type func,
                      vector<SC>& funcParameter
                      )
 {
    ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblemLaplace.xml");
    UN FEloc = checkFE(dim,FEType);
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    Teuchos::ArrayRCP< SC > valuesRhs = a->getDataNonConst(0);
    int parameters;
    double x;
    vec2D_dbl_Type nodes;
    std::vector<double> valueFunc(dim);
    SC* paras = &(funcParameter[0]);
    
    func( &x, &valueFunc[0], paras );
    SC value;
    // Tupel
    tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
    tuple_ssii_Type vel ("RHS","P2",3,6);   // numnodes, geomertry Type string continuus lagrage (Ansatzraum typ, klassisch Lagrange), als statisches Objekt
    problemDisk->push_back(vel);
    if(assemblyFEElements_.size() != elements->numberElements())
         TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );
    for (UN T=0; T<assemblyFEElements_.size(); T++) {
        assemblyFEElements_[T]->addRHSFunc(func);
        assemblyFEElements_[T]->assembleRHS();
        vec_dbl_Type elementVec =  assemblyFEElements_[T]->getRHS();
        addFeMv(a, elementVec, elements->getElement(T),1); // if they are both multivectors its actually super simple! Import entries and add
    }
}
/*!
 \brief Inserting element stiffness matrices into global stiffness matrix
@todo column indices pre determine
@param[in] &A Global Matrix
@param[in] elementMatrix Stiffness matrix of one element
@param[in] element Corresponding finite element
@param[in] map Map that corresponds to repeated nodes
*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::addFeMatrix(MatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type map, int dofs){      
        int numNodes = elementMatrix->size()/dofs;
        Teuchos::Array<SC> value( numNodes, 0. );
        Teuchos::Array<GO> columnIndices( numNodes, 0 );
        //Teuchos::ArrayView<const LO> indices;
        //Teuchos::ArrayView<const SC> values;
        for (UN i=0; i < numNodes ; i++) {
            for(int d=0; d<dofs; d++){
                GO row =GO (dofs* map->getGlobalElement( element.getNode(i) )+d);
                for(int k=0;k<dofs;k++){
                    for (UN j=0; j < columnIndices.size(); j++){
                        columnIndices[j] = GO ( dofs * map->getGlobalElement( element.getNode(j) )+k);
                        value[j] = (*elementMatrix)[dofs*i+d][dofs*j+k];
                    }
                    A->insertGlobalValues( row, columnIndices(), value() ); // Automatically adds entries if a value already exists                                   
                }           
            }
        }
}


/*!
 \brief Inserting element rhs vector into global rhs vector
 
@param[in] res BlockMultiVector of residual vec; Repeated distribution; 2 blocks
@param[in] rhsVec sorted the same way as residual vec
@param[in] element of block1

*/

template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::addFeMv(MultiVectorPtr_Type &res, vec_dbl_Type rhsVec, FiniteElement elementBlock, int dofs){

    Teuchos::ArrayRCP<SC>  resArray = res->getDataNonConst(0);

	vec_LO_Type nodeList = elementBlock.getVectorNodeList();

	for(int i=0; i< nodeList.size() ; i++){
		for(int d=0; d<dofs; d++)
			resArray[nodeList[i]*dofs+d] += rhsVec[i*dofs+d];
	}

}
// ----------------------------------------------------------
// Helper Functions from FE Class 
/*! \brief Matches Finite Element Discretization to domain. Useful 
    for AssemblyFE class.
*/
template <class SC, class LO, class GO, class NO>
int FE_Test<SC,LO,GO,NO>::checkFE(int dim, string FEType){
    int FEloc;
    std::vector<int> matches;
    for (int i = 0; i < domainVec_.size(); i++) {
        if (domainVec_.at(i)->getDimension() == dim)
            matches.push_back(i);
    }
    bool found = false;
    for (int i = 0; i < matches.size();i++) {
        if (domainVec_.at( matches.at(i) )->getFEType() == FEType) {
            FEloc = matches.at(i);
            found = true;
        }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error   ,"Combination of dimenson(2/3) and FE Type(P1/P2) not defined yet. Use addFE(domain)");
    return FEloc;
}
/*!
 \brief Returns coordinates of local node ids
@param[in] localIDs
@param[in] points
@param[out] coordinates 
*/
template <class SC, class LO, class GO, class NO>
vec2D_dbl_Type FE_Test<SC,LO,GO,NO>::getCoordinates(vec_LO_Type localIDs, vec2D_dbl_ptr_Type points){
    vec2D_dbl_Type coordinates(0,vec_dbl_Type( points->at(0).size()));
    for(int i=0; i < localIDs.size() ; i++){
        coordinates.push_back(points->at(localIDs[i]));
    }
    return coordinates;
}
/*!
 \brief Returns entries of u of element
@param[in] localIDs
@param[in] points
@param[out] coordinates 
*/
template <class SC, class LO, class GO, class NO>
vec_dbl_Type FE_Test<SC,LO,GO,NO>::getSolution(vec_LO_Type localIDs, MultiVectorPtr_Type u_rep, int dofsVelocity){
    Teuchos::ArrayRCP<SC>  uArray = u_rep->getDataNonConst(0);
    
    vec_dbl_Type solution(0);
    for(int i=0; i < localIDs.size() ; i++){
        for(int d=0; d<dofsVelocity; d++){
            solution.push_back(uArray[localIDs[i]*dofsVelocity+d]);
        }
    }
    return solution;
}
};
#endif // FE_TEST_DEF_hpp
