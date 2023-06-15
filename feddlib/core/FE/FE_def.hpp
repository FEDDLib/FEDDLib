#ifndef FE_DEF_hpp
#define FE_DEF_hpp

#ifdef FEDD_HAVE_ACEGENINTERFACE
#include "aceinterface.hpp"
#endif

#include "FE_decl.hpp"

/*!
 Definition of FE

 @brief  FE
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


int MMAInitialisationCode[]={
    0,0
};


using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::outArg;

namespace FEDD {
DataElement::DataElement():
ht_(1,0.),
hp_(1,0.)
{
    
}

DataElement::DataElement(int size):
ht_(size,0.),
hp_(size,0.)
{
    
}

std::vector<double> DataElement::getHp()
{
    return hp_;
}

std::vector<double> DataElement::getHt()
{
    return ht_;
}

void DataElement::setHp( double* ht )
{
    for (int i=0; i<hp_.size(); i++)
        hp_[i] = ht[i];
}


template <class SC, class LO, class GO, class NO>
FE<SC,LO,GO,NO>::FE(bool saveAssembly):
domainVec_(0),
es_(),
setZeros_(false),
myeps_(),
ed_(0),
saveAssembly_(saveAssembly)
{
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::addFE(DomainConstPtr_Type domain){
    
    if (saveAssembly_){
        DomainPtr_Type domainNC = Teuchos::rcp_const_cast<Domain_Type>( domain );
        domainNC->initializeFEData();
    }
    domainVec_.push_back(domain);

}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::doSetZeros(double eps){

    setZeros_ = true;
    myeps_ = eps;

}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::applyBTinv( vec3D_dbl_ptr_Type& dPhiIn,
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

/*!

 \brief Assembly of Jacobian 
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 

*/

// Check the order of chemistry and solid in system matrix
/*template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::globalAssembly(string ProblemType,
                        int dim,                    
                        int degree,                       
                        MultiVectorPtr_Type sol_rep,
                        BlockMatrixPtr_Type &A,
                        BlockMultiVectorPtr_Type &resVec,
                        ParameterListPtr_Type params,
                        string assembleMode,
                        bool callFillComplete,
                        int FELocExternal){

    int problemSize = domainVec_.size(); // Problem size, as in number of subproblems

    // Depending on problem size we extract necessary information 

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numChem=3;
    if(FETypeChem == "P2"){
        numChem=6;
    }    
	if(dim==3){
		numChem=4;
        if(FETypeChem == "P2")
            numChem=10;
	}
    int numSolid=3;
    if(FETypeSolid == "P2")
        numSolid=6;
        
	if(dim==3){
		numSolid=4;
        if(FETypeSolid == "P2")
            numSolid=10;
	}

	BlockMultiVectorPtr_Type resVecRep = Teuchos::rcp( new BlockMultiVector_Type( 2) );

	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
    for(int i=0; i< problemSize; i++){
        int numNodes=0.;
        if(domainVec.at(i)->getFEType() == "P2"){
            numNodes=6;
        }    
        if(dim==3){
            numNodes=4;
            if(domainVec.at(i)->getFEType() == "P2")
                numNodes=10;
        }
        tuple_ssii_Type infoTuple (domainVec.at(i)->getPhysic(),domainVec.at(i)->getFEType(),domainVec.at(i)->getDofs(),numChem);
        problemDisk->push_back(infoTuple);

        MultiVectorPtr_Type resVec;
        if(domainVec.at(i)->getDofs() == 1)
            resVec = Teuchos::rcp( new MultiVector_Type( domainVec_.at(i)->getMapRepeated(), 1 ) );
        else
            resVec = Teuchos::rcp( new MultiVector_Type( domainVec_.at(i)->getMapVecFieldRepeated(), 1 ) );

        resVecRep->addBlock(resVec,i);


    }
	

	if(assemblyFEElements_.size()== 0){
       	initAssembleFEElements(ProblemType,problemDisk,domainVec_.at(0)->getElementsC(), params,domainVec_.at(0)->getPointsRepeated(),domainVec_.at(0)->getElementMap());
    }
	else if(assemblyFEElements_.size() != domainVec_.at(0)->getElementsC()->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );

    vec_dbl_Type solution_tmp;
    for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);
        vec_FE_Type elements;
        for(int i=0; i< problemSize; i++){
		    solution_tmp = getSolution(domainVec_.at(i)->getElementsC()->getElement(T).getVectorNodeList(), sol_rep->getBlock(i),domainVec.at(i)->getDofs());      
            solution.insert( solution.end(), solution_tmp.begin(), solution_tmp.end() );


        }   
        
  		assemblyFEElements_[T]->updateSolution(solution);

 		SmallMatrixPtr_Type elementMatrix;
	    vec_dbl_ptr_Type rhsVec;

     	if(assembleMode == "Jacobian"){
			assemblyFEElements_[T]->assembleJacobian();

            elementMatrix = assemblyFEElements_[T]->getJacobian(); 
            //elementMatrix->print();
			assemblyFEElements_[T]->advanceNewtonStep(); // n genereal non linear solver step
			
			addFeBlockMatrix(A, elementMatrix, domainVec_.at(0)->getElementsC()->getElement(T), problemDisk);
    
		}
		if(assembleMode == "Rhs"){
		    assemblyFEElements_[T]->assembleRHS();
		    rhsVec = assemblyFEElements_[T]->getRHS(); 
			addFeBlockMv(resVecRep, rhsVec, elementsSolid->getElement(T),elementsChem->getElement(T), dofsSolid,dofsChem);

		}
      		
	}
	if ( assembleMode == "Jacobian"){
        for(int probRow = 0; probRow <  problemSize; probRow++){
            for(int probCol = 0; probCol <  problemSize; probCol++){
                MapConstPtr_Type rowMap;
                MapConstPtr_Type domainMap;

                if(domainVec.at(probRow)->getDofs() == 1)
                    rowMap = domainVec_.at(FElocRow)->getMapUnique();
                else
                    rowMap = domainVec_.at(probRow)->getMapVecFieldUnique();

                if(domainVec.at(probCol)->getDofs() == 1)
                    domainMap = domainVec_.at(probCol)->getMapUnique()
                else
                    domainMap = domainVec_.at(probCol)->getMapVecFieldUnique()

                A->getBlock(probRow,probCol)->fillComplete(domainMap,rowMap);
            }
        }
	}
    else if(assembleMode == "Rhs"){

        for(int probRow = 0; probRow <  problemSize; probRow++){
            MapConstPtr_Type rowMap;
             if(domainVec.at(probRow)->getDofs() == 1)
                rowMap = domainVec_.at(FElocRow)->getMapUnique();
            else
                rowMap = domainVec_.at(probRow)->getMapVecFieldUnique();

		    MultiVectorPtr_Type resVecUnique = Teuchos::rcp( new MultiVector_Type( rowMap, 1 ) );

            resVecUnique->putScalar(0.);

            resVecUnique->exportFromVector( resVecRep, true, "Add" );

            resVec->addBlock(resVecUnique,probRow);
        }
	}

}*/

/*!

 \brief Inserting element stiffness matrices into global stiffness matrix


@param[in] &A Global Block Matrix
@param[in] elementMatrix Stiffness matrix of one element
@param[in] element Corresponding finite element
@param[in] map Map that corresponds to repeated nodes of first block
@param[in] map Map that corresponds to repeated nodes of second block

*/
/*template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement_vec element, tuple_disk_vec_ptr_Type problemDisk){
		
    int numDisk = problemDisk->size();

    for(int probRow = 0; probRow < numDisk; probRow++){
        for(int probCol = 0; probCol < numDisk; probCol++){
            int dofs1 = std::get<2>(problemDisk->at(probRow));
            int dofs2 = std::get<2>(problemDisk->at(probCol));

            int numNodes1 = std::get<3>(problemDisk->at(probRow));
            int numNodes2=std::get<3>(problemDisk->at(probCol));

            int offsetRow= numNodes1*dofs1*probRow;
            int offsetCol= numNodes1*dofs1*probCol;

            mapRow = domainVec.at(probRow)->getMapRepeated();
            mapCol = domainVec.at(probCol)->getMapRepeated();

            Teuchos::Array<SC> value2( numNodes2, 0. );
            Teuchos::Array<GO> columnIndices2( numNodes2, 0 );
            for (UN i=0; i < numNodes1; i++){
                for(int di=0; di<dofs1; di++){				
                    GO row =GO (dofs1* mapRow->getGlobalElement( element[probRow].getNode(i) )+di);
                    for(int d=0; d<dofs2; d++){				
                        for (UN j=0; j < numNodes2 ; j++) {
                            value2[j] = (*elementMatrix)[offsetRow+i*dofs1+di][offsetCol+j*dofs2+d];			    				    		
                            columnIndices2[j] =GO (dofs2* mapCol->getGlobalElement( element.getNode(j) )+d);
                        }
                        A->getBlock(probRow,probCol)->insertGlobalValues( row, columnIndices2(), value2() ); // Automatically adds entries if a value already exists                            
                    }
                }      
            }

        }
    }
}

*/

/*!

 \brief Inserting local rhsVec into global residual Mv;


@param[in] res BlockMultiVector of residual vec; Repeated distribution; 2 blocks
@param[in] rhsVec sorted the same way as residual vec
@param[in] element of block1

*/
/*
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec,tuple_disk_vec_ptr_Type problemDisk, FiniteElement_vec element)
    int numDisk = problemDisk->size();

    for(int probRow = 0; probRow < numDisk; probRow++){
        Teuchos::ArrayRCP<SC>  resArray_block = res->getBlockNonConst(probRow)->getDataNonConst(0);
      	vec_LO_Type nodeList_block = element[probRow].getVectorNodeList();
        int dofs = std::get<2>(problemDisk->at(probRow));
        int numNodes = std::get<3>(problemDisk->at(probRow));

        int offset = numNodes*dofs; 
        for(int i=0; i< nodeList_block.size() ; i++){
            for(int d=0; d<dofs; d++){
                resArray_block[nodeList_block[i]*dofs+d] += (*rhsVec)[i*dofs+d+offset];
            }
        }
    }
}
*/

/*!

 \brief Assembly of Jacobian for Linear Elasticity
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 

*/

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLinearElasticity(int dim,
	                                    string FEType,
	                                    int degree,
										int dofs,
										MultiVectorPtr_Type d_rep,
	                                    BlockMatrixPtr_Type &A,
										BlockMultiVectorPtr_Type &resVec,
 										ParameterListPtr_Type params,
 										bool reAssemble,
 										string assembleMode,
	                                    bool callFillComplete,
	                                    int FELocExternal){
	
	ElementsPtr_Type elements = domainVec_.at(0)->getElementsC();

	int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(0)->getPointsRepeated();

	MapConstPtr_Type mapVel = domainVec_.at(0)->getMapRepeated();

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
	 	initAssembleFEElements("LinearElasticity",problemDisk,elements, params,pointsRep,domainVec_.at(0)->getElementMap());
	else if(assemblyFEElements_.size() != elements->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );

	MultiVectorPtr_Type resVec_d = Teuchos::rcp( new MultiVector_Type( domainVec_.at(0)->getMapVecFieldRepeated(), 1 ) );
	
	BlockMultiVectorPtr_Type resVecRep = Teuchos::rcp( new BlockMultiVector_Type( 1) );
	resVecRep->addBlock(resVec_d,0);

 	SmallMatrixPtr_Type elementMatrix;
	for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);

		solution_d = getSolution(elements->getElement(T).getVectorNodeList(), d_rep,dofs);

		solution.insert( solution.end(), solution_d.begin(), solution_d.end() );

		assemblyFEElements_[T]->updateSolution(solution);
 
		assemblyFEElements_[T]->assembleJacobian();

		elementMatrix = assemblyFEElements_[T]->getJacobian(); 
			
		assemblyFEElements_[T]->advanceNewtonStep();


		addFeBlock(A, elementMatrix, elements->getElement(T), mapVel, 0, 0, problemDisk);
			
	}
	if (callFillComplete)
	    A->getBlock(0,0)->fillComplete( domainVec_.at(0)->getMapVecFieldUnique(),domainVec_.at(0)->getMapVecFieldUnique());
	



}

/*!

 \brief Assembly of Jacobian for nonlinear Elasticity
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 

*/

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyNonLinearElasticity(int dim,
	                                    string FEType,
	                                    int degree,
										int dofs,
										MultiVectorPtr_Type d_rep,
	                                    BlockMatrixPtr_Type &A,
										BlockMultiVectorPtr_Type &resVec,
 										ParameterListPtr_Type params, 									
	                                    bool callFillComplete,
	                                    int FELocExternal){
	
	ElementsPtr_Type elements = domainVec_.at(0)->getElementsC();

	int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(0)->getPointsRepeated();

	MapConstPtr_Type map = domainVec_.at(0)->getMapRepeated();

	vec_dbl_Type solution(0);
	vec_dbl_Type solution_d;

	vec_dbl_ptr_Type rhsVec;

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numNodes=6;
	if(dim==3){
		numNodes=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type displacement ("Displacement",FEType,dofs,numNodes);
	problemDisk->push_back(displacement);

    int neoHookeNum = params->sublist("Parameter").get("Neo-Hooke Modell",1);

    string nonLinElasModell = "NonLinearElasticity2";
    if(neoHookeNum == 1)
        nonLinElasModell = "NonLinearElasticity";

    //cout << " ######## Assembly Modell: " << nonLinElasModell << " ############ " <<  endl;


	if(assemblyFEElements_.size()== 0)
	 	initAssembleFEElements(nonLinElasModell,problemDisk,elements, params,pointsRep,domainVec_.at(0)->getElementMap());
	else if(assemblyFEElements_.size() != elements->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );
	

 	SmallMatrixPtr_Type elementMatrix;
	for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);

		solution_d = getSolution(elements->getElement(T).getVectorNodeList(), d_rep,dofs);

		solution.insert( solution.end(), solution_d.begin(), solution_d.end() );

		assemblyFEElements_[T]->updateSolution(solution);

        assemblyFEElements_[T]->assembleJacobian();
        elementMatrix = assemblyFEElements_[T]->getJacobian();              
        addFeBlock(A, elementMatrix, elements->getElement(T), map, 0, 0, problemDisk);

        assemblyFEElements_[T]->assembleRHS();
        rhsVec = assemblyFEElements_[T]->getRHS(); 
        addFeBlockMv(resVec, rhsVec, elements->getElement(T),  dofs);

        assemblyFEElements_[T]->advanceNewtonStep();


	}
	if (callFillComplete)
	    A->getBlock(0,0)->fillComplete( domainVec_.at(0)->getMapVecFieldUnique(),domainVec_.at(0)->getMapVecFieldUnique());
	
}

/*!

 \brief Assembly of Jacobian for nonlinear Elasticity
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] A Resulting matrix
@param[in] callFillComplete If Matrix A should be completely filled at end of function
@param[in] FELocExternal 

*/				

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyNonLinearElasticity(int dim,
	                                    string FEType,
	                                    int degree,
										int dofs,
										MultiVectorPtr_Type d_rep,
	                                    BlockMatrixPtr_Type &A,
										BlockMultiVectorPtr_Type &resVec,
 										ParameterListPtr_Type params, 									
                                        DomainConstPtr_Type domain,
                                        MultiVectorPtr_Type eModVec,
	                                    bool callFillComplete,
	                                    int FELocExternal){
	
	ElementsPtr_Type elements = domain->getElementsC();

	int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domain->getPointsRepeated();

	MapConstPtr_Type map = domain->getMapRepeated();

	vec_dbl_Type solution(0);
	vec_dbl_Type solution_d;

	vec_dbl_ptr_Type rhsVec;

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numNodes=6;
	if(dim==3){
		numNodes=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type displacement ("Displacement",FEType,dofs,numNodes);
	problemDisk->push_back(displacement);

    int neoHookeNum = params->sublist("Parameter").get("Neo-Hooke Modell",1);

    string nonLinElasModell = "NonLinearElasticity2";
    if(neoHookeNum == 1)
        nonLinElasModell = "NonLinearElasticity";

    //cout << " ######## Assembly Modell: " << nonLinElasModell << " ############ " <<  endl;

	if(assemblyFEElements_.size()== 0)
	 	initAssembleFEElements(nonLinElasModell,problemDisk,elements, params,pointsRep,domain->getElementMap());
	else if(assemblyFEElements_.size() != elements->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );
	
    Teuchos::ArrayRCP<SC>  eModVecA = eModVec->getDataNonConst(0);

 	SmallMatrixPtr_Type elementMatrix;
	for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);

		solution_d = getSolution(elements->getElement(T).getVectorNodeList(), d_rep,dofs);

		solution.insert( solution.end(), solution_d.begin(), solution_d.end() );

		assemblyFEElements_[T]->updateSolution(solution);
        assemblyFEElements_[T]->updateParameter("E",eModVecA[T]);
        assemblyFEElements_[T]->assembleJacobian();
        elementMatrix = assemblyFEElements_[T]->getJacobian();              
        addFeBlock(A, elementMatrix, elements->getElement(T), map, 0, 0, problemDisk);

        assemblyFEElements_[T]->assembleRHS();
        rhsVec = assemblyFEElements_[T]->getRHS(); 
        addFeBlockMv(resVec, rhsVec, elements->getElement(T),  dofs);

        assemblyFEElements_[T]->advanceNewtonStep();


	}
	if (callFillComplete)
	    A->getBlock(0,0)->fillComplete( domainVec_.at(0)->getMapVecFieldUnique(),domainVec_.at(0)->getMapVecFieldUnique());
	
}



/*!

 \brief Inserting local rhsVec into global residual Mv;


@param[in] res BlockMultiVector of residual vec; Repeated distribution; 2 blocks
@param[in] rhsVec sorted the same way as residual vec
@param[in] element of block1

*/

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock, int dofs){

    Teuchos::ArrayRCP<SC>  resArray_block = res->getBlockNonConst(0)->getDataNonConst(0);

	vec_LO_Type nodeList_block = elementBlock.getVectorNodeList();

	for(int i=0; i< nodeList_block.size() ; i++){
		for(int d=0; d<dofs; d++)
			resArray_block[nodeList_block[i]*dofs+d] += (*rhsVec)[i*dofs+d];
	}
}

// Check the order of chemistry and solid in system matrix
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyAceDeformDiffu(int dim,
                        string FETypeChem,
                        string FETypeSolid,
                        int degree,
                        int dofsChem,
                        int dofsSolid,
                        MultiVectorPtr_Type c_rep,
                        MultiVectorPtr_Type d_rep,
                        BlockMatrixPtr_Type &A,
                        BlockMultiVectorPtr_Type &resVec,
                        ParameterListPtr_Type params,
                        string assembleMode,
                        bool callFillComplete,
                        int FELocExternal){

    if((FETypeChem != "P2") || (FETypeSolid != "P2") || dim != 3)
    	TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "No AceGen Implementation available for Discretization and Dimension." );


    UN FElocChem = 1; //checkFE(dim,FETypeChem); // Checks for different domains which belongs to a certain fetype
    UN FElocSolid = 0; //checkFE(dim,FETypeSolid); // Checks for different domains which belongs to a certain fetype

	ElementsPtr_Type elementsChem= domainVec_.at(FElocChem)->getElementsC();

	ElementsPtr_Type elementsSolid = domainVec_.at(FElocSolid)->getElementsC();

    //this->domainVec_.at(FElocChem)->info();
    //this->domainVec_.at(FElocSolid)->info();
	//int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FElocSolid)->getPointsRepeated();

	MapConstPtr_Type mapChem = domainVec_.at(FElocChem)->getMapRepeated();

	MapConstPtr_Type mapSolid = domainVec_.at(FElocSolid)->getMapRepeated();

	vec_dbl_Type solution_c;
	vec_dbl_Type solution_d;

	vec_dbl_ptr_Type rhsVec;

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numChem=3;
    if(FETypeChem == "P2"){
        numChem=6;
    }    
	if(dim==3){
		numChem=4;
        if(FETypeChem == "P2")
            numChem=10;
	}
    int numSolid=3;
    if(FETypeSolid == "P2")
        numSolid=6;
        
	if(dim==3){
		numSolid=4;
        if(FETypeSolid == "P2")
            numSolid=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type chem ("Chemistry",FETypeChem,dofsChem,numChem);
	tuple_ssii_Type solid ("Solid",FETypeSolid,dofsSolid,numSolid);
	problemDisk->push_back(solid);
	problemDisk->push_back(chem);

	tuple_disk_vec_ptr_Type problemDiskChem = Teuchos::rcp(new tuple_disk_vec_Type(0));
    problemDiskChem->push_back(chem);

	string SCIModel = params->sublist("Parameter").get("Structure Model","SCI_simple");

	if(assemblyFEElements_.size()== 0){
       	initAssembleFEElements(SCIModel,problemDisk,elementsChem, params,pointsRep,domainVec_.at(FElocSolid)->getElementMap());
    }
	else if(assemblyFEElements_.size() != elementsChem->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );

	//SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement));

	MultiVectorPtr_Type resVec_c = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocChem)->getMapRepeated(), 1 ) );
    MultiVectorPtr_Type resVec_d = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocSolid)->getMapVecFieldRepeated(), 1 ) );
	
	BlockMultiVectorPtr_Type resVecRep = Teuchos::rcp( new BlockMultiVector_Type( 2) );
    resVecRep->addBlock(resVec_d,0);
    resVecRep->addBlock(resVec_c,1);
   
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);


    for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);

		solution_c = getSolution(elementsChem->getElement(T).getVectorNodeList(), c_rep,dofsChem);
		solution_d = getSolution(elementsSolid->getElement(T).getVectorNodeList(), d_rep,dofsSolid);

        // First Solid, then Chemistry
		solution.insert( solution.end(), solution_d.begin(), solution_d.end() );
		solution.insert( solution.end(), solution_c.begin(), solution_c.end() );
        
  		assemblyFEElements_[T]->updateSolution(solution);

 		SmallMatrixPtr_Type elementMatrix;

        // ------------------------
        /*buildTransformation(elementsSolid->getElement(T).getVectorNodeList(), pointsRep, B, FETypeSolid);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);
        cout << " Determinante " << detB << endl;*/
        // ------------------------




		if(assembleMode == "Jacobian"){
			assemblyFEElements_[T]->assembleJacobian();

            elementMatrix = assemblyFEElements_[T]->getJacobian(); 
           // elementMatrix->print();
			assemblyFEElements_[T]->advanceNewtonStep(); // n genereal non linear solver step
			
			addFeBlockMatrix(A, elementMatrix, elementsSolid->getElement(T), elementsSolid->getElement(T), mapSolid, mapChem, problemDisk);

          

		}
		if(assembleMode == "Rhs"){
		    assemblyFEElements_[T]->assembleRHS();
		    rhsVec = assemblyFEElements_[T]->getRHS(); 
			addFeBlockMv(resVecRep, rhsVec, elementsSolid->getElement(T),elementsChem->getElement(T), dofsSolid,dofsChem);

		}
        if(assembleMode=="MassMatrix"){
            assemblyFEElements_[T]->assembleJacobian();

            AssembleFE_SCI_SMC_Active_Growth_Reorientation_Ptr_Type elTmp = Teuchos::rcp_dynamic_cast<AssembleFE_SCI_SMC_Active_Growth_Reorientation_Type>(assemblyFEElements_[T] );
            elTmp->getMassMatrix(elementMatrix);
            //elementMatrix->print();
   			addFeBlock(A, elementMatrix, elementsChem->getElement(T), mapChem, 0, 0, problemDiskChem);


        }

        

			
	}
	if ( assembleMode == "Jacobian"){
		A->getBlock(0,0)->fillComplete();
	    A->getBlock(1,0)->fillComplete(domainVec_.at(FElocSolid)->getMapVecFieldUnique(),domainVec_.at(FElocChem)->getMapUnique());
	    A->getBlock(0,1)->fillComplete(domainVec_.at(FElocChem)->getMapUnique(),domainVec_.at(FElocSolid)->getMapVecFieldUnique());
	    A->getBlock(1,1)->fillComplete();
	}
    else if(assembleMode == "Rhs"){

		MultiVectorPtr_Type resVecUnique_d = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocSolid)->getMapVecFieldUnique(), 1 ) );
		MultiVectorPtr_Type resVecUnique_c = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocChem)->getMapUnique(), 1 ) );

		resVecUnique_d->putScalar(0.);
		resVecUnique_c->putScalar(0.);

		resVecUnique_d->exportFromVector( resVec_d, true, "Add" );
		resVecUnique_c->exportFromVector( resVec_c, true, "Add" );

		resVec->addBlock(resVecUnique_d,0);
		resVec->addBlock(resVecUnique_c,1);
	}
    else if(assembleMode == "MassMatrix"){
   		A->getBlock(0,0)->fillComplete();
    }

}

// Check the order of chemistry and solid in system matrix
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyAceDeformDiffuBlock(int dim,
                        string FETypeChem,
                        string FETypeSolid,
                        int degree,
                        int dofsChem,
                        int dofsSolid,
                        MultiVectorPtr_Type c_rep,
                        MultiVectorPtr_Type d_rep,
                        BlockMatrixPtr_Type &A,
                        int blockRow,
                        int blockCol,
                        BlockMultiVectorPtr_Type &resVec,
                        int block,
                        ParameterListPtr_Type params,
                        string assembleMode,
                        bool callFillComplete,
                        int FELocExternal){

    if((FETypeChem != "P2") || (FETypeSolid != "P2") || dim != 3)
    	TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "No AceGen Implementation available for Discretization and Dimension." );
    if((blockRow != blockCol) && blockRow != 0)
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Block assemblyDeformDiffu AceGEN Version only implemented for 0,0 block right now" );

    //UN FElocChem = 1; //checkFE(dim,FETypeChem); // Checks for different domains which belongs to a certain fetype
    UN FElocSolid = 0; //checkFE(dim,FETypeSolid); // Checks for different domains which belongs to a certain fetype

	ElementsPtr_Type elementsChem= domainVec_.at(FElocSolid)->getElementsC();

	ElementsPtr_Type elementsSolid = domainVec_.at(FElocSolid)->getElementsC();

    //this->domainVec_.at(FElocChem)->info();
    //this->domainVec_.at(FElocSolid)->info();
	//int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FElocSolid)->getPointsRepeated();

	//MapConstPtr_Type mapChem = domainVec_.at(FElocChem)->getMapRepeated();

	MapConstPtr_Type mapSolid = domainVec_.at(FElocSolid)->getMapRepeated();

	vec_dbl_Type solution_c;
	vec_dbl_Type solution_d;

	vec_dbl_ptr_Type rhsVec;

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numChem=3;
    if(FETypeChem == "P2"){
        numChem=6;
    }    
	if(dim==3){
		numChem=4;
        if(FETypeChem == "P2")
            numChem=10;
	}
    int numSolid=3;
    if(FETypeSolid == "P2")
        numSolid=6;
        
	if(dim==3){
		numSolid=4;
        if(FETypeSolid == "P2")
            numSolid=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type chem ("Chemistry",FETypeChem,dofsChem,numChem);
	tuple_ssii_Type solid ("Solid",FETypeSolid,dofsSolid,numSolid);
	problemDisk->push_back(solid);
	problemDisk->push_back(chem);
	
	string SCIModel = params->sublist("Parameter").get("Structure Model","SCI_simple");

	if(assemblyFEElements_.size()== 0){
       	initAssembleFEElements(SCIModel,problemDisk,elementsChem, params,pointsRep,domainVec_.at(FElocSolid)->getElementMap());
    }
	else if(assemblyFEElements_.size() != elementsChem->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );

	//SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement));

	//MultiVectorPtr_Type resVec_c = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocChem)->getMapRepeated(), 1 ) );
    MultiVectorPtr_Type resVec_d = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocSolid)->getMapVecFieldRepeated(), 1 ) );
	
	BlockMultiVectorPtr_Type resVecRep = Teuchos::rcp( new BlockMultiVector_Type( 1) );
    resVecRep->addBlock(resVec_d,0);
    //resVecRep->addBlock(resVec_c,1);

    for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);

		solution_c = getSolution(elementsChem->getElement(T).getVectorNodeList(), c_rep,dofsChem);
		solution_d = getSolution(elementsSolid->getElement(T).getVectorNodeList(), d_rep,dofsSolid);

        // First Solid, then Chemistry
		solution.insert( solution.end(), solution_d.begin(), solution_d.end() );
		solution.insert( solution.end(), solution_c.begin(), solution_c.end() );
        
  		assemblyFEElements_[T]->updateSolution(solution);

 		SmallMatrixPtr_Type elementMatrix;

		if(assembleMode == "Jacobian"){
			assemblyFEElements_[T]->assembleJacobian();

            elementMatrix = assemblyFEElements_[T]->getJacobian(); 
           // elementMatrix->print();
			assemblyFEElements_[T]->advanceNewtonStep(); // n genereal non linear solver step
			addFeBlock(A, elementMatrix, elementsSolid->getElement(T), mapSolid, 0, 0, problemDisk);
			//addFeBlockMatrix(A, elementMatrix, elementsSolid->getElement(T),  mapSolid, mapChem, problemDisk);
		}
		if(assembleMode == "Rhs"){
		    assemblyFEElements_[T]->assembleRHS();
		    rhsVec = assemblyFEElements_[T]->getRHS(); 
			addFeBlockMv(resVecRep, rhsVec, elementsSolid->getElement(T), dofsSolid);
		}
        //if(assembleMode=="compute")
        //    assemblyFEElements_[T]->compute();

        

			
	}
	if ( assembleMode != "Rhs"){
		A->getBlock(0,0)->fillComplete();
	    //A->getBlock(1,0)->fillComplete(domainVec_.at(FElocSolid)->getMapVecFieldUnique(),domainVec_.at(FElocChem)->getMapUnique());
	    //A->getBlock(0,1)->fillComplete(domainVec_.at(FElocChem)->getMapUnique(),domainVec_.at(FElocSolid)->getMapVecFieldUnique());
	    //A->getBlock(1,1)->fillComplete();
	}

    if(assembleMode == "Rhs"){

		MultiVectorPtr_Type resVecUnique_d = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocSolid)->getMapVecFieldUnique(), 1 ) );
		//MultiVectorPtr_Type resVecUnique_c = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocChem)->getMapUnique(), 1 ) );

		resVecUnique_d->putScalar(0.);
		//resVecUnique_c->putScalar(0.);

		resVecUnique_d->exportFromVector( resVec_d, true, "Add" );
		//resVecUnique_c->exportFromVector( resVec_c, true, "Add" );

		resVec->addBlock(resVecUnique_d,0);
		//resVec->addBlock(resVecUnique_c,1);
	}


}
/*!

 \brief Inserting element stiffness matrices into global stiffness matrix


@param[in] &A Global Block Matrix
@param[in] elementMatrix Stiffness matrix of one element
@param[in] element Corresponding finite element
@param[in] map Map that corresponds to repeated nodes of first block
@param[in] map Map that corresponds to repeated nodes of second block

*/
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element1, FiniteElement element2, MapConstPtr_Type mapFirstRow,MapConstPtr_Type mapSecondRow, tuple_disk_vec_ptr_Type problemDisk){
		
		int numDisk = problemDisk->size();

		int dofs1 = std::get<2>(problemDisk->at(0));
		int dofs2 = std::get<2>(problemDisk->at(1));

		int numNodes1 = std::get<3>(problemDisk->at(0));
		int numNodes2=std::get<3>(problemDisk->at(1));

		int dofsBlock1 = dofs1*numNodes1;
		int dofsBlock2  = dofs2*numNodes2;

		Teuchos::Array<SC> value1( numNodes1, 0. );
        Teuchos::Array<GO> columnIndices1( numNodes1, 0 );

        Teuchos::Array<SC> value2( numNodes2, 0. );
        Teuchos::Array<GO> columnIndices2( numNodes2, 0 );

		for (UN i=0; i < numNodes1 ; i++) {
			for(int di=0; di<dofs1; di++){
				GO row =GO (dofs1* mapFirstRow->getGlobalElement( element1.getNode(i) )+di);
				for(int d=0; d<dofs1; d++){
					for (UN j=0; j < columnIndices1.size(); j++){
		                columnIndices1[j] = GO ( dofs1 * mapFirstRow->getGlobalElement( element1.getNode(j) ) + d );
						value1[j] = (*elementMatrix)[dofs1*i+di][dofs1*j+d];	
					}
			  		A->getBlock(0,0)->insertGlobalValues( row, columnIndices1(), value1() ); // Automatically adds entries if a value already exists 
				}          
            }
		}
		int offset= numNodes1*dofs1;


        Teuchos::Array<SC> value( 1, 0. );
        Teuchos::Array<GO> columnIndex( 1, 0 );
        for (UN i=0; i < numNodes2 ; i++) {
            for(int di=0; di<dofs2; di++){
                GO row =GO (dofs2* mapSecondRow->getGlobalElement( element2.getNode(i) )+di);
                for(int d=0; d<dofs2; d++){
                    for (UN j=0; j < columnIndices2.size(); j++){
                        double tmpValue =  (*elementMatrix)[offset+dofs2*i+di][offset+dofs2*j+d];
                        if(std::fabs(tmpValue) > 1.e-13){
                            columnIndex[0] = GO ( dofs2 * mapSecondRow->getGlobalElement( element2.getNode(j) ) + d );
                            value[0] = tmpValue;
                            A->getBlock(1,1)->insertGlobalValues( row, columnIndex(), value() ); // Automatically adds entries if a value already exists 

                        }
                    }
                }          
            }
        }
        
		for (UN i=0; i < numNodes1; i++){
			for(int di=0; di<dofs1; di++){				
                GO row =GO (dofs1* mapFirstRow->getGlobalElement( element1.getNode(i) )+di);
                for(int d=0; d<dofs2; d++){				
                    for (UN j=0; j < numNodes2 ; j++) {
                        value2[j] = (*elementMatrix)[i*dofs1+di][offset+j*dofs2+d];			    				    		
                        columnIndices2[j] =GO (dofs2* mapSecondRow->getGlobalElement( element2.getNode(j) )+d);
                    }
                    A->getBlock(0,1)->insertGlobalValues( row, columnIndices2(), value2() ); // Automatically adds entries if a value already exists                            
                }
            }      
		}

        for (UN j=0; j < numNodes2; j++){
            for(int di=0; di<dofs2; di++){	
                GO row = GO (dofs2* mapSecondRow->getGlobalElement( element2.getNode(j) ) +di );
                for(int d=0; d<dofs1; d++){				
                    for (UN i=0; i < numNodes1 ; i++) {
                        value1[i] = (*elementMatrix)[offset+j*dofs2+di][dofs1*i+d];			    				    		
                        columnIndices1[i] =GO (dofs1* mapFirstRow->getGlobalElement( element1.getNode(i) )+d);
                    }
                    A->getBlock(1,0)->insertGlobalValues( row, columnIndices1(), value1() ); // Automatically adds entries if a value already exists                       
                }    
            }  
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
void FE<SC,LO,GO,NO>::assemblyNavierStokes(int dim,
	                                    string FETypeVelocity,
	                                    string FETypePressure,
	                                    int degree,
										int dofsVelocity,
										int dofsPressure,
										MultiVectorPtr_Type u_rep,
										MultiVectorPtr_Type p_rep,
	                                    BlockMatrixPtr_Type &A,
										BlockMultiVectorPtr_Type &resVec,
										SmallMatrix_Type coeff,
 										ParameterListPtr_Type params,
 										bool reAssemble,
 										string assembleMode,
	                                    bool callFillComplete,
	                                    int FELocExternal){
	

    UN FElocVel = checkFE(dim,FETypeVelocity); // Checks for different domains which belongs to a certain fetype
    UN FElocPres = checkFE(dim,FETypePressure); // Checks for different domains which belongs to a certain fetype

	ElementsPtr_Type elements = domainVec_.at(FElocVel)->getElementsC();

	ElementsPtr_Type elementsPres = domainVec_.at(FElocPres)->getElementsC();

	int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FElocVel)->getPointsRepeated();

	MapConstPtr_Type mapVel = domainVec_.at(FElocVel)->getMapRepeated();

	MapConstPtr_Type mapPres = domainVec_.at(FElocPres)->getMapRepeated();

	vec_dbl_Type solution(0);
	vec_dbl_Type solution_u;
	vec_dbl_Type solution_p;

	vec_dbl_ptr_Type rhsVec;

	/// Tupel construction follows follwing pattern:
	/// string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numVelo=3;
    if(FETypeVelocity == "P2")
        numVelo=6;
        
	if(dim==3){
		numVelo=4;
        if(FETypeVelocity == "P2")
            numVelo=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type vel ("Velocity",FETypeVelocity,dofsVelocity,numVelo);
	tuple_ssii_Type pres ("Pressure",FETypePressure,dofsPressure,dim+1);
	problemDisk->push_back(vel);
	problemDisk->push_back(pres);

	if(assemblyFEElements_.size()== 0){
        if(params->sublist("Parameter").get("Newtonian",true) == false)
	 	    initAssembleFEElements("NavierStokesNonNewtonian",problemDisk,elements, params,pointsRep,domainVec_.at(FElocVel)->getElementMap()); // In cas of non Newtonian Fluid
        else
        	initAssembleFEElements("NavierStokes",problemDisk,elements, params,pointsRep,domainVec_.at(FElocVel)->getElementMap());
    }
	else if(assemblyFEElements_.size() != elements->numberElements())
	     TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );

	//SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement));

	MultiVectorPtr_Type resVec_u = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocVel)->getMapVecFieldRepeated(), 1 ) );
    MultiVectorPtr_Type resVec_p = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocPres)->getMapRepeated(), 1 ) );
	
	BlockMultiVectorPtr_Type resVecRep = Teuchos::rcp( new BlockMultiVector_Type( 2) );
	resVecRep->addBlock(resVec_u,0);
	resVecRep->addBlock(resVec_p,1);


	for (UN T=0; T<assemblyFEElements_.size(); T++) {
		vec_dbl_Type solution(0);

		solution_u = getSolution(elements->getElement(T).getVectorNodeList(), u_rep,dofsVelocity);
		solution_p = getSolution(elementsPres->getElement(T).getVectorNodeList(), p_rep,dofsPressure);

		solution.insert( solution.end(), solution_u.begin(), solution_u.end() );
		solution.insert( solution.end(), solution_p.begin(), solution_p.end() );

		assemblyFEElements_[T]->updateSolution(solution);
 
 		SmallMatrixPtr_Type elementMatrix;

		if(assembleMode == "Jacobian"){
			assemblyFEElements_[T]->assembleJacobian();
		    
            elementMatrix = assemblyFEElements_[T]->getJacobian(); 

			assemblyFEElements_[T]->advanceNewtonStep(); // n genereal non linear solver step

			if(reAssemble)
				addFeBlock(A, elementMatrix, elements->getElement(T), mapVel, 0, 0, problemDisk);
			else
				addFeBlockMatrix(A, elementMatrix, elements->getElement(T), elementsPres->getElement(T), mapVel, mapPres, problemDisk);
		}
   		if(assembleMode == "FixedPoint"){

            AssembleFENavierStokesPtr_Type elTmp = Teuchos::rcp_dynamic_cast<AssembleFENavierStokes_Type>(assemblyFEElements_[T] );
            
            elTmp->assembleFixedPoint();
		    
            elementMatrix = elTmp->getFixedPointMatrix(); 

			assemblyFEElements_[T]->advanceNewtonStep(); // n genereal non linear solver step

			if(reAssemble)
				addFeBlock(A, elementMatrix, elements->getElement(T), mapVel, 0, 0, problemDisk);
			else
				addFeBlockMatrix(A, elementMatrix, elements->getElement(T), elementsPres->getElement(T),mapVel, mapPres, problemDisk);

        }
		if(assembleMode == "Rhs"){
			AssembleFENavierStokesPtr_Type elTmp = Teuchos::rcp_dynamic_cast<AssembleFENavierStokes_Type>(assemblyFEElements_[T] );
			elTmp->setCoeff(coeff);// Coeffs from time discretization. Right now default [1][1] // [0][0]
		    assemblyFEElements_[T]->assembleRHS();
		    rhsVec = assemblyFEElements_[T]->getRHS(); 
			addFeBlockMv(resVecRep, rhsVec, elements->getElement(T),elementsPres->getElement(T), dofsVelocity,dofsPressure);
		}

			
	}
	if (callFillComplete && reAssemble && assembleMode != "Rhs" )
	    A->getBlock(0,0)->fillComplete( domainVec_.at(FElocVel)->getMapVecFieldUnique(),domainVec_.at(FElocVel)->getMapVecFieldUnique());
	else if(callFillComplete && !reAssemble && assembleMode != "Rhs"){
		A->getBlock(0,0)->fillComplete();
	    A->getBlock(1,0)->fillComplete(domainVec_.at(FElocVel)->getMapVecFieldUnique(),domainVec_.at(FElocPres)->getMapUnique());
	    A->getBlock(0,1)->fillComplete(domainVec_.at(FElocPres)->getMapUnique(),domainVec_.at(FElocVel)->getMapVecFieldUnique());
	    A->getBlock(1,1)->fillComplete();
	}

	if(assembleMode == "Rhs"){

		MultiVectorPtr_Type resVecUnique_u = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocVel)->getMapVecFieldUnique(), 1 ) );
		MultiVectorPtr_Type resVecUnique_p = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocPres)->getMapUnique(), 1 ) );

		resVecUnique_u->putScalar(0.);
		resVecUnique_p->putScalar(0.);

		resVecUnique_u->exportFromVector( resVec_u, true, "Add" );
		resVecUnique_p->exportFromVector( resVec_p, true, "Add" );

		resVec->addBlock(resVecUnique_u,0);
		resVec->addBlock(resVecUnique_p,1);
	}


}

/*!

 \brief Inserting local rhsVec into global residual Mv;


@param[in] res BlockMultiVector of residual vec; Repeated distribution; 2 blocks
@param[in] rhsVec sorted the same way as residual vec
@param[in] element of block1

*/

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock1,FiniteElement elementBlock2, int dofs1, int dofs2 ){

    Teuchos::ArrayRCP<SC>  resArray_block1 = res->getBlockNonConst(0)->getDataNonConst(0);

    Teuchos::ArrayRCP<SC>  resArray_block2 = res->getBlockNonConst(1)->getDataNonConst(0);

	vec_LO_Type nodeList_block1 = elementBlock1.getVectorNodeList();

	vec_LO_Type nodeList_block2 = elementBlock2.getVectorNodeList();

	for(int i=0; i< nodeList_block1.size() ; i++){
		for(int d=0; d<dofs1; d++){
			resArray_block1[nodeList_block1[i]*dofs1+d] += (*rhsVec)[i*dofs1+d];
        }
	}
	int offset = nodeList_block1.size()*dofs1;

	for(int i=0; i < nodeList_block2.size(); i++){
		for(int d=0; d<dofs2; d++)
			resArray_block2[nodeList_block2[i]*dofs2+d] += (*rhsVec)[i*dofs2+d+offset];
	}

}


/*!

 \brief Adding FEBlock (row,column) to FE Blockmatrix
@todo column indices pre determine

@param[in] &A Global Block Matrix
@param[in] elementMatrix Stiffness matrix of one element
@param[in] element Corresponding finite element
@param[in] map Map that corresponds to repeated nodes of first block
@param[in] map Map that corresponds to repeated nodes of second block

*/
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::addFeBlock(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type mapRow, int row, int column, tuple_disk_vec_ptr_Type problemDisk){
		
		int dofs1 = std::get<2>(problemDisk->at(row));

		int numNodes1 = std::get<3>(problemDisk->at(row));

		int dofsBlock1 = dofs1*numNodes1;

		Teuchos::Array<SC> value( numNodes1, 0. );
        Teuchos::Array<GO> columnIndices( numNodes1, 0 );

		for (UN i=0; i < numNodes1 ; i++) {
			for(int di=0; di<dofs1; di++){
				GO rowID =GO (dofs1* mapRow->getGlobalElement( element.getNode(i) )+di);
				for(int d=0; d<dofs1; d++){
					for (UN j=0; j < columnIndices.size(); j++){
		                columnIndices[j] = GO ( dofs1 * mapRow->getGlobalElement( element.getNode(j) ) + d );
						value[j] = (*elementMatrix)[dofs1*i+di][dofs1*j+d];	
					}
			  		A->getBlock(row,column)->insertGlobalValues( rowID, columnIndices(), value() ); // Automatically adds entries if a value already exists 
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
void FE<SC,LO,GO,NO>::initAssembleFEElements(string elementType,tuple_disk_vec_ptr_Type problemDisk,ElementsPtr_Type elements, ParameterListPtr_Type params,vec2D_dbl_ptr_Type pointsRep, MapConstPtr_Type elementMap){
    
	vec2D_dbl_Type nodes;
	for (UN T=0; T<elements->numberElements(); T++) {
		
		nodes = getCoordinates(elements->getElement(T).getVectorNodeList(), pointsRep);

		AssembleFEFactory<SC,LO,GO,NO> assembleFEFactory;

		AssembleFEPtr_Type assemblyFE = assembleFEFactory.build(elementType,elements->getElement(T).getFlag(),nodes, params,problemDisk);
        //  
        assemblyFE->setGlobalElementID(elementMap->getGlobalElement(T));

		assemblyFEElements_.push_back(assemblyFE);

	}

}

/*!

 \brief Returns coordinates of local node ids

@param[in] localIDs
@param[in] points
@param[out] coordinates 

*/

template <class SC, class LO, class GO, class NO>
vec2D_dbl_Type FE<SC,LO,GO,NO>::getCoordinates(vec_LO_Type localIDs, vec2D_dbl_ptr_Type points){

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
vec_dbl_Type FE<SC,LO,GO,NO>::getSolution(vec_LO_Type localIDs, MultiVectorPtr_Type u_rep, int dofsVelocity){

    Teuchos::ArrayRCP<SC>  uArray = u_rep->getDataNonConst(0);
	
	vec_dbl_Type solution(0);
	for(int i=0; i < localIDs.size() ; i++){
		for(int d=0; d<dofsVelocity; d++)
			solution.push_back(uArray[localIDs[i]*dofsVelocity+d]);
	}

    return solution;
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::applyBTinv( vec3D_dbl_ptr_Type& dPhiIn,
                                    vec3D_dbl_Type& dPhiOut,
                                    const SmallMatrix<SC>& Binv){
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
    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyEmptyMatrix(MatrixPtr_Type &A){
    A->fillComplete();
}
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyIdentity(MatrixPtr_Type &A){
    Teuchos::Array<SC> value(1, Teuchos::ScalarTraits<SC>::one() );
    Teuchos::Array<GO> index(1);
    MapConstPtr_Type map = A->getMap();
    for (int i=0; i<A->getNodeNumRows(); i++) {
        index[0] = map->getGlobalElement( i );
        A->insertGlobalValues( index[0], index(), value() );
    }
    A->fillComplete();
}

// Assembling the nonlinear reaction part of Reaction-Diffusion equation
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyReactionTerm(int dim,
                                           std::string FEType,
                                           MatrixPtr_Type &A,
                                           MultiVectorPtr_Type u,
                                           bool callFillComplete,
                     					   std::vector<SC>& funcParameter,
										   RhsFunc_Type reactionFunc){

    //TEUCHOS_TEST_FOR_EXCEPTION( u->getNumVectors()>1, std::logic_error, "Implement for numberMV > 1 ." );
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    
    UN FEloc = checkFE(dim,FEType);
    
	ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

	MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

	vec2D_dbl_ptr_Type     phi;
	vec_dbl_ptr_Type    weights = Teuchos::rcp(new vec_dbl_Type(0));

	UN extraDeg = Helper::determineDegree( dim, FEType, Std); //Elementwise assembly of grad u

	UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Std, extraDeg);

	Helper::getPhi(phi, weights, dim, FEType, deg);
	
    // We have a scalar value of concentration in each point
	vec_dbl_Type uLoc( weights->size() , -1. );
	Teuchos::ArrayRCP< const SC > uArray = u->getData(0);

    std::vector<double> valueFunc(1);

    SC* paras = &(funcParameter[0]);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);

	for (UN T=0; T<elements->numberElements(); T++) {
        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        // Building u
        for (int w=0; w<phi->size(); w++){ //quadpoints
            uLoc[w] = 0.;
            for (int i=0; i < phi->at(0).size(); i++) { // points of element
                LO index = elements->getElement(T).getNode(i);
                uLoc[w] += uArray[index] * phi->at(w).at(i);
            }            
        }

        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<SC> value( phi->at(0).size(), 0. );
            Teuchos::Array<GO> indices( phi->at(0).size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<phi->size(); w++) {
                    value[j] += weights->at(w) * uLoc[w] * (*phi)[w][i] ;                                         
                }
                reactionFunc(&value[j], &valueFunc[0] ,paras);

                value[j] *= valueFunc[0] * absDetB;
                if (setZeros_ && std::fabs(value[j]) < myeps_) {
                    value[j] = 0.;
                }
                indices[j] = GO( map->getGlobalElement( elements->getElement(T).getNode(j) ));

            }
            GO row = GO ( map->getGlobalElement( elements->getElement(T).getNode(i) ) );
            /*if( domainVec_.at(FEloc)->getComm()->getRank() == 1){
                cout << " Inserting in Row " << row ;
                for(int j=0 ; j< indices.size() ; j++)
                    cout << " with indices " << indices[j] << " " ;
                cout << "on rank " <<domainVec_.at(FEloc)->getComm()->getRank() << endl;
            }*/

            A->insertGlobalValues( row, indices(), value() );     
            if( domainVec_.at(FEloc)->getComm()->getRank() == 1)
                cout << " Inserted " << row << endl;      
        }
    }
    
    if (callFillComplete)
        A->fillComplete();
}


// Assembling the nonlinear reaction part of Reaction-Diffusion equation
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLinearReactionTerm(int dim,
                                           std::string FEType,
                                           MatrixPtr_Type &A,
                                           bool callFillComplete,
                     					   std::vector<SC>& funcParameter,
										   RhsFunc_Type reactionFunc){

    //TEUCHOS_TEST_FOR_EXCEPTION( u->getNumVectors()>1, std::logic_error, "Implement for numberMV > 1 ." );
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    
    UN FEloc = checkFE(dim,FEType);
    
	ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

	MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

	vec2D_dbl_ptr_Type     phi;
	vec_dbl_ptr_Type    weights = Teuchos::rcp(new vec_dbl_Type(0));

	UN extraDeg = Helper::determineDegree( dim, FEType, Std); //Elementwise assembly of grad u

	UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Std, extraDeg);

	Helper::getPhi(phi, weights, dim, FEType, deg);
	
    std::vector<double> valueFunc(1);

    SC* paras = &(funcParameter[0]);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);

	for (UN T=0; T<elements->numberElements(); T++) {
        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<SC> value( phi->at(0).size(), 0. );
            Teuchos::Array<GO> indices( phi->at(0).size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<phi->size(); w++) {
                    value[j] += weights->at(w) * (*phi)[w][j] * (*phi)[w][i] ;                                         
                }
                reactionFunc(&value[j], &valueFunc[0] ,paras);

                value[j] *= valueFunc[0] * absDetB;
                if (setZeros_ && std::fabs(value[j]) < myeps_) {
                    value[j] = 0.;
                }
                indices[j] = GO( map->getGlobalElement( elements->getElement(T).getNode(j) ));

            }
            GO row = GO ( map->getGlobalElement( elements->getElement(T).getNode(i) ) );
          
            A->insertGlobalValues( row, indices(), value() );     
            
        }
    }
    
    if (callFillComplete)
        A->fillComplete();
}

// Assembling the nonlinear reaction part of Reaction-Diffusion equation
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyDReactionTerm(int dim,
                                           std::string FEType,
                                           MatrixPtr_Type &A,
                                           MultiVectorPtr_Type u,
                                           bool callFillComplete,
                     					   std::vector<SC>& funcParameter,
										   RhsFunc_Type reactionFunc){

    //TEUCHOS_TEST_FOR_EXCEPTION( u->getNumVectors()>1, std::logic_error, "Implement for numberMV > 1 ." );
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    
    UN FEloc = checkFE(dim,FEType);
    
	ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

	MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

	vec2D_dbl_ptr_Type     phi;
    vec3D_dbl_ptr_Type 	dPhi;
	vec_dbl_ptr_Type    weights = Teuchos::rcp(new vec_dbl_Type(0));

	UN extraDeg = Helper::determineDegree( dim, FEType, Std); //Elementwise assembly of grad u

	UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Std, extraDeg);

	Helper::getPhi(phi, weights, dim, FEType, deg);
	
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);

    // We have a scalar value of concentration in each point
	vec2D_dbl_Type duLoc( weights->size() ,vec_dbl_Type(dim ,-1. ));
	Teuchos::ArrayRCP< const SC > uArray = u->getData(0);

    std::vector<double> valueFunc(1);

    SC* paras = &(funcParameter[0]);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);

	for (UN T=0; T<elements->numberElements(); T++) {

        
        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyBTinv( dPhi, dPhiTrans, Binv );

        for (int w=0; w<dPhiTrans.size(); w++){ //quads points
            for (int i=0; i < dPhiTrans[0].size(); i++) {
                LO index = elements->getElement(T).getNode(i) ;
                for (int d2=0; d2<dim; d2++)
                    duLoc[w][d2] += uArray[index] * dPhiTrans[w][i][d2];
            }
            
        }

        
        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<SC> value( phi->at(0).size(), 0. );
            Teuchos::Array<GO> indices( phi->at(0).size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN d2=0; d2<dim; d2++){
                    for (UN w=0; w<phi->size(); w++) {
                        value[j] += weights->at(w) * duLoc[w][d2] * (*phi)[w][i] ;                                         
                    }
                }
                reactionFunc(&value[j], &valueFunc[0] ,paras);

                value[j] *= valueFunc[0] * absDetB;
                if (setZeros_ && std::fabs(value[j]) < myeps_) {
                    value[j] = 0.;
                }
                indices[j] = GO( map->getGlobalElement( elements->getElement(T).getNode(j) ));

            }
            GO row = GO ( map->getGlobalElement( elements->getElement(T).getNode(i) ) );
            A->insertGlobalValues( row, indices(), value() );           
        }
    }
    
    if (callFillComplete)
        A->fillComplete();
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
void FE<SC,LO,GO,NO>::assemblyLaplaceAssFE(int dim,
                                        string FEType,
                                        int degree,
                                        int dofs,
                                        BlockMatrixPtr_Type &A,
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
        initAssembleFEElements("Laplace",problemDisk,elements, params,pointsRep,domainVec_.at(0)->getElementMap());
    else if(assemblyFEElements_.size() != elements->numberElements())
         TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Number Elements not the same as number assembleFE elements." );
    for (UN T=0; T<elements->numberElements(); T++) {
        assemblyFEElements_[T]->assembleJacobian();
        SmallMatrixPtr_Type elementMatrix = assemblyFEElements_[T]->getJacobian(); 
      	addFeBlock(A, elementMatrix, elements->getElement(T), map, 0, 0, problemDisk);

    }
    if(callFillComplete)
        A->getBlock(0,0)->fillComplete();
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLaplaceDiffusion(int dim,
                                        std::string FEType,
                                        int degree,
                                        MatrixPtr_Type &A,
										vec2D_dbl_Type diffusionTensor,
                                        bool callFillComplete,
                                        int FELocExternal){
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    UN FEloc;
    if (FELocExternal<0)
        FEloc = checkFE(dim,FEType);
    else
        FEloc = FELocExternal;
    
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,FEType,Grad,Grad);//+1;
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);


 	SmallMatrix<SC> diffusionT(dim);
	// Linear Diffusion Tensor
	if(diffusionTensor.size()==0 || diffusionTensor.size() < dim ){
		vec2D_dbl_Type diffusionTensor(3,vec_dbl_Type(3,0));
		for(int i=0; i< dim; i++){
			diffusionTensor[i][i]=1.;
		}
	}

	for(int i=0; i< dim; i++){
		for(int j=0; j<dim; j++){
			diffusionT[i][j]=diffusionTensor[i][j];
		}
	}
	//Teuchos::ArrayRCP< SC >  linearDiff = diffusionTensor->getDataNonConst( 0 );
	//cout << "Assembly Info " << "num Elements " <<  elements->numberElements() << " num Nodes " << pointsRep->size()  << endl;
    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyBTinv( dPhi, dPhiTrans, Binv );

        vec3D_dbl_Type dPhiTransDiff( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyDiff( dPhiTrans, dPhiTransDiff, diffusionT );

        for (UN i=0; i < dPhiTrans[0].size(); i++) {
            Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
            Teuchos::Array<GO> indices( dPhiTrans[0].size(), 0 );

            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<dPhiTrans.size(); w++) {
                    for (UN d=0; d<dim; d++){
                        value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTransDiff[w][j][d];
                    }
                }
                value[j] *= absDetB;
                indices[j] = map->getGlobalElement( elements->getElement(T).getNode(j) );
                if (setZeros_ && std::fabs(value[j]) < myeps_) {
                    value[j] = 0.;
                }
            }
            GO row = map->getGlobalElement( elements->getElement(T).getNode(i) );

            A->insertGlobalValues( row, indices(), value() );
        }


    }
    if (callFillComplete)
        A->fillComplete();

}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::applyDiff( vec3D_dbl_Type& dPhiIn,
                                    vec3D_dbl_Type& dPhiOut,
                                    SmallMatrix<SC>& diffT){
    UN dim = diffT.size();
    for (UN w=0; w<dPhiIn.size(); w++){
        for (UN i=0; i < dPhiIn[w].size(); i++) {
            for (UN d1=0; d1<dim; d1++) {
                for (UN d2=0; d2<dim; d2++) {
                    dPhiOut[w][i][d1] += dPhiIn[w][i][d2]* diffT[d2][d1];
                }
            }
        }
    }
}
    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyAceGenTPM(    MatrixPtr_Type &A00,
                                            MatrixPtr_Type &A01,
                                            MatrixPtr_Type &A10,
                                            MatrixPtr_Type &A11,
                                            MultiVectorPtr_Type &F0,
                                            MultiVectorPtr_Type &F1,
                                            MapPtr_Type &mapRepeated1,
                                            MapPtr_Type &mapRepeated2,
                                            ParameterListPtr_Type parameterList,
                                            MultiVectorPtr_Type u_repeatedNewton,
                                            MultiVectorPtr_Type p_repeatedNewton,
                                            MultiVectorPtr_Type u_repeatedTime,
                                            MultiVectorPtr_Type p_repeatedTime,
                                            bool update,
                                            bool updateHistory)
{
    

    std::string tpmType = parameterList->sublist("Parameter").get("TPM Type","Biot");
    
    int dim = domainVec_[0]->getDimension();
    int idata = 1; //= we should init this
    int ic = -1; int ng = -1;
    
    //ed.hp:history previous (timestep); previous solution (velocity and acceleration)
    //ed.ht:? same length as hp
    ElementsPtr_Type elements1 = domainVec_[0]->getElementsC();
    ElementsPtr_Type elements2 = domainVec_[1]->getElementsC();
    
    int sizeED = 24; /* 2D case for P2 elements:
                      12 velocities, 12 accelerations (2 dof per P2 node)
                      */
    if (dim==3)
        sizeED = 60;/* 3D case for P2 elements:
                       30 velocities, 30 accelerations (3 dof per P2 node)
                    */
    if (ed_.size()==0){
        for (UN T=0; T<elements1->numberElements(); T++)
            ed_.push_back( Teuchos::rcp(new DataElement( sizeED )) );
    }
    
    std::vector<ElementSpec> es_vec( parameterList->sublist("Parameter").get("Number of materials",1) , ElementSpec());
    vec2D_dbl_Type dataVec( parameterList->sublist("Parameter").get("Number of materials",1), vec_dbl_Type(6,0.) );
    
    for (int i=0; i<dataVec.size(); i++) {
        if (tpmType == "Biot") {
            if (dim==2) {
                dataVec[i][0] = parameterList->sublist("Timestepping Parameter").get("Newmark gamma",0.5);
                dataVec[i][1] = parameterList->sublist("Timestepping Parameter").get("Newmark beta",0.25);
                dataVec[i][2] = parameterList->sublist("Parameter").get("initial volume fraction solid material"+std::to_string(i+1),0.5); //do we need this?
                dataVec[i][3] = parameterList->sublist("Parameter").get("Darcy parameter material"+std::to_string(i+1),1.e-2);
                dataVec[i][4] = parameterList->sublist("Parameter").get("Youngs modulus material"+std::to_string(i+1),60.e6);
                dataVec[i][5] = parameterList->sublist("Parameter").get("Poisson ratio material"+std::to_string(i+1),0.3);
            }
            else if (dim==3) {
                dataVec[i].resize(12);
                dataVec[i][0] = parameterList->sublist("Parameter").get("Youngs modulus material"+std::to_string(i+1),2.e5);
                dataVec[i][1] = parameterList->sublist("Parameter").get("Poisson ratio material"+std::to_string(i+1),0.3);
                dataVec[i][2] = 0.; //body force x
                dataVec[i][3] = 0.; //body force y
                dataVec[i][4] = parameterList->sublist("Parameter").get("body force z"+std::to_string(i+1),0.);; //body force z
                dataVec[i][5] = parameterList->sublist("Parameter").get("initial volume fraction solid material"+std::to_string(i+1),0.67);
                dataVec[i][6] = parameterList->sublist("Parameter").get("Darcy parameter material"+std::to_string(i+1),0.01);
                dataVec[i][7] = 2000.; //effective density solid
                dataVec[i][8] = 1000.; //effective density fluid?
                dataVec[i][9] = 9.81;  // gravity
                dataVec[i][10] = parameterList->sublist("Timestepping Parameter").get("Newmark gamma",0.5);
                dataVec[i][11] = parameterList->sublist("Timestepping Parameter").get("Newmark beta",0.25);
            }
        }
        
        
        else if (tpmType == "Biot-StVK") {
            dataVec[i][0] = parameterList->sublist("Parameter").get("Youngs modulus material"+std::to_string(i+1),60.e6);
            dataVec[i][1] = parameterList->sublist("Parameter").get("Poisson ratio material"+std::to_string(i+1),0.3);
            dataVec[i][2] = parameterList->sublist("Parameter").get("initial volume fraction solid material"+std::to_string(i+1),0.5); //do we need this?
            dataVec[i][3] = parameterList->sublist("Parameter").get("Darcy parameter material"+std::to_string(i+1),1.e-2);
            dataVec[i][4] = parameterList->sublist("Timestepping Parameter").get("Newmark gamma",0.5);
            dataVec[i][5] = parameterList->sublist("Timestepping Parameter").get("Newmark beta",0.25);
        }
    }
    
    for (int i=0; i<es_vec.size(); i++){
        if(tpmType == "Biot"){
            if (dim==2)
                this->SMTSetElSpecBiot(&es_vec[i] ,&idata, ic, ng, dataVec[i]);
            else if(dim==3)
                this->SMTSetElSpecBiot3D(&es_vec[i] ,&idata, ic, ng, dataVec[i]);
        }
        else if(tpmType == "Biot-StVK")
            this->SMTSetElSpecBiotStVK(&es_vec[i] ,&idata, ic, ng, dataVec[i]);
    }
    LO elementSizePhase = elements1->nodesPerElement();
    LO sizePhase = dim * elementSizePhase;
    LO sizePressure = elements2->nodesPerElement();
    GO sizePhaseGlobal = A00->getMap()->getMaxAllGlobalIndex()+1;
    int workingVectorSize;
    if(tpmType == "Biot"){
        if (dim==2)
            workingVectorSize = 5523;
        else if(dim==3)
            workingVectorSize = 1817;
    }
    else if(tpmType == "Biot-StVK")
        workingVectorSize = 5223;
    
    double* v = new double [workingVectorSize];

    // nd sind Nodalwerte, Anzahl an structs in nd sollte den Knoten entsprechen, bei P2-P1 in 2D also 9
    // In X stehen die Koordinaten, X[0] ist x-Koordinate, X[1] ist y-Koordinate, etc.
    // nd->X[0]
    // at ist die Loesung im letzten Newtonschritt.
    // nd[0]->at[0];
    // ap ist die Loesung im letzten Zeitschritt.
    // nd[0]->ap[0]
    // rdata ist die Zeitschrittweite, RD_TimeIncrement wird in sms.h definiert, entsprechend wird auch die Laenge von rdata dort definiert. Standard 400, aber auch nicht gesetzt. Wert muss selber initialisiert werden; eventuell kuerzer moeglich.

    std::vector<double> rdata(RD_TimeIncrement+1, 0.);

    rdata[RD_TimeIncrement] = parameterList->sublist("Timestepping Parameter").get("dt",0.01);
    
    NodeSpec *ns=NULL;//dummy not need in SKR
        
    NodeData** nd = new NodeData*[ elementSizePhase + sizePressure ];

    for (int i=0; i<elementSizePhase + sizePressure; i++){
        nd[i] = new NodeData();
    }

    int numNodes = elementSizePhase + sizePressure;
    
    vec2D_dbl_Type xFull( numNodes, vec_dbl_Type(dim,0.) );
    vec2D_dbl_Type atFull( numNodes, vec_dbl_Type(dim,0.) );
    vec2D_dbl_Type apFull( numNodes, vec_dbl_Type(dim,0.) );
    
    for (int i=0; i<elementSizePhase + sizePressure; i++) {
        nd[i]->X = &(xFull[i][0]);
        nd[i]->at = &(atFull[i][0]);
        nd[i]->ap = &(apFull[i][0]);
    }
    
    GO offsetMap1 = dim * mapRepeated1->getMaxAllGlobalIndex()+1;
    vec2D_dbl_ptr_Type pointsRepU = domainVec_.at(0)->getPointsRepeated();
    vec2D_dbl_ptr_Type pointsRepP = domainVec_.at(1)->getPointsRepeated();
    
    Teuchos::ArrayRCP< const SC > uArrayNewton = u_repeatedNewton->getData(0);
    Teuchos::ArrayRCP< const SC > pArrayNewton = p_repeatedNewton->getData(0);
    Teuchos::ArrayRCP< const SC > uArrayTime = u_repeatedTime->getData(0);
    Teuchos::ArrayRCP< const SC > pArrayTime = p_repeatedTime->getData(0);

    double** mat = new double*[sizePhase+sizePressure];
    for (int i=0; i<sizePhase+sizePressure; i++){
        mat[i] = new double[sizePhase+sizePressure];
    }
    
    Teuchos::ArrayRCP<SC> fValues0 = F0->getDataNonConst(0);
    Teuchos::ArrayRCP<SC> fValues1 = F1->getDataNonConst(0);
    
    // Element loop

    ElementData ed = ElementData();
    for (UN T=0; T<elements1->numberElements(); T++) {
        
        std::vector<double> tmpHp = ed_[T]->getHp(); // Dies sind die alten Daten
        std::vector<double> tmpHt = ed_[T]->getHt(); // Dies sind die neuen Daten nachdem das Element aufgerufen wurde, wir hier eigentlich nicht als Variable in ed_ benoetigt.
        ed.hp = &tmpHp[0];
        ed.ht = &tmpHt[0];
        
        int materialFlag = elements1->getElement(T).getFlag();
        TEUCHOS_TEST_FOR_EXCEPTION( materialFlag>es_vec.size()-1, std::runtime_error, "There are not enought material parameters initialized." ) ;
        int counter=0;
        //Newtonloesung at und Zeitschrittloesung ap
        for (int j=0; j<elementSizePhase; j++) {
            for (int d=0; d<dim; d++) {
                LO index = dim * elements1->getElement(T).getNode(j)+d;//dim * elements1->at(T).at( j ) + d;
                atFull[j][d] = uArrayNewton[index];
                apFull[j][d] = uArrayTime[index];
            }
        }
        for (int j=0; j<sizePressure; j++) {
            LO index = elements2->getElement(T).getNode(j);//elements2->at(T).at( j );
            atFull[elementSizePhase+j][0] = pArrayNewton[index];
            apFull[elementSizePhase+j][0] = pArrayTime[index];
        }
        
        //Nodes
        for (int j=0; j<elementSizePhase; j++ ) {
            LO index = elements1->getElement(T).getNode(j);
            for (int d=0; d<dim; d++) {
                xFull[j][d] = (*pointsRepU)[index][d];
            }
        }
        for (int j=0; j<sizePressure; j++ ) {
            LO index = elements2->getElement(T).getNode(j);
            for (int d=0; d<dim; d++) {
                xFull[elementSizePhase+j][d] = (*pointsRepP)[index][d];
            }
        }
        vec_dbl_Type p( sizePhase+sizePressure , 0. );

        for (int i=0; i<sizePhase+sizePressure; i++){
            for (int j=0; j<sizePhase+sizePressure; j++)
                mat[i][j] = 0.;
        }
        // element assembly
        if(tpmType == "Biot"){
            if(dim==2)
                this->SKR_Biot( v, &es_vec[materialFlag], &ed, &ns, nd , &rdata[0], &idata, &p[0], mat  );
            else if (dim==3)
                this->SKR_Biot3D( v, &es_vec[materialFlag], &ed, &ns, nd , &rdata[0], &idata, &p[0], mat  );
        }
        else if(tpmType == "Biot-StVK")
            this->SKR_Biot_StVK( v, &es_vec[materialFlag], &ed, &ns, nd , &rdata[0], &idata, &p[0], mat  );
        
        if (updateHistory)
            ed_[T]->setHp( ed.ht );
        
        if (update) {
                    
            // A00 & A01
            for (UN i=0; i < sizePhase; i++) {
                Teuchos::Array<SC> value00( sizePhase, 0. );
                Teuchos::Array<GO> indices00( sizePhase, 0 );
                for (UN j=0; j < value00.size(); j++) {
                    
                    value00[j] = mat[i][j];
                    
                    LO tmpJ = j/dim;
                    LO index = elements1->getElement(T).getNode(tmpJ);
                    if (j%dim==0)
                        indices00[j] = dim * mapRepeated1->getGlobalElement( index );
                    else if (j%dim==1)
                        indices00[j] = dim * mapRepeated1->getGlobalElement( index ) + 1;
                    else if (j%dim==2)
                        indices00[j] = dim * mapRepeated1->getGlobalElement( index ) + 2;
                }
                
                Teuchos::Array<SC> value01( sizePressure, 0. );
                Teuchos::Array<GO> indices01( sizePressure, 0 );

                for (UN j=0; j < value01.size(); j++) {
                    value01[j] = mat[i][sizePhase+j];
                    LO index = elements2->getElement(T).getNode(j);
                    indices01[j] = mapRepeated2->getGlobalElement( index );
                }
                
                GO row;
                LO tmpI = i/dim;
                LO index = elements1->getElement(T).getNode(tmpI);
                if (i%dim==0)
                    row = dim * mapRepeated1->getGlobalElement( index );
                else if (i%dim==1)
                    row = dim * mapRepeated1->getGlobalElement( index ) + 1;
                else if (i%dim==2)
                    row = dim * mapRepeated1->getGlobalElement( index ) + 2;
                
                A00->insertGlobalValues( row, indices00(), value00() );
                A01->insertGlobalValues( row, indices01(), value01() );
                
                if (i%dim==0)
                    fValues0[ dim*index ] += p[ i ];
                else if (i%dim==1)
                    fValues0[ dim*index+1 ] += p[ i ];
                else if (i%dim==2)
                    fValues0[ dim*index+2 ] += p[ i ];
            }
            // A10 & A11
            for (UN i=0; i < sizePressure; i++) {
                Teuchos::Array<SC> value10( sizePhase   , 0. );
                Teuchos::Array<GO> indices10( sizePhase   , 0 );
                for (UN j=0; j < value10.size(); j++) {
                    value10[j] = mat[sizePhase+i][j];
                    
                    LO tmpJ = j/dim;
                    LO index = elements1->getElement(T).getNode(tmpJ);
                    if (j%dim==0)
                        indices10[j] = dim * mapRepeated1->getGlobalElement( index );
                    else if (j%dim==1)
                        indices10[j] = dim * mapRepeated1->getGlobalElement( index ) + 1;
                    else if (j%dim==2)
                        indices10[j] = dim * mapRepeated1->getGlobalElement( index ) + 2;
                }
                
                Teuchos::Array<SC> value11( sizePressure, 0. );
                Teuchos::Array<GO> indices11( sizePressure, 0 );
                for (UN j=0; j < value11.size(); j++) {
                    value11[j] = mat[sizePhase+i][sizePhase+j];
                    
                    LO index = elements2->getElement(T).getNode(j);
                    indices11[j] = mapRepeated2->getGlobalElement( index );
                }

                
                LO index2 = elements2->getElement(T).getNode(i);
                GO row = mapRepeated2->getGlobalElement( index2 );
                A10->insertGlobalValues( row, indices10(), value10() );
                A11->insertGlobalValues( row, indices11(), value11() );
                
                fValues1[ index2 ] += p[ sizePhase + i ];
            }
        }
    }
    
    for (int i=0; i<sizePhase+sizePressure; i++)
        delete [] mat[i];
    delete [] mat;
    
    delete [] v;
    
    for (int i=0; i<elementSizePhase+sizePressure; i++)
        delete nd[i];
    
    delete [] nd;
    
    
    A00->fillComplete( A00->getMap("row"), A00->getMap("row") );
    A01->fillComplete( A10->getMap("row"), A00->getMap("row") );
    A10->fillComplete( A00->getMap("row"), A10->getMap("row") );
    A11->fillComplete( A10->getMap("row"), A10->getMap("row") );
    
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyMass(int dim,
                                     std::string FEType,
                                     std::string fieldType,
                                     MatrixPtr_Type &A,
                                     bool callFillComplete){

    TEUCHOS_TEST_FOR_EXCEPTION( FEType == "P0", std::logic_error, "Not implemented for P0" );
    UN FEloc = checkFE(dim,FEType);
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree(dim,FEType,FEType,Std,Std);

    Helper::getPhi( phi, weights, dim, FEType, deg );

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
        detB = B.computeDet( );
        absDetB = std::fabs(detB);

        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<SC> value( phi->at(0).size(), 0. );
            Teuchos::Array<GO> indices( phi->at(0).size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<phi->size(); w++) {
                    value[j] += weights->at(w) * (*phi)[w][i] * (*phi)[w][j];

                }
                value[j] *= absDetB;
                if (!fieldType.compare("Scalar")) {
                    indices[j] = map->getGlobalElement( elements->getElement(T).getNode(j) );
                }

            }
            if (!fieldType.compare("Scalar")) {
                GO row = map->getGlobalElement( elements->getElement(T).getNode(i) );
                A->insertGlobalValues( row, indices(), value() );
            }
            else if (!fieldType.compare("Vector")) {
                for (UN d=0; d<dim; d++) {
                    for (int j=0; j<indices.size(); j++) {
                        indices[j] = (GO) ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d );
                    }
                    GO row = (GO) ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d );
                    A->insertGlobalValues( row, indices(), value() );
                }
            }
            else
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Specify valid vieldType for assembly of mass matrix.");
        }

    }

    if (callFillComplete)
        A->fillComplete();
}


// Ueberladung der Assemblierung der Massematrix fuer FSI, da
// checkFE sonst auch fuer das Strukturproblem FEloc = 1 liefert (= Fluid)
// und somit die welche domain und Map in der Assemblierung genutzt wird.
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyMass(int dim,
                                     std::string FEType,
                                     std::string fieldType,
                                     MatrixPtr_Type &A,
                                     int FEloc, // 0 = Fluid, 2 = Struktur
                                     bool callFillComplete){

    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type	weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree(dim,FEType,FEType,Std,Std);

    Helper::getPhi( phi, weights, dim, FEType, deg );

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
        detB = B.computeDet( );
        absDetB = std::fabs(detB);

        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<SC> value( phi->at(0).size(), 0. );
            Teuchos::Array<GO> indices( phi->at(0).size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<phi->size(); w++) {
                    value[j] += weights->at(w) * (*phi)[w][i] * (*phi)[w][j];
                }
                value[j] *= absDetB;
                if (!fieldType.compare("Scalar")) {
                    indices[j] = map->getGlobalElement( elements->getElement(T).getNode(j) );
                }

            }
            if (!fieldType.compare("Scalar")) {
                GO row = map->getGlobalElement( elements->getElement(T).getNode(i) );
                A->insertGlobalValues( row, indices(), value() );
            }
            else if (!fieldType.compare("Vector")) {
                for (UN d=0; d<dim; d++) {
                    for (int j=0; j<indices.size(); j++) {
                        indices[j] = (GO) ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d );
                    }
                    GO row = (GO) ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d );
                    A->insertGlobalValues( row, indices(), value() );
                }
            }
            else
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Specify valid vieldType for assembly of mass matrix.");
        }


    }
    if (callFillComplete)
        A->fillComplete();
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLaplace(int dim,
                                        std::string FEType,
                                        int degree,
                                        MatrixPtr_Type &A,
                                        bool callFillComplete,
                                        int FELocExternal){
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    UN FEloc;
    if (FELocExternal<0)
        FEloc = checkFE(dim,FEType);
    else
        FEloc = FELocExternal;
    
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,FEType,Grad,Grad);
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyBTinv( dPhi, dPhiTrans, Binv );
        for (UN i=0; i < dPhiTrans[0].size(); i++) {
            Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
            Teuchos::Array<GO> indices( dPhiTrans[0].size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<dPhiTrans.size(); w++) {
                    for (UN d=0; d<dim; d++){
                        value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                    }
                }
                value[j] *= absDetB;
                indices[j] = map->getGlobalElement( elements->getElement(T).getNode(j) );
            }
            GO row = map->getGlobalElement( elements->getElement(T).getNode(i) );

            A->insertGlobalValues( row, indices(), value() );
        }


    }
    if (callFillComplete)
        A->fillComplete();

}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLaplaceVecField(int dim,
                                                std::string FEType,
                                                int degree,
                                                MatrixPtr_Type &A,
                                                bool callFillComplete){

    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P1-disc" || FEType == "P0",std::logic_error, "Not implemented for P0 or P1-disc");
    UN FEloc = checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree(dim,FEType,FEType,Grad,Grad);

    Helper::getDPhi(dPhi, weights, dim, FEType, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);


    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyBTinv( dPhi, dPhiTrans, Binv );

        for (UN i=0; i < dPhiTrans[0].size(); i++) {
            Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
            Teuchos::Array<GO> indices( dPhiTrans[0].size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<dPhiTrans.size(); w++) {
                    for (UN d=0; d<dim; d++)
                        value[j] += weights->at(w) * dPhiTrans[w][i][d] * dPhiTrans[w][j][d];
                }
                value[j] *= absDetB;
                if (setZeros_ && std::fabs(value[j]) < myeps_) {
                    value[j] = 0.;
                }
            }
            for (UN d=0; d<dim; d++) {
                for (UN j=0; j < indices.size(); j++)
                    indices[j] = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d );

                GO row = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d );
                A->insertGlobalValues( row, indices(), value() );
            }
        }
    }
    if (callFillComplete)
        A->fillComplete();
}
//this assembly used blas matrix-matrix multiplications. It determines the local stiffness matrix at once, but has some overhead due to zero off-diagonal blocks which are computed.
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLaplaceVecFieldV2(int dim,
                                                std::string FEType,
                                                int degree,
                                                MatrixPtr_Type &A,
                                                bool callFillComplete){

    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    UN FEloc = checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree(dim,FEType,FEType,Grad,Grad);

    Helper::getDPhi(dPhi, weights, dim, FEType, deg);

    Teuchos::BLAS<int, SC> teuchosBLAS;

    int nmbQuadPoints = dPhi->size();
    int nmbScalarDPhi = dPhi->at(0).size();
    int nmbAllDPhi = nmbScalarDPhi * dim;
    int nmbAllDPhiAllQaud = nmbQuadPoints * nmbAllDPhi;
    int sizeLocStiff = dim*dim;
    Teuchos::Array<SmallMatrix<double> > dPhiMat( nmbAllDPhiAllQaud, SmallMatrix<double>(dim) );
    this->buildFullDPhi( dPhi, dPhiMat ); //builds matrix from gradient of scalar phi


    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        Teuchos::Array<SmallMatrix<double> > allDPhiMatTrans( dPhiMat.size(), SmallMatrix<double>() );

        for (int i=0; i<allDPhiMatTrans.size(); i++) {
            SmallMatrix<double> res = dPhiMat[i] * Binv;
            allDPhiMatTrans[i] = res;
        }

        SmallMatrix<double> locStiffMat( nmbAllDPhi, 0. );

        for (int p=0; p<nmbQuadPoints; p++){

            double* allDPhiBlas = new double[ nmbAllDPhi * sizeLocStiff ];

            int offset = p * nmbAllDPhi;
            int offsetInArray = 0;
            for (int i=0; i<nmbAllDPhi; i++) {
                fillMatrixArray( allDPhiMatTrans[ offset + i ], allDPhiBlas, "rows",offsetInArray );
                offsetInArray += sizeLocStiff;
            }

            double* locStiffMatBlas = new double[ nmbAllDPhi * nmbAllDPhi ];

            teuchosBLAS.GEMM (Teuchos::TRANS, Teuchos::NO_TRANS, nmbAllDPhi, nmbAllDPhi, sizeLocStiff, 1., allDPhiBlas, sizeLocStiff/*lda of A not trans(A)! Otherwise result is wrong*/, allDPhiBlas, sizeLocStiff, 0., locStiffMatBlas, nmbAllDPhi);

            for (int i=0; i<nmbAllDPhi; i++) {
                for (int j=0; j<nmbAllDPhi; j++) {
                    locStiffMat[i][j] += weights->at(p) * locStiffMatBlas[ j * nmbAllDPhi + i ];
                }
            }

            delete [] allDPhiBlas;
            delete [] locStiffMatBlas;

        }

        for (UN i=0; i < nmbScalarDPhi; i++) {
            Teuchos::Array<SC> value( nmbAllDPhi, 0. );
            Teuchos::Array<GO> indices( nmbAllDPhi, 0 );
            for (UN d=0; d<dim; d++) {
                for (UN j=0; j < nmbScalarDPhi; j++){
                    value[ j * dim + d ] = absDetB * locStiffMat[dim * i + d][j];
                    indices[ j * dim + d ] = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d );
                }
                GO row = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d );
                A->insertGlobalValues( row, indices(), value() );
            }
        }
    }
    if (callFillComplete)
        A->fillComplete();
}

/*!
 \brief Assembly of Jacobian for nonlinear Laplace example
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] u_rep The current solution
@param[in] A Resulting matrix
@param[in] resVec Resulting residual
@param[in] params Params needed by the problem. Placeholder for now.
@param[in] assembleMode What should be assembled i.e. Rhs (residual) or the
Jacobian
@param[in] callFillComplete If Matrix A should be redistributed across MPI procs
at end of function
@param[in] FELocExternal ?
*/

template <class SC, class LO, class GO, class NO>
void FE<SC, LO, GO, NO>::assemblyNonlinearLaplace(
    int dim, std::string FEType, int degree, MultiVectorPtr_Type u_rep,
    BlockMatrixPtr_Type &A, BlockMultiVectorPtr_Type &resVec,
    ParameterListPtr_Type params, string assembleMode, bool callFillComplete,
    int FELocExternal) {

    ElementsPtr_Type elements = this->domainVec_.at(0)->getElementsC();

    // Only scalar laplace
    int dofs = 1;

    vec2D_dbl_ptr_Type pointsRep = this->domainVec_.at(0)->getPointsRepeated();
    MapConstPtr_Type map = this->domainVec_.at(0)->getMapRepeated();

    vec_dbl_Type solution_u;
    vec_dbl_ptr_Type rhsVec;

    int numNodes = 3;
    if (FEType == "P2") {
        numNodes = 6;
    }
    if (dim == 3) {
        numNodes = 4;
        if (FEType == "P2") {
            numNodes = 10;
        }
    }

    // Tupel construction follows follwing pattern:
    // string: Physical Entity (i.e. Velocity) , string: Discretisation (i.e.
    // "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per
    // element)
    tuple_disk_vec_ptr_Type problemDisk =
        Teuchos::rcp(new tuple_disk_vec_Type(0));
    tuple_ssii_Type temp("Solution", FEType, dofs, numNodes);
    problemDisk->push_back(temp);

    // Construct an assembler for each element if not already done
    if (assemblyFEElements_.size() == 0) {
        initAssembleFEElements("NonLinearLaplace", problemDisk, elements, params, pointsRep, domainVec_.at(0)->getElementMap());
    } else if (assemblyFEElements_.size() != elements->numberElements()) {
        TEUCHOS_TEST_FOR_EXCEPTION(
            true, std::logic_error,
            "Number Elements not the same as number assembleFE elements.");
    }

    MultiVectorPtr_Type resVec_u;
    BlockMultiVectorPtr_Type resVecRep;

    if (assembleMode != "Rhs") {
        // add new or overwrite existing block (0,0) of system matrix
        // This is done in specific problem class for most other problems
        // Placing it here instead as more fitting
        auto A_block_zero_zero = Teuchos::rcp(
            new Matrix_Type(this->domainVec_.at(0)->getMapUnique(), this->domainVec_.at(0)->getApproxEntriesPerRow()));

        A->addBlock(A_block_zero_zero, 0, 0);
    } else {
        // Or same for the residual vector
        resVec_u = Teuchos::rcp(new MultiVector_Type(map, 1));
        resVecRep = Teuchos::rcp(new BlockMultiVector_Type(1));
        resVecRep->addBlock(resVec_u, 0);
    }
    // Call assembly routines on each element
    for (UN T = 0; T < assemblyFEElements_.size(); T++) {
        vec_dbl_Type solution(0);

        // Update solution on the element
        solution_u = getSolution(elements->getElement(T).getVectorNodeList(),
                                 u_rep, dofs);
        solution.insert(solution.end(), solution_u.begin(), solution_u.end());
        assemblyFEElements_[T]->updateSolution(solution);

        if (assembleMode == "Jacobian") {
            SmallMatrixPtr_Type elementMatrix;
            assemblyFEElements_[T]->assembleJacobian();
            elementMatrix = assemblyFEElements_[T]->getJacobian();
            // Insert (additive) the local element Jacobian into the global
            // matrix
            assemblyFEElements_[T]
                ->advanceNewtonStep(); // n genereal non linear solver step
            addFeBlock(A, elementMatrix, elements->getElement(T), map, 0, 0,
                       problemDisk);
        }

        if (assembleMode == "Rhs") {
            assemblyFEElements_[T]->assembleRHS();
            rhsVec = assemblyFEElements_[T]->getRHS();
            // Name RHS comes from solving linear systems
            // For nonlinear systems RHS synonymous to residual
            // Insert (additive) the updated residual into the global residual
            // vector
            addFeBlockMv(resVecRep, rhsVec, elements->getElement(T), dofs);
        }
    }
    if (callFillComplete && assembleMode != "Rhs") {
        // Signal that editing A has finished. This causes the entries of A to
        // be redistributed across the MPI ranks
        A->getBlock(0, 0)->fillComplete(domainVec_.at(0)->getMapUnique(),
                                        domainVec_.at(0)->getMapUnique());
    }
    if (assembleMode == "Rhs") {
        // Export from overlapping residual to unique residual
        MultiVectorPtr_Type resVecUnique = Teuchos::rcp(
            new MultiVector_Type(domainVec_.at(0)->getMapUnique(), 1));
        resVecUnique->putScalar(0.);
        resVecUnique->exportFromVector(resVec_u, true, "Add");
        resVec->addBlock(resVecUnique, 0);
    }
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyElasticityJacobianAndStressAceFEM(int dim,
                                                                std::string FEType,
                                                                MatrixPtr_Type &A,
                                                                MultiVectorPtr_Type &f,
                                                                MultiVectorPtr_Type u,
                                                                ParameterListPtr_Type pList,
                                                                double C,
                                                                bool callFillComplete){
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::runtime_error, "Not implemented for P0");
    UN FEloc = checkFE(dim,FEType);
    
    
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,FEType,Grad,Grad);
    
    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    
    Teuchos::BLAS<int,SC> teuchosBLAS;
    
    int nmbQuadPoints = dPhi->size();
    int nmbScalarDPhi = dPhi->at(0).size();
    int nmbAllDPhi = nmbScalarDPhi * dim;
    int nmbAllDPhiAllQaud = nmbQuadPoints * nmbAllDPhi;
    int sizeLocStiff = dim*dim;
    Teuchos::Array<SmallMatrix<SC> > dPhiMat( nmbAllDPhiAllQaud, SmallMatrix<SC>(dim) );
    
    this->buildFullDPhi( dPhi, dPhiMat ); //builds matrix from gradient of scalar phi
    
    std::string material_model = pList->sublist("Parameter").get("Material model","Neo-Hooke");
    
    double poissonRatio = pList->sublist("Parameter").get("Poisson Ratio",0.4);
    double mue = pList->sublist("Parameter").get("Mu",2.0e+6);
    double mue1 = pList->sublist("Parameter").get("Mu1",2.0e+6);
    double mue2 = pList->sublist("Parameter").get("Mu2",2.0e+6);
    // Berechne daraus nun E (Youngsches Modul) und die erste Lamé-Konstante \lambda
    double E = pList->sublist("Parameter").get("E",3.0e+6); // For StVK mue_*2.*(1. + poissonRatio_);
    double E1 = pList->sublist("Parameter").get("E1",3.0e+6); // For StVK mue_*2.*(1. + poissonRatio_);
    double E2 = pList->sublist("Parameter").get("E2",3.0e+6); // For StVK mue_*2.*(1. + poissonRatio_);
    
    if (material_model=="Saint Venant-Kirchhoff") {
        E = mue*2.*(1. + poissonRatio);
        E1 = mue1*2.*(1. + poissonRatio);
        E2 = mue2*2.*(1. + poissonRatio);
    }
    
    // For StVK (poissonRatio*E)/((1 + poissonRatio)*(1 - 2*poissonRatio));
    double lambda = (poissonRatio*E)/((1 + poissonRatio)*(1 - 2*poissonRatio));
    double lambda1 = (poissonRatio*E1)/((1 + poissonRatio)*(1 - 2*poissonRatio));
    double lambda2 = (poissonRatio*E2)/((1 + poissonRatio)*(1 - 2*poissonRatio));
    
    if (dim == 2){
        double* v;
        if(!material_model.compare("Saint Venant-Kirchhoff"))
            v = new double[154];
        else
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only Saint Venant-Kirchhoff in 2D.");
        
        double** Pmat = new double*[2];
        for (int i=0; i<2; i++)
            Pmat[i] = new double[2];
        
        double** F = new double*[2];
        for (int i=0; i<2; i++)
            F[i] = new double[2];
        
        double**** Amat = new double***[2];
        for (int i=0; i<2; i++){
            Amat[i] = new double**[2];
            for (int j=0; j<2; j++) {
                Amat[i][j] = new double*[2];
                for (int k=0; k<2; k++)
                    Amat[i][j][k] = new double[2];
            }
        }
        
        Teuchos::ArrayRCP< const SC > uArray = u->getData(0);
        
        Teuchos::ArrayRCP<SC> fValues = f->getDataNonConst(0);
        
        Teuchos::Array<int> indices(2);
        for (int T=0; T<elements->numberElements(); T++) {
            
            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);
            
            Teuchos::Array<SmallMatrix<SC> > all_dPhiMat_Binv( dPhiMat.size(), SmallMatrix<SC>() );
            
            for (int i=0; i<all_dPhiMat_Binv.size(); i++) {
                SmallMatrix<SC> res = dPhiMat[i] * Binv;
                all_dPhiMat_Binv[i] = res;
            }
            
            SmallMatrix<SC> locStiffMat( nmbAllDPhi, 0. );
            std::vector<SC> locStresses( nmbAllDPhi, 0. );
            int elementFlag = 0;
            for (int p=0; p<nmbQuadPoints; p++){
                
                SmallMatrix<SC> Fmat( dim, 0. );
                SmallMatrix<SC> tmpForScaling( dim, 0. );
                Fmat[0][0] = 1.; Fmat[1][1] = 1.;
                
                for (int i=0; i<nmbScalarDPhi; i++) {
                    indices.at(0) = dim * elements->getElement(T).getNode(i);
                    indices.at(1) = dim * elements->getElement(T).getNode(i) + 1;
                    
                    for (int j=0; j<dim; j++) {
                        tmpForScaling = all_dPhiMat_Binv[ p * nmbAllDPhi + dim * i + j ]; //we should not copy here
                        SC v = uArray[indices.at(j)];
                        tmpForScaling.scale( v );
                        Fmat += tmpForScaling;
                    }
                }
                
                for (int i=0; i<Fmat.size(); i++) {
                    for (int j=0; j<Fmat.size(); j++) {
                        F[i][j] = Fmat[i][j]; //fix so we dont need to copy.
                    }
                }
                                
                elementFlag = elements->getElement(T).getFlag();
                if (elementFlag == 1){
                    lambda = lambda1;
                    mue = mue1;
                    E = E1;
                }
                else if (elementFlag == 2){
                    lambda = lambda2;
                    mue = mue2;
                    E = E2;
                }
                
                if ( !material_model.compare("Saint Venant-Kirchhoff") )
                    stvk2d(v, &lambda, &mue, F, Pmat, Amat);
                
                SmallMatrix<SC> Aloc(dim*dim);
                for (int i=0; i<2; i++) {
                    for (int j=0; j<2; j++) {
                        for (int k=0; k<2; k++) {
                            for (int l=0; l<2; l++) {
                                Aloc[ 2 * i + j ][ 2 * k + l ] = Amat[i][j][k][l];
                            }
                        }
                    }
                }
                
                double* aceFEMFunc = new double[ sizeLocStiff * sizeLocStiff ];
                double* allDPhiBlas = new double[ nmbAllDPhi * sizeLocStiff ];
                
                //jacobian
                double* resTmp = new double[ nmbAllDPhi * sizeLocStiff ];
                // all_dPhiMat_Binv: quadpoints -> basisfunction vector field
                fillMatrixArray(Aloc, aceFEMFunc, "cols"); //blas uses column-major
                
                int offset = p * nmbAllDPhi;
                int offsetInArray = 0;
                for (int i=0; i<nmbAllDPhi; i++) {
                    fillMatrixArray( all_dPhiMat_Binv[ offset + i ], allDPhiBlas, "rows",offsetInArray );
                    offsetInArray += sizeLocStiff;
                }
                
                teuchosBLAS.GEMM (Teuchos::NO_TRANS, Teuchos::NO_TRANS, sizeLocStiff, nmbAllDPhi, sizeLocStiff, 1., aceFEMFunc, sizeLocStiff, allDPhiBlas, sizeLocStiff, 0., resTmp, sizeLocStiff);
                
                
                double* locStiffMatBlas = new double[ nmbAllDPhi * nmbAllDPhi ];
                
                teuchosBLAS.GEMM (Teuchos::TRANS, Teuchos::NO_TRANS, nmbAllDPhi, nmbAllDPhi, sizeLocStiff, 1., allDPhiBlas, sizeLocStiff/*lda of A not trans(A)! Otherwise result is wrong*/, resTmp, sizeLocStiff, 0., locStiffMatBlas, nmbAllDPhi);
                
                for (int i=0; i<nmbAllDPhi; i++) {
                    for (int j=0; j<nmbAllDPhi; j++)
                        locStiffMat[i][j] += weights->at(p) * locStiffMatBlas[ j * nmbAllDPhi + i ];
                }
                
                delete [] resTmp;
                delete [] locStiffMatBlas;
                
                
                //stress
                double* fArray = new double[ sizeLocStiff ];
                for (int i=0; i<dim; i++) {
                    for (int j=0; j<dim; j++) {
                        fArray[i * dim + j] = Pmat[i][j]; //is this correct?
                    }
                }
                
                double* res = new double[ nmbAllDPhi ];
                teuchosBLAS.GEMV(Teuchos::TRANS, sizeLocStiff, nmbAllDPhi, 1., allDPhiBlas, sizeLocStiff, fArray, 1, 0., res, 1);
                for (int i=0; i<locStresses.size(); i++) {
                    locStresses[i] += weights->at(p) * res[i];
                }
                
                delete [] res;
                delete [] aceFEMFunc;
                delete [] allDPhiBlas;
                delete [] fArray;
            }
            
            for (int i=0; i<nmbScalarDPhi; i++) {
                for (int d1=0; d1<dim; d1++) {
                    
                    LO rowLO = dim * elements->getElement(T).getNode(i) + d1;
                    SC v = absDetB * locStresses[ dim * i + d1 ];
                    fValues[rowLO] += v;
                    
                    Teuchos::Array<SC> value( nmbAllDPhi, 0. );
                    Teuchos::Array<GO> indices( nmbAllDPhi, 0 );
                    LO counter = 0;
                    for (UN j=0; j < nmbScalarDPhi; j++){
                        for (UN d2=0; d2<dim; d2++) {
                            indices[counter] = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d2 );
                            value[counter] = absDetB * locStiffMat[dim*i+d1][dim*j+d2];
                            counter++;
                        }
                    }
                    GO row = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d1 );
                    A->insertGlobalValues( row, indices(), value() );
                }
            }
        }
        
        delete [] v;
        for (int i=0; i<2; i++)
            delete [] Pmat[i];
        delete [] Pmat;
        for (int i=0; i<2; i++)
            delete [] F[i];
        delete [] F;
        
        for (int i=0; i<2; i++){
            for (int j=0; j<2; j++) {
                for (int k=0; k<2; k++)
                    delete [] Amat[i][j][k];
                delete [] Amat[i][j];
            }
            delete [] Amat[i];
        }
        delete [] Amat;
        
        
    }
    else if (dim == 3) {
        double* v;
        if (!material_model.compare("Neo-Hooke"))
            v = new double[466];
        else if(!material_model.compare("Mooney-Rivlin"))
            v = new double[476];
        else if(!material_model.compare("Saint Venant-Kirchhoff"))
            v = new double[279];
        else{
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only Neo-Hooke, Mooney-Rivlin and Saint Venant-Kirchhoff.");
        }
        
        double** Pmat = new double*[3];
        for (int i=0; i<3; i++)
            Pmat[i] = new double[3];
        
        double** F = new double*[3];
        for (int i=0; i<3; i++)
            F[i] = new double[3];
        
        double**** Amat = new double***[3];
        for (int i=0; i<3; i++){
            Amat[i] = new double**[3];
            for (int j=0; j<3; j++) {
                Amat[i][j] = new double*[3];
                for (int k=0; k<3; k++)
                    Amat[i][j][k] = new double[3];
            }
        }
        
        Teuchos::ArrayRCP< const SC > uArray = u->getData(0);

        Teuchos::ArrayRCP<SC> fValues = f->getDataNonConst(0);

        Teuchos::Array<int> indices(3);
        for (int T=0; T<elements->numberElements(); T++) {
            
            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);
            
            Teuchos::Array<SmallMatrix<SC> > all_dPhiMat_Binv( dPhiMat.size(), SmallMatrix<SC>() );
            
            for (int i=0; i<all_dPhiMat_Binv.size(); i++) {
                SmallMatrix<SC> res = dPhiMat[i] * Binv;
                all_dPhiMat_Binv[i] = res;
            }
            
            SmallMatrix<SC> locStiffMat( nmbAllDPhi, 0. );
            std::vector<SC> locStresses( nmbAllDPhi, 0. );
            int elementFlag = 0;
            for (int p=0; p<nmbQuadPoints; p++){
                
                SmallMatrix<SC> Fmat( dim, 0. );
                SmallMatrix<SC> tmpForScaling( dim, 0. );
                Fmat[0][0] = 1.; Fmat[1][1] = 1.; Fmat[2][2] = 1.;
                
                for (int i=0; i<nmbScalarDPhi; i++) {
                    indices.at(0) = dim * elements->getElement(T).getNode(i);
                    indices.at(1) = dim * elements->getElement(T).getNode(i) + 1;
                    indices.at(2) = dim * elements->getElement(T).getNode(i) + 2;
                    
                    for (int j=0; j<dim; j++) {
                        tmpForScaling = all_dPhiMat_Binv[ p * nmbAllDPhi + dim * i + j ]; //we should not copy here
                        SC v = uArray[indices.at(j)];
                        tmpForScaling.scale( v );
                        Fmat += tmpForScaling;
                    }
                }
                
                for (int i=0; i<Fmat.size(); i++) {
                    for (int j=0; j<Fmat.size(); j++) {
                        F[i][j] = Fmat[i][j]; //fix so we dont need to copy.
                    }
                }
                
                elementFlag = elements->getElement(T).getFlag();
                if (elementFlag == 1){
                    lambda = lambda1;
                    mue = mue1;
                    E = E1;
                }
                else if (elementFlag == 2){
                    lambda = lambda2;
                    mue = mue2;
                    E = E2;
                }
                
                if ( !material_model.compare("Neo-Hooke") )
                    nh3d(v, &E, &poissonRatio, F, Pmat, Amat);
                else if ( !material_model.compare("Mooney-Rivlin") )
                    mr3d(v, &E, &poissonRatio, &C, F, Pmat, Amat);
                else if ( !material_model.compare("Saint Venant-Kirchhoff") )
                    stvk3d(v, &lambda, &mue, F, Pmat, Amat);
                                
                SmallMatrix<SC> Aloc(dim*dim);
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        for (int k=0; k<3; k++) {
                            for (int l=0; l<3; l++) {
                                Aloc[ 3 * i + j ][ 3 * k + l ] = Amat[i][j][k][l];
                            }
                        }
                    }
                }
                
                double* aceFEMFunc = new double[ sizeLocStiff * sizeLocStiff ];
                double* allDPhiBlas = new double[ nmbAllDPhi * sizeLocStiff ];
                
                //jacobian
                double* resTmp = new double[ nmbAllDPhi * sizeLocStiff ];
                // all_dPhiMat_Binv: quadpoints -> basisfunction vector field
                fillMatrixArray(Aloc, aceFEMFunc, "cols"); //blas uses column-major
                
                int offset = p * nmbAllDPhi;
                int offsetInArray = 0;
                for (int i=0; i<nmbAllDPhi; i++) {
                    fillMatrixArray( all_dPhiMat_Binv[ offset + i ], allDPhiBlas, "rows",offsetInArray );
                    offsetInArray += sizeLocStiff;
                }
                
                teuchosBLAS.GEMM (Teuchos::NO_TRANS, Teuchos::NO_TRANS, sizeLocStiff, nmbAllDPhi, sizeLocStiff, 1., aceFEMFunc, sizeLocStiff, allDPhiBlas, sizeLocStiff, 0., resTmp, sizeLocStiff);

                
                double* locStiffMatBlas = new double[ nmbAllDPhi * nmbAllDPhi ];
                
                teuchosBLAS.GEMM (Teuchos::TRANS, Teuchos::NO_TRANS, nmbAllDPhi, nmbAllDPhi, sizeLocStiff, 1., allDPhiBlas, sizeLocStiff/*lda of A not trans(A)! Otherwise result is wrong*/, resTmp, sizeLocStiff, 0., locStiffMatBlas, nmbAllDPhi);
                
                for (int i=0; i<nmbAllDPhi; i++) {
                    for (int j=0; j<nmbAllDPhi; j++)
                        locStiffMat[i][j] += weights->at(p) * locStiffMatBlas[ j * nmbAllDPhi + i ];
                }
                
                delete [] resTmp;
                delete [] locStiffMatBlas;
                
                
                //stress
                double* fArray = new double[ sizeLocStiff ];
                for (int i=0; i<dim; i++) {
                    for (int j=0; j<dim; j++) {
                        fArray[i * dim + j] = Pmat[i][j]; //is this correct?
                    }
                }
                
                double* res = new double[ nmbAllDPhi ];
                teuchosBLAS.GEMV(Teuchos::TRANS, sizeLocStiff, nmbAllDPhi, 1., allDPhiBlas, sizeLocStiff, fArray, 1, 0., res, 1);
                for (int i=0; i<locStresses.size(); i++) {
                    locStresses[i] += weights->at(p) * res[i];
                }
                
                delete [] res;
                delete [] aceFEMFunc;
                delete [] allDPhiBlas;
                delete [] fArray;
            }
            
            for (int i=0; i<nmbScalarDPhi; i++) {
                for (int d1=0; d1<dim; d1++) {
                    
                    LO rowLO = dim * elements->getElement(T).getNode(i) + d1;
                    SC v = absDetB * locStresses[ dim * i + d1 ];
                    fValues[rowLO] += v;

                    Teuchos::Array<SC> value( nmbAllDPhi, 0. );
                    Teuchos::Array<GO> indices( nmbAllDPhi, 0 );
                    LO counter = 0;
                    for (UN j=0; j < nmbScalarDPhi; j++){
                        for (UN d2=0; d2<dim; d2++) {
                            indices[counter] = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d2 );
                            value[counter] = absDetB * locStiffMat[dim*i+d1][dim*j+d2];
                                                    
                            counter++;
                        }
                    }
                    GO row = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d1 );
                    A->insertGlobalValues( row, indices(), value() );
                }
            }
        }
        
        delete [] v;
        for (int i=0; i<3; i++)
            delete [] Pmat[i];
        delete [] Pmat;
        for (int i=0; i<3; i++)
            delete [] F[i];
        delete [] F;
        
        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++) {
                for (int k=0; k<3; k++)
                    delete [] Amat[i][j][k];
                delete [] Amat[i][j];
            }
            delete [] Amat[i];
        }
        delete [] Amat;
        
    }
    if (callFillComplete)
        A->fillComplete();
    
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyElasticityJacobianAceFEM(int dim,
                                                  std::string FEType,
                                                  MatrixPtr_Type &A,
                                                  MultiVectorPtr_Type u,
                                                  std::string material_model,
                                                  double E,
                                                  double nu,
                                                  double C,
                                                  bool callFillComplete){
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    UN FEloc = checkFE(dim,FEType);

    vec2D_int_ptr_Type elements = domainVec_.at(FEloc)->getElements();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree(dim,FEType,FEType,Grad,Grad);

    Helper::getDPhi(dPhi, weights, dim, FEType, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);

    Teuchos::BLAS<int, SC> teuchosBLAS;

    int nmbQuadPoints = dPhi->size();
    int nmbScalarDPhi = dPhi->at(0).size();
    int nmbAllDPhi = nmbScalarDPhi * dim;
    int nmbAllDPhiAllQaud = nmbQuadPoints * nmbAllDPhi;
    int sizeLocStiff = dim*dim;
    Teuchos::Array<SmallMatrix<SC> > dPhiMat( nmbAllDPhiAllQaud, SmallMatrix<SC>(dim) );

    this->buildFullDPhi( dPhi, dPhiMat ); //builds matrix from gradient of scalar phi

    if (dim == 2){
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only for 3D.");
    }
    else if (dim == 3) {

        double* v;
        if (!material_model.compare("Neo-Hooke"))
            v = new double[466];
        else if(!material_model.compare("Mooney-Rivlin"))
            v = new double[476];
        else{
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only Neo-Hooke and Mooney-Rivlin.");
        }


        double** Pmat = new double*[3];
        for (int i=0; i<3; i++)
            Pmat[i] = new double[3];

        double** F = new double*[3];
        for (int i=0; i<3; i++)
            F[i] = new double[3];

        double**** Amat = new double***[3];
        for (int i=0; i<3; i++){
            Amat[i] = new double**[3];
            for (int j=0; j<3; j++) {
                Amat[i][j] = new double*[3];
                for (int k=0; k<3; k++)
                    Amat[i][j][k] = new double[3];
            }
        }

        Teuchos::ArrayRCP< const SC > uArray = u->getData(0);

        Teuchos::Array<int> indices(3);
        for (int T=0; T<elements->size(); T++) {

            Helper::buildTransformation(elements->at(T), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            Teuchos::Array<SmallMatrix<SC> > all_dPhiMat_Binv( dPhiMat.size(), SmallMatrix<SC>() );

            for (int i=0; i<all_dPhiMat_Binv.size(); i++) {
                SmallMatrix<SC> res = dPhiMat[i] * Binv;
                all_dPhiMat_Binv[i] = res;
            }

            SmallMatrix<SC> locStiffMat( nmbAllDPhi, 0. );

            for (int p=0; p<nmbQuadPoints; p++){

                SmallMatrix<SC> Fmat( dim, 0. );
                SmallMatrix<SC> tmpForScaling( dim, 0. );
                Fmat[0][0] = 1.; Fmat[1][1] = 1.; Fmat[2][2] = 1.;

                for (int i=0; i<nmbScalarDPhi; i++) {
                    indices.at(0) = dim * elements->at(T).at(i);
                    indices.at(1) = dim * elements->at(T).at(i) + 1;
                    indices.at(2) = dim * elements->at(T).at(i) + 2;

                    for (int j=0; j<dim; j++) {
                        tmpForScaling = all_dPhiMat_Binv[ p * nmbAllDPhi + dim * i + j ]; //we should not copy here
                        SC v = uArray[indices.at(j)];
                        tmpForScaling.scale( v );
                        Fmat += tmpForScaling;
                    }
                }

                for (int i=0; i<Fmat.size(); i++) {
                    for (int j=0; j<Fmat.size(); j++) {
                        F[i][j] = Fmat[i][j]; //fix so we dont need to copy.
                    }
                }
                if ( !material_model.compare("Neo-Hooke") )
                    nh3d(v, &E, &nu, F, Pmat, Amat);
                else if ( !material_model.compare("Mooney-Rivlin") )
                    mr3d(v, &E, &nu, &C, F, Pmat, Amat);

                SmallMatrix<SC> Aloc(dim*dim);
                for (int i=0; i<3; i++) {
                    for (int j=0; j<3; j++) {
                        for (int k=0; k<3; k++) {
                            for (int l=0; l<3; l++) {
                                Aloc[ 3 * i + j ][ 3 * k + l ] = Amat[i][j][k][l];
                            }
                        }
                    }
                }

                double* aceFEMFunc = new double[ sizeLocStiff * sizeLocStiff ];
                double* allDPhiBlas = new double[ nmbAllDPhi * sizeLocStiff ];
                double* resTmp = new double[ nmbAllDPhi * sizeLocStiff ];
                // all_dPhiMat_Binv: quadpoints -> basisfunction vector field
                fillMatrixArray(Aloc, aceFEMFunc, "cols"); //blas uses column-major

                int offset = p * nmbAllDPhi;
                int offsetInArray = 0;
                for (int i=0; i<nmbAllDPhi; i++) {
                    fillMatrixArray( all_dPhiMat_Binv[ offset + i ], allDPhiBlas, "rows",offsetInArray );
                    offsetInArray += sizeLocStiff;
                }

                teuchosBLAS.GEMM (Teuchos::NO_TRANS, Teuchos::NO_TRANS, sizeLocStiff, nmbAllDPhi, sizeLocStiff, 1., aceFEMFunc, sizeLocStiff, allDPhiBlas, sizeLocStiff, 0., resTmp, sizeLocStiff);

                double* locStiffMatBlas = new double[ nmbAllDPhi * nmbAllDPhi ];

                teuchosBLAS.GEMM (Teuchos::TRANS, Teuchos::NO_TRANS, nmbAllDPhi, nmbAllDPhi, sizeLocStiff, 1., allDPhiBlas, sizeLocStiff/*lda of A not trans(A)! Otherwise result is wrong*/, resTmp, sizeLocStiff, 0., locStiffMatBlas, nmbAllDPhi);

                for (int i=0; i<nmbAllDPhi; i++) {
                    for (int j=0; j<nmbAllDPhi; j++)
                        locStiffMat[i][j] += weights->at(p) * locStiffMatBlas[ j * nmbAllDPhi + i ];
                }

                delete [] aceFEMFunc;
                delete [] allDPhiBlas;
                delete [] resTmp;
                delete [] locStiffMatBlas;

            }
            for (int i=0; i<nmbScalarDPhi; i++) {
                for (int d1=0; d1<dim; d1++) {
                    Teuchos::Array<SC> value( nmbAllDPhi, 0. );
                    Teuchos::Array<GO> indices( nmbAllDPhi, 0 );
                    LO counter = 0;
                    for (UN j=0; j < nmbScalarDPhi; j++){
                        for (UN d2=0; d2<dim; d2++) {
                            indices[counter] = GO ( dim * map->getGlobalElement( elements->at(T).at(j) ) + d2 );
                            value[counter] = absDetB * locStiffMat[dim*i+d1][dim*j+d2];
                            counter++;
                        }
                    }
                    GO row = GO ( dim * map->getGlobalElement( elements->at(T).at(i) ) + d1 );
                    A->insertGlobalValues( row, indices(), value() );
                }
            }
        }

        delete [] v;
        for (int i=0; i<3; i++)
            delete [] Pmat[i];
        delete [] Pmat;
        for (int i=0; i<3; i++)
            delete [] F[i];
        delete [] F;

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++) {
                for (int k=0; k<3; k++)
                    delete [] Amat[i][j][k];
                delete [] Amat[i][j];
            }
            delete [] Amat[i];
        }
        delete [] Amat;

    }
    if (callFillComplete)
        A->fillComplete();

}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyElasticityStressesAceFEM(int dim,
                                                             std::string FEType,
                                                             MultiVectorPtr_Type &f,
                                                             MultiVectorPtr_Type u,
                                                             std::string material_model,
                                                             double E,
                                                             double nu,
                                                             double C,
                                                             bool callFillComplete){
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    UN FEloc = checkFE(dim,FEType);

    vec2D_int_ptr_Type elements = domainVec_.at(FEloc)->getElements();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 	dPhi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree(dim,FEType,FEType,Grad,Grad);

    Helper::getDPhi(dPhi, weights, dim, FEType, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);

    Teuchos::BLAS<int, SC> teuchosBLAS;

    int nmbQuadPoints = dPhi->size();
    int nmbScalarDPhi = dPhi->at(0).size();
    int nmbAllDPhi = nmbScalarDPhi * dim;
    int nmbAllDPhiAllQaud = nmbQuadPoints * nmbAllDPhi;
    int sizeLocStiff = dim*dim;
    Teuchos::Array<SmallMatrix<SC> > dPhiMat( nmbAllDPhiAllQaud, SmallMatrix<SC>(dim) );

    this->buildFullDPhi( dPhi, dPhiMat ); //builds matrix from gradient of scalar phi

    if (dim == 2){
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only for 3D.");
    }
    else if (dim == 3) {

        double* v;
        if (!material_model.compare("Neo-Hooke"))
            v = new double[466];
        else if(!material_model.compare("Mooney-Rivlin"))
            v = new double[476];
        else{
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only Neo-Hooke and Mooney-Rivlin.");
        }


        double** Pmat = new double*[3];
        for (int i=0; i<3; i++)
            Pmat[i] = new double[3];

        double** F = new double*[3];
        for (int i=0; i<3; i++)
            F[i] = new double[3];

        double**** Amat = new double***[3];
        for (int i=0; i<3; i++){
            Amat[i] = new double**[3];
            for (int j=0; j<3; j++) {
                Amat[i][j] = new double*[3];
                for (int k=0; k<3; k++)
                    Amat[i][j][k] = new double[3];
            }
        }

        Teuchos::ArrayRCP< const SC > uArray = u->getData(0);
        
        Teuchos::ArrayRCP<SC> fValues = f->getDataNonConst(0);
        
        Teuchos::Array<int> indices(3);
        for (int T=0; T<elements->size(); T++) {

            Helper::buildTransformation(elements->at(T), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            Teuchos::Array<SmallMatrix<SC> > all_dPhiMat_Binv( dPhiMat.size(), SmallMatrix<SC>() );

            for (int i=0; i<all_dPhiMat_Binv.size(); i++) {
                SmallMatrix<SC> res = dPhiMat[i] * Binv;
                all_dPhiMat_Binv[i] = res;
            }
            std::vector<double> locStresses( nmbAllDPhi, 0. );

            for (int p=0; p<nmbQuadPoints; p++){

                SmallMatrix<SC> Fmat( dim, 0. );
                SmallMatrix<SC> tmpForScaling( dim, 0. );
                Fmat[0][0] = 1.; Fmat[1][1] = 1.; Fmat[2][2] = 1.;

                for (int i=0; i<nmbScalarDPhi; i++) {
                    indices.at(0) = dim * elements->at(T).at(i);
                    indices.at(1) = dim * elements->at(T).at(i) + 1;
                    indices.at(2) = dim * elements->at(T).at(i) + 2;

                    for (int j=0; j<dim; j++) {
                        tmpForScaling = all_dPhiMat_Binv[ p * nmbAllDPhi + dim * i + j ]; //we should not copy here
                        SC v = uArray[indices.at(j)];
                        tmpForScaling.scale( v );
                        Fmat += tmpForScaling;
                    }
                }

                for (int i=0; i<Fmat.size(); i++) {
                    for (int j=0; j<Fmat.size(); j++) {
                        F[i][j] = Fmat[i][j]; //fix so we dont need to copy.
                    }
                }
                if ( !material_model.compare("Neo-Hooke") )
                    nh3d(v, &E, &nu, F, Pmat, Amat);
                else if ( !material_model.compare("Mooney-Rivlin") )
                    mr3d(v, &E, &nu, &C, F, Pmat, Amat);

                double* aceFEMFunc = new double[ sizeLocStiff * sizeLocStiff ];
                double* allDPhiBlas = new double[ nmbAllDPhi * sizeLocStiff ];

                int offset = p * nmbAllDPhi;
                int offsetInArray = 0;
                for (int i=0; i<nmbAllDPhi; i++) {
                    fillMatrixArray( all_dPhiMat_Binv[ offset + i ], allDPhiBlas, "rows",offsetInArray );
                    offsetInArray += sizeLocStiff;
                }


                double* fArray = new double[ sizeLocStiff ];
                for (int i=0; i<dim; i++) {
                    for (int j=0; j<dim; j++) {
                        fArray[i * dim + j] = Pmat[i][j]; //is this correct?
                    }
                }

                double* res = new double[ nmbAllDPhi ];
                teuchosBLAS.GEMV(Teuchos::TRANS, sizeLocStiff, nmbAllDPhi, 1., allDPhiBlas, sizeLocStiff, fArray, 1, 0., res, 1);
                for (int i=0; i<locStresses.size(); i++) {
                    locStresses[i] += weights->at(p) * res[i];
                }

                delete [] aceFEMFunc;
                delete [] allDPhiBlas;
                delete [] fArray;
            }

            
            
            for (int i=0; i<nmbScalarDPhi; i++) {
                for (int d1=0; d1<dim; d1++) {
                    LO row = dim * elements->at(T).at(i) + d1;
                    SC v = absDetB * locStresses[ dim * i + d1 ];
                    fValues[row] = v;
                }
            }
        }


        delete [] v;
        for (int i=0; i<3; i++)
            delete [] Pmat[i];
        delete [] Pmat;
        for (int i=0; i<3; i++)
            delete [] F[i];
        delete [] F;

        for (int i=0; i<3; i++){
            for (int j=0; j<3; j++) {
                for (int k=0; k<3; k++)
                    delete [] Amat[i][j][k];
                delete [] Amat[i][j];
            }
            delete [] Amat[i];
        }
        delete [] Amat;

    }
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyAdvectionVecField(int dim,
                                                  std::string FEType,
                                                  MatrixPtr_Type &A,
                                                  MultiVectorPtr_Type u,
                                                  bool callFillComplete){

    TEUCHOS_TEST_FOR_EXCEPTION( u->getNumVectors()>1, std::logic_error, "Implement for numberMV > 1 ." );
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    
    UN FEloc = checkFE(dim,FEType);
    
    if (saveAssembly_) {
//        ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsNew();
//        vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
//        MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
//
//        vec3D_dbl_ptr_Type     dPhi;
//        vec2D_dbl_ptr_Type     phi;
//        vec_dbl_ptr_Type    weights = Teuchos::rcp(new vec_dbl_Type(0));
//
//        UN extraDeg = determineDegree( dim, FEType, Std); //Elementwise assembly of grad u
//
//        UN deg = determineDegree( dim, FEType, FEType, Grad, Std, extraDeg);
//
//        getDPhi(dPhi, weights, dim, FEType, deg);
//        getPhi(phi, weights, dim, FEType, deg);
//        GO glob_i, glob_j;
//        vec_dbl_Type v_i(dim);
//        vec_dbl_Type v_j(dim);
//
//        vec2D_dbl_Type uLoc( dim, vec_dbl_Type( weights->size() , -1. ) );
//        Teuchos::ArrayRCP< const SC > uArray = u->getData(0);
//
//        for (UN T=0; T<elements->numberElements(); T++) {
//
//            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
//            applyBTinv( dPhi, dPhiTrans, elements->getBTinv( T ) );
//
//            for (int w=0; w<phi->size(); w++){ //quads points
//                for (int d=0; d<dim; d++) {
//                    uLoc[d][w] = 0.;
//                    for (int i=0; i < phi->at(0).size(); i++) {
//                        LO index = dim * elements->getElement( T ).getNode( i ) + d;
//                        uLoc[d][w] += uArray[index] * phi->at(w).at(i);
//                    }
//                }
//
//            }
//
//            for (UN i=0; i < phi->at(0).size(); i++) {
//                Teuchos::Array<SC> value( dPhiTrans[0].size(), 0. );
//                Teuchos::Array<GO> indices( dPhiTrans[0].size(), 0 );
//                for (UN j=0; j < value.size(); j++) {
//                    for (UN w=0; w<dPhiTrans.size(); w++) {
//                        for (UN d=0; d<dim; d++)
//                            value[j] += weights->at(w) * uLoc[d][w] * (*phi)[w][i] * dPhiTrans[w][j][d];
//                    }
//                    value[j] *= elements->getDetBTinv(T);
//                    if (setZeros_ && std::fabs(value[j]) < myeps_) {
//                        value[j] = 0.;
//                    }
//
//                    GO row = GO ( dim * map->getGlobalElement( elements->at(T).at(i) )  );
//                    GO glob_j = GO ( dim * map->getGlobalElement( elements->at(T).at(j) )  );
//                }
//                for (UN d=0; d<dim; d++) {
//                    for (UN j=0; j < indices.size(); j++)
//                        indices[j] = GO ( dim * map->getGlobalElement( elements->at(T).at(j) ) + d );
//
//                    GO row = GO ( dim * map->getGlobalElement( elements->at(T).at(i) ) + d );
//                    A->insertGlobalValues( row, indices(), value() );
//                }
//            }
//        }
    } else {
        ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

        vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

        MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

        vec3D_dbl_ptr_Type     dPhi;
        vec2D_dbl_ptr_Type     phi;
        vec_dbl_ptr_Type    weights = Teuchos::rcp(new vec_dbl_Type(0));

        UN extraDeg = Helper::determineDegree( dim, FEType, Std); //Elementwise assembly of grad u

        UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Std, extraDeg);

        Helper::getDPhi(dPhi, weights, dim, FEType, deg);
        Helper::getPhi(phi, weights, dim, FEType, deg);
        SC detB;
        SC absDetB;
        SmallMatrix<SC> B(dim);
        SmallMatrix<SC> Binv(dim);
        GO glob_i, glob_j;
        vec_dbl_Type v_i(dim);
        vec_dbl_Type v_j(dim);

        vec2D_dbl_Type uLoc( dim, vec_dbl_Type( weights->size() , -1. ) );
        Teuchos::ArrayRCP< const SC > uArray = u->getData(0);

        for (UN T=0; T<elements->numberElements(); T++) {

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv );

            for (int w=0; w<phi->size(); w++){ //quads points
                for (int d=0; d<dim; d++) {
                    uLoc[d][w] = 0.;
                    for (int i=0; i < phi->at(0).size(); i++) {
                        LO index = dim * elements->getElement(T).getNode(i) + d;
                        uLoc[d][w] += uArray[index] * phi->at(w).at(i);
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
                    if (setZeros_ && std::fabs(value[j]) < myeps_) {
                        value[j] = 0.;
                    }

                    GO row = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) )  );
                    GO glob_j = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) )  );
                }
                for (UN d=0; d<dim; d++) {
                    for (UN j=0; j < indices.size(); j++)
                        indices[j] = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d );

                    GO row = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d );
                    A->insertGlobalValues( row, indices(), value() );
                }
            }
        }
    }
    
    if (callFillComplete)
        A->fillComplete();
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyAdvectionInUVecField(int dim,
                                                  std::string FEType,
                                                  MatrixPtr_Type &A,
                                                  MultiVectorPtr_Type u,
                                                  bool callFillComplete){

    TEUCHOS_TEST_FOR_EXCEPTION( u->getNumVectors()>1, std::logic_error, "Implement for numberMV > 1 ." );
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    UN FEloc = checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN extraDeg = Helper::determineDegree( dim, FEType, Grad); //Elementwise assembly of u

    UN deg = Helper::determineDegree( dim, FEType, FEType, Std, Std, extraDeg);

    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    Helper::getPhi(phi, weights, dim, FEType, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    Teuchos::ArrayRCP< const SC > uArray = u->getData(0);

    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyBTinv( dPhi, dPhiTrans, Binv );

        std::vector<SmallMatrix<SC> > duLoc( weights->size(), SmallMatrix<SC>(dim) ); //for all quad points p_i each matrix is [u_x * grad Phi(p_i), u_y * grad Phi(p_i), u_z * grad Phi(p_i) (if 3D) ], duLoc[w] = [[phixx;phixy],[phiyx;phiyy]] (2D)

        for (int w=0; w<dPhiTrans.size(); w++){ //quads points
            for (int d1=0; d1<dim; d1++) {
                for (int i=0; i < dPhiTrans[0].size(); i++) {
                    LO index = dim * elements->getElement(T).getNode(i) + d1;
                    for (int d2=0; d2<dim; d2++)
                        duLoc[w][d2][d1] += uArray[index] * dPhiTrans[w][i][d2];
                }
            }
        }

        for (UN i=0; i < phi->at(0).size(); i++) {
            for (UN d1=0; d1<dim; d1++) {
                Teuchos::Array<SC> value( dim*phi->at(0).size(), 0. ); //These are value (W_ix,W_iy,W_iz)
                Teuchos::Array<GO> indices( dim*phi->at(0).size(), 0 );
                for (UN j=0; j < phi->at(0).size(); j++) {
                    for (UN d2=0; d2<dim; d2++){
                        for (UN w=0; w<phi->size(); w++) {
                            value[ dim * j + d2 ] += weights->at(w) * duLoc[w][d2][d1] * (*phi)[w][i] * (*phi)[w][j];
                        }
                        value[ dim * j + d2 ] *= absDetB;

                        if (setZeros_ && std::fabs(value[ dim * j + d2 ]) < myeps_) {
                            value[ dim * j + d2 ] = 0.;
                        }
                    }
                }
                for (UN j=0; j < phi->at(0).size(); j++){
                    for (UN d2=0; d2<dim; d2++){
                        indices[ dim * j + d2 ] = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(j) ) + d2 );
                    }
                }

                GO row = GO ( dim * map->getGlobalElement( elements->getElement(T).getNode(i) ) + d1 );
                A->insertGlobalValues( row, indices(), value() );
            }
        }
    }
    if (callFillComplete)
        A->fillComplete();
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyDivAndDivT( int dim,
                                            std::string FEType1,
                                            std::string FEType2,
                                            int degree,
                                            MatrixPtr_Type &Bmat,
                                            MatrixPtr_Type &BTmat,
                                            MapConstPtr_Type map1,
                                            MapConstPtr_Type map2,
                                            bool callFillComplete) {


    UN FEloc1 = checkFE(dim,FEType1);
    UN FEloc2 = checkFE(dim,FEType2);

    ElementsPtr_Type elements1 = domainVec_.at(FEloc1)->getElementsC();
    ElementsPtr_Type elements2 = domainVec_.at(FEloc2)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep1 = domainVec_.at(FEloc1)->getPointsRepeated();

    MapConstPtr_Type mapping1 = domainVec_.at(FEloc1)->getMapRepeated();
    MapConstPtr_Type mapping2;

    if (FEType2 == "P0")
        mapping2 = domainVec_.at(FEloc2)->getElementMap();
    else
        mapping2 = domainVec_.at(FEloc2)->getMapRepeated();

    vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree( dim, FEType1, FEType2, Grad, Std);

    Helper::getDPhi(dPhi, weights, dim, FEType1, deg);

    if (FEType2=="P1-disc-global")
        Helper::getPhiGlobal(phi, weights, dim, FEType2, deg);
    if (FEType2=="P1-disc" && FEType1=="Q2" )
        Helper::getPhi(phi, weights, dim, FEType2, deg, FEType1);
    else
        Helper::getPhi(phi, weights, dim, FEType2, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    for (UN T=0; T<elements1->numberElements(); T++) {

        Helper::buildTransformation(elements1->getElement(T).getVectorNodeList(), pointsRep1, B, FEType1);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);

        vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyBTinv( dPhi, dPhiTrans, Binv );

        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<Teuchos::Array<SC> >valueVec( dim, Teuchos::Array<SC>( dPhiTrans[0].size(), 0. ) );
            Teuchos::Array<GO> indices( dPhiTrans[0].size(), 0 );

            for (UN j=0; j < valueVec[0].size(); j++) {
                for (UN w=0; w<dPhiTrans.size(); w++) {
                    for (UN d=0; d<dim; d++)
                        valueVec[d][j] += weights->at(w) * phi->at(w)[i] * dPhiTrans[w][j][d];
                }
                for (UN d=0; d<dim; d++){
                    valueVec[d][j] *= absDetB;
                    if (setZeros_ && std::fabs(valueVec[d][j]) < myeps_) {
                        valueVec[d][j] = 0.;
                    }
                }
            }
            for (UN d=0; d<dim; d++) {
                for (UN j=0; j < indices.size(); j++)
                    indices[j] = GO ( dim * mapping1->getGlobalElement( elements1->getElement(T).getNode(j) ) + d );

                GO row;
                if (FEType2=="P0")
                    row = GO ( mapping2->getGlobalElement( T ) );
                else
                    row = GO ( mapping2->getGlobalElement( elements2->getElement(T).getNode(i) ) );
                Bmat->insertGlobalValues( row, indices(), valueVec[d]() );
            }
        }

        // We compute value twice, maybe we should change this
        for (UN i=0; i < dPhiTrans[0].size(); i++) {

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

    }
    if (callFillComplete) {
        Bmat->fillComplete( map1, map2 );
        BTmat->fillComplete( map2, map1 );
    }

}

    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyDivAndDivTFast( int dim,
                                             std::string FEType1,
                                             std::string FEType2,
                                             int degree,
                                             MatrixPtr_Type &Bmat,
                                             MatrixPtr_Type &BTmat,
                                             MapConstPtr_Type map1,
                                             MapConstPtr_Type map2,
                                             bool callFillComplete) {
    
    
    UN FEloc1 = checkFE(dim,FEType1);
    UN FEloc2 = checkFE(dim,FEType2);
    
    ElementsPtr_Type elements1 = domainVec_.at(FEloc1)->getElementsC();
    ElementsPtr_Type elements2 = domainVec_.at(FEloc2)->getElementsC();
    
    vec2D_dbl_ptr_Type pointsRep1 = domainVec_.at(FEloc1)->getPointsRepeated();
    
    MapConstPtr_Type mapping1 = domainVec_.at(FEloc1)->getMapRepeated();
    MapConstPtr_Type mapping2;
    
    if (FEType2 == "P0")
        mapping2 = domainVec_.at(FEloc2)->getElementMap();
    else
        mapping2 = domainVec_.at(FEloc2)->getMapRepeated();
    
    vec3D_dbl_ptr_Type 	dPhi;
    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree( dim, FEType1, FEType2, Grad, Std);
    
    Helper::getDPhi(dPhi, weights, dim, FEType1, deg);
    
    if (FEType2=="P1-disc-global")
        Helper::getPhiGlobal(phi, weights, dim, FEType2, deg);
    if (FEType2=="P1-disc" && FEType1=="Q2" )
        Helper::getPhi(phi, weights, dim, FEType2, deg, FEType1);
    else
        Helper::getPhi(phi, weights, dim, FEType2, deg);
    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);
    
    Teuchos::Array<GO> colIndex( 1, 0 );
    Teuchos::Array<GO> rowIndex( 1, 0 );
    Teuchos::Array<SC> value(1, 0.);

    for (UN T=0; T<elements1->numberElements(); T++) {
        
        Helper::buildTransformation(elements1->getElement(T).getVectorNodeList(), pointsRep1, B, FEType1);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);
        
        vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
        applyBTinv( dPhi, dPhiTrans, Binv );
        
        for (UN i=0; i < phi->at(0).size(); i++) {
            if (FEType2=="P0")
                rowIndex[0] = GO ( mapping2->getGlobalElement( T ) );
            else
                rowIndex[0] = GO ( mapping2->getGlobalElement( elements2->getElement(T).getNode(i) ) );

            for (UN j=0; j < dPhiTrans[0].size(); j++) {
                for (UN d=0; d<dim; d++){
                    value[0] = 0.;
                    for (UN w=0; w<dPhiTrans.size(); w++)
                        value[0] += weights->at(w) * phi->at(w)[i] * dPhiTrans[w][j][d];
                    value[0] *= absDetB;
                    colIndex[0] = GO ( dim * mapping1->getGlobalElement( elements1->getElement(T).getNode(j) ) + d );
                    Bmat->insertGlobalValues( rowIndex[0], colIndex(), value() );
                    BTmat->insertGlobalValues( colIndex[0], rowIndex(), value() );	
					
                }
            }
		}
	    	
    }
    if (callFillComplete) {
        Bmat->fillComplete( map1, map2 );
        BTmat->fillComplete( map2, map1 );
    }
    
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyBDStabilization(int dim,
                                              std::string FEType,
                                              MatrixPtr_Type &A,
                                              bool callFillComplete){
     
    TEUCHOS_TEST_FOR_EXCEPTION(FEType != "P1",std::logic_error, "Only implemented for P1. Q1 is equivalent but we need to adjust scaling for the reference element.");
    UN FEloc = checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec2D_dbl_ptr_Type 	phi;

    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    
    UN deg = Helper::determineDegree(dim,FEType,FEType,Std,Std);

    Helper::getPhi( phi, weights, dim, FEType, deg );

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    SC refElementSize;
    SC refElementScale;
    if (dim==2) {
        refElementSize = 0.5;
        refElementScale = 1./9.;
    }
    else if(dim==3){
        refElementSize = 1./6.;
        refElementScale = 1./16.;
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Only implemented for 2D and 3D.");

    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
        detB = B.computeDet( );
        absDetB = std::fabs(detB);

        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<SC> value( phi->at(0).size(), 0. );
            Teuchos::Array<GO> indices( phi->at(0).size(), 0 );
            for (UN j=0; j < value.size(); j++) {
                for (UN w=0; w<phi->size(); w++) {
                    value[j] += weights->at(w) * (*phi)[w][i] * (*phi)[w][j];
                }
                value[j] *= absDetB;
                value[j] -= refElementSize * absDetB * refElementScale;

                indices[j] = map->getGlobalElement( elements->getElement(T).getNode(j) );
            }

            GO row = map->getGlobalElement( elements->getElement(T).getNode(i) );
            A->insertGlobalValues( row, indices(), value() );
        }

    }

    if (callFillComplete)
        A->fillComplete();
}



template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLaplaceXDim(int dim,
                            std::string FEType,
                            MatrixPtr_Type &A,
                            CoeffFuncDbl_Type func,
                            double* parameters,
                            bool callFillComplete)
{
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    int FEloc = this->checkFE(dim,FEType);

    DomainConstPtr_Type domain = domainVec_.at(FEloc);
    ElementsPtr_Type elements = domain->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domain->getPointsRepeated();
    MapConstPtr_Type map = domain->getMapRepeated();

    vec3D_dbl_ptr_Type 			dPhi;
    vec_dbl_ptr_Type			weightsDPhi = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type			quadPts;

    // double val, value1_j, value2_j , value1_i, value2_i;

    UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Grad);

    Helper::getDPhi(dPhi, weightsDPhi, dim, FEType, deg);
    Helper::getQuadratureValues(dim, deg, quadPts, weightsDPhi, FEType);

    // SC = double, GO = long, UN = int
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;


    vec_dbl_ptr_Type dist = domain->getDistancesToInterface();
    if (dim == 2)
    {
        double val, value1_j, value2_j , value1_i, value2_i;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0);

        double distance1, distance2, distance3;
        vec_dbl_Type distance_mean(1); // Durchschnittliche Distanz des elements T
        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));

            distance1 = dist->at(elements->getElement(T).getNode(0));
            distance2 = dist->at(elements->getElement(T).getNode(1));
            distance3 = dist->at(elements->getElement(T).getNode(2));

            distance_mean.at(0) = (distance1 + distance2 + distance3)/3.0; // Mittelwert
            double funcvalue = func(&distance_mean.at(0),parameters);

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    val = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {

                        value1_j = dPhiTrans.at(k).at(j).at(0);
                        value2_j = dPhiTrans.at(k).at(j).at(1);

                        value1_i = dPhiTrans.at(k).at(i).at(0);
                        value2_i = dPhiTrans.at(k).at(i).at(1);

                        val = val + funcvalue * weightsDPhi->at(k) * ( value1_j*value1_i + value2_j*value2_i );
                    }
                    val = absDetB * val;
                    value[0] = val;

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;

                                        
                    A->insertGlobalValues(glob_i, indices(), value());
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i+1, indices(), value());
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
    else if(dim == 3)
    {
        double val, value1_j, value2_j ,value3_j, value1_i, value2_i ,value3_i;

        long long glob_i, glob_j;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0), p4(3,0.0);

        double distance1, distance2, distance3, distance4;
        vec_dbl_Type distance_mean(1); // Durchschnittliche Distanz des elements T
        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));
            p4 = pointsRep->at(elements->getElement(T).getNode(3));

            distance1 = dist->at(elements->getElement(T).getNode(0));
            distance2 = dist->at(elements->getElement(T).getNode(1));
            distance3 = dist->at(elements->getElement(T).getNode(2));
            distance4 = dist->at(elements->getElement(T).getNode(3));

            distance_mean.at(0) = (distance1 + distance2 + distance3 + distance4)/4.0; //Mittelwert
            double funcvalue = func(&distance_mean.at(0),parameters);

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    val = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {
                        value1_j = dPhiTrans.at(k).at(j).at(0);
                        value2_j = dPhiTrans.at(k).at(j).at(1);
                        value3_j = dPhiTrans.at(k).at(j).at(2);

                        value1_i = dPhiTrans.at(k).at(i).at(0);
                        value2_i = dPhiTrans.at(k).at(i).at(1);
                        value3_i = dPhiTrans.at(k).at(i).at(2);

                        val = val + funcvalue * weightsDPhi->at(k) * (value1_j*value1_i + value2_j*value2_i + value3_j*value3_i);
                    }
                    val = absDetB * val;
                    value[0] = val;

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value());
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i+1, indices(), value());
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i+2, indices(), value());

                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }

}



template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyStress(int dim,
                                     std::string FEType,
                                     MatrixPtr_Type &A,
                                     CoeffFunc_Type func,
                                     int* parameters,
                                     bool callFillComplete)
{
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    int FEloc = this->checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 			dPhi;
    vec_dbl_ptr_Type			weightsDPhi = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type			quadPts;

    // double value, value1_j, value2_j , value1_i, value2_i;

    UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Grad);
    Helper::getDPhi(dPhi, weightsDPhi, dim, FEType, deg);
    Helper::getQuadratureValues(dim, deg, quadPts, weightsDPhi,FEType);

    // SC = double, GO = long, UN = int
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;

    if (dim == 2)
    {
        double v11, v12, v21, v22, value1_j, value2_j , value1_i, value2_i;
        double e_11_j_1,e_12_j_1,e_21_j_1,e_22_j_1;
        double e_11_j_2,e_12_j_2,e_21_j_2,e_22_j_2;
        double e_11_i_1,e_12_i_1,e_21_i_1,e_22_i_1;
        double e_11_i_2,e_12_i_2,e_21_i_2,e_22_i_2;

        SmallMatrix<double> tmpRes1(dim);
        SmallMatrix<double> tmpRes2(dim);
        SmallMatrix<double> e1i(dim);
        SmallMatrix<double> e2i(dim);
        SmallMatrix<double> e1j(dim);
        SmallMatrix<double> e2j(dim);

        long long glob_i, glob_j;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0);

        vec_dbl_Type xy(2);
        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1)
            // Also \hat{grad_phi}.
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j=0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0;v12 = 0.0;v21 = 0.0;v22 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {
                        // Mappen der Gausspunkte (definiert auf T_ref) auf T (bzw. \Omega)
                        // xy = F(quadPts) = B*quadPts + b, mit b = p1 (affin lineare Transformation)
                        xy[0]=0.; xy[1]=0.;
                        for (int r=0; r<2; r++) {
                            xy[0] += B[0][r]*quadPts->at(k).at(r);
                            xy[1] += B[1][r]*quadPts->at(k).at(r);
                        }
                        xy[0] += p1[0];
                        xy[1] += p1[1];

                        value1_j = dPhiTrans.at(k).at(j).at(0);
                        value2_j = dPhiTrans.at(k).at(j).at(1);

                        value1_i = dPhiTrans.at(k).at(i).at(0);
                        value2_i = dPhiTrans.at(k).at(i).at(1);

                        tmpRes1[0][0] = value1_j;
                        tmpRes1[0][1] = value2_j;
                        tmpRes1[1][0] = 0.;
                        tmpRes1[1][1] = 0.;

                        tmpRes2[0][0] = value1_j;
                        tmpRes2[0][1] = 0.;
                        tmpRes2[1][0] = value2_j;
                        tmpRes2[1][1] = 0.;

                        tmpRes1.add(tmpRes2,e1j/*result*/);

                        e1i[0][0] = value1_i;
                        e1i[0][1] = value2_i;


                        tmpRes1[0][0] = 0.;
                        tmpRes1[0][1] = 0.;
                        tmpRes1[1][0] = value1_j;
                        tmpRes1[1][1] = value2_j;

                        tmpRes2[0][0] = 0.;
                        tmpRes2[0][1] = value1_j;
                        tmpRes2[1][0] = 0.;
                        tmpRes2[1][1] = value2_j;

                        tmpRes1.add(tmpRes2,e2j/*result*/);

                        e2i[1][0] = value1_i;
                        e2i[1][1] = value2_i;

                        double funcvalue = func(&xy.at(0),parameters);
                        v11 = v11 + funcvalue * weightsDPhi->at(k) * e1i.innerProduct(e1j);
                        v12 = v12 + funcvalue * weightsDPhi->at(k) * e1i.innerProduct(e2j);
                        v21 = v21 + funcvalue * weightsDPhi->at(k) * e2i.innerProduct(e1j);
                        v22 = v22 + funcvalue * weightsDPhi->at(k) * e2i.innerProduct(e2j);

                    }

                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;

                    value11[0] = v11;
                    value12[0] = v12;
                    value21[0] = v21;
                    value22[0] = v22;

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    
                    A->insertGlobalValues(glob_i, indices(), value11());
                    A->insertGlobalValues(glob_i+1, indices(), value21());
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value12());
                    A->insertGlobalValues(glob_i+1, indices(), value22());
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
    else if(dim == 3)
    {
        double v11, v12, v13, v21, v22, v23, v31, v32, v33, value1_j, value2_j, value3_j , value1_i, value2_i, value3_i;

        SmallMatrix<double> e1i(dim);
        SmallMatrix<double> e2i(dim);
        SmallMatrix<double> e3i(dim);
        SmallMatrix<double> e1j(dim);
        SmallMatrix<double> e2j(dim);
        SmallMatrix<double> e3j(dim);

        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0), p4(3,0.0);
        vec_dbl_Type xyz(3);

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));
            p4 = pointsRep->at(elements->getElement(T).getNode(3));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value13( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<SC> value23( 1, 0. );
                Teuchos::Array<SC> value31( 1, 0. );
                Teuchos::Array<SC> value32( 1, 0. );
                Teuchos::Array<SC> value33( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0;v12 = 0.0;v13 = 0.0;v21 = 0.0;v22 = 0.0;v23 = 0.0;v31 = 0.0;v32 = 0.0;v33 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {

                        xyz[0]=0.; xyz[1]=0.; xyz[2]=0.;
                        for (int r = 0; r < 3; r++)
                        {
                            xyz[0] += B[0][r]*quadPts->at(k).at(r);
                            xyz[1] += B[1][r]*quadPts->at(k).at(r);
                            xyz[2] += B[2][r]*quadPts->at(k).at(r);
                        }
                        xyz[0] += p1[0];
                        xyz[1] += p1[1];
                        xyz[2] += p1[2];



                        value1_j = dPhiTrans.at(k).at(j).at(0);
                        value2_j = dPhiTrans.at(k).at(j).at(1);
                        value3_j = dPhiTrans.at(k).at(j).at(2);


                        value1_i = dPhiTrans.at(k).at(i).at(0);
                        value2_i = dPhiTrans.at(k).at(i).at(1);
                        value3_i = dPhiTrans.at(k).at(i).at(2);


                        e1j[0][0] = 2.*value1_j;
                        e1j[0][1] = value2_j;
                        e1j[0][2] = value3_j;
                        e1j[1][0] = value2_j;
                        e1j[2][0] = value3_j;

                        e1i[0][0] = value1_i;
                        e1i[0][1] = value2_i;
                        e1i[0][2] = value3_i;


                        e2j[1][0] = value1_j;
                        e2j[1][1] = 2.*value2_j;
                        e2j[1][2] = value3_j;
                        e2j[0][1] = value1_j;
                        e2j[2][1] = value3_j;

                        e2i[1][0] = value1_i;
                        e2i[1][1] = value2_i;
                        e2i[1][2] = value3_i;


                        e3j[2][0] = value1_j;
                        e3j[2][1] = value2_j;
                        e3j[2][2] = 2.*value3_j;
                        e3j[0][2] = value1_j;
                        e3j[1][2] = value2_j;

                        e3i[2][0] = value1_i;
                        e3i[2][1] = value2_i;
                        e3i[2][2] = value3_i;

                        double funcvalue = func(&xyz.at(0),parameters);

                        v11 = v11 + funcvalue * weightsDPhi->at(k) * e1i.innerProduct(e1j);
                        v12 = v12 + funcvalue * weightsDPhi->at(k) * e1i.innerProduct(e2j);
                        v13 = v13 + funcvalue * weightsDPhi->at(k) * e1i.innerProduct(e3j);

                        v21 = v21 + funcvalue * weightsDPhi->at(k) * e2i.innerProduct(e1j);
                        v22 = v22 + funcvalue * weightsDPhi->at(k) * e2i.innerProduct(e2j);
                        v23 = v23 + funcvalue * weightsDPhi->at(k) * e2i.innerProduct(e3j);

                        v31 = v31 + funcvalue * weightsDPhi->at(k) * e3i.innerProduct(e1j);
                        v32 = v32 + funcvalue * weightsDPhi->at(k) * e3i.innerProduct(e2j);
                        v33 = v33 + funcvalue * weightsDPhi->at(k) * e3i.innerProduct(e3j);


                    }
                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v13 = absDetB * v13;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;
                    v23 = absDetB * v23;
                    v31 = absDetB * v31;
                    v32 = absDetB * v32;
                    v33 = absDetB * v33;

                    value11[0] = v11;
                    value12[0] = v12;
                    value13[0] = v13;
                    value21[0] = v21;
                    value22[0] = v22;
                    value23[0] = v23;
                    value31[0] = v31;
                    value32[0] = v32;
                    value33[0] = v33;

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value11());
                    A->insertGlobalValues(glob_i+1, indices(), value21());
                    A->insertGlobalValues(glob_i+2, indices(), value31());

                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value12());
                    A->insertGlobalValues(glob_i+1, indices(), value22());
                    A->insertGlobalValues(glob_i+2, indices(), value32());

                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value13());
                    A->insertGlobalValues(glob_i+1, indices(), value23());
                    A->insertGlobalValues(glob_i+2, indices(), value33());
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }

}



template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLinElasXDim(int dim,
                                          std::string FEType,
                                          MatrixPtr_Type &A,
                                          double lambda,
                                          double mu,
                                          bool callFillComplete)
{
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    int FEloc = this->checkFE(dim,FEType);

    // Hole Elemente und Knotenliste
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 			dPhi;
    vec_dbl_ptr_Type			weightsDPhi = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type			quadPts;

    UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Grad);

    // Hole die grad_phi, hier DPhi
    Helper::getDPhi(dPhi, weightsDPhi, dim, FEType, deg);
    Helper::getQuadratureValues(dim, deg, quadPts, weightsDPhi,FEType);

    // Definiere die Basisfunktion \phi_i bzw. \phi_j
    // vec_dbl_Type basisValues_i(dim,0.), basisValues_j(dim,0.);

    // SC = double, GO = long, UN = int
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;

    // Fuer Zwischenergebniss
    SC res;

    // Fuer die Berechnung der Spur
    double res_trace_i, res_trace_j;    
    
    if (dim == 2)
    {

        double v11, v12, v21, v22;
        // Setzte Vektoren der Groesse 2 und initialisiere mit 0.0 (double)
        vec_dbl_Type p1(2,0.0), p2(2,0.0), p3(2,0.0);

        // Matrizen der Groesse (2x2) in denen die einzelnen Epsilon-Tensoren berechnet werden.
        // Siehe unten fuer mehr.
        SmallMatrix<double> epsilonValuesMat1_i(dim), epsilonValuesMat2_i(dim),
        epsilonValuesMat1_j(dim), epsilonValuesMat2_j(dim);

        for (int T = 0; T < elements->numberElements(); T++)
        {
            // Hole die Eckknoten des Dreiecks
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));

            // Berechne die Transormationsmatrix B fuer das jeweilige Element (2D)
            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1).
            // Also \hat{grad_phi}.
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0; v12 = 0.0; v21 = 0.0; v22 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {
                        // In epsilonValuesMat1_i (2x2 Matrix) steht fuer die Ansatzfunktion i bzw. \phi_i
                        // der epsilonTensor fuer eine skalare Ansatzfunktion fuer die Richtung 1 (vgl. Mat1).
                        // Also in Mat1_i wird dann also phi_i = (phi_scalar_i, 0) gesetzt und davon \eps berechnet.

                        // Stelle \hat{grad_phi_i} = basisValues_i auf, also B^(-T)*grad_phi_i
                        // GradPhiOnRef( dPhi->at(k).at(i), b_T_inv, basisValues_i );

                        // \eps(v) = \eps(phi_i)
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat1_i, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat2_i, 1); // y-Richtung

                        // Siehe oben, nur fuer j
                        // GradPhiOnRef( DPhi->at(k).at(j), b_T_inv, basisValues_j );

                        // \eps(u) = \eps(phi_j)
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat1_j, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat2_j, 1); // y-Richtung

                        // Nun berechnen wir \eps(u):\eps(v) = \eps(phi_j):\eps(phi_i).
                        // Das Ergebniss steht in res.
                        // Berechne zudem noch die Spur der Epsilon-Tensoren tr(\eps(u)) (j) und tr(\eps(v)) (i)
                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat1_j, res); // x-x
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v11 = v11 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat2_j, res); // x-y
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v12 = v12 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat1_j, res); // y-x
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v21 = v21 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat2_j, res); // y-y
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v22 = v22 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);


                    }
                    // Noch mit der abs(det(B)) skalieren
                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;

                    value11[0] = v11;
                    value12[0] = v12;
                    value21[0] = v21;
                    value22[0] = v22;

                    // Hole die globale Zeile und Spalte in der die Eintraege hingeschrieben werden sollen
                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value11()); // x-x
                    A->insertGlobalValues(glob_i+1, indices(), value21()); // y-x
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value12()); // x-y
                    A->insertGlobalValues(glob_i+1, indices(), value22()); // y-y
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
    else if(dim == 3)
    {

        double v11, v12, v13, v21, v22, v23, v31, v32, v33;

        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0), p4(3,0.0);
        SmallMatrix<double> epsilonValuesMat1_i(dim), epsilonValuesMat2_i(dim), epsilonValuesMat3_i(dim),
        epsilonValuesMat1_j(dim), epsilonValuesMat2_j(dim), epsilonValuesMat3_j(dim);

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));
            p4 = pointsRep->at(elements->getElement(T).getNode(3));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value13( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<SC> value23( 1, 0. );
                Teuchos::Array<SC> value31( 1, 0. );
                Teuchos::Array<SC> value32( 1, 0. );
                Teuchos::Array<SC> value33( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0; v12 = 0.0; v13 = 0.0; v21 = 0.0; v22 = 0.0; v23 = 0.0; v31 = 0.0; v32 = 0.0; v33 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {

                        // GradPhiOnRef( DPhi->at(k).at(i), b_T_inv, basisValues_i );

                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat1_i, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat2_i, 1); // y-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat3_i, 2); // z-Richtung


                        // GradPhiOnRef( DPhi->at(k).at(j), b_T_inv, basisValues_j );

                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat1_j, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat2_j, 1); // y-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat3_j, 2); // z-Richtung

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat1_j, res); // x-x
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v11 = v11 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat2_j, res); // x-y
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v12 = v12 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat3_j, res); // x-z
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v13 = v13 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat1_j, res); // y-x
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v21 = v21 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat2_j, res); // y-y
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v22 = v22 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat3_j, res); // y-z
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v23 = v23 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat1_j, res); // z-x
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v31 = v31 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat2_j, res); // z-y
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v32 = v32 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat3_j, res); // z-z
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v33 = v33 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                    }
                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v13 = absDetB * v13;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;
                    v23 = absDetB * v23;
                    v31 = absDetB * v31;
                    v32 = absDetB * v32;
                    v33 = absDetB * v33;

                    value11[0] = v11;
                    value12[0] = v12;
                    value13[0] = v13;
                    value21[0] = v21;
                    value22[0] = v22;
                    value23[0] = v23;
                    value31[0] = v31;
                    value32[0] = v32;
                    value33[0] = v33;

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value11()); // x-x
                    A->insertGlobalValues(glob_i+1, indices(), value21()); // y-x
                    A->insertGlobalValues(glob_i+2, indices(), value31()); // z-x
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value12()); // x-y
                    A->insertGlobalValues(glob_i+1, indices(), value22()); // y-y
                    A->insertGlobalValues(glob_i+2, indices(), value32()); // z-y
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value13()); // x-z
                    A->insertGlobalValues(glob_i+1, indices(), value23()); // y-z
                    A->insertGlobalValues(glob_i+2, indices(), value33()); // z-z
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
}

// Determine the change of emodule depending on concentration
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::determineEMod(std::string FEType, MultiVectorPtr_Type solution,MultiVectorPtr_Type &eModVec, DomainConstPtr_Type domain, 	ParameterListPtr_Type params){


    ElementsPtr_Type elements = domain->getElementsC();

    int dim = domain->getDimension();
    vec2D_dbl_ptr_Type pointsRep = domain->getPointsRepeated();

    //MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec2D_dbl_ptr_Type 	phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));

    UN deg = Helper::determineDegree(dim,FEType,FEType,Std,Std);

    Helper::getPhi( phi, weights, dim, FEType, deg );

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);

    
    Teuchos::ArrayRCP< const SC > uArray = solution->getData(0);
    Teuchos::ArrayRCP< SC > eModVecA = eModVec->getDataNonConst(0);

    double E0 = params->sublist("Parameter Solid").get("E",3.0e+6);
    double E1 = params->sublist("Parameter Solid").get("E1",3.0e+5);
    double c1 = params->sublist("Parameter Solid").get("c1",1.0);
    double eModMax =E1;
    double eModMin = E0;

    int nodesElement = elements->getElement(0).getVectorNodeList().size();
    for (UN T=0; T<elements->numberElements(); T++) {
   
        /*buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeInverse(Binv);
        absDetB = std::fabs(detB);*/
        
        double uLoc = 0.;

       //for (int w=0; w<phi->size(); w++){ //quads points
       //     for (int i=0; i < phi->at(0).size(); i++) {
            for(int i=0; i< nodesElement;i++){
                LO index = elements->getElement(T).getNode(i) ;
                uLoc += 1./nodesElement*uArray[index];
            }
       //     } 
       // }
        //uLoc = uLoc*absDetB;           

        eModVecA[T] = E0-(E0-E1)*(uLoc/(uLoc+c1));
        if(eModVecA[T] > eModMax )
            eModMax = eModVecA[T];
        if(eModVecA[T] < eModMin)
            eModMin = eModVecA[T];
    }
    Teuchos::reduceAll<int, double> (*(domain->getComm()), Teuchos::REDUCE_MIN, eModMin, Teuchos::outArg (eModMin));
    Teuchos::reduceAll<int, double> (*(domain->getComm()), Teuchos::REDUCE_MAX, eModMax, Teuchos::outArg (eModMax));

    if(domain->getComm()->getRank()==0)
        cout << " #################  eMOD Min: " << eModMin << " \t eModMax: " << eModMax<< " ############# " <<endl;


}


/// \brief Same as assemblyLinElasXDim except for changing E Module Value
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyLinElasXDimE(int dim,
                                          std::string FEType,
                                          MatrixPtr_Type &A,
                                          MultiVectorPtr_Type eModVec,
                                          double nu,
                                          bool callFillComplete)
{
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    int FEloc = this->checkFE(dim,FEType);
    // Hole Elemente und Knotenliste
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    vec3D_dbl_ptr_Type 			dPhi;
    vec_dbl_ptr_Type			weightsDPhi = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type			quadPts;

    UN deg = Helper::determineDegree( dim, FEType, FEType, Grad, Grad);

    // Hole die grad_phi, hier DPhi
    Helper::getDPhi(dPhi, weightsDPhi, dim, FEType, deg);
    Helper::getQuadratureValues(dim, deg, quadPts, weightsDPhi,FEType);

    // Definiere die Basisfunktion \phi_i bzw. \phi_j
    // vec_dbl_Type basisValues_i(dim,0.), basisValues_j(dim,0.);

    // SC = double, GO = long, UN = int
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;

    // Fuer Zwischenergebniss
    SC res;

    // Fuer die Berechnung der Spur
    double res_trace_i, res_trace_j;    
    
    Teuchos::ArrayRCP< const SC > E = eModVec->getData(0);
    double lambda;
    double mu ;

  
    if (dim == 2)
    {

        double v11, v12, v21, v22;
        // Setzte Vektoren der Groesse 2 und initialisiere mit 0.0 (double)
        vec_dbl_Type p1(2,0.0), p2(2,0.0), p3(2,0.0);

        // Matrizen der Groesse (2x2) in denen die einzelnen Epsilon-Tensoren berechnet werden.
        // Siehe unten fuer mehr.
        SmallMatrix<double> epsilonValuesMat1_i(dim), epsilonValuesMat2_i(dim),
        epsilonValuesMat1_j(dim), epsilonValuesMat2_j(dim);

        for (int T = 0; T < elements->numberElements(); T++)
        {
            /// \lambda = E(T)* \nu / ( (1+\nu))*(1-2*nu))
            lambda = E[T]* nu / ((1.+nu)*(1.-2.*nu));
            mu = E[T] / (2.*(1.+nu));

            // Hole die Eckknoten des Dreiecks
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));

            // Berechne die Transormationsmatrix B fuer das jeweilige Element (2D)
            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also B^(-T) * \grad_phi bzw. \grad_phi^T * B^(-1).
            // Also \hat{grad_phi}.
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0; v12 = 0.0; v21 = 0.0; v22 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {
                        // In epsilonValuesMat1_i (2x2 Matrix) steht fuer die Ansatzfunktion i bzw. \phi_i
                        // der epsilonTensor fuer eine skalare Ansatzfunktion fuer die Richtung 1 (vgl. Mat1).
                        // Also in Mat1_i wird dann also phi_i = (phi_scalar_i, 0) gesetzt und davon \eps berechnet.

                        // Stelle \hat{grad_phi_i} = basisValues_i auf, also B^(-T)*grad_phi_i
                        // GradPhiOnRef( dPhi->at(k).at(i), b_T_inv, basisValues_i );

                        // \eps(v) = \eps(phi_i)
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat1_i, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat2_i, 1); // y-Richtung

                        // Siehe oben, nur fuer j
                        // GradPhiOnRef( DPhi->at(k).at(j), b_T_inv, basisValues_j );

                        // \eps(u) = \eps(phi_j)
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat1_j, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat2_j, 1); // y-Richtung

                        // Nun berechnen wir \eps(u):\eps(v) = \eps(phi_j):\eps(phi_i).
                        // Das Ergebniss steht in res.
                        // Berechne zudem noch die Spur der Epsilon-Tensoren tr(\eps(u)) (j) und tr(\eps(v)) (i)
                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat1_j, res); // x-x
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v11 = v11 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat2_j, res); // x-y
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v12 = v12 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat1_j, res); // y-x
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v21 = v21 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat2_j, res); // y-y
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v22 = v22 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);


                    }
                    // Noch mit der abs(det(B)) skalieren
                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;

                    value11[0] = v11;
                    value12[0] = v12;
                    value21[0] = v21;
                    value22[0] = v22;

                    // Hole die globale Zeile und Spalte in der die Eintraege hingeschrieben werden sollen
                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value11()); // x-x
                    A->insertGlobalValues(glob_i+1, indices(), value21()); // y-x
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value12()); // x-y
                    A->insertGlobalValues(glob_i+1, indices(), value22()); // y-y
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
    else if(dim == 3)
    {

        double v11, v12, v13, v21, v22, v23, v31, v32, v33;

        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0), p4(3,0.0);
        SmallMatrix<double> epsilonValuesMat1_i(dim), epsilonValuesMat2_i(dim), epsilonValuesMat3_i(dim),
        epsilonValuesMat1_j(dim), epsilonValuesMat2_j(dim), epsilonValuesMat3_j(dim);

        for (int T = 0; T < elements->numberElements(); T++)
        {
            lambda = E[T]* nu / ((1.+nu)*(1.-2.*nu));
            mu = E[T] / (2.*(1.+nu));

            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));
            p4 = pointsRep->at(elements->getElement(T).getNode(3));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                
                Teuchos::Array<SC> value11( 1, 0. );
                Teuchos::Array<SC> value12( 1, 0. );
                Teuchos::Array<SC> value13( 1, 0. );
                Teuchos::Array<SC> value21( 1, 0. );
                Teuchos::Array<SC> value22( 1, 0. );
                Teuchos::Array<SC> value23( 1, 0. );
                Teuchos::Array<SC> value31( 1, 0. );
                Teuchos::Array<SC> value32( 1, 0. );
                Teuchos::Array<SC> value33( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    v11 = 0.0; v12 = 0.0; v13 = 0.0; v21 = 0.0; v22 = 0.0; v23 = 0.0; v31 = 0.0; v32 = 0.0; v33 = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {

                        // GradPhiOnRef( DPhi->at(k).at(i), b_T_inv, basisValues_i );

                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat1_i, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat2_i, 1); // y-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(i), epsilonValuesMat3_i, 2); // z-Richtung


                        // GradPhiOnRef( DPhi->at(k).at(j), b_T_inv, basisValues_j );

                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat1_j, 0); // x-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat2_j, 1); // y-Richtung
                        epsilonTensor( dPhiTrans.at(k).at(j), epsilonValuesMat3_j, 2); // z-Richtung

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat1_j, res); // x-x
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v11 = v11 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat2_j, res); // x-y
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v12 = v12 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat1_i.innerProduct(epsilonValuesMat3_j, res); // x-z
                        epsilonValuesMat1_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v13 = v13 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat1_j, res); // y-x
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v21 = v21 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat2_j, res); // y-y
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v22 = v22 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat2_i.innerProduct(epsilonValuesMat3_j, res); // y-z
                        epsilonValuesMat2_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v23 = v23 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat1_j, res); // z-x
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat1_j.trace(res_trace_j);
                        v31 = v31 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat2_j, res); // z-y
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat2_j.trace(res_trace_j);
                        v32 = v32 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                        epsilonValuesMat3_i.innerProduct(epsilonValuesMat3_j, res); // z-z
                        epsilonValuesMat3_i.trace(res_trace_i);
                        epsilonValuesMat3_j.trace(res_trace_j);
                        v33 = v33 + weightsDPhi->at(k)*(2*mu*res + lambda*res_trace_j*res_trace_i);

                    }
                    v11 = absDetB * v11;
                    v12 = absDetB * v12;
                    v13 = absDetB * v13;
                    v21 = absDetB * v21;
                    v22 = absDetB * v22;
                    v23 = absDetB * v23;
                    v31 = absDetB * v31;
                    v32 = absDetB * v32;
                    v33 = absDetB * v33;

                    value11[0] = v11;
                    value12[0] = v12;
                    value13[0] = v13;
                    value21[0] = v21;
                    value22[0] = v22;
                    value23[0] = v23;
                    value31[0] = v31;
                    value32[0] = v32;
                    value33[0] = v33;

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value11()); // x-x
                    A->insertGlobalValues(glob_i+1, indices(), value21()); // y-x
                    A->insertGlobalValues(glob_i+2, indices(), value31()); // z-x
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value12()); // x-y
                    A->insertGlobalValues(glob_i+1, indices(), value22()); // y-y
                    A->insertGlobalValues(glob_i+2, indices(), value32()); // z-y
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value13()); // x-z
                    A->insertGlobalValues(glob_i+1, indices(), value23()); // y-z
                    A->insertGlobalValues(glob_i+2, indices(), value33()); // z-z
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyAdditionalConvection(int dim,
                                  std::string FEType,
                                  MatrixPtr_Type &A,
                                  MultiVectorPtr_Type w,
                                  bool callFillComplete)
{

    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    int FEloc = this->checkFE(dim,FEType);

    DomainConstPtr_Type domain = domainVec_.at(FEloc);
    ElementsPtr_Type elements = domain->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domain->getPointsRepeated();
    MapConstPtr_Type map = domain->getMapRepeated();

    vec3D_dbl_ptr_Type 			dPhi;
    vec2D_dbl_ptr_Type 	        phi;
    vec_dbl_ptr_Type			weights = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type			quadPts;

    UN extraDeg = Helper::determineDegree( dim, FEType, Grad); // Fuer diskretes (\grad \cdot w) in den Gausspuntken
    UN deg = Helper::determineDegree( dim, FEType, FEType, Std, Std, extraDeg);

    Helper::getDPhi(dPhi, weights, dim, FEType, deg);
    Helper::getPhi(phi, weights, dim, FEType, deg);
    Helper::getQuadratureValues(dim, deg, quadPts, weights,FEType);

    // SC = double, GO = long, UN = int
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;

    // Der nichtlineare Teil als Array
    Teuchos::ArrayRCP< const SC > wArray = w->getData(0);

    if (dim == 2)
    {
        double val;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0);

        vec2D_dbl_Type w11(1, vec_dbl_Type(weights->size(), -1.)); // diskretes w_11. Siehe unten.
        vec2D_dbl_Type w22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type divergenz(1, vec_dbl_Type(weights->size(), -1.));

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            // Diskretes \div(w) = (\grad \cdot w) = w_11 + w_22 berechnen,
            // wobei w_ij = \frac{\partial w_i}{\partial x_j} ist.
            for(int k = 0; k < dPhiTrans.size(); k++) // Quadraturpunkte
            {
                w11[0][k] = 0.0;
                w22[0][k] = 0.0;
                for(int i = 0; i < dPhiTrans[0].size(); i++)
                {
                    LO index1 = dim * elements->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements->getElement(T).getNode(i) + 1; // y
                    w11[0][k] += wArray[index1] * dPhiTrans[k][i][0];
                    w22[0][k] += wArray[index2] * dPhiTrans[k][i][1];

                    // TEST
                    // LO indexTest1 = dim * i + 0;
                    // LO indexTest2 = dim * i + 1;
                    // w11[0][k] += wTest[indexTest1] * dPhiTrans[k][i][0];
                    // w22[0][k] += wTest[indexTest2] * dPhiTrans[k][i][1];
                }
            }

            for(int k = 0; k < dPhiTrans.size(); k++) // Quadraturpunkte
            {
                divergenz[0][k] = w11[0][k] + w22[0][k];
                // if(T == 0)
                // {
                //     std::cout << "k: " << k << " Divergenz: " << divergenz[0][k] << '\n';
                // }
            }


            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    val = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {
                        val = val + divergenz[0][k] * weights->at(k) * (*phi)[k][i] * (*phi)[k][j];
                    }
                    val = absDetB * val;
                    value[0] = val;

                    // if(T == 0)
                    // {
                    //     std::cout << "i: " << i << " j: " << j << " val: " << val << '\n';
                    // }

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;

                    A->insertGlobalValues(glob_i, indices(), value());
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i+1, indices(), value());
                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
    else if(dim == 3)
    {
        double val;

        // long long glob_i, glob_j;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0), p4(3,0.0);

        vec2D_dbl_Type w11(1, vec_dbl_Type(weights->size(), -1.)); // diskretes w_11. Siehe unten.
        vec2D_dbl_Type w22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w33(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type divergenz(1, vec_dbl_Type(weights->size(), -1.));

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));
            p4 = pointsRep->at(elements->getElement(T).getNode(3));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTrans( dPhi->size(), vec2D_dbl_Type( dPhi->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhi, dPhiTrans, Binv ); //dPhiTrans berechnen

            // Diskretes \div(w) = (\grad \cdot w) = w_11 + w_22 + w33 berechnen,
            // wobei w_ij = \frac{\partial w_i}{\partial x_j} ist.
            for(int k = 0; k < dPhiTrans.size(); k++) // Quadraturpunkte
            {
                w11[0][k] = 0.0;
                w22[0][k] = 0.0;
                w33[0][k] = 0.0;
                for(int i = 0; i < dPhiTrans[0].size(); i++)
                {
                    LO index1 = dim * elements->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements->getElement(T).getNode(i) + 1; // y
                    LO index3 = dim * elements->getElement(T).getNode(i) + 2; // z
                    w11[0][k] += wArray[index1] * dPhiTrans[k][i][0];
                    w22[0][k] += wArray[index2] * dPhiTrans[k][i][1];
                    w33[0][k] += wArray[index3] * dPhiTrans[k][i][2];
                }
            }

            for(int k = 0; k < dPhiTrans.size(); k++) // Quadraturpunkte
            {
                divergenz[0][k] = w11[0][k] + w22[0][k] + w33[0][k];
            }

            for (int i = 0; i < dPhi->at(0).size(); i++)
            {
                Teuchos::Array<SC> value( 1, 0. );
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhi->at(0).size(); j++)
                {
                    val = 0.0;
                    for (int k = 0; k < dPhi->size(); k++)
                    {
                        val = val + divergenz[0][k] * weights->at(k) * (*phi)[k][i] * (*phi)[k][j];
                    }
                    val = absDetB * val;
                    value[0] = val;

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i, indices(), value());
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i+1, indices(), value());
                    glob_j++;
                    indices[0] = glob_j;
                    A->insertGlobalValues(glob_i+2, indices(), value());

                }
            }
        }
        if (callFillComplete)
        {
            A->fillComplete();
        }
    }
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyDummyCoupling(int dim,
                                          std::string FEType,
                                          MatrixPtr_Type &C,
                                          int FEloc, // 0 = Fluid, 2 = Struktur
                                          bool callFillComplete)
{
    DomainConstPtr_Type domain = domainVec_.at(FEloc);
    
    MapConstPtr_Type mapInterfaceVecField = domain->getInterfaceMapVecFieldUnique(); // Interface-Map in der Interface-Nummerierung
    MapConstPtr_Type mapGlobalInterfaceVecField = domain->getGlobalInterfaceMapVecFieldUnique(); // Interface-Map in der globalen Nummerierung
    
    MapConstPtr_Type mapFieldPartial = domain->getGlobalInterfaceMapVecFieldPartial();
    
    Teuchos::Array<SC> value( 1, 0. );
    value[0] = 1.0; // da Einheitsmatrix
    Teuchos::Array<GO> indices( 1, 0 );
    
    GO dofGlobal, dofLocal;
    
    for(int k = 0; k < mapGlobalInterfaceVecField->getNodeNumElements(); k++)
    {
        dofGlobal = mapGlobalInterfaceVecField->getGlobalElement(k);
        if ( mapFieldPartial->getLocalElement( dofGlobal ) == Teuchos::OrdinalTraits<LO>::invalid() ) {
            // Globale ID des Interface-Knotens bzgl. der Globalen oder Interface-Nummerierung
            dofGlobal = mapInterfaceVecField->getGlobalElement( k );
            indices[0] = dofGlobal;
            C->insertGlobalValues(dofGlobal, indices(), value());
        }
    }
    
    if (callFillComplete)
        C->fillComplete(mapInterfaceVecField, mapInterfaceVecField);

}
    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyFSICoupling(int dim,
                         std::string FEType,
                         MatrixPtr_Type &C,
                         MatrixPtr_Type &C_T,
                         int FEloc1, // 0 = Fluid, 2 = Struktur
                         int FEloc2, // 0 = Fluid, 2 = Struktur
                         MapConstPtr_Type map1, // DomainMap: InterfaceMapVecFieldUnique = Spalten von C_T
                         MapConstPtr_Type map2, // RangeMap: this->getDomain(0)->getMapVecFieldUnique() = Zeilen von C_T
                         bool callFillComplete)
{
    // int FEloc = this->checkFE(dim,FEType);

    DomainConstPtr_Type domain1 = domainVec_.at(FEloc1);

    MapConstPtr_Type mapInterfaceVecField = domain1->getInterfaceMapVecFieldUnique(); // Interface-Map in der Interface-Nummerierung

    MapConstPtr_Type mapGlobalInterfaceVecField;
    MapConstPtr_Type mapFieldPartial;
    if (FEloc1!=FEloc2){
        mapFieldPartial = domain1->getOtherGlobalInterfaceMapVecFieldPartial();
        mapGlobalInterfaceVecField = domain1->getOtherGlobalInterfaceMapVecFieldUnique();
    }
    else{
        mapFieldPartial = domain1->getGlobalInterfaceMapVecFieldPartial();
        mapGlobalInterfaceVecField = domain1->getGlobalInterfaceMapVecFieldUnique();
    }
    
    Teuchos::Array<SC> value( 1, 0. );
    value[0] = 1.0; // da Einheitsmatrix
    Teuchos::Array<GO> indices( 1, 0 );

    GO dofGlobal, dofLocal;
    if (mapFieldPartial.is_null()) {
        for(int k = 0; k < mapGlobalInterfaceVecField->getNodeNumElements(); k++)
        {
            // Globale ID des Interface-Knotens bzgl. der Globalen oder Interface-Nummerierung
            dofGlobal = mapGlobalInterfaceVecField->getGlobalElement(k);
            dofLocal = mapInterfaceVecField->getGlobalElement(k);
            
            indices[0] = dofLocal;
            C_T->insertGlobalValues(dofGlobal, indices(), value());
            indices[0] = dofGlobal;
            C->insertGlobalValues(dofLocal, indices(), value());
            
        }
    }
    else{
        for(int k = 0; k < mapGlobalInterfaceVecField->getNodeNumElements(); k++) {
            dofGlobal = mapGlobalInterfaceVecField->getGlobalElement(k);
            if ( mapFieldPartial->getLocalElement( dofGlobal ) != Teuchos::OrdinalTraits<LO>::invalid() ) {

                dofLocal = mapInterfaceVecField->getGlobalElement(k);
                
                indices[0] = dofLocal;
                C_T->insertGlobalValues(dofGlobal, indices(), value());
                indices[0] = dofGlobal;
                C->insertGlobalValues(dofLocal, indices(), value());
            }
        }
    }

    if (callFillComplete)
    {
        // Erstes Argument: Domain (=Spalten von C_T)
        // Zweites Arguement: Range (=Zeilen von C_T)
        C_T->fillComplete(map1, map2);
        C->fillComplete(map2, map1);
    }
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyGeometryCoupling(int dim,
                              std::string FEType,
                              MatrixPtr_Type &C,
                              int FEloc, // 0 = Fluid, 2 = Struktur, 4 = 0 = Geometrie
                              MapConstPtr_Type map1, // Fluid-Interface-Map
                              MapConstPtr_Type map2, // DomainMap: this->getDomain(2)->getMapVecFieldUnique() = Spalten von C
                              MapConstPtr_Type map3, // RangeMap: this->getDomain(4)->getMapVecFieldUnique() = Zeilen von C
                              bool callFillComplete)
{

    DomainConstPtr_Type domain = domainVec_.at(FEloc);

    MapConstPtr_Type mapInt = domain->getGlobalInterfaceMapVecFieldUnique(); // Interface-Map in der globalen Nummerierung
    MapConstPtr_Type mapOtherInt = domain->getOtherGlobalInterfaceMapVecFieldUnique(); // Interface-Map in der globalen Nummerierung von other. For FELoc=0 or =4, otherInterface has solid dofs
    MapConstPtr_Type mapPartInt = domain->getGlobalInterfaceMapVecFieldPartial();
    MapConstPtr_Type mapOtherPartInt = domain->getOtherGlobalInterfaceMapVecFieldPartial();
    Teuchos::Array<SC> value( 1, 0. );
    value[0] = 1.0; // da Einheitsmatrix
    Teuchos::Array<GO> indices( 1, 0 );

    GO dofRow;
    if (mapPartInt.is_null()) {
        for(int k = 0; k < mapInt->getNodeNumElements(); k++){
            dofRow = mapInt->getGlobalElement(k);
            indices[0] = mapOtherInt->getGlobalElement(k);
            C->insertGlobalValues(dofRow, indices(), value());
        }
    }
    else{
        for(int k = 0; k < mapPartInt->getNodeNumElements(); k++){
            dofRow = mapPartInt->getGlobalElement(k);
            indices[0] = mapOtherPartInt->getGlobalElement(k);
            C->insertGlobalValues(dofRow, indices(), value());
        }
    }
    if (callFillComplete)
    {
        // (Domain, Range)
        C->fillComplete(map2, map3);
    }
}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyShapeDerivativeVelocity(int dim,
                                    std::string FEType1, // P2
                                    std::string FEType2, // P1
                                    MatrixPtr_Type &D,
                                    int FEloc, // 0 = Fluid (Velocity)
                                    MultiVectorPtr_Type u, // Geschwindigkeit
                                    MultiVectorPtr_Type w, // Beschleunigung Gitter
                                    MultiVectorPtr_Type p, // Druck
                                    double dt, // Zeitschrittweite
                                    double rho, // Dichte vom Fluid
                                    double nu, // Viskositaet vom Fluid
                                    bool callFillComplete)
{
    // int FEloc = this->checkFE(dim,FEType1);

    DomainConstPtr_Type domain = domainVec_.at(FEloc);
    ElementsPtr_Type elements = domain->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domain->getPointsRepeated();
    MapConstPtr_Type map = domain->getMapRepeated();

    vec3D_dbl_ptr_Type 			dPhiU;
    vec2D_dbl_ptr_Type 	        phiU;
    vec2D_dbl_ptr_Type 	        phiP;
    vec_dbl_ptr_Type			weights = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type			quadPts;

    // Hoechste Quadraturordnung angeben (= Zusaetzlicher Term wg. non-conservativ); bei P2/P1 hier Ordnung 6
    UN extraDeg = Helper::determineDegree( dim, FEType1, Grad) + Helper::determineDegree( dim, FEType1, Grad);
    UN deg = Helper::determineDegree( dim, FEType1, FEType1, Std, Std, extraDeg);

    Helper::getDPhi(dPhiU, weights, dim, FEType1, deg);
    Helper::getPhi(phiU, weights, dim, FEType1, deg);
    Helper::getPhi(phiP, weights, dim, FEType2, deg);
    Helper::getQuadratureValues(dim, deg, quadPts, weights,FEType1);

    // SC = double, GO = long, UN = int
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;

    // Der nichtlineare Teil als Array
    Teuchos::ArrayRCP< const SC > uArray = u->getData(0);
    Teuchos::ArrayRCP< const SC > wArray = w->getData(0);
    Teuchos::ArrayRCP< const SC > pArray = p->getData(0);

    if (dim == 2)
    {
        double val11, val12, val21, val22;
        double valDK1_11, valDK1_12, valDK1_21, valDK1_22;
        double valDK2_11, valDK2_12, valDK2_21, valDK2_22;
        double valDN_11, valDN_12, valDN_21, valDN_22;
        double valDW_11, valDW_12, valDW_21, valDW_22;
        double valDP_11, valDP_12, valDP_21, valDP_22;
        double valDM_11, valDM_12, valDM_21, valDM_22;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0);

        // Alle diskreten Vektoren aufstellen, dabei bezeichnet Xij = X_ij,
        // also i-te Komponenten von X nach der j-ten Variablen abgeleitet.
        // Der Gradient ist bei mir wie folgt definiert: \grad(u) = [u11, u12; u21 u22] = [grad(u_1)^T; grad(u_2)^T]
        vec2D_dbl_Type u1Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u2Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w1Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w2Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type pLoc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma22(1, vec_dbl_Type(weights->size(), -1.));

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTransU( dPhiU->size(), vec2D_dbl_Type( dPhiU->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhiU, dPhiTransU, Binv ); //dPhiTrans berechnen

            // Diskrete Vektoren u1, u2, w1 und w2 berechnen
            for(int k = 0; k < phiU->size(); k++) // Quadraturpunkte
            {
                u1Loc[0][k] = 0.0;
                u2Loc[0][k] = 0.0;
                w1Loc[0][k] = 0.0;
                w2Loc[0][k] = 0.0;
                for(int i = 0; i < phiU->at(0).size(); i++)
                {
                    LO index1 = dim * elements->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements->getElement(T).getNode(i) + 1; // y
                    u1Loc[0][k] += uArray[index1] * phiU->at(k).at(i);
                    u2Loc[0][k] += uArray[index2] * phiU->at(k).at(i);
                    w1Loc[0][k] += wArray[index1] * phiU->at(k).at(i);
                    w2Loc[0][k] += wArray[index2] * phiU->at(k).at(i);
                   
                }
            }

            // Diskreten Vektor p berechnen
            // Beachte: phiP->size() = phiU->size()
            for(int k = 0; k < phiP->size(); k++) // Quadraturpunkte
            {
                pLoc[0][k] = 0.0;
                for(int i = 0; i < phiP->at(0).size(); i++)
                {
                    // Die ersten Eintraege in der Elementliste sind P1
                    // Alternativ elements2 holen
                    LO index = elements->getElement(T).getNode(i) + 0;
                    pLoc[0][k] += pArray[index] * phiP->at(k).at(i);
                    
                }
            }

            // Diskrete Grad-Vektoren berechnen,
            // wobei z.B. w_ij = \frac{\partial w_i}{\partial x_j} ist.
            for(int k = 0; k < dPhiTransU.size(); k++) // Quadraturpunkte
            {
                u11[0][k] = 0.0;
                u12[0][k] = 0.0;
                u21[0][k] = 0.0;
                u22[0][k] = 0.0;
                w11[0][k] = 0.0;
                w12[0][k] = 0.0;
                w21[0][k] = 0.0;
                w22[0][k] = 0.0;
                for(int i = 0; i < dPhiTransU[0].size(); i++)
                {
                    LO index1 = dim * elements->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements->getElement(T).getNode(i) + 1; // y
                    u11[0][k] += uArray[index1] * dPhiTransU[k][i][0];
                    u12[0][k] += uArray[index1] * dPhiTransU[k][i][1];
                    u21[0][k] += uArray[index2] * dPhiTransU[k][i][0];
                    u22[0][k] += uArray[index2] * dPhiTransU[k][i][1];
                    w11[0][k] += wArray[index1] * dPhiTransU[k][i][0];
                    w12[0][k] += wArray[index1] * dPhiTransU[k][i][1];
                    w21[0][k] += wArray[index2] * dPhiTransU[k][i][0];
                    w22[0][k] += wArray[index2] * dPhiTransU[k][i][1];

                }
            }

            // Diskretes \sigma = \rho * \nu * ( grad u + (grad u)^T ) - pI berechnen
            // Beachte: phiP->size() = phiU->size()
            for(int k = 0; k < dPhiTransU.size(); k++) // Quadraturpunkte
            {
                sigma11[0][k] = rho * nu * (u11[0][k] + u11[0][k]) - pLoc[0][k];
                sigma12[0][k] = rho * nu * (u12[0][k] + u21[0][k]);
                sigma21[0][k] = rho * nu * (u21[0][k] + u12[0][k]);
                sigma22[0][k] = rho * nu * (u22[0][k] + u22[0][k]) - pLoc[0][k];
            }


            for (int i = 0; i < dPhiU->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. ); // x-x
                Teuchos::Array<SC> value12( 1, 0. ); // x-y
                Teuchos::Array<SC> value21( 1, 0. ); // y-x
                Teuchos::Array<SC> value22( 1, 0. ); // y-y
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhiU->at(0).size(); j++)
                {
                    // DK1
                    valDK1_11 = 0.0;
                    valDK1_12 = 0.0;
                    valDK1_21 = 0.0;
                    valDK1_22 = 0.0;

                    // DK2
                    valDK2_11 = 0.0;
                    valDK2_12 = 0.0;
                    valDK2_21 = 0.0;
                    valDK2_22 = 0.0;

                    // DN
                    valDN_11 = 0.0;
                    valDN_12 = 0.0;
                    valDN_21 = 0.0;
                    valDN_22 = 0.0;

                    // DW
                    valDW_11 = 0.0;
                    valDW_12 = 0.0;
                    valDW_21 = 0.0;
                    valDW_22 = 0.0;

                    // DP
                    valDP_11 = 0.0;
                    valDP_12 = 0.0;
                    valDP_21 = 0.0;
                    valDP_22 = 0.0;

                    // DM
                    valDM_11 = 0.0;
                    valDM_12 = 0.0;
                    valDM_21 = 0.0;
                    valDM_22 = 0.0;

                    for (int k = 0; k < dPhiU->size(); k++)
                    {
                        // DK1
                        valDK1_11 = valDK1_11 +  weights->at(k) *
                                    ( 2 * u11[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][0] +
                                    u11[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][1] +
                                    u21[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] );
                        valDK1_12 = valDK1_12 +  weights->at(k) *
                                    ( 2 * u12[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][0] +
                                    u12[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][1] +
                                    u22[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] );
                        valDK1_21 = valDK1_21 +  weights->at(k) *
                                    ( u11[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    u21[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][0] +
                                    2 * u21[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][1] );
                        valDK1_22 = valDK1_22 +  weights->at(k) *
                                    ( u12[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    u22[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][0] +
                                    2 * u22[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][1] );

                        // DK2
                        valDK2_11 = valDK2_11 +  weights->at(k) *
                                    ( -sigma12[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    sigma12[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] );
                        valDK2_12 = valDK2_12 +  weights->at(k) *
                                    ( sigma11[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    -sigma11[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] );
                        valDK2_21 = valDK2_21 +  weights->at(k) *
                                    ( -sigma22[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    sigma22[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] );
                        valDK2_22 = valDK2_22 +  weights->at(k) *
                                    ( sigma21[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    -sigma21[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] );

                        // DN
                        valDN_11 = valDN_11 +  weights->at(k) *
                                    ( -(u2Loc[0][k] - w2Loc[0][k]) * dPhiTransU[k][j][1] * u11[0][k] * phiU->at(k).at(i) +
                                    (u2Loc[0][k] - w2Loc[0][k]) * dPhiTransU[k][j][0] * u12[0][k] * phiU->at(k).at(i) );
                        valDN_12 = valDN_12 +  weights->at(k) *
                                    ( (u1Loc[0][k] - w1Loc[0][k]) * dPhiTransU[k][j][1] * u11[0][k] * phiU->at(k).at(i) -
                                    (u1Loc[0][k] - w1Loc[0][k]) * dPhiTransU[k][j][0] * u12[0][k] * phiU->at(k).at(i) );
                        valDN_21 = valDN_21 +  weights->at(k) *
                                    ( -(u2Loc[0][k] - w2Loc[0][k]) * dPhiTransU[k][j][1] * u21[0][k] * phiU->at(k).at(i) +
                                    (u2Loc[0][k] - w2Loc[0][k]) * dPhiTransU[k][j][0] * u22[0][k] * phiU->at(k).at(i) );
                        valDN_22 = valDN_22 +  weights->at(k) *
                                    ( (u1Loc[0][k] - w1Loc[0][k]) * dPhiTransU[k][j][1] * u21[0][k] * phiU->at(k).at(i) -
                                    (u1Loc[0][k] - w1Loc[0][k]) * dPhiTransU[k][j][0] * u22[0][k] * phiU->at(k).at(i) );

                        // DW
                        valDW_11 = valDW_11 +  weights->at(k) *
                                    ( u11[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_12 = valDW_12 +  weights->at(k) *
                                    ( u12[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_21 = valDW_21 +  weights->at(k) *
                                    ( u21[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_22 = valDW_22 +  weights->at(k) *
                                    ( u22[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );

                        // DP
                        valDP_11 = valDP_11 +  weights->at(k) *
                                    ( ( -w21[0][k] * dPhiTransU[k][j][1] + w22[0][k] * dPhiTransU[k][j][0] ) * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDP_12 = valDP_12 +  weights->at(k) *
                                    ( ( w11[0][k] * dPhiTransU[k][j][1] - w12[0][k] * dPhiTransU[k][j][0] ) * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDP_21 = valDP_21 +  weights->at(k) *
                                    ( ( -w21[0][k] * dPhiTransU[k][j][1] + w22[0][k] * dPhiTransU[k][j][0] ) * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDP_22 = valDP_22 +  weights->at(k) *
                                    ( ( w11[0][k] * dPhiTransU[k][j][1] - w12[0][k] * dPhiTransU[k][j][0] ) * u2Loc[0][k] * phiU->at(k).at(i) );

                        // DM
                        valDM_11 = valDM_11 +  weights->at(k) *
                                    ( dPhiTransU[k][j][0] * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDM_12 = valDM_12 +  weights->at(k) *
                                    ( dPhiTransU[k][j][1] * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDM_21 = valDM_21 +  weights->at(k) *
                                    ( dPhiTransU[k][j][0] * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDM_22 = valDM_22 +  weights->at(k) *
                                    ( dPhiTransU[k][j][1] * u2Loc[0][k] * phiU->at(k).at(i) );
                    }

                    val11 = -rho*nu*valDK1_11 + valDK2_11 + rho*valDN_11 - rho*valDP_11 - (1.0/dt)*rho*valDW_11 + (0.5/dt)*rho*valDM_11;
                    val12 = -rho*nu*valDK1_12 + valDK2_12 + rho*valDN_12 - rho*valDP_12 - (1.0/dt)*rho*valDW_12 + (0.5/dt)*rho*valDM_12;
                    val21 = -rho*nu*valDK1_21 + valDK2_21 + rho*valDN_21 - rho*valDP_21 - (1.0/dt)*rho*valDW_21 + (0.5/dt)*rho*valDM_21;
                    val22 = -rho*nu*valDK1_22 + valDK2_22 + rho*valDN_22 - rho*valDP_22 - (1.0/dt)*rho*valDW_22 + (0.5/dt)*rho*valDM_22;

                    val11 = absDetB * val11;
                    val12 = absDetB * val12;
                    val21 = absDetB * val21;
                    val22 = absDetB * val22;

                    value11[0] = val11; // x-x
                    value12[0] = val12; // x-y
                    value21[0] = val21; // y-x
                    value22[0] = val22; // y-y

                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;

                    D->insertGlobalValues(glob_i, indices(), value11()); // x-x
                    D->insertGlobalValues(glob_i+1, indices(), value21()); // y-x
                    glob_j++;
                    indices[0] = glob_j;
                    D->insertGlobalValues(glob_i, indices(), value12()); // x-y
                    D->insertGlobalValues(glob_i+1, indices(), value22()); // y-y
                }
            }
        }
        if (callFillComplete)
        {
            D->fillComplete();
        }
    }
    else if(dim == 3)
    {
        double val11, val12, val13, val21, val22, val23, val31, val32, val33;
        double valDK1_11, valDK1_12, valDK1_13, valDK1_21, valDK1_22, valDK1_23, valDK1_31, valDK1_32, valDK1_33;
        double valDK2_11, valDK2_12, valDK2_13, valDK2_21, valDK2_22, valDK2_23, valDK2_31, valDK2_32, valDK2_33;
        double valDN_11, valDN_12, valDN_13, valDN_21, valDN_22, valDN_23, valDN_31, valDN_32, valDN_33;
        double valDW_11, valDW_12, valDW_13, valDW_21, valDW_22, valDW_23, valDW_31, valDW_32, valDW_33;
        double valDP_11, valDP_12, valDP_13, valDP_21, valDP_22, valDP_23, valDP_31, valDP_32, valDP_33;
        double valDM_11, valDM_12, valDM_13, valDM_21, valDM_22, valDM_23, valDM_31, valDM_32, valDM_33;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0), p4(3,0.0);

        // Alle diskreten Vektoren aufstellen, dabei bezeichnet Xij = X_ij,
        // also i-te Komponenten von X nach der j-ten Variablen abgeleitet.
        // Der Gradient ist bei mir wie folgt definiert: \grad(u) = [u11, u12; u21 u22] = [grad(u_1)^T; grad(u_2)^T]
        vec2D_dbl_Type u1Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u2Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u3Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w1Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w2Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w3Loc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type pLoc(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u13(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u23(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u31(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u32(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u33(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w13(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w23(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w31(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w32(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type w33(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma13(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma23(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma31(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma32(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type sigma33(1, vec_dbl_Type(weights->size(), -1.));

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));
            p4 = pointsRep->at(elements->getElement(T).getNode(3));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTransU( dPhiU->size(), vec2D_dbl_Type( dPhiU->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhiU, dPhiTransU, Binv ); //dPhiTrans berechnen

            // Diskrete Vektoren u1, u2, w1 und w2 berechnen
            for(int k = 0; k < phiU->size(); k++) // Quadraturpunkte
            {
                u1Loc[0][k] = 0.0;
                u2Loc[0][k] = 0.0;
                u3Loc[0][k] = 0.0;
                w1Loc[0][k] = 0.0;
                w2Loc[0][k] = 0.0;
                w3Loc[0][k] = 0.0;
                for(int i = 0; i < phiU->at(0).size(); i++)
                {
                    LO index1 = dim * elements->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements->getElement(T).getNode(i) + 1; // y
                    LO index3 = dim * elements->getElement(T).getNode(i) + 2; // z
                    u1Loc[0][k] += uArray[index1] * phiU->at(k).at(i);
                    u2Loc[0][k] += uArray[index2] * phiU->at(k).at(i);
                    u3Loc[0][k] += uArray[index3] * phiU->at(k).at(i);
                    w1Loc[0][k] += wArray[index1] * phiU->at(k).at(i);
                    w2Loc[0][k] += wArray[index2] * phiU->at(k).at(i);
                    w3Loc[0][k] += wArray[index3] * phiU->at(k).at(i);
                }
            }

            // Diskreten Vektor p berechnen
            // Beachte: phiP->size() = phiU->size()
            for(int k = 0; k < phiP->size(); k++) // Quadraturpunkte
            {
                pLoc[0][k] = 0.0;
                for(int i = 0; i < phiP->at(0).size(); i++)
                {
                    // Die ersten Eintraege in der Elementliste sind P1
                    LO index = elements->getElement(T).getNode(i) + 0;
                    pLoc[0][k] += pArray[index] * phiP->at(k).at(i);
                }
            }

            // Diskrete Grad-Vektoren berechnen,
            // wobei z.B. w_ij = \frac{\partial w_i}{\partial x_j} ist.
            for(int k = 0; k < dPhiTransU.size(); k++) // Quadraturpunkte
            {
                u11[0][k] = 0.0;
                u12[0][k] = 0.0;
                u13[0][k] = 0.0;
                u21[0][k] = 0.0;
                u22[0][k] = 0.0;
                u23[0][k] = 0.0;
                u31[0][k] = 0.0;
                u32[0][k] = 0.0;
                u33[0][k] = 0.0;
                w11[0][k] = 0.0;
                w12[0][k] = 0.0;
                w13[0][k] = 0.0;
                w21[0][k] = 0.0;
                w22[0][k] = 0.0;
                w23[0][k] = 0.0;
                w31[0][k] = 0.0;
                w32[0][k] = 0.0;
                w33[0][k] = 0.0;
                for(int i = 0; i < dPhiTransU[0].size(); i++)
                {
                    LO index1 = dim * elements->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements->getElement(T).getNode(i) + 1; // y
                    LO index3 = dim * elements->getElement(T).getNode(i) + 2; // z
                    u11[0][k] += uArray[index1] * dPhiTransU[k][i][0];
                    u12[0][k] += uArray[index1] * dPhiTransU[k][i][1];
                    u13[0][k] += uArray[index1] * dPhiTransU[k][i][2];
                    u21[0][k] += uArray[index2] * dPhiTransU[k][i][0];
                    u22[0][k] += uArray[index2] * dPhiTransU[k][i][1];
                    u23[0][k] += uArray[index2] * dPhiTransU[k][i][2];
                    u31[0][k] += uArray[index3] * dPhiTransU[k][i][0];
                    u32[0][k] += uArray[index3] * dPhiTransU[k][i][1];
                    u33[0][k] += uArray[index3] * dPhiTransU[k][i][2];
                    w11[0][k] += wArray[index1] * dPhiTransU[k][i][0];
                    w12[0][k] += wArray[index1] * dPhiTransU[k][i][1];
                    w13[0][k] += wArray[index1] * dPhiTransU[k][i][2];
                    w21[0][k] += wArray[index2] * dPhiTransU[k][i][0];
                    w22[0][k] += wArray[index2] * dPhiTransU[k][i][1];
                    w23[0][k] += wArray[index2] * dPhiTransU[k][i][2];
                    w31[0][k] += wArray[index3] * dPhiTransU[k][i][0];
                    w32[0][k] += wArray[index3] * dPhiTransU[k][i][1];
                    w33[0][k] += wArray[index3] * dPhiTransU[k][i][2];
                }
            }

            // Diskretes \sigma = \rho * \nu * ( grad u + (grad u)^T ) - pI berechnen
            // Beachte: phiP->size() = phiU->size()
            for(int k = 0; k < dPhiTransU.size(); k++) // Quadraturpunkte
            {
                sigma11[0][k] = rho * nu * (u11[0][k] + u11[0][k]) - pLoc[0][k];
                sigma12[0][k] = rho * nu * (u12[0][k] + u21[0][k]);
                sigma13[0][k] = rho * nu * (u13[0][k] + u31[0][k]);
                sigma21[0][k] = rho * nu * (u21[0][k] + u12[0][k]);
                sigma22[0][k] = rho * nu * (u22[0][k] + u22[0][k]) - pLoc[0][k];
                sigma23[0][k] = rho * nu * (u23[0][k] + u32[0][k]);
                sigma31[0][k] = rho * nu * (u31[0][k] + u13[0][k]);
                sigma32[0][k] = rho * nu * (u32[0][k] + u23[0][k]);
                sigma33[0][k] = rho * nu * (u33[0][k] + u33[0][k]) - pLoc[0][k];
            }


            for (int i = 0; i < dPhiU->at(0).size(); i++)
            {
                Teuchos::Array<SC> value11( 1, 0. ); // x-x
                Teuchos::Array<SC> value12( 1, 0. ); // x-y
                Teuchos::Array<SC> value13( 1, 0. ); // x-z
                Teuchos::Array<SC> value21( 1, 0. ); // y-x
                Teuchos::Array<SC> value22( 1, 0. ); // y-y
                Teuchos::Array<SC> value23( 1, 0. ); // y-z
                Teuchos::Array<SC> value31( 1, 0. ); // z-x
                Teuchos::Array<SC> value32( 1, 0. ); // z-y
                Teuchos::Array<SC> value33( 1, 0. ); // z-z
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhiU->at(0).size(); j++)
                {
                    // DK1
                    valDK1_11 = 0.0;
                    valDK1_12 = 0.0;
                    valDK1_13 = 0.0;
                    valDK1_21 = 0.0;
                    valDK1_22 = 0.0;
                    valDK1_23 = 0.0;
                    valDK1_31 = 0.0;
                    valDK1_32 = 0.0;
                    valDK1_33 = 0.0;

                    // DK2
                    valDK2_11 = 0.0;
                    valDK2_12 = 0.0;
                    valDK2_13 = 0.0;
                    valDK2_21 = 0.0;
                    valDK2_22 = 0.0;
                    valDK2_23 = 0.0;
                    valDK2_31 = 0.0;
                    valDK2_32 = 0.0;
                    valDK2_33 = 0.0;

                    // DN
                    valDN_11 = 0.0;
                    valDN_12 = 0.0;
                    valDN_13 = 0.0;
                    valDN_21 = 0.0;
                    valDN_22 = 0.0;
                    valDN_23 = 0.0;
                    valDN_31 = 0.0;
                    valDN_32 = 0.0;
                    valDN_33 = 0.0;

                    // DW
                    valDW_11 = 0.0;
                    valDW_12 = 0.0;
                    valDW_13 = 0.0;
                    valDW_21 = 0.0;
                    valDW_22 = 0.0;
                    valDW_23 = 0.0;
                    valDW_31 = 0.0;
                    valDW_32 = 0.0;
                    valDW_33 = 0.0;

                    // DP
                    valDP_11 = 0.0;
                    valDP_12 = 0.0;
                    valDP_13 = 0.0;
                    valDP_21 = 0.0;
                    valDP_22 = 0.0;
                    valDP_23 = 0.0;
                    valDP_31 = 0.0;
                    valDP_32 = 0.0;
                    valDP_33 = 0.0;

                    // DM
                    valDM_11 = 0.0;
                    valDM_12 = 0.0;
                    valDM_13 = 0.0;
                    valDM_21 = 0.0;
                    valDM_22 = 0.0;
                    valDM_23 = 0.0;
                    valDM_31 = 0.0;
                    valDM_32 = 0.0;
                    valDM_33 = 0.0;

                    for (int k = 0; k < dPhiU->size(); k++)
                    {
                        // DK1
                        valDK1_11 = valDK1_11 +  weights->at(k) *
                                    ( 2 * u11[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][0] +
                                    ( u11[0][k] * dPhiTransU[k][j][1] + u21[0][k] * dPhiTransU[k][j][0] ) * dPhiTransU[k][i][1] +
                                    ( u11[0][k] * dPhiTransU[k][j][2] + u31[0][k] * dPhiTransU[k][j][0] ) * dPhiTransU[k][i][2] );
                        valDK1_12 = valDK1_12 +  weights->at(k) *
                                    ( 2 * u12[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][0] +
                                    ( u12[0][k] * dPhiTransU[k][j][1] + u22[0][k] * dPhiTransU[k][j][0] ) * dPhiTransU[k][i][1] +
                                    ( u12[0][k] * dPhiTransU[k][j][2] + u32[0][k] * dPhiTransU[k][j][0] ) * dPhiTransU[k][i][2] );
                        valDK1_13 = valDK1_13 +  weights->at(k) *
                                    ( 2 * u13[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][0] +
                                    ( u13[0][k] * dPhiTransU[k][j][1] + u23[0][k] * dPhiTransU[k][j][0] ) * dPhiTransU[k][i][1] +
                                    ( u13[0][k] * dPhiTransU[k][j][2] + u33[0][k] * dPhiTransU[k][j][0] ) * dPhiTransU[k][i][2] );
                        valDK1_21 = valDK1_21 +  weights->at(k) *
                                    ( ( u21[0][k] * dPhiTransU[k][j][0] + u11[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][0] +
                                    2 * u21[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][1] +
                                    ( u21[0][k] * dPhiTransU[k][j][2] + u31[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][2] );
                        valDK1_22 = valDK1_22 +  weights->at(k) *
                                    ( ( u22[0][k] * dPhiTransU[k][j][0] + u12[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][0] +
                                    2 * u22[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][1] +
                                    ( u22[0][k] * dPhiTransU[k][j][2] + u32[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][2] );
                        valDK1_23 = valDK1_23 +  weights->at(k) *
                                    ( ( u23[0][k] * dPhiTransU[k][j][0] + u13[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][0] +
                                    2 * u23[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][1] +
                                    ( u23[0][k] * dPhiTransU[k][j][2] + u33[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][2] );
                        valDK1_31 = valDK1_31 +  weights->at(k) *
                                    ( ( u31[0][k] * dPhiTransU[k][j][0] + u11[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][0] +
                                    ( u31[0][k] * dPhiTransU[k][j][1] + u21[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][1] ) +
                                    2 * u31[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][2];
                        valDK1_32 = valDK1_32 +  weights->at(k) *
                                    ( ( u32[0][k] * dPhiTransU[k][j][0] + u12[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][0] +
                                    ( u32[0][k] * dPhiTransU[k][j][1] + u22[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][1] ) +
                                    2 * u32[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][2];
                        valDK1_33 = valDK1_33 +  weights->at(k) *
                                    ( ( u33[0][k] * dPhiTransU[k][j][0] + u13[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][0] +
                                    ( u33[0][k] * dPhiTransU[k][j][1] + u23[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][1] ) +
                                    2 * u33[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][2];

                        // DK2
                        valDK2_11 = valDK2_11 +  weights->at(k) *
                                    ( ( -sigma12[0][k] * dPhiTransU[k][j][1] - sigma13[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][0] +
                                    sigma12[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] +
                                    sigma13[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][2] );
                        valDK2_12 = valDK2_12 +  weights->at(k) *
                                    ( sigma11[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    ( -sigma11[0][k] * dPhiTransU[k][j][0] - sigma13[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][1] +
                                    sigma13[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][2] );
                        valDK2_13 = valDK2_13 +  weights->at(k) *
                                    ( sigma11[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][0] +
                                    sigma12[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][1] +
                                    ( -sigma11[0][k] * dPhiTransU[k][j][0] - sigma12[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][2] );
                        valDK2_21 = valDK2_21 +  weights->at(k) *
                                    ( ( -sigma22[0][k] * dPhiTransU[k][j][1] - sigma23[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][0] +
                                    sigma22[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] +
                                    sigma23[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][2] );
                        valDK2_22 = valDK2_22 +  weights->at(k) *
                                    ( sigma21[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    ( -sigma21[0][k] * dPhiTransU[k][j][0] - sigma23[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][1] +
                                    sigma23[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][2] );
                        valDK2_23 = valDK2_23 +  weights->at(k) *
                                    ( sigma21[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][0] +
                                    sigma22[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][1] +
                                    ( -sigma21[0][k] * dPhiTransU[k][j][0] - sigma22[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][2] );
                        valDK2_31 = valDK2_31 +  weights->at(k) *
                                    ( ( -sigma32[0][k] * dPhiTransU[k][j][1] - sigma33[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][0] +
                                    sigma32[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][1] +
                                    sigma33[0][k] * dPhiTransU[k][j][0] * dPhiTransU[k][i][2] );
                        valDK2_32 = valDK2_32 +  weights->at(k) *
                                    ( sigma31[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][0] +
                                    ( -sigma31[0][k] * dPhiTransU[k][j][0] - sigma33[0][k] * dPhiTransU[k][j][2] ) * dPhiTransU[k][i][1] +
                                    sigma33[0][k] * dPhiTransU[k][j][1] * dPhiTransU[k][i][2] );
                        valDK2_33 = valDK2_33 +  weights->at(k) *
                                    ( sigma31[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][0] +
                                    sigma32[0][k] * dPhiTransU[k][j][2] * dPhiTransU[k][i][1] +
                                    ( -sigma31[0][k] * dPhiTransU[k][j][0] - sigma32[0][k] * dPhiTransU[k][j][1] ) * dPhiTransU[k][i][2] );

                        // DN
                        double ZN_11; // Die Z_i fuer das DN wie in der Masterarbeit definiert
                        double ZN_12;
                        double ZN_13;
                        double ZN_21;
                        double ZN_22;
                        double ZN_23;
                        double ZN_31;
                        double ZN_32;
                        double ZN_33;
                        ZN_11 = - ( u2Loc[0][k] - w2Loc[0][k] ) * dPhiTransU[k][j][1] - ( u3Loc[0][k] - w3Loc[0][k] ) * dPhiTransU[k][j][2];
                        ZN_12 = ( u2Loc[0][k] - w2Loc[0][k] ) * dPhiTransU[k][j][0];
                        ZN_13 = ( u3Loc[0][k] - w3Loc[0][k] ) * dPhiTransU[k][j][0];
                        ZN_21 = ( u1Loc[0][k] - w1Loc[0][k] ) * dPhiTransU[k][j][1];
                        ZN_22 = - ( u1Loc[0][k] - w1Loc[0][k] ) * dPhiTransU[k][j][0] - ( u3Loc[0][k] - w3Loc[0][k] ) * dPhiTransU[k][j][2];
                        ZN_23 = ( u3Loc[0][k] - w3Loc[0][k] ) * dPhiTransU[k][j][1];
                        ZN_31 = ( u1Loc[0][k] - w1Loc[0][k] ) * dPhiTransU[k][j][2];
                        ZN_32 = ( u2Loc[0][k] - w2Loc[0][k] ) * dPhiTransU[k][j][2];
                        ZN_33 = - ( u1Loc[0][k] - w1Loc[0][k] ) * dPhiTransU[k][j][0] - ( u2Loc[0][k] - w2Loc[0][k] ) * dPhiTransU[k][j][1];

                        valDN_11 = valDN_11 +  weights->at(k) *
                                    ( ZN_11 * u11[0][k] * phiU->at(k).at(i) +
                                    ZN_12 * u12[0][k] * phiU->at(k).at(i) +
                                    ZN_13 * u13[0][k] * phiU->at(k).at(i) );
                        valDN_12 = valDN_12 +  weights->at(k) *
                                    ( ZN_21 * u11[0][k] * phiU->at(k).at(i) +
                                    ZN_22 * u12[0][k] * phiU->at(k).at(i) +
                                    ZN_23 * u13[0][k] * phiU->at(k).at(i) );
                        valDN_13 = valDN_13 +  weights->at(k) *
                                    ( ZN_31 * u11[0][k] * phiU->at(k).at(i) +
                                    ZN_32 * u12[0][k] * phiU->at(k).at(i) +
                                    ZN_33 * u13[0][k] * phiU->at(k).at(i) );
                        valDN_21 = valDN_21 +  weights->at(k) *
                                    ( ZN_11 * u21[0][k] * phiU->at(k).at(i) +
                                    ZN_12 * u22[0][k] * phiU->at(k).at(i) +
                                    ZN_13 * u23[0][k] * phiU->at(k).at(i) );
                        valDN_22 = valDN_22 +  weights->at(k) *
                                    ( ZN_21 * u21[0][k] * phiU->at(k).at(i) +
                                    ZN_22 * u22[0][k] * phiU->at(k).at(i) +
                                    ZN_23 * u23[0][k] * phiU->at(k).at(i) );
                        valDN_23 = valDN_23 +  weights->at(k) *
                                    ( ZN_31 * u21[0][k] * phiU->at(k).at(i) +
                                    ZN_32 * u22[0][k] * phiU->at(k).at(i) +
                                    ZN_33 * u23[0][k] * phiU->at(k).at(i) );
                        valDN_31 = valDN_31 +  weights->at(k) *
                                    ( ZN_11 * u31[0][k] * phiU->at(k).at(i) +
                                    ZN_12 * u32[0][k] * phiU->at(k).at(i) +
                                    ZN_13 * u33[0][k] * phiU->at(k).at(i) );
                        valDN_32 = valDN_32 +  weights->at(k) *
                                    ( ZN_21 * u31[0][k] * phiU->at(k).at(i) +
                                    ZN_22 * u32[0][k] * phiU->at(k).at(i) +
                                    ZN_23 * u33[0][k] * phiU->at(k).at(i) );
                        valDN_33 = valDN_33 +  weights->at(k) *
                                    ( ZN_31 * u31[0][k] * phiU->at(k).at(i) +
                                    ZN_32 * u32[0][k] * phiU->at(k).at(i) +
                                    ZN_33 * u33[0][k] * phiU->at(k).at(i) );

                        // DW
                        valDW_11 = valDW_11 +  weights->at(k) *
                                    ( u11[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_12 = valDW_12 +  weights->at(k) *
                                    ( u12[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_13 = valDW_13 +  weights->at(k) *
                                    ( u13[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_21 = valDW_21 +  weights->at(k) *
                                    ( u21[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_22 = valDW_22 +  weights->at(k) *
                                    ( u22[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_23 = valDW_23 +  weights->at(k) *
                                    ( u23[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_31 = valDW_31 +  weights->at(k) *
                                    ( u31[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_32 = valDW_32 +  weights->at(k) *
                                    ( u32[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );
                        valDW_33 = valDW_33 +  weights->at(k) *
                                    ( u33[0][k] * phiU->at(k).at(j) * phiU->at(k).at(i) );

                        // DP
                        double ZP_1; // Die Z_i fuer das DP wie in der Masterarbeit definiert
                        double ZP_2;
                        double ZP_3;
                        ZP_1 = -w21[0][k] * dPhiTransU[k][j][1] + w22[0][k] * dPhiTransU[k][j][0] -
                                w31[0][k] * dPhiTransU[k][j][2] + w33[0][k] * dPhiTransU[k][j][0];
                        ZP_2 = w11[0][k] * dPhiTransU[k][j][1] - w12[0][k] * dPhiTransU[k][j][0] -
                                w32[0][k] * dPhiTransU[k][j][2] + w33[0][k] * dPhiTransU[k][j][1];
                        ZP_3 = w11[0][k] * dPhiTransU[k][j][2] - w13[0][k] * dPhiTransU[k][j][0] +
                                w22[0][k] * dPhiTransU[k][j][2] - w23[0][k] * dPhiTransU[k][j][1];

                        valDP_11 = valDP_11 +  weights->at(k) *
                                    ( ZP_1 * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDP_12 = valDP_12 +  weights->at(k) *
                                    ( ZP_2 * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDP_13 = valDP_13 +  weights->at(k) *
                                    ( ZP_3 * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDP_21 = valDP_21 +  weights->at(k) *
                                    ( ZP_1 * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDP_22 = valDP_22 +  weights->at(k) *
                                    ( ZP_2 * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDP_23 = valDP_23 +  weights->at(k) *
                                    ( ZP_3 * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDP_31 = valDP_31 +  weights->at(k) *
                                    ( ZP_1 * u3Loc[0][k] * phiU->at(k).at(i) );
                        valDP_32 = valDP_32 +  weights->at(k) *
                                    ( ZP_2 * u3Loc[0][k] * phiU->at(k).at(i) );
                        valDP_33 = valDP_33 +  weights->at(k) *
                                    ( ZP_3 * u3Loc[0][k] * phiU->at(k).at(i) );

                        // DM
                        valDM_11 = valDM_11 +  weights->at(k) *
                                    ( dPhiTransU[k][j][0] * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDM_12 = valDM_12 +  weights->at(k) *
                                    ( dPhiTransU[k][j][1] * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDM_13 = valDM_13 +  weights->at(k) *
                                    ( dPhiTransU[k][j][2] * u1Loc[0][k] * phiU->at(k).at(i) );
                        valDM_21 = valDM_21 +  weights->at(k) *
                                    ( dPhiTransU[k][j][0] * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDM_22 = valDM_22 +  weights->at(k) *
                                    ( dPhiTransU[k][j][1] * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDM_23 = valDM_23 +  weights->at(k) *
                                    ( dPhiTransU[k][j][2] * u2Loc[0][k] * phiU->at(k).at(i) );
                        valDM_31 = valDM_31 +  weights->at(k) *
                                    ( dPhiTransU[k][j][0] * u3Loc[0][k] * phiU->at(k).at(i) );
                        valDM_32 = valDM_32 +  weights->at(k) *
                                    ( dPhiTransU[k][j][1] * u3Loc[0][k] * phiU->at(k).at(i) );
                        valDM_33 = valDM_33 +  weights->at(k) *
                                    ( dPhiTransU[k][j][2] * u3Loc[0][k] * phiU->at(k).at(i) );
                    }

                    val11 = -rho*nu*valDK1_11 + valDK2_11 + rho*valDN_11 - rho*valDP_11 - (1.0/dt)*rho*valDW_11 + (0.5/dt)*rho*valDM_11;
                    val12 = -rho*nu*valDK1_12 + valDK2_12 + rho*valDN_12 - rho*valDP_12 - (1.0/dt)*rho*valDW_12 + (0.5/dt)*rho*valDM_12;
                    val13 = -rho*nu*valDK1_13 + valDK2_13 + rho*valDN_13 - rho*valDP_13 - (1.0/dt)*rho*valDW_13 + (0.5/dt)*rho*valDM_13;
                    val21 = -rho*nu*valDK1_21 + valDK2_21 + rho*valDN_21 - rho*valDP_21 - (1.0/dt)*rho*valDW_21 + (0.5/dt)*rho*valDM_21;
                    val22 = -rho*nu*valDK1_22 + valDK2_22 + rho*valDN_22 - rho*valDP_22 - (1.0/dt)*rho*valDW_22 + (0.5/dt)*rho*valDM_22;
                    val23 = -rho*nu*valDK1_23 + valDK2_23 + rho*valDN_23 - rho*valDP_23 - (1.0/dt)*rho*valDW_23 + (0.5/dt)*rho*valDM_23;
                    val31 = -rho*nu*valDK1_31 + valDK2_31 + rho*valDN_31 - rho*valDP_31 - (1.0/dt)*rho*valDW_31 + (0.5/dt)*rho*valDM_31;
                    val32 = -rho*nu*valDK1_32 + valDK2_32 + rho*valDN_32 - rho*valDP_32 - (1.0/dt)*rho*valDW_32 + (0.5/dt)*rho*valDM_32;
                    val33 = -rho*nu*valDK1_33 + valDK2_33 + rho*valDN_33 - rho*valDP_33 - (1.0/dt)*rho*valDW_33 + (0.5/dt)*rho*valDM_33;

                    val11 = absDetB * val11;
                    val12 = absDetB * val12;
                    val13 = absDetB * val13;
                    val21 = absDetB * val21;
                    val22 = absDetB * val22;
                    val23 = absDetB * val23;
                    val31 = absDetB * val31;
                    val32 = absDetB * val32;
                    val33 = absDetB * val33;

                    value11[0] = val11; // x-x
                    value12[0] = val12; // x-y
                    value13[0] = val13; // x-z
                    value21[0] = val21; // y-x
                    value22[0] = val22; // y-y
                    value23[0] = val23; // y-z
                    value31[0] = val31; // z-x
                    value32[0] = val32; // z-y
                    value33[0] = val33; // z-z


                    glob_j = dim * map->getGlobalElement(elements->getElement(T).getNode(j));
                    glob_i = dim * map->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;

                    D->insertGlobalValues(glob_i, indices(), value11()); // x-x
                    D->insertGlobalValues(glob_i+1, indices(), value21()); // y-x
                    D->insertGlobalValues(glob_i+2, indices(), value31()); // z-x
                    glob_j++;
                    indices[0] = glob_j;
                    D->insertGlobalValues(glob_i, indices(), value12()); // x-y
                    D->insertGlobalValues(glob_i+1, indices(), value22()); // y-y
                    D->insertGlobalValues(glob_i+2, indices(), value32()); // z-y
                    glob_j++;
                    indices[0] = glob_j;
                    D->insertGlobalValues(glob_i, indices(), value13()); // x-z
                    D->insertGlobalValues(glob_i+1, indices(), value23()); // y-z
                    D->insertGlobalValues(glob_i+2, indices(), value33()); // z-z
                }
            }
        }
        if (callFillComplete)
        {
            D->fillComplete();
        }
    }

}


template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyShapeDerivativeDivergence(int dim,
                                       std::string FEType1,
                                       std::string FEType2,
                                       MatrixPtr_Type &DB,
                                       int FEloc1, // 1 = Fluid-Pressure
                                       int FEloc2, // 0 = Fluid-Velocity
                                       MapConstPtr_Type map1_unique, // Pressure-Map
                                       MapConstPtr_Type map2_unique, // Velocity-Map unique als VecField
                                       MultiVectorPtr_Type u, // Geschwindigkeit
                                       bool callFillComplete)
{
    DomainConstPtr_Type domain1 = domainVec_.at(FEloc1);
    ElementsPtr_Type elements = domain1->getElementsC();
    vec2D_dbl_ptr_Type pointsRep = domain1->getPointsRepeated();
    MapConstPtr_Type map1_rep = domain1->getMapRepeated();

    // Fuer die Fluid-Velocity-Map
    DomainConstPtr_Type domain2 = domainVec_.at(FEloc2);
    MapConstPtr_Type map2_rep = domain2->getMapRepeated();
    ElementsPtr_Type elements2 = domain2->getElementsC();

    vec3D_dbl_ptr_Type 			dPhiU;
    vec2D_dbl_ptr_Type 	        phiU;
    vec2D_dbl_ptr_Type 	        phiP;
    vec_dbl_ptr_Type			weights = Teuchos::rcp(new vec_dbl_Type(0));
    vec2D_dbl_ptr_Type			quadPts;

    UN extraDeg = Helper::determineDegree( dim, FEType1, Grad);
    UN deg = Helper::determineDegree( dim, FEType1, FEType1, Std, Grad, extraDeg);

    Helper::getDPhi(dPhiU, weights, dim, FEType1, deg);
    Helper::getPhi(phiU, weights, dim, FEType1, deg);
    Helper::getPhi(phiP, weights, dim, FEType2, deg);
    Helper::getQuadratureValues(dim, deg, quadPts, weights, FEType1);

    // SC = double, GO = long, UN = int
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    SmallMatrix<SC> Binv(dim);
    GO glob_i, glob_j;

    // Der nichtlineare Teil als Array
    Teuchos::ArrayRCP< const SC > uArray = u->getData(0);
   
    if (dim == 2)
    {
        double val1, val2;
        double valDB_1, valDB_2;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0);

        // Alle diskreten Vektoren aufstellen, dabei bezeichnet Xij = X_ij,
        // also i-te Komponenten von X nach der j-ten Variablen abgeleitet.
        // Der Gradient ist bei mir wie folgt definiert: \grad(u) = [u11, u12; u21 u22] = [grad(u_1)^T; grad(u_2)^T]
        vec2D_dbl_Type u11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u22(1, vec_dbl_Type(weights->size(), -1.));

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTransU( dPhiU->size(), vec2D_dbl_Type( dPhiU->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhiU, dPhiTransU, Binv ); //dPhiTrans berechnen

            // Diskrete Grad-Vektoren berechnen,
            // wobei z.B. w_ij = \frac{\partial w_i}{\partial x_j} ist.
            for(int k = 0; k < dPhiTransU.size(); k++) // Quadraturpunkte
            {
                u11[0][k] = 0.0;
                u12[0][k] = 0.0;
                u21[0][k] = 0.0;
                u22[0][k] = 0.0;
                for(int i = 0; i < dPhiTransU[0].size(); i++)
                {
                    LO index1 = dim * elements2->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements2->getElement(T).getNode(i) + 1; // y
                    u11[0][k] += uArray[index1] * dPhiTransU[k][i][0];
                    u12[0][k] += uArray[index1] * dPhiTransU[k][i][1];
                    u21[0][k] += uArray[index2] * dPhiTransU[k][i][0];
                    u22[0][k] += uArray[index2] * dPhiTransU[k][i][1];

                }
            }

            for (int i = 0; i < phiP->at(0).size(); i++)
            {
                Teuchos::Array<SC> value1( 1, 0. ); // p-x
                Teuchos::Array<SC> value2( 1, 0. ); // p-y
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhiU->at(0).size(); j++)
                {
                    valDB_1 = 0.0;
                    valDB_2 = 0.0;

                    for (int k = 0; k < dPhiU->size(); k++)
                    {
                        // DB
                        valDB_1 = valDB_1 +  weights->at(k) *
                                    ( phiP->at(k).at(i) * ( -u21[0][k] * dPhiTransU[k][j][1] + u22[0][k] * dPhiTransU[k][j][0] ) );
                        valDB_2 = valDB_2 +  weights->at(k) *
                                    ( phiP->at(k).at(i) * ( u11[0][k] * dPhiTransU[k][j][1] - u12[0][k] * dPhiTransU[k][j][0] ) );
                    }

                    val1 = valDB_1;
                    val2 = valDB_2;

                    val1 = absDetB * val1;
                    val2 = absDetB * val2;

                    value1[0] = val1; // p-x
                    value2[0] = val2; // p-y

                    glob_j = dim * map2_rep->getGlobalElement(elements2->getElement(T).getNode(j));
                    glob_i = map1_rep->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;

                    DB->insertGlobalValues(glob_i, indices(), value1()); // p-x
                    glob_j++;
                    indices[0] = glob_j;
                    DB->insertGlobalValues(glob_i, indices(), value2()); // p-y
                }
            }
        }
        if (callFillComplete)
        {
            DB->fillComplete(map2_unique, map1_unique);
        }
    }
    else if(dim == 3)
    {
        double val1, val2, val3;
        double valDB_1, valDB_2, valDB_3;
        vec_dbl_Type p1(3,0.0), p2(3,0.0), p3(3,0.0), p4(3,0.0);

        // Alle diskreten Vektoren aufstellen, dabei bezeichnet Xij = X_ij,
        // also i-te Komponenten von X nach der j-ten Variablen abgeleitet.
        // Der Gradient ist bei mir wie folgt definiert: \grad(u) = [u11, u12; u21 u22] = [grad(u_1)^T; grad(u_2)^T]
        vec2D_dbl_Type u11(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u12(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u13(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u21(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u22(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u23(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u31(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u32(1, vec_dbl_Type(weights->size(), -1.));
        vec2D_dbl_Type u33(1, vec_dbl_Type(weights->size(), -1.));

        for (int T = 0; T < elements->numberElements(); T++)
        {
            p1 = pointsRep->at(elements->getElement(T).getNode(0));
            p2 = pointsRep->at(elements->getElement(T).getNode(1));
            p3 = pointsRep->at(elements->getElement(T).getNode(2));
            p4 = pointsRep->at(elements->getElement(T).getNode(3));

            Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B);
            detB = B.computeInverse(Binv);
            absDetB = std::fabs(detB);

            // dPhiTrans sind die transformierten Basifunktionen, also \grad_phi * B^(-T)
            vec3D_dbl_Type dPhiTransU( dPhiU->size(), vec2D_dbl_Type( dPhiU->at(0).size(), vec_dbl_Type(dim,0.) ) );
            applyBTinv( dPhiU, dPhiTransU, Binv ); //dPhiTrans berechnen

            // Diskrete Grad-Vektoren berechnen,
            // wobei z.B. w_ij = \frac{\partial w_i}{\partial x_j} ist.
            for(int k = 0; k < dPhiTransU.size(); k++) // Quadraturpunkte
            {
                u11[0][k] = 0.0;
                u12[0][k] = 0.0;
                u13[0][k] = 0.0;
                u21[0][k] = 0.0;
                u22[0][k] = 0.0;
                u23[0][k] = 0.0;
                u31[0][k] = 0.0;
                u32[0][k] = 0.0;
                u33[0][k] = 0.0;

                for(int i = 0; i < dPhiTransU[0].size(); i++)
                {
                    LO index1 = dim * elements2->getElement(T).getNode(i) + 0; // x
                    LO index2 = dim * elements2->getElement(T).getNode(i) + 1; // y
                    LO index3 = dim * elements2->getElement(T).getNode(i) + 2; // z
                    u11[0][k] += uArray[index1] * dPhiTransU[k][i][0];
                    u12[0][k] += uArray[index1] * dPhiTransU[k][i][1];
                    u13[0][k] += uArray[index1] * dPhiTransU[k][i][2];
                    u21[0][k] += uArray[index2] * dPhiTransU[k][i][0];
                    u22[0][k] += uArray[index2] * dPhiTransU[k][i][1];
                    u23[0][k] += uArray[index2] * dPhiTransU[k][i][2];
                    u31[0][k] += uArray[index3] * dPhiTransU[k][i][0];
                    u32[0][k] += uArray[index3] * dPhiTransU[k][i][1];
                    u33[0][k] += uArray[index3] * dPhiTransU[k][i][2];
                }
            }

            for (int i = 0; i < phiP->at(0).size(); i++)
            {
                Teuchos::Array<SC> value1( 1, 0. ); // p-x
                Teuchos::Array<SC> value2( 1, 0. ); // p-y
                Teuchos::Array<SC> value3( 1, 0. ); // p-z
                Teuchos::Array<GO> indices( 1, 0 );

                for (int j = 0; j < dPhiU->at(0).size(); j++)
                {
                    valDB_1 = 0.0;
                    valDB_2 = 0.0;
                    valDB_3 = 0.0;

                    for (int k = 0; k < dPhiU->size(); k++)
                    {
                        // DB
                        valDB_1 = valDB_1 +  weights->at(k) *
                                    ( phiP->at(k).at(i) * ( -u21[0][k] * dPhiTransU[k][j][1] + u22[0][k] * dPhiTransU[k][j][0] -
                                                            u31[0][k] * dPhiTransU[k][j][2] + u33[0][k] * dPhiTransU[k][j][0] ) );
                        valDB_2 = valDB_2 +  weights->at(k) *
                                    ( phiP->at(k).at(i) * ( u11[0][k] * dPhiTransU[k][j][1] - u12[0][k] * dPhiTransU[k][j][0] -
                                                            u32[0][k] * dPhiTransU[k][j][2] + u33[0][k] * dPhiTransU[k][j][1] ) );
                        valDB_3 = valDB_3 +  weights->at(k) *
                                    ( phiP->at(k).at(i) * ( u11[0][k] * dPhiTransU[k][j][2] - u13[0][k] * dPhiTransU[k][j][0] +
                                                            u22[0][k] * dPhiTransU[k][j][2] - u23[0][k] * dPhiTransU[k][j][1] ) );
                    }

                    val1 = valDB_1;
                    val2 = valDB_2;
                    val3 = valDB_3;

                    val1 = absDetB * val1;
                    val2 = absDetB * val2;
                    val3 = absDetB * val3;

                    value1[0] = val1; // p-x
                    value2[0] = val2; // p-y
                    value3[0] = val3; // p-z

                    glob_j = dim * map2_rep->getGlobalElement(elements2->getElement(T).getNode(j));
                    glob_i = map1_rep->getGlobalElement(elements->getElement(T).getNode(i));
                    indices[0] = glob_j;

                    DB->insertGlobalValues(glob_i, indices(), value1()); // p-x
                    glob_j++;
                    indices[0] = glob_j;
                    DB->insertGlobalValues(glob_i, indices(), value2()); // p-y
                    glob_j++;
                    indices[0] = glob_j;
                    DB->insertGlobalValues(glob_i, indices(), value3()); // p-z
                }
            }
        }
        if (callFillComplete)
        {
            DB->fillComplete(map2_unique, map1_unique);
        }
    }

}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblySurfaceIntegralExternal(int dim,
                                              std::string FEType,
                                              MultiVectorPtr_Type f,
                                              MultiVectorPtr_Type d_rep,
                                              std::vector<SC>& funcParameter,
                                              RhsFunc_Type func,
                                              ParameterListPtr_Type params,
                                              int FEloc) {
    
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    
    SC elScaling;
    SmallMatrix<SC> B(dim);
    vec_dbl_Type b(dim);
    f->putScalar(0.);
    Teuchos::ArrayRCP< SC > valuesF = f->getDataNonConst(0);
    
    int flagSurface = params->sublist("Parameter Solid").get("Flag Surface",5); 
    
    std::vector<double> valueFunc(dim);

    SC* paramsFunc = &(funcParameter[0]);

    // The second last entry is a placeholder for the surface element flag. It will be set below
    for (UN T=0; T<elements->numberElements(); T++) {
        FiniteElement fe = elements->getElement( T );
        ElementsPtr_Type subEl = fe.getSubElements(); // might be null
        for (int surface=0; surface<fe.numSubElements(); surface++) {
            FiniteElement feSub = subEl->getElement( surface  );
            if(subEl->getDimension() == dim-1 ){
               
                vec_int_Type nodeList = feSub.getVectorNodeListNonConst ();

            
		        vec_dbl_Type solution_d = getSolution(nodeList, d_rep,dim);
                vec2D_dbl_Type nodes;
		        nodes = getCoordinates(nodeList, pointsRep);


                double positions[18];
                int count =0;
                for(int i=0;i<6;i++)
                    for(int j=0;j<3;j++){
                        positions[count] = nodes[i][j];
                        count++;

                    }
               		
               
                paramsFunc[ funcParameter.size() - 1 ] = feSub.getFlag();
                vec_dbl_Type p1 = {0.,0.,0.}; // Dummy vector
                func( &p1[0], &valueFunc[0], paramsFunc);
  
                if(valueFunc[0] != 0.){

                    double *residuumVector;
                    #ifdef FEDD_HAVE_ACEGENINTERFACE
                    
                    AceGenInterface::PressureTriangle3D6 pt(valueFunc[0], 1., 35, &positions[0], &solution_d[0]);
                    pt.computeTangentResidual();
                    residuumVector = pt.getResiduum();
                    #endif
                   
                    for(int i=0; i< nodeList.size() ; i++){
                            for(int d=0; d<dim; d++)
                                valuesF[nodeList[i]*dim+d] += residuumVector[i*dim+d];
                    }

                    // free(residuumVector);
                }
                    
            }
        }
    }
    //f->scale(-1.);

}
    
/// @brief  assemblyNonlinearSurfaceIntegralExternal -
/// @brief This force is assembled in AceGEN as deformation-dependent load. This force is applied as Pressure boundary in opposite direction of surface normal.

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyNonlinearSurfaceIntegralExternal(int dim,
                                              std::string FEType,
                                              MultiVectorPtr_Type f,
                                              MultiVectorPtr_Type d_rep,
                                              MatrixPtr_Type &Kext,
                                              std::vector<SC>& funcParameter,
                                              RhsFunc_Type func,
                                              ParameterListPtr_Type params,
                                              int FEloc) {
    
    // degree of function funcParameter[0]

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    SC elScaling;
    SmallMatrix<SC> B(dim);
    vec_dbl_Type b(dim);
    f->putScalar(0.);
    Teuchos::ArrayRCP< SC > valuesF = f->getDataNonConst(0);
        
    std::vector<double> valueFunc(dim);

    SC* paramsFunc = &(funcParameter[0]);

    // The second last entry is a placeholder for the surface element flag. It will be set below
    for (UN T=0; T<elements->numberElements(); T++) {
        FiniteElement fe = elements->getElement( T );
        ElementsPtr_Type subEl = fe.getSubElements(); // might be null
        for (int surface=0; surface<fe.numSubElements(); surface++) {
            FiniteElement feSub = subEl->getElement( surface  );
            if(subEl->getDimension() == dim-1 ){
                vec_int_Type nodeList = feSub.getVectorNodeListNonConst ();

		        vec_dbl_Type solution_d = getSolution(nodeList, d_rep,dim);
                vec2D_dbl_Type nodes;
		        nodes = getCoordinates(nodeList, pointsRep);


                double positions[18];
                int count =0;
                for(int i=0;i<6;i++){
                    for(int j=0;j<3;j++){
                        positions[count] = nodes[i][j];
                        count++;

                    }
                }

                vec_dbl_Type p1 = {0.,0.,0.}; // Dummy vector
                paramsFunc[ funcParameter.size() - 1 ] = feSub.getFlag();          
                func( &p1[0], &valueFunc[0], paramsFunc);
  
                if(valueFunc[0] != 0.){
                    
                    double *residuumVector;
                    double **stiffMat;

                    #ifdef FEDD_HAVE_ACEGENINTERFACE
                    AceGenInterface::PressureTriangle3D6 pt(valueFunc[0], 1.0, 35, &positions[0], &solution_d[0]);
                    pt.computeTangentResidual();

                    residuumVector = pt.getResiduum();
                    stiffMat = pt.getStiffnessMatrix();
                    #endif

                 

                    int dofs1 = dim;
                    int numNodes1 =nodeList.size();

                    SmallMatrix_Type elementMatrixPrint(18,0.);
                    for(int i=0; i< 18 ; i++){
                        for(int j=0; j< 18; j++){
                           if(fabs(stiffMat[i][j]) >1e-13)
                                elementMatrixPrint[i][j] = stiffMat[i][j];

                        }
                    }

                    SmallMatrix_Type elementMatrixWrite(18,0.);

                    SmallMatrix_Type elementMatrixIDsRow(18,0.);
                    SmallMatrix_Type elementMatrixIDsCol(18,0.);


                    for (UN i=0; i < numNodes1 ; i++) {
                        for(int di=0; di<dim; di++){
                            Teuchos::Array<SC> value1( numNodes1*dim, 0. );
                            Teuchos::Array<GO> columnIndices1( numNodes1*dim, 0 );
                            GO row =GO (dim* map->getGlobalElement( nodeList[i] )+di);
                            LO rowLO = dim*i+di;
                            // Zeilenweise werden die Einträge global assembliert
                            for (UN j=0; j <numNodes1; j++){
                                for(int d=0; d<dim; d++){
                                    columnIndices1[dim*j+d] = GO ( dim * map->getGlobalElement( nodeList[j] ) + d );
                                    value1[dim*j+d] = stiffMat[rowLO][dim*j+d];	
                                }
                            }  
                            Kext->insertGlobalValues( row, columnIndices1(), value1() ); // Automatically adds entries if a value already exists 
                        }
                    }
            

                                        
                    for(int i=0; i< nodeList.size() ; i++){
                        for(int d=0; d<dim; d++){
                            valuesF[nodeList[i]*dim+d] += residuumVector[i*dim+d];
                        }
                    }

                
                }
                
                    
            }
        }
    }
    //f->scale(-1.);
    Kext->fillComplete(domainVec_.at(FEloc)->getMapVecFieldUnique(),domainVec_.at(FEloc)->getMapVecFieldUnique());
    // Kext->writeMM("K_ext1");
}

/// Compute Surface Normal based on surface nodes,
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::computeSurfaceNormal(int dim,
                                            vec2D_dbl_ptr_Type pointsRep,
                                            vec_int_Type nodeList,
                                            vec_dbl_Type &v_E,
                                            double &norm_v_E)
{

    vec_dbl_Type p1(dim),p2(dim);

    if(dim==2){
        v_E[0] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[1]).at(1);
        v_E[1] = -(pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[1]).at(0));
        norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2));	
        
    }
    else if(dim==3){

        p1[0] = pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[1]).at(0);
        p1[1] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[1]).at(1);
        p1[2] = pointsRep->at(nodeList[0]).at(2) - pointsRep->at(nodeList[1]).at(2);

        p2[0] = pointsRep->at(nodeList[0]).at(0) - pointsRep->at(nodeList[2]).at(0);
        p2[1] = pointsRep->at(nodeList[0]).at(1) - pointsRep->at(nodeList[2]).at(1);
        p2[2] = pointsRep->at(nodeList[0]).at(2) - pointsRep->at(nodeList[2]).at(2);

        v_E[0] = p1[1]*p2[2] - p1[2]*p2[1];
        v_E[1] = p1[2]*p2[0] - p1[0]*p2[2];
        v_E[2] = p1[0]*p2[1] - p1[1]*p2[0];
        
        norm_v_E = sqrt(pow(v_E[0],2)+pow(v_E[1],2)+pow(v_E[2],2));
        
    }

}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblySurfaceIntegral(int dim,
                                              std::string FEType,
                                              MultiVectorPtr_Type f,
                                              std::string fieldType,
                                              RhsFunc_Type func,
                                              std::vector<SC>& funcParameter) {
    
    // degree of function funcParameter[0]
    //TEUCHOS_TEST_FOR_EXCEPTION( funcParameter[funcParameter.size()-1] > 0., std::logic_error, "We only support constant functions for now.");
    UN FEloc = checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec2D_dbl_ptr_Type phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    UN degFunc = funcParameter[funcParameter.size()-1] + 1.e-14;
    UN deg = Helper::determineDegree( dim-1, FEType, Std);// + 1.0;

    Helper::getPhi(phi, weights, dim-1, FEType, deg);

    vec2D_dbl_ptr_Type quadPoints;
    vec_dbl_ptr_Type w = Teuchos::rcp(new vec_dbl_Type(0));
    Helper::getQuadratureValues(dim-1, deg, quadPoints, w, FEType);
    w.reset();

    SC elScaling;
    SmallMatrix<SC> B(dim);
    vec_dbl_Type b(dim);
    f->putScalar(0.);
    Teuchos::ArrayRCP< SC > valuesF = f->getDataNonConst(0);
    int parameters;
    
    
    std::vector<double> valueFunc(dim);
    // The second last entry is a placeholder for the surface element flag. It will be set below
    SC* params = &(funcParameter[0]);
    for (UN T=0; T<elements->numberElements(); T++) {
        FiniteElement fe = elements->getElement( T );
        ElementsPtr_Type subEl = fe.getSubElements(); // might be null
        for (int surface=0; surface<fe.numSubElements(); surface++) {
            FiniteElement feSub = subEl->getElement( surface  );
            if(subEl->getDimension() == dim-1){
                // Setting flag to the placeholder (second last entry). The last entry at (funcParameter.size() - 1) should always be the degree of the surface function
                params[ funcParameter.size() - 1 ] = feSub.getFlag();
               
                vec_int_Type nodeList = feSub.getVectorNodeListNonConst ();

                vec_dbl_Type v_E(dim,1.);
                double norm_v_E=1.;

                Helper::computeSurfaceNormal(dim, pointsRep,nodeList,v_E,norm_v_E);

		        Helper::buildTransformationSurface( nodeList, pointsRep, B, b, FEType);
                elScaling = B.computeScaling( );
                // loop over basis functions
                for (UN i=0; i < phi->at(0).size(); i++) {
                    Teuchos::Array<SC> value(0);
                    if ( fieldType == "Scalar" )
                        value.resize( 1, 0. );
                    else if ( fieldType == "Vector" )
                        value.resize( dim, 0. );
                    // loop over basis functions quadrature points
                    for (UN w=0; w<phi->size(); w++) {
                        vec_dbl_Type x(dim,0.); //coordinates
                        for (int k=0; k<dim; k++) {// transform quad points to global coordinates
                            for (int l=0; l<dim-1; l++){
                                x[ k ] += B[k][l] * (*quadPoints)[ w ][ l ];
                            }   
                            x[k] += b[k];
                        }
                       
                        func( &x[0], &valueFunc[0], params);
                        if ( fieldType == "Scalar" )
                            value[0] += weights->at(w) * valueFunc[0] * (*phi)[w][i];
                        else if ( fieldType == "Vector" ){
                            for (int j=0; j<value.size(); j++){
                                value[j] += weights->at(w) * valueFunc[j]*v_E[j]/norm_v_E * (*phi)[w][i];
                            }
                        }
                    }

                    for (int j=0; j<value.size(); j++)
                        value[j] *= elScaling;
                    
                    if ( fieldType== "Scalar" )
                        valuesF[ nodeList[ i ] ] += value[0];


                    else if ( fieldType== "Vector" ){
                        for (int j=0; j<value.size(); j++)
                            valuesF[ dim * nodeList[ i ] + j ] += value[j];
                    }
                }
            }
        }
    }
    f->scale(-1.); // Generally the pressure is definied in opposite direction of the normal
}
    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblySurfaceIntegralFlag(int dim,
                                              std::string FEType,
                                              MultiVectorPtr_Type  f,
                                              std::string fieldType,
                                              BC_func_Type func,
                                              std::vector<SC>& funcParameter) {
    
// degree of function funcParameter[0]
    TEUCHOS_TEST_FOR_EXCEPTION(funcParameter[0]!=0,std::logic_error, "We only support constant functions for now.");
    
    UN FEloc = checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec2D_dbl_ptr_Type phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    UN degFunc = funcParameter[0] + 1.e-14;
    UN deg = Helper::determineDegree( dim-1, FEType, Std);// + degFunc;

    Helper::getPhi(phi, weights, dim-1, FEType, deg);

    vec2D_dbl_ptr_Type quadPoints;
    vec_dbl_ptr_Type w = Teuchos::rcp(new vec_dbl_Type(0));
    Helper::getQuadratureValues(dim-1, deg, quadPoints, w, FEType);
    w.reset();

    SC elScaling;
    SmallMatrix<SC> B(dim);
    vec_dbl_Type b(dim);
    f->putScalar(0.);
    Teuchos::ArrayRCP< SC > valuesF = f->getDataNonConst(0);
    int parameters;

    std::vector<double> valueFunc(dim);
    SC* params = &(funcParameter[1]);
    for (UN T=0; T<elements->numberElements(); T++) {
        FiniteElement fe = elements->getElement( T );
        ElementsPtr_Type subEl = fe.getSubElements(); // might be null
        for (int surface=0; surface<fe.numSubElements(); surface++) {
            FiniteElement feSub = subEl->getElement( surface  );
            if (params[1] == feSub.getFlag()){
                FiniteElement feSub = subEl->getElement( surface  );
                vec_int_Type nodeList = feSub.getVectorNodeListNonConst ();
                Helper::buildTransformationSurface( nodeList, pointsRep, B, b, FEType);
                elScaling = B.computeScaling( );
                // loop over basis functions
                for (UN i=0; i < phi->at(0).size(); i++) {
                    Teuchos::Array<SC> value(0);
                    if ( fieldType == "Scalar" )
                        value.resize( 1, 0. );
                    else if ( fieldType == "Vector" )
                        value.resize( dim, 0. );
                    // loop over basis functions quadrature points
                    for (UN w=0; w<phi->size(); w++) {
                        vec_dbl_Type x(dim,0.); //coordinates
                        for (int k=0; k<dim; k++) {// transform quad points to global coordinates
                            for (int l=0; l<dim-1; l++)
                                x[ k ] += B[k][l] * (*quadPoints)[ w ][ l ];
                        	x[k] += b[k];
                        }
                        
                        func( &x[0], &valueFunc[0], params[0], params);
//                        func( &x[0], &valueFunc[0], params);
                        if ( fieldType == "Scalar" )
                            value[0] += weights->at(w) * valueFunc[0] * (*phi)[w][i];
                        else if ( fieldType == "Vector" ){
                            for (int j=0; j<value.size(); j++){
                                value[j] += weights->at(w) * valueFunc[j] * (*phi)[w][i];
                            }
                        }
                    }

                    for (int j=0; j<value.size(); j++)
                        value[j] *= elScaling;

                    if ( fieldType== "Scalar" )
                        valuesF[ nodeList[ i ] ] += value[0];


                    else if ( fieldType== "Vector" ){
                        for (int j=0; j<value.size(); j++)
                            valuesF[ dim * nodeList[ i ] + j ] += value[j];
                    }
                }
            }
        }

    }
}
    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyRHS( int dim,
                                   std::string FEType,
                                   MultiVectorPtr_Type  a,
                                   std::string fieldType,
                                   RhsFunc_Type func,
                                   std::vector<SC>& funcParameter
                                  ) {

    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");

    TEUCHOS_TEST_FOR_EXCEPTION( a.is_null(), std::runtime_error, "MultiVector in assemblyConstRHS is null." );
    TEUCHOS_TEST_FOR_EXCEPTION( a->getNumVectors()>1, std::logic_error, "Implement for numberMV > 1 ." );
    UN FEloc;
    FEloc = checkFE(dim,FEType);

    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec2D_dbl_ptr_Type phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    // last parameter should alwayss be the degree
    UN degFunc = funcParameter[funcParameter.size()-1] + 1.e-14;
    UN deg = Helper::determineDegree( dim, FEType, Std);// + degFunc;

    vec2D_dbl_ptr_Type quadPoints;
    Helper::getQuadratureValues(dim, deg, quadPoints, weights, FEType); // quad points for rhs values


    Helper::getPhi(phi, weights, dim, FEType, deg);

    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);

    Teuchos::ArrayRCP< SC > valuesRhs = a->getDataNonConst(0);
    int parameters;
    double x;
    //for now just const!
    std::vector<double> valueFunc(dim);
    SC* paras = &(funcParameter[0]);
    
    func( &x, &valueFunc[0], paras );
    SC value;

    for (UN T=0; T<elements->numberElements(); T++) {

        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, FEType);
        detB = B.computeDet( );
        absDetB = std::fabs(detB);

		vec2D_dbl_Type quadPointsTrans(weights->size(),vec_dbl_Type(dim));
		for(int i=0; i< weights->size(); i++){
			 for(int j=0; j< dim ; j++){
				for(int k=0; k< dim; k++){
		 			quadPointsTrans[i][j] += B[j][k]* quadPoints->at(i).at(k) ; 
				}
				quadPointsTrans[i][j] += pointsRep->at(elements->getElement(T).getNode(0)).at(j); 
			 }
		}
        for (UN i=0; i < phi->at(0).size(); i++) {
            if ( !fieldType.compare("Scalar") ) {
		  	    value = Teuchos::ScalarTraits<SC>::zero();
 				for (UN w=0; w<weights->size(); w++){
					func(&quadPointsTrans[w][0], &valueFunc[0] ,paras);
	           		value += weights->at(w) * phi->at(w).at(i)*valueFunc[0];
				}
	            value *= absDetB;
	            LO row = (LO) elements->getElement(T).getNode(i);
	            valuesRhs[row] += value;
            }
            else if( !fieldType.compare("Vector") ) {
                for (UN d=0; d<dim; d++) {
		    		value = Teuchos::ScalarTraits<SC>::zero();
 					for (UN w=0; w<weights->size(); w++){
						func(&quadPointsTrans[w][0], &valueFunc[0] ,paras);
               			value += weights->at(w) * phi->at(w).at(i)*valueFunc[d];
					}
              		value *= absDetB;
                    SC v_i = value;
                    LO row = (LO) ( dim * elements->getElement(T).getNode(i)  + d );
                    valuesRhs[row] += v_i;
                }
            }
            else
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Invalid field type." );
        }
    }
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::assemblyRHSDegTest( int dim,
                                          std::string FEType,
                                          MultiVectorPtr_Type  a,
                                          std::string fieldType,
                                          RhsFunc_Type func,
                                          std::vector<SC>& funcParameter,
                                          int degree) {
    
    TEUCHOS_TEST_FOR_EXCEPTION(FEType == "P0",std::logic_error, "Not implemented for P0");
    
    TEUCHOS_TEST_FOR_EXCEPTION( a.is_null(), std::runtime_error, "MultiVector in assemblyConstRHS is null." );
    TEUCHOS_TEST_FOR_EXCEPTION( a->getNumVectors()>1, std::logic_error, "Implement for numberMV > 1 ." );
    
    UN FEloc = checkFE(dim,FEType);
    
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();
    
    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();
    
    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();
    vec2D_dbl_ptr_Type phi;
    vec_dbl_ptr_Type weights = Teuchos::rcp(new vec_dbl_Type(0));
    UN degFunc = funcParameter[0] + 1.e-14;
    UN deg = Helper::determineDegree( dim, FEType, Std) + degFunc;
    Helper::getPhi(phi, weights, dim, FEType, degree);
    
    vec2D_dbl_ptr_Type quadPoints;
    vec_dbl_ptr_Type w = Teuchos::rcp(new vec_dbl_Type(0));
    Helper::getQuadratureValues(dim, degree, quadPoints, w, FEType);
    w.reset();

    
    SC detB;
    SC absDetB;
    SmallMatrix<SC> B(dim);
    vec_dbl_Type b(dim);
    GO glob_i, glob_j;
    vec_dbl_Type v_i(dim);
    vec_dbl_Type v_j(dim);
    
    Teuchos::ArrayRCP< SC > valuesRhs = a->getDataNonConst(0);
    int parameters;
    double x;
    //for now just const!
    std::vector<double> valueFunc(dim);
    SC* params = &(funcParameter[1]);
    for (UN T=0; T<elements->numberElements(); T++) {
        
        Helper::buildTransformation(elements->getElement(T).getVectorNodeList(), pointsRep, B, b, FEType);
        detB = B.computeDet( );
        absDetB = std::fabs(detB);
        
        for (UN i=0; i < phi->at(0).size(); i++) {
            Teuchos::Array<SC> value(1);
            for (UN w=0; w<weights->size(); w++){
                vec_dbl_Type x(dim,0.); //coordinates
                for (int k=0; k<dim; k++) {// transform quad points to global coordinates
                    for (int l=0; l<dim; l++)
                        x[ k ] += B[k][l] * (*quadPoints)[ w ][ l ] + b[k];
                }
                
                func( &x[0], &valueFunc[0], params);
                if ( !fieldType.compare("Scalar") ) {
                    value[0] += weights->at(w) * valueFunc[0] * (*phi)[w][i];
                }
                else if( !fieldType.compare("Vector") ) {
                    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "No test for field type Vector." );
                }
                else
                    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error, "Invalid field type." );

            }
            value[0] *= absDetB;
            LO row = (LO) elements->getElement(T).getNode(i);
            valuesRhs[row] += value[0];

        }
    }
}
    

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::buildFullDPhi(vec3D_dbl_ptr_Type dPhi, Teuchos::Array<SmallMatrix<double> >& dPhiMat){

    TEUCHOS_TEST_FOR_EXCEPTION(dPhi->size()*dPhi->at(0).size()*dPhi->at(0).at(0).size() != dPhiMat.size(), std::logic_error, "Wrong sizes for dPhi and dPhiMat.");

    int dim = dPhi->at(0).at(0).size();
    int nmbBasisFunc = dPhi->at(0).size();
    int nmbTotalBasisFunc = nmbBasisFunc * dim;
    if (dim==2) {
        for (int p=0; p<dPhi->size(); p++) { //loop over quad points
            for (int i=0; i<nmbBasisFunc; i++) { //loop over basis functions
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i ][0][0] = dPhi->at(p).at(i).at(0);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i ][0][1] = dPhi->at(p).at(i).at(1);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i ][1][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i ][1][1] = 0.;

                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][0][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][0][1] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][1][0] = dPhi->at(p).at(i).at(0);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][1][1] = dPhi->at(p).at(i).at(1);
            }
        }
    }
    else if(dim==3){
        for (int p=0; p<dPhi->size(); p++) { //loop over quad points
            for (int i=0; i<nmbBasisFunc; i++) { //loop over basis functions
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][0][0] = dPhi->at(p).at(i).at(0);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][0][1] = dPhi->at(p).at(i).at(1);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][0][2] = dPhi->at(p).at(i).at(2);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][1][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][1][1] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][1][2] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][2][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][2][1] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i  ][2][2] = 0.;

                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][0][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][0][1] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][0][2] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][1][0] = dPhi->at(p).at(i).at(0);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][1][1] = dPhi->at(p).at(i).at(1);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][1][2] = dPhi->at(p).at(i).at(2);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][2][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][2][1] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 1 ][2][2] = 0.;

                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][0][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][0][1] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][0][2] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][1][0] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][1][1] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][1][2] = 0.;
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][2][0] = dPhi->at(p).at(i).at(0);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][2][1] = dPhi->at(p).at(i).at(1);
                dPhiMat[ p * nmbTotalBasisFunc  + dim*i + 2 ][2][2] = dPhi->at(p).at(i).at(2);
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::fillMatrixArray(SmallMatrix<double> &matIn, double* matArrayOut, std::string order,int offset){
    if (!order.compare("cols")) {
        for (int j=0; j<matIn.size(); j++) {
            for (int i=0; i<matIn.size(); i++) {
                matArrayOut[ j * matIn.size() + i + offset ] = matIn[i][j]; //Spalten der Matrix werden hintereinander in array geschrieben
            }
        }
    }
    else if(!order.compare("rows")) {
        for (int i=0; i<matIn.size(); i++) {
            for (int j=0; j<matIn.size(); j++) {
                matArrayOut[ i * matIn.size() + j + offset ] = matIn[i][j]; //Zeilen der Matrix werden hintereinander in array geschrieben
            }
        }
    }
    else
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Unknown ordering for matrix to array conversion. Choose rows or cols.");
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::epsilonTensor(vec_dbl_Type &basisValues, SmallMatrix<SC> &epsilonValues, int activeDof){

    for (int i=0; i<epsilonValues.size(); i++) {
        for (int j=0; j<epsilonValues.size(); j++) {
            epsilonValues[i][j] = 0.;
            if (i==activeDof) {
                epsilonValues[i][j] += 0.5*basisValues.at(j);
            }
            if (j==activeDof) {
                epsilonValues[i][j] += 0.5*basisValues.at(i);
            }
        }
    }
}

/*template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::phi(int dim,
                          int intFE,
                          int i,
                          vec_dbl_Type &p,
                          double* value){
    
    if (dim==1) {
        switch (intFE) {
            case 0: //P0
                switch (i) {
                    case 0:
                        *value = 1.;
                        break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        *value = ( 1. - p.at(0) );
                        break;
                    case 1:
                        *value = p.at(0);
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        *value = ( 1. - 3. * p[0] + 2. * p[0] *  p[0] );
                        break;
                    case 1:
                        *value = ( - p[0] + 2. * p[0] *  p[0] );
                        break;
                    case 2:
                        *value = ( 4. * p[0] - 4. * p[0] *  p[0] );
                        break;
                        
                }
                break;
            default:
                TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "Only P0,P1,P2 1D basis functions available." );
                break;
        }
    }
    else if (dim==2) {
        switch (intFE) {
            case 0://P0
                switch (i) {
                    case 0:
                        *value = 1.;
                        break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        *value = (1. - p.at(0)-p.at(1));
                        break;
                    case 1:
                        *value = p.at(0);
                        break;
                    case 2:
                        *value = p.at(1);
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        *value = -(1. - p.at(0)-p.at(1)) * (1 - 2.*(1-p.at(0) - p.at(1)));
                        break;
                    case 1:
                        *value = -p.at(0) *  (1 - 2*p.at(0));
                        break;
                    case 2:
                        *value = -p.at(1) *  (1 - 2*p.at(1));
                        break;
                    case 3:
                        *value = 4*p.at(0) * (1 - p.at(0)-p.at(1));
                        break;
                    case 4:
                        *value = 4*p.at(0)*p.at(1);
                        break;
                    case 5:
                        *value = 4*p.at(1) * (1 - p.at(0)-p.at(1));
                        break;
                }
                break;
        }
    }
    else if(dim==3){
        switch (intFE) {
            case 1://P1
                switch (i) {
                    case 0:
                        *value = (1. - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 1:
                        *value = p.at(0);
                        break;
                    case 2:
                        *value = p.at(1);
                        break;
                    case 3:
                        *value = p.at(2);
                        break;
                }
                break;
            case 2: //P2
                switch (i) {
                    case 0:
                        *value = (1. - p.at(0)-p.at(1)-p.at(2)) * (1 - 2*p.at(0) - 2*p.at(1) - 2*p.at(2));
                        break;
                    case 1:
                        *value = p.at(0) * (2*p.at(0) - 1);
                        break;
                    case 2:
                        *value = p.at(1) * (2*p.at(1) - 1);
                        break;
                    case 3:
                        *value = p.at(2) * (2*p.at(2) - 1);
                        break;
                    case 4:
                        *value = 4*p.at(0) * (1 - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 5:
                        *value = 4*p.at(0)*p.at(1);
                        break;
                    case 6:
                        *value = 4*p.at(1) * (1 - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 7:
                        *value = 4*p.at(2) * (1 - p.at(0)-p.at(1)-p.at(2));
                        break;
                    case 8:
                        *value = 4*p.at(0)*p.at(2);
                        break;
                    case 9:
                        *value = 4*p.at(1)*p.at(2);
                        break;
                }
                break;
            case 3: //Q1
            {
                double a = 1./8;
                switch (i) {
                    case 0:
                        *value = a * ( 1 - p[0] ) * ( 1 - p[1] ) * ( 1 - p[2] );
                        break;
                    case 1:
                        *value = a * ( 1 + p[0] ) * ( 1 - p[1] ) * ( 1 - p[2] );
                        break;
                    case 2:
                        *value = a * ( 1 + p[0] ) * ( 1 + p[1] ) * ( 1 - p[2] );
                        break;
                    case 3:
                        *value = a * ( 1 - p[0] ) * ( 1 + p[1] ) * ( 1 - p[2] );
                        break;
                    case 4:
                        *value = a * ( 1 - p[0] ) * ( 1 - p[1] ) * ( 1 + p[2] );
                        break;
                    case 5:
                        *value = a * ( 1 + p[0] ) * ( 1 - p[1] ) * ( 1 + p[2] );
                        break;
                    case 6:
                        *value = a * ( 1 + p[0] ) * ( 1 + p[1] ) * ( 1 + p[2] );
                        break;
                    case 7:
                        *value = a * ( 1 - p[0] ) * ( 1 + p[1] ) * ( 1 + p[2] );
                        break;
                    default:
                        break;
                }
                break;
            }
            case 4: //Q2
            {
                double a = 1./8;
                double b = 1./4;
                double c = 1./2;
                switch (i) {
                    case 0:
                        *value = 0.125*p[0]*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 1:
                        *value = 0.125*p[0]*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 2:
                        *value = -0.125*p[0]*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 3:
                        *value = -0.125*p[0]*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 4:
                        *value = -0.125*p[0]*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 5:
                        *value = -0.125*p[0]*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 6:
                        *value = 0.125*p[0]*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 7:
                        *value = 0.125*p[0]*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[1]*p[2] + 0.125*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 8:
                        *value = 0.250*p[1]*p[2] + -0.250*p[1]*p[1]*p[2] + -0.250*p[2]*p[2]*p[1] + -0.250*p[0]*p[0]*p[1]*p[2] + 0.250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 9:
                        *value = -0.250*p[0]*p[2] + -0.250*p[0]*p[0]*p[2] + 0.250*p[2]*p[2]*p[0] + 0.250*p[1]*p[1]*p[0]*p[2] + 0.250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 10:
                        *value = -0.250*p[1]*p[2] + -0.250*p[1]*p[1]*p[2] + 0.250*p[2]*p[2]*p[1] + 0.250*p[0]*p[0]*p[1]*p[2] + 0.250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 11:
                        *value = 0.250*p[0]*p[2] + -0.250*p[0]*p[0]*p[2] + -0.250*p[2]*p[2]*p[0] + -0.250*p[1]*p[1]*p[0]*p[2] + 0.250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 12:
                        *value = -0.250*p[1]*p[2] + 0.250*p[1]*p[1]*p[2] + -0.250*p[2]*p[2]*p[1] + 0.250*p[0]*p[0]*p[1]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 13:
                        *value = 0.250*p[0]*p[2] + 0.250*p[0]*p[0]*p[2] + 0.250*p[2]*p[2]*p[0] + -0.250*p[1]*p[1]*p[0]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 14:
                        *value = 0.250*p[1]*p[2] + 0.250*p[1]*p[1]*p[2] + 0.250*p[2]*p[2]*p[1] + -0.250*p[0]*p[0]*p[1]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 15:
                        *value = -0.250*p[0]*p[2] + 0.250*p[0]*p[0]*p[2] + -0.250*p[2]*p[2]*p[0] + 0.250*p[1]*p[1]*p[0]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[2]*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 16:
                        *value = 0.250*p[0]*p[1] + -0.250*p[0]*p[0]*p[1] + -0.250*p[1]*p[1]*p[0] + -0.250*p[2]*p[2]*p[0]*p[1] + 0.250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 17:
                        *value = -0.250*p[0]*p[1] + -0.250*p[0]*p[0]*p[1] + 0.250*p[1]*p[1]*p[0] + 0.250*p[2]*p[2]*p[0]*p[1] + 0.250*p[0]*p[0]*p[2]*p[2]*p[1] + -0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 18:
                        *value =0.250*p[0]*p[1] + 0.250*p[0]*p[0]*p[1] + 0.250*p[1]*p[1]*p[0] + -0.250*p[2]*p[2]*p[0]*p[1] + -0.250*p[0]*p[0]*p[2]*p[2]*p[1] + -0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 19:
                        *value = -0.250*p[0]*p[1] + 0.250*p[0]*p[0]*p[1] + -0.250*p[1]*p[1]*p[0] + 0.250*p[2]*p[2]*p[0]*p[1] + -0.250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 20:
                        *value = -0.500*p[1] + 0.500*p[1]*p[1] + 0.500*p[0]*p[0]*p[1] + 0.500*p[2]*p[2]*p[1] + -0.500*p[0]*p[0]*p[2]*p[2]*p[1] + -0.500*p[0]*p[0]*p[1]*p[1] + -0.500*p[1]*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 21:
                        *value = 0.500*p[0] + 0.500*p[0]*p[0] + -0.500*p[1]*p[1]*p[0] + -0.500*p[2]*p[2]*p[0] + 0.500*p[1]*p[1]*p[2]*p[2]*p[0] + -0.500*p[0]*p[0]*p[1]*p[1] + -0.500*p[0]*p[0]*p[2]*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 22:
                        *value = 0.500*p[1] + 0.500*p[1]*p[1] + -0.500*p[0]*p[0]*p[1] + -0.500*p[2]*p[2]*p[1] + 0.500*p[0]*p[0]*p[2]*p[2]*p[1] + -0.500*p[0]*p[0]*p[1]*p[1] + -0.500*p[1]*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 23:
                        *value = -0.500*p[0] + 0.500*p[0]*p[0] + 0.500*p[1]*p[1]*p[0] + 0.500*p[2]*p[2]*p[0] + -0.500*p[1]*p[1]*p[2]*p[2]*p[0] + -0.500*p[0]*p[0]*p[1]*p[1] + -0.500*p[0]*p[0]*p[2]*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 24:
                        *value = -0.500*p[2] + 0.500*p[2]*p[2] + 0.500*p[0]*p[0]*p[2] + 0.500*p[1]*p[1]*p[2] + -0.500*p[0]*p[0]*p[1]*p[1]*p[2] + -0.500*p[0]*p[0]*p[2]*p[2] + -0.500*p[1]*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 25:
                        *value = 1.000 + -1.000*p[0]*p[0] + -1.000*p[1]*p[1] + -1.000*p[2]*p[2] + 1.000*p[0]*p[0]*p[1]*p[1] + 1.000*p[0]*p[0]*p[2]*p[2] + 1.000*p[1]*p[1]*p[2]*p[2] + -1.000*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    case 26:
                        *value = 0.500*p[2] + 0.500*p[2]*p[2] + -0.500*p[0]*p[0]*p[2] + -0.500*p[1]*p[1]*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*p[2] + -0.500*p[0]*p[0]*p[2]*p[2] + -0.500*p[1]*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*p[2]*p[2];
                        break;
                    default:
                        break;
                    }
                    break;
                }
            case 5: //Q2-20
            {
                switch (i) {
                    case 0:
                        *value = -0.0625 + 0.1250*p[0]*p[0]*p[1]*p[2] + 0.1250*p[1]*p[1]*p[0]*p[2] + 0.1250*p[2]*p[2]*p[0]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + -0.1250*p[0]*p[1]*p[2];
                        
                        break;
                    case 1:
                        *value = -0.0625 + 0.1250*p[0]*p[0]*p[1]*p[2] + -0.1250*p[1]*p[1]*p[0]*p[2] + -0.1250*p[2]*p[2]*p[0]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + 0.1250*p[0]*p[1]*p[2];
                        
                        
                        
                        break;
                    case 2:
                        *value = -0.0625 + -0.1250*p[0]*p[0]*p[1]*p[2] + -0.1250*p[1]*p[1]*p[0]*p[2] + 0.1250*p[2]*p[2]*p[0]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + -0.1250*p[0]*p[1]*p[2];
                        
                        
                        
                        break;
                    case 3:
                        *value = -0.0625 + -0.1250*p[0]*p[0]*p[1]*p[2] + 0.1250*p[1]*p[1]*p[0]*p[2] + -0.1250*p[2]*p[2]*p[0]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + 0.1250*p[0]*p[1]*p[2];
                        
                        
                        
                        break;
                    case 4:
                        *value = -0.0625 + -0.1250*p[0]*p[0]*p[1]*p[2] + -0.1250*p[1]*p[1]*p[0]*p[2] + 0.1250*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + 0.1250*p[0]*p[1]*p[2];
                        
                        
                        
                        break;
                    case 5:
                        *value = -0.0625 + -0.1250*p[0]*p[0]*p[1]*p[2] + 0.1250*p[1]*p[1]*p[0]*p[2] + -0.1250*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + -0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + -0.1250*p[0]*p[1]*p[2];
                        
                        
                        break;
                    case 6:
                        *value = -0.0625 + 0.1250*p[0]*p[0]*p[1]*p[2] + 0.1250*p[1]*p[1]*p[0]*p[2] + 0.1250*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + 0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + 0.1250*p[0]*p[1]*p[2];
                        
                        
                        
                        break;
                    case 7:
                        *value = -0.0625 + 0.1250*p[0]*p[0]*p[1]*p[2] + -0.1250*p[1]*p[1]*p[0]*p[2] + -0.1250*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1]*p[2] + 0.1250*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[1]*p[1]*p[2]*p[2]*p[0] + 0.0625*p[0]*p[0]*p[1]*p[1] + 0.0625*p[0]*p[0]*p[2]*p[2] + 0.0625*p[1]*p[1]*p[2]*p[2] + -0.1250*p[0]*p[1]*p[2];
                        
                        
                        break;
                    case 8:
                        *value = 0.1250 + 0.2500*p[1]*p[2] + -0.2500*p[1]*p[1]*p[2] + -0.2500*p[2]*p[2]*p[1] + -0.2500*p[0]*p[0]*p[1]*p[2] + 0.2500*p[0]*p[0]*p[1]*p[1]*p[2] + 0.2500*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + 0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 9:
                        *value = 0.1250 + -0.2500*p[2] + -0.2500*p[0]*p[2] + 0.2500*p[1]*p[1]*p[2] + 0.2500*p[2]*p[2]*p[0] + 0.2500*p[1]*p[1]*p[0]*p[2] + -0.2500*p[1]*p[1]*p[2]*p[2]*p[0] + -0.1250*p[0]*p[0]*p[1]*p[1] + 0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 10:
                        *value = 0.1250 + -0.2500*p[1]*p[2] + -0.2500*p[1]*p[1]*p[2] + 0.2500*p[2]*p[2]*p[1] + 0.2500*p[0]*p[0]*p[1]*p[2] + 0.2500*p[0]*p[0]*p[1]*p[1]*p[2] + -0.2500*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + 0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 11:
                        *value = 0.1250 + -0.2500*p[2] + 0.2500*p[0]*p[2] + 0.2500*p[1]*p[1]*p[2] + -0.2500*p[2]*p[2]*p[0] + -0.2500*p[1]*p[1]*p[0]*p[2] + 0.2500*p[1]*p[1]*p[2]*p[2]*p[0] + -0.1250*p[0]*p[0]*p[1]*p[1] + 0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 12:
                        *value = 0.1250 + -0.2500*p[1]*p[2] + 0.2500*p[1]*p[1]*p[2] + -0.2500*p[2]*p[2]*p[1] + 0.2500*p[0]*p[0]*p[1]*p[2] + -0.2500*p[0]*p[0]*p[1]*p[1]*p[2] + 0.2500*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + 0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 13:
                        *value = 0.1250 + 0.2500*p[2] + 0.2500*p[0]*p[2] + -0.2500*p[1]*p[1]*p[2] + 0.2500*p[2]*p[2]*p[0] + -0.2500*p[1]*p[1]*p[0]*p[2] + -0.0000*p[2]*p[2]*p[0]*p[1] + -0.2500*p[1]*p[1]*p[2]*p[2]*p[0] + -0.1250*p[0]*p[0]*p[1]*p[1] + 0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 14:
                        *value = 0.1250 + 0.2500*p[1]*p[2] + 0.2500*p[1]*p[1]*p[2] + 0.2500*p[2]*p[2]*p[1] + -0.2500*p[0]*p[0]*p[1]*p[2] + 0.0000*p[1]*p[1]*p[0]*p[2] + -0.2500*p[0]*p[0]*p[1]*p[1]*p[2] + -0.2500*p[0]*p[0]*p[2]*p[2]*p[1] + -0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + 0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 15:
                        *value = 0.1250 + 0.2500*p[2] + -0.2500*p[0]*p[2] + -0.2500*p[1]*p[1]*p[2] + -0.2500*p[2]*p[2]*p[0] + -0.0000*p[0]*p[0]*p[1]*p[2] + 0.2500*p[1]*p[1]*p[0]*p[2] + 0.2500*p[1]*p[1]*p[2]*p[2]*p[0] + -0.1250*p[0]*p[0]*p[1]*p[1] + 0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 16:
                        *value = 0.1250 + -0.2500*p[0] + -0.2500*p[1] + 0.2500*p[0]*p[1] + 0.2500*p[2]*p[2]*p[0] + 0.2500*p[2]*p[2]*p[1] + -0.2500*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 17:
                        *value = 0.1250 + 0.2500*p[0] + -0.2500*p[1] + -0.2500*p[0]*p[1] + -0.2500*p[2]*p[2]*p[0] + 0.2500*p[2]*p[2]*p[1] + 0.2500*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        
                        break;
                    case 18:
                        *value = 0.1250 + 0.2500*p[0] + 0.2500*p[1] + 0.2500*p[0]*p[1] + -0.2500*p[2]*p[2]*p[0] + -0.2500*p[2]*p[2]*p[1] + -0.2500*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        break;
                    case 19:
                        *value = 0.1250 + -0.2500*p[0] + 0.2500*p[1] + -0.2500*p[0]*p[1] + 0.2500*p[2]*p[2]*p[0] + -0.2500*p[2]*p[2]*p[1] + 0.2500*p[2]*p[2]*p[0]*p[1] + 0.1250*p[0]*p[0]*p[1]*p[1] + -0.1250*p[0]*p[0]*p[2]*p[2] + -0.1250*p[1]*p[1]*p[2]*p[2];
                        
                        
                        break;
                    default:
                        break;
                }
                break;
            }
                
        }

    }
}
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::buildTransformation(const vec_int_Type& element,
                                          vec2D_dbl_ptr_Type pointsRep,
                                          SmallMatrix<SC>& B,
                                          std::string FEType){

    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = element.at(0);
    if (FEType[0]=='P') {
        for (UN j=0; j<B.size(); j++) {
            index = element.at(j+1);
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = pointsRep->at(index).at(i) - pointsRep->at(index0).at(i);
            }
        }
    }
    else if (FEType[0]=='Q'){
        TEUCHOS_TEST_FOR_EXCEPTION( B.size()!=3, std::logic_error, "Transformation for quadrilateral elements only in 3D.");
        std::vector<int> indexVec(3);
        indexVec[0] = element[1]; indexVec[1] = element[3]; indexVec[2] = element[4];
        for (UN j=0; j<B.size(); j++) {
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = ( pointsRep->at( indexVec[j] ).at(i) - pointsRep->at( index0 ).at(i) ) / 2.;
            }
        }
    }
}
    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::buildTransformation(const vec_int_Type& element,
                                          vec2D_dbl_ptr_Type pointsRep,
                                          SmallMatrix<SC>& B,
                                          vec_dbl_Type& b,
                                          std::string FEType){
    
    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = element.at(0);
    if (FEType[0]=='P') {
        for (UN j=0; j<B.size(); j++) {
            index = element.at(j+1);
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = pointsRep->at(index).at(i) - pointsRep->at(index0).at(i);
            }
        }
        for (UN i=0; i<B.size(); i++)
            b[i] = pointsRep->at(index0).at(i);
    }
    else if (FEType[0]=='Q'){
        TEUCHOS_TEST_FOR_EXCEPTION( B.size()!=3, std::logic_error, "Transformation for quadrilateral elements only in 3D.");
        std::vector<int> indexVec(3);
        indexVec[0] = element[1]; indexVec[1] = element[3]; indexVec[2] = element[4];
        for (UN j=0; j<B.size(); j++) {
            for (UN i=0; i<B.size(); i++) {
                B[i][j] = ( pointsRep->at( indexVec[j] ).at(i) - pointsRep->at( index0 ).at(i) ) / 2.;
            }
        }
        for (UN i=0; i<B.size(); i++)
            b[i] = pointsRep->at(index0).at(i);

    }
}
    
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::buildTransformationSurface(const vec_int_Type& element,
                                                 vec2D_dbl_ptr_Type pointsRep,
                                                 SmallMatrix<SC>& B,
                                                 vec_dbl_Type& b,
                                                 std::string FEType){
    // small matrix always square
    TEUCHOS_TEST_FOR_EXCEPTION( (B.size()<2 || B.size()>3), std::logic_error, "Initialize SmallMatrix for transformation.");
    UN index;
    UN index0 = element.at(0);
    if (FEType[0]=='P') {
        for (UN j=0; j<B.size()-1; j++) {
            index = element.at(j+1);
            for (UN i=0; i<B.size(); i++) { // dimension
                B[i][j] = pointsRep->at(index).at(i) - pointsRep->at(index0).at(i);
            }
        }
        for (UN i=0; i<B.size(); i++)
            b[i] = pointsRep->at(index0).at(i);
    }
    else if (FEType[0]=='Q'){
        TEUCHOS_TEST_FOR_EXCEPTION( B.size()!=3, std::logic_error, "No Transformation for surface integrals.");
    }
}

template <class SC, class LO, class GO, class NO>
UN FE<SC,LO,GO,NO>::determineDegree(UN dim, std::string FEType1, std::string FEType2, VarType type1,VarType type2, UN extraDeg){

    TEUCHOS_TEST_FOR_EXCEPTION( dim==2 && ( FEType1=="P2-CR" || FEType2=="P2-CR"), std::runtime_error, "P2-CR should be only available in 3D.");
    UN deg1, deg2;
    if (!FEType1.compare("P0")) {
        deg1 = 0;
    }
    else if ( !FEType1.compare("P1") || !FEType1.compare("P1-disc") ) {
        if (type1==Std)
            deg1 = 1;
        else if (type1==Grad)
            deg1 = 0;
    }
    else if (!FEType1.compare("P2")) {
        if (type1==Std)
            deg1 = 2;
        else if (type1==Grad)
            deg1 = 1;
    }
    else if (!FEType1.compare("P2-CR")) {
        if (type1==Std)
            deg1 = 4;
        else if (type1==Grad)
            deg1 = 3;
    }

    else if (!FEType1.compare("Q2")) {
        if (type1==Std)
            deg1 = 2;
        else if (type1==Grad)
            deg1 = 2;
    }
    else if (!FEType1.compare("Q2-20")) {
        if (type1==Std)
            deg1 = 2;
        else if (type1==Grad)
            deg1 = 2;
    }
    


    if (!FEType2.compare("P0")) {
        deg2 = 0;
    }
    else if ( !FEType2.compare("P1") || !FEType2.compare("P1-disc") ) {
        if (type2==Std)
            deg2 = 1;
        else if (type2==Grad)
            deg2 = 0;
    }
    else if (!FEType2.compare("P2")) {
        if (type2==Std)
            deg2 = 2;
        else if (type2==Grad)
            deg2 = 1;
    }
    else if (!FEType2.compare("P2-CR")) {
        if (type2==Std)
            deg2 = 4;
        else if (type2==Grad)
            deg2 = 3;
    }

    else if (!FEType2.compare("Q2")) {
        if (type2==Std)
            deg2 = 2;
        else if (type2==Grad)
            deg2 = 2;
    }
    else if (!FEType2.compare("Q2-20")) {
        if (type2==Std)
            deg2 = 2;
        else if (type2==Grad)
            deg2 = 2;
    }

    UN deg = deg1+deg2+extraDeg;
    if (deg==0)
        deg = 1;
    
    return deg;
}


template <class SC, class LO, class GO, class NO>
UN FE<SC,LO,GO,NO>::determineDegree(UN dim, std::string FEType, VarType type){
    UN deg;
    if (!FEType.compare("P0")) {
        deg = 0;
    }
    else if (!FEType.compare("P1")) {
        if (type==Std)
            deg = 1;
        else if (type==Grad)
            deg = 0;
    }
    else if (!FEType.compare("P2")) {
        if (type==Std)
            deg = 2;
        else if (type==Grad)
            deg = 1;
    }
    else if (!FEType.compare("Q2")) {
        if (type==Std)
            deg = 2;
        else if (type==Grad)
            deg = 2;
    }
    
    if (deg==0)
        deg = 1;
    return deg;
}

template <class SC, class LO, class GO, class NO>
UN FE<SC,LO,GO,NO>::determineDegree(UN dim, std::string FEType, UN degFunc){
    UN deg;
    if (!FEType.compare("P0"))
        deg = 0;
    else if (!FEType.compare("P1"))
        deg = 1;
    else if (!FEType.compare("P2"))
        deg = 2;
    else if (!FEType.compare("Q2"))
        deg = 2;
    
    deg += degFunc;

    if (deg==0)
        deg = 1;
    return deg;
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::gradPhi(int dim,
                int intFE,
                int i,
                vec_dbl_Type &p,
                vec_dbl_ptr_Type &value){
    if (dim==2) {
        switch (intFE) {
            case 0://P0
                switch (i) {
                    case 0:
                        value->at(0)= 0.;
                        value->at(1)= 0.;
                        break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        value->at(0)= -1.;
                        value->at(1)= -1.;
                        break;
                    case 1:
                        value->at(0)= 1.;
                        value->at(1)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 1.;
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        value->at(0)= 1. - 4.*(1 - p[0] - p[1]);
                        value->at(1)= 1. - 4.*(1 - p[0] - p[1]);
                        break;
                    case 1:
                        value->at(0)= 4.*p[0] - 1;
                        value->at(1)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 4.*p[1] - 1;
                        break;
                    case 3:
                        value->at(0)= 4 * (1. - 2*p[0] - p[1]);
                        value->at(1)= -4 * p[0];
                        break;
                    case 4:
                        value->at(0)= 4.*p[1];
                        value->at(1)= 4.*p[0];
                        break;
                    case 5:
                        value->at(0)= - 4.*p[1];
                        value->at(1)= 4 * (1. - p[0] - 2*p[1]);
                        break;
                }
                break;
        }
    }
    else if(dim==3) {
        switch (intFE) {
            case 0://P0
                switch (i) {
                    case 0:
                    value->at(0)= 0.;
                    value->at(1)= 0.;
                    value->at(2)= 0.;
                    break;
                }
                break;
            case 1://P1
                switch (i) {
                    case 0:
                        value->at(0)= -1.;
                        value->at(1)= -1.;
                        value->at(2)= -1.;
                        break;
                    case 1:
                        value->at(0)= 1.;
                        value->at(1)= 0.;
                        value->at(2)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 1.;
                        value->at(2)= 0.;
                        break;
                    case 3:
                        value->at(0)= 0.;
                        value->at(1)= 0.;
                        value->at(2)= 1.;
                        break;
                }
                break;
            case 2://P2
                switch (i) {
                    case 0:
                        value->at(0)= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
                        value->at(1)= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
                        value->at(2)= -3. + 4.*p[0] + 4.*p[1] + 4.*p[2];
                        break;
                    case 1:
                        value->at(0)= 4.*p[0] - 1;
                        value->at(1)= 0.;
                        value->at(2)= 0.;
                        break;
                    case 2:
                        value->at(0)= 0.;
                        value->at(1)= 4.*p[1] - 1;
                        value->at(2)= 0.;
                        break;
                    case 3:
                        value->at(0)= 0.;
                        value->at(1)= 0.;
                        value->at(2)= 4.*p[2] - 1;
                        break;
                    case 4:
                        value->at(0)= 4. - 8.*p[0] - 4.*p[1] - 4.*p[2];
                        value->at(1)= - 4.*p[0];
                        value->at(2)= - 4.*p[0];
                        break;
                    case 5:
                        value->at(0)= 4.*p[1];
                        value->at(1)= 4.*p[0];
                        value->at(2)= 0.;
                        break;
                    case 6:
                        value->at(0)= - 4.*p[1];
                        value->at(1)= 4. - 4.*p[0] - 8.*p[1] - 4.*p[2];
                        value->at(2)= - 4.*p[1];
                        break;
                    case 7:
                        value->at(0)= - 4.*p[2];
                        value->at(1)= - 4.*p[2];
                        value->at(2)= 4. - 4.*p[0] - 4.*p[1] - 8.*p[2];
                        break;
                    case 8:
                        value->at(0)= 4.*p[2];
                        value->at(1)= 0.;
                        value->at(2)= 4.*p[0];
                        break;
                    case 9:
                        value->at(0)= 0.;
                        value->at(1)= 4.*p[2];
                        value->at(2)= 4.*p[1];
                        break;
                }
                break;
            case 3: //Q1
            {
                double a = 1./8;
                switch (i) {
                    case 0:
                        value->at(0) = - a * ( 1 - p[1] ) * ( 1 - p[2] );
                        value->at(1) = - a * ( 1 - p[0] ) * ( 1 - p[2] );
                        value->at(2) = - a * ( 1 - p[0] ) * ( 1 - p[1] );
                        break;
                    case 1:
                        value->at(0) = a * ( 1 - p[1] ) * ( 1 - p[2] );
                        value->at(1) = - a * ( 1 + p[0] ) * ( 1 - p[2] );
                        value->at(2) = - a * ( 1 + p[0] ) * ( 1 - p[1] );
                        break;
                    case 2:
                        value->at(0) = a * ( 1 + p[1] ) * ( 1 - p[2] );
                        value->at(1) = a * ( 1 + p[0] ) * ( 1 - p[2] );
                        value->at(2) = - a * ( 1 + p[0] ) * ( 1 + p[1] );
                        break;
                    case 3:
                        value->at(0) = - a * ( 1 + p[1] ) * ( 1 - p[2] );
                        value->at(1) = a * ( 1 - p[0] ) * ( 1 - p[2] );
                        value->at(2) = - a * ( 1 - p[0] ) * ( 1 + p[1] );
                        break;
                    case 4:
                        value->at(0) = - a * ( 1 - p[1] ) * ( 1 + p[2] );
                        value->at(1) = - a * ( 1 - p[0] ) * ( 1 + p[2] );
                        value->at(2) = a * ( 1 - p[0] ) * ( 1 - p[1] );
                        break;
                    case 5:
                        value->at(0) = a * ( 1 - p[1] ) * ( 1 + p[2] );
                        value->at(1) = - a * ( 1 + p[0] ) * ( 1 + p[2] );
                        value->at(2) = a * ( 1 + p[0] ) * ( 1 - p[1] );
                        break;
                    case 6:
                        value->at(0) = a * ( 1 + p[1] ) * ( 1 + p[2] );
                        value->at(1) = a * ( 1 + p[0] ) * ( 1 + p[2] );
                        value->at(2) = a * ( 1 + p[0] ) * ( 1 + p[1] );
                        break;
                    case 7:
                        value->at(0) = - a * ( 1 + p[1] ) * ( 1 + p[2] );
                        value->at(1) = a * ( 1 - p[0] ) * ( 1 + p[2] );
                        value->at(2) = a * ( 1 - p[0] ) * ( 1 + p[1] );
                        break;
                    default:
                        break;
                }
                break;
            }
            case 4: //Q2
            {
                switch (i) {
                    case 0:
                        value->at(0) = 0.125*2*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[2] + 0.125*p[2]*p[2]*p[1] + -0.125*2*p[0]*p[1]*p[1]*p[2] + -0.125*2*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2] + -0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.125*p[0]*p[0]*p[2] + 0.125*2*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0] + -0.125*p[0]*p[0]*2*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2] + -0.125*2*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.125*p[0]*p[0]*p[1] + 0.125*p[1]*p[1]*p[0] + 0.125*2*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1] + -0.125*p[0]*p[0]*2*p[2]*p[1] + -0.125*p[1]*p[1]*2*p[2]*p[0] + -0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 1:
                        value->at(0) = 0.125*2*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[2] + -0.125*p[2]*p[2]*p[1] + -0.125*2*p[0]*p[1]*p[1]*p[2] + -0.125*2*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2] + 0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.125*p[0]*p[0]*p[2] + -0.125*2*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0] + -0.125*p[0]*p[0]*2*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2] + 0.125*2*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.125*p[0]*p[0]*p[1] + -0.125*p[1]*p[1]*p[0] + -0.125*2*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1] + -0.125*p[0]*p[0]*2*p[2]*p[1] + 0.125*p[1]*p[1]*2*p[2]*p[0] + 0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 2:
                        value->at(0) = -0.125*2*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[2] + 0.125*p[2]*p[2]*p[1] + -0.125*2*p[0]*p[1]*p[1]*p[2] + 0.125*2*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2] + -0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.125*p[0]*p[0]*p[2] + -0.125*2*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0] + -0.125*p[0]*p[0]*2*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2] + 0.125*2*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.125*p[0]*p[0]*p[1] + -0.125*p[1]*p[1]*p[0] + 0.125*2*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1] + 0.125*p[0]*p[0]*2*p[2]*p[1] + 0.125*p[1]*p[1]*2*p[2]*p[0] + -0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 3:
                        value->at(0) = -0.125*2*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[2] + -0.125*p[2]*p[2]*p[1] + -0.125*2*p[0]*p[1]*p[1]*p[2] + 0.125*2*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2] + 0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.125*p[0]*p[0]*p[2] + 0.125*2*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0] + -0.125*p[0]*p[0]*2*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2] + -0.125*2*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.125*p[0]*p[0]*p[1] + 0.125*p[1]*p[1]*p[0] + -0.125*2*p[2]*p[0]*p[1] + -0.125*p[0]*p[0]*p[1]*p[1] + 0.125*p[0]*p[0]*2*p[2]*p[1] + -0.125*p[1]*p[1]*2*p[2]*p[0] + 0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 4:
                        value->at(0) = -0.125*2*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[2] + 0.125*p[2]*p[2]*p[1] + 0.125*2*p[0]*p[1]*p[1]*p[2] + -0.125*2*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2] + 0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.125*p[0]*p[0]*p[2] + -0.125*2*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0] + 0.125*p[0]*p[0]*2*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2] + -0.125*2*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.125*p[0]*p[0]*p[1] + -0.125*p[1]*p[1]*p[0] + 0.125*2*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1] + -0.125*p[0]*p[0]*2*p[2]*p[1] + -0.125*p[1]*p[1]*2*p[2]*p[0] + 0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 5:
                        value->at(0) = -0.125*2*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[2] + -0.125*p[2]*p[2]*p[1] + 0.125*2*p[0]*p[1]*p[1]*p[2] + -0.125*2*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2] + -0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.125*p[0]*p[0]*p[2] + 0.125*2*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0] + 0.125*p[0]*p[0]*2*p[1]*p[2] + -0.125*p[0]*p[0]*p[2]*p[2] + 0.125*2*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.125*p[0]*p[0]*p[1] + 0.125*p[1]*p[1]*p[0] + -0.125*2*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1] + -0.125*p[0]*p[0]*2*p[2]*p[1] + 0.125*p[1]*p[1]*2*p[2]*p[0] + -0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 6:
                        value->at(0) = 0.125*2*p[0]*p[1]*p[2] + 0.125*p[1]*p[1]*p[2] + 0.125*p[2]*p[2]*p[1] + 0.125*2*p[0]*p[1]*p[1]*p[2] + 0.125*2*p[0]*p[2]*p[2]*p[1] + 0.125*p[1]*p[1]*p[2]*p[2] + 0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.125*p[0]*p[0]*p[2] + 0.125*2*p[1]*p[0]*p[2] + 0.125*p[2]*p[2]*p[0] + 0.125*p[0]*p[0]*2*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2] + 0.125*2*p[1]*p[2]*p[2]*p[0] + 0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.125*p[0]*p[0]*p[1] + 0.125*p[1]*p[1]*p[0] + 0.125*2*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1] + 0.125*p[0]*p[0]*2*p[2]*p[1] + 0.125*p[1]*p[1]*2*p[2]*p[0] + 0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 7:
                        value->at(0) = 0.125*2*p[0]*p[1]*p[2] + -0.125*p[1]*p[1]*p[2] + -0.125*p[2]*p[2]*p[1] + 0.125*2*p[0]*p[1]*p[1]*p[2] + 0.125*2*p[0]*p[2]*p[2]*p[1] + -0.125*p[1]*p[1]*p[2]*p[2] + -0.125*p[1]*p[2] + 0.125*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.125*p[0]*p[0]*p[2] + -0.125*2*p[1]*p[0]*p[2] + -0.125*p[2]*p[2]*p[0] + 0.125*p[0]*p[0]*2*p[1]*p[2] + 0.125*p[0]*p[0]*p[2]*p[2] + -0.125*2*p[1]*p[2]*p[2]*p[0] + -0.125*p[0]*p[2] + 0.125*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.125*p[0]*p[0]*p[1] + -0.125*p[1]*p[1]*p[0] + -0.125*2*p[2]*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1] + 0.125*p[0]*p[0]*2*p[2]*p[1] + -0.125*p[1]*p[1]*2*p[2]*p[0] + -0.125*p[0]*p[1] + 0.125*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 8:
                        value->at(0) = -0.250*2*p[0]*p[1]*p[2] + 0.250*2*p[0]*p[1]*p[1]*p[2] + 0.250*2*p[0]*p[2]*p[2]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.250*p[2] + -0.250*2*p[1]*p[2] + -0.250*p[2]*p[2] + -0.250*p[0]*p[0]*p[2] + 0.250*p[0]*p[0]*2*p[1]*p[2] + 0.250*p[0]*p[0]*p[2]*p[2] + 0.250*2*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.250*p[1] + -0.250*p[1]*p[1] + -0.250*2*p[2]*p[1] + -0.250*p[0]*p[0]*p[1] + 0.250*p[0]*p[0]*p[1]*p[1] + 0.250*p[0]*p[0]*2*p[2]*p[1] + 0.250*p[1]*p[1]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 9:
                        value->at(0) = -0.250*p[2] + -0.250*2*p[0]*p[2] + 0.250*p[2]*p[2] + 0.250*p[1]*p[1]*p[2] + 0.250*2*p[0]*p[1]*p[1]*p[2] + -0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[2]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.250*p[0]*p[0]*p[2] + 0.250*2*p[1]*p[0]*p[2] + 0.250*p[0]*p[0]*2*p[1]*p[2] + -0.250*2*p[1]*p[2]*p[2]*p[0] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.250*p[0] + -0.250*p[0]*p[0] + 0.250*2*p[2]*p[0] + 0.250*p[1]*p[1]*p[0] + 0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[1]*p[1]*2*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 10:
                        value->at(0) = 0.250*2*p[0]*p[1]*p[2] + 0.250*2*p[0]*p[1]*p[1]*p[2] + -0.250*2*p[0]*p[2]*p[2]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.250*p[2] + -0.250*2*p[1]*p[2] + 0.250*p[2]*p[2] + 0.250*p[0]*p[0]*p[2] + 0.250*p[0]*p[0]*2*p[1]*p[2] + -0.250*p[0]*p[0]*p[2]*p[2] + 0.250*2*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.250*p[1] + -0.250*p[1]*p[1] + 0.250*2*p[2]*p[1] + 0.250*p[0]*p[0]*p[1] + 0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[0]*p[0]*2*p[2]*p[1] + 0.250*p[1]*p[1]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 11:
                        value->at(0) = 0.250*p[2] + -0.250*2*p[0]*p[2] + -0.250*p[2]*p[2] + -0.250*p[1]*p[1]*p[2] + 0.250*2*p[0]*p[1]*p[1]*p[2] + 0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[2]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.250*p[0]*p[0]*p[2] + -0.250*2*p[1]*p[0]*p[2] + 0.250*p[0]*p[0]*2*p[1]*p[2] + 0.250*2*p[1]*p[2]*p[2]*p[0] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.250*p[0] + -0.250*p[0]*p[0] + -0.250*2*p[2]*p[0] + -0.250*p[1]*p[1]*p[0] + 0.250*p[0]*p[0]*p[1]*p[1] + 0.250*p[1]*p[1]*2*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 12:
                        value->at(0) = 0.250*2*p[0]*p[1]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2] + 0.250*2*p[0]*p[2]*p[2]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.250*p[2] + 0.250*2*p[1]*p[2] + -0.250*p[2]*p[2] + 0.250*p[0]*p[0]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2] + 0.250*p[0]*p[0]*p[2]*p[2] + 0.250*2*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.250*p[1] + 0.250*p[1]*p[1] + -0.250*2*p[2]*p[1] + 0.250*p[0]*p[0]*p[1] + -0.250*p[0]*p[0]*p[1]*p[1] + 0.250*p[0]*p[0]*2*p[2]*p[1] + 0.250*p[1]*p[1]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 13:
                        value->at(0) = 0.250*p[2] + 0.250*2*p[0]*p[2] + 0.250*p[2]*p[2] + -0.250*p[1]*p[1]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2] + -0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[2]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.250*p[0]*p[0]*p[2] + -0.250*2*p[1]*p[0]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2] + -0.250*2*p[1]*p[2]*p[2]*p[0] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.250*p[0] + 0.250*p[0]*p[0] + 0.250*2*p[2]*p[0] + -0.250*p[1]*p[1]*p[0] + -0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[1]*p[1]*2*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 14:
                        value->at(0) = -0.250*2*p[0]*p[1]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2] + -0.250*2*p[0]*p[2]*p[2]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.250*p[2] + 0.250*2*p[1]*p[2] + 0.250*p[2]*p[2] + -0.250*p[0]*p[0]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2] + -0.250*p[0]*p[0]*p[2]*p[2] + 0.250*2*p[1]*p[2]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.250*p[1] + 0.250*p[1]*p[1] + 0.250*2*p[2]*p[1] + -0.250*p[0]*p[0]*p[1] + -0.250*p[0]*p[0]*p[1]*p[1] + -0.250*p[0]*p[0]*2*p[2]*p[1] + 0.250*p[1]*p[1]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 15:
                        value->at(0) = -0.250*p[2] + 0.250*2*p[0]*p[2] + -0.250*p[2]*p[2] + 0.250*p[1]*p[1]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2] + 0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[2]*p[2] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.250*p[0]*p[0]*p[2] + 0.250*2*p[1]*p[0]*p[2] + -0.250*p[0]*p[0]*2*p[1]*p[2] + 0.250*2*p[1]*p[2]*p[2]*p[0] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.250*p[0] + 0.250*p[0]*p[0] + -0.250*2*p[2]*p[0] + 0.250*p[1]*p[1]*p[0] + -0.250*p[0]*p[0]*p[1]*p[1] + 0.250*p[1]*p[1]*2*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[2] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 16:
                        value->at(0) = 0.250*p[1] + -0.250*2*p[0]*p[1] + -0.250*p[1]*p[1] + -0.250*p[2]*p[2]*p[1] + 0.250*2*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[1]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.250*p[0] + -0.250*p[0]*p[0] + -0.250*2*p[1]*p[0] + -0.250*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[2]*p[2] + 0.250*2*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[1] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.250*2*p[2]*p[0]*p[1] + 0.250*p[0]*p[0]*2*p[2]*p[1] + 0.250*p[1]*p[1]*2*p[2]*p[0] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 17:
                        value->at(0) = -0.250*p[1] + -0.250*2*p[0]*p[1] + 0.250*p[1]*p[1] + 0.250*p[2]*p[2]*p[1] + 0.250*2*p[0]*p[2]*p[2]*p[1] + -0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[1]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.250*p[0] + -0.250*p[0]*p[0] + 0.250*2*p[1]*p[0] + 0.250*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*p[2]*p[2] + -0.250*2*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[1] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.250*2*p[2]*p[0]*p[1] + 0.250*p[0]*p[0]*2*p[2]*p[1] + -0.250*p[1]*p[1]*2*p[2]*p[0] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 18:
                        value->at(0) = 0.250*p[1] + 0.250*2*p[0]*p[1] + 0.250*p[1]*p[1] + -0.250*p[2]*p[2]*p[1] + -0.250*2*p[0]*p[2]*p[2]*p[1] + -0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[1]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.250*p[0] + 0.250*p[0]*p[0] + 0.250*2*p[1]*p[0] + -0.250*p[2]*p[2]*p[0] + -0.250*p[0]*p[0]*p[2]*p[2] + -0.250*2*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[1] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.250*2*p[2]*p[0]*p[1] + -0.250*p[0]*p[0]*2*p[2]*p[1] + -0.250*p[1]*p[1]*2*p[2]*p[0] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 19:
                        value->at(0) = -0.250*p[1] + 0.250*2*p[0]*p[1] + -0.250*p[1]*p[1] + 0.250*p[2]*p[2]*p[1] + -0.250*2*p[0]*p[2]*p[2]*p[1] + 0.250*p[1]*p[1]*p[2]*p[2] + 0.250*2*p[0]*p[1]*p[1] + -0.250*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.250*p[0] + 0.250*p[0]*p[0] + -0.250*2*p[1]*p[0] + 0.250*p[2]*p[2]*p[0] + -0.250*p[0]*p[0]*p[2]*p[2] + 0.250*2*p[1]*p[2]*p[2]*p[0] + 0.250*p[0]*p[0]*2*p[1] + -0.250*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.250*2*p[2]*p[0]*p[1] + -0.250*p[0]*p[0]*2*p[2]*p[1] + 0.250*p[1]*p[1]*2*p[2]*p[0] + -0.250*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 20:
                        value->at(0) = 0.500*2*p[0]*p[1] + -0.500*2*p[0]*p[2]*p[2]*p[1] + -0.500*2*p[0]*p[1]*p[1] + 0.500*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.500 + 0.500*2*p[1] + 0.500*p[0]*p[0] + 0.500*p[2]*p[2] + -0.500*p[0]*p[0]*p[2]*p[2] + -0.500*p[0]*p[0]*2*p[1] + -0.500*2*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.500*2*p[2]*p[1] + -0.500*p[0]*p[0]*2*p[2]*p[1] + -0.500*p[1]*p[1]*2*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 21:
                        value->at(0) = 0.500 + 0.500*2*p[0] + -0.500*p[1]*p[1] + -0.500*p[2]*p[2] + 0.500*p[1]*p[1]*p[2]*p[2] + -0.500*2*p[0]*p[1]*p[1] + -0.500*2*p[0]*p[2]*p[2] + 0.500*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.500*2*p[1]*p[0] + 0.500*2*p[1]*p[2]*p[2]*p[0] + -0.500*p[0]*p[0]*2*p[1] + 0.500*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.500*2*p[2]*p[0] + 0.500*p[1]*p[1]*2*p[2]*p[0] + -0.500*p[0]*p[0]*2*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 22:
                        value->at(0) = -0.500*2*p[0]*p[1] + 0.500*2*p[0]*p[2]*p[2]*p[1] + -0.500*2*p[0]*p[1]*p[1] + 0.500*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.500 + 0.500*2*p[1] + -0.500*p[0]*p[0] + -0.500*p[2]*p[2] + 0.500*p[0]*p[0]*p[2]*p[2] + -0.500*p[0]*p[0]*2*p[1] + -0.500*2*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.500*2*p[2]*p[1] + 0.500*p[0]*p[0]*2*p[2]*p[1] + -0.500*p[1]*p[1]*2*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 23:
                        value->at(0) = -0.500 + 0.500*2*p[0] + 0.500*p[1]*p[1] + 0.500*p[2]*p[2] + -0.500*p[1]*p[1]*p[2]*p[2] + -0.500*2*p[0]*p[1]*p[1] + -0.500*2*p[0]*p[2]*p[2] + 0.500*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.500*2*p[1]*p[0] + -0.500*2*p[1]*p[2]*p[2]*p[0] + -0.500*p[0]*p[0]*2*p[1] + 0.500*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.500*2*p[2]*p[0] + -0.500*p[1]*p[1]*2*p[2]*p[0] + -0.500*p[0]*p[0]*2*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 24:
                        value->at(0) = 0.500*2*p[0]*p[2] + -0.500*2*p[0]*p[1]*p[1]*p[2] + -0.500*2*p[0]*p[2]*p[2] + 0.500*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = 0.500*p[0]*p[0]*p[2] + 0.500*2*p[1]*p[2] + -0.500*p[0]*p[0]*2*p[1]*p[2] + -0.500*2*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -0.500 + 0.500*2*p[2] + 0.500*p[0]*p[0] + 0.500*p[1]*p[1] + -0.500*p[0]*p[0]*p[1]*p[1] + -0.500*p[0]*p[0]*2*p[2] + -0.500*p[1]*p[1]*2*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 25:
                        value->at(0) = -1.000*2*p[0] + 1.000*2*p[0]*p[1]*p[1] + 1.000*2*p[0]*p[2]*p[2] + -1.000*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -1.000*2*p[1] + 1.000*p[0]*p[0]*2*p[1] + 1.000*2*p[1]*p[2]*p[2] + -1.000*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = -1.000*2*p[2] + 1.000*p[0]*p[0]*2*p[2] + 1.000*p[1]*p[1]*2*p[2] + -1.000*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    case 26:
                        value->at(0) = -0.500*2*p[0]*p[2] + 0.500*2*p[0]*p[1]*p[1]*p[2] + -0.500*2*p[0]*p[2]*p[2] + 0.500*2*p[0]*p[1]*p[1]*p[2]*p[2];
                        value->at(1) = -0.500*p[0]*p[0]*p[2] + -0.500*2*p[1]*p[2] + 0.500*p[0]*p[0]*2*p[1]*p[2] + -0.500*2*p[1]*p[2]*p[2] + 0.500*p[0]*p[0]*2*p[1]*p[2]*p[2];
                        value->at(2) = 0.500 + 0.500*2*p[2] + -0.500*p[0]*p[0] + -0.500*p[1]*p[1] + 0.500*p[0]*p[0]*p[1]*p[1] + -0.500*p[0]*p[0]*2*p[2] + -0.500*p[1]*p[1]*2*p[2] + 0.500*p[0]*p[0]*p[1]*p[1]*2*p[2];
                        break;
                    default:
                        break;
                }
                break;
            }
            case 5: //Q2-20
            {
                std::cout << "Warning! Q2-20 not working correct!" << std::endl;
                double a = 1./8;
                double b = 1./4;
                switch (i) {
                    case 0:
                        (*value)[0] = a * (p[1] - 1) * (p[2] - 1) * (2 * p[0] + p[1] + p[2] + 1);
                        (*value)[1] = a * (p[0] - 1) * (p[2] - 1) * (p[0] + 2 * p[1] + p[2] + 1);
                        (*value)[2] = a * (p[0] - 1) * (p[1] - 1) * (p[0] + p[1] + 2 * p[2] + 1);
                        break;
                    case 1:
                        (*value)[0] = - a * ( p[1] - 1 ) * ( p[2] - 1 ) * ( -2*p[0] + p[1] + p[2] + 1);
                        (*value)[1] = a * ( p[0] + 1 ) * ( p[2] - 1 ) * ( p[0] - 2*p[1] - p[2] - 1);
                        (*value)[2] = a * ( p[0] + 1 ) * ( p[1] - 1 ) * ( p[0] - p[1] - 2*p[2] - 1);
                        break;
                    case 2:
                        (*value)[0] = a * ( p[1] + 1 ) * ( p[2] - 1 ) * ( -2*p[0] + p[1] - p[2] - 1);
                        (*value)[1] = - a * ( p[0] - 1 ) * ( p[2] - 1 ) * ( p[0] - 2*p[1] + p[2] + 1);
                        (*value)[2] = - a * ( p[0] + 1 ) * ( p[1] + 1 ) * ( p[0] + p[1] - 2*p[2] - 1);
                        break;
                    case 3:
                        (*value)[0] = a * (1 + p[1]) * (-1 - 2 * p[0] + p[1] - p[2]) * (-1 + p[2]);
                        (*value)[1] = - a * ((-1 + p[0]) * (-1 + p[2]) * (1 + p[0] - 2 * p[1] + p[2]));
                        (*value)[2] = - a * ((-1 + p[0]) * (1 + p[1]) * (1 + p[0] - p[1] + 2 * p[2]));
                        break;
                    case 4:
                        (*value)[0] = -a*(p[1] - 1) * (p[2] + 1) * (2 * p[0] + p[1] - p[2] + 1);
                        (*value)[1] = -a*(p[0] - 1) * (p[2] + 1) * (p[0] + 2 * p[1] - p[2] + 1);
                        (*value)[2] = -a*(p[0] - 1) * (p[1] - 1) * (p[0] + p[1] - 2 * p[2] + 1);
                        break;
                    case 5:
                        (*value)[0] = a*(p[1] - 1) * (p[2] + 1) * (-2 * p[0] + p[1] - p[2] + 1);
                        (*value)[1] = -a*(p[0] + 1) * (p[2] + 1) * (p[0] - 2 * p[1] + p[2] - 1);
                        (*value)[2] = -a*(p[0] + 1) * (p[1] - 1) * (p[0] - p[1] + 2 * p[2] - 1);
                        break;
                    case 6:
                        (*value)[0] = a*(p[1] + 1) * (p[2] + 1) * (2 * p[0] + p[1] + p[2] - 1);
                        (*value)[1] = a*(p[0] + 1) * (p[2] + 1) * (p[0] + 2 * p[1] + p[2] - 1);
                        (*value)[2] = a*(p[0] + 1) * (p[1] + 1) * (p[0] + p[1] + 2 * p[2] - 1);
                        break;
                    case 7:
                        (*value)[0] = -a*(p[1] + 1) * (p[2] + 1) * (-2 * p[0] + p[1] + p[2] - 1);
                        (*value)[1] = a*(p[0] - 1) * (p[2] + 1) * (p[0] - 2 * p[1] - p[2] + 1);
                        (*value)[2] = a*(p[0] - 1) * (p[1] + 1) * (p[0] - p[1] - 2 * p[2] + 1);
                        break;
                    case 8:
                        (*value)[0] = -2 * b * p[0] * (p[1] - 1) * (p[2] - 1);
                        (*value)[1] = -b * (p[0]*p[0] - 1) * (p[2] - 1);
                        (*value)[2] = -b * (p[0]*p[0] - 1) * (p[1] - 1);
                        break;
                    case 9:
                        (*value)[0] = b * (p[1]*p[1] - 1) * (p[2] - 1);
                        (*value)[1] = 2 * b * (p[0] + 1) * p[1] * (p[2] - 1);
                        (*value)[2] = b * (p[0] + 1) * (p[1]*p[1] - 1);
                        break;
                    case 10:
                        (*value)[0] = 2 * b * p[0] * (p[1] + 1) * (p[2] - 1);
                        (*value)[1] = b * (p[0]*p[0] - 1) * (p[2] - 1);
                        (*value)[2] = b * (p[0]*p[0] - 1) * (p[1] + 1);
                        break;
                    case 11:
                        (*value)[0] = -b * (p[1]*p[1] - 1) * (p[2] - 1);
                        (*value)[1] = -2 * b * (p[0] - 1) * p[1] * (p[2] - 1);
                        (*value)[2] = -b * (p[0] - 1) * (p[1]*p[1] - 1);
                        break;
                    case 12:
                        (*value)[0] = 2* b * p[0] * (p[1] - 1) * (p[2] + 1);
                        (*value)[1] = b * (p[0]*p[0] - 1) * (p[2] + 1);
                        (*value)[2] = b * (p[0]*p[0] - 1) * (p[1] - 1);
                        break;
                    case 13:
                        (*value)[0] = -b * (p[1]*p[1] - 1) * (p[2] + 1);
                        (*value)[1] = -2 * b * (p[0] + 1) * p[1] * (p[2] + 1);
                        (*value)[2] = -b * (p[0] + 1) * (p[1]*p[1] - 1);
                        break;
                    case 14:
                        (*value)[0] = -2 * b * p[0] * (p[1] + 1) * (p[2] + 1);
                        (*value)[1] = -b * (p[0]*p[0] - 1) * (p[2] + 1);
                        (*value)[2] = -b *(p[0]*p[0] - 1) * (p[1] + 1);
                        break;
                    case 15:
                        (*value)[0] = b * (p[1]*p[1] - 1) * (p[2] + 1);
                        (*value)[1] = 2 * b * (p[0] - 1) * p[1] * (p[2] + 1);
                        (*value)[2] = b * (p[0] - 1) * (p[1]*p[1] - 1);
                        break;
                    case 16:
                        (*value)[0] = -b * (p[1] - 1) * (p[2]*p[2] - 1);
                        (*value)[1] = -b * (p[0] - 1) * (p[2]*p[2] - 1);
                        (*value)[2] = -2 * b * (p[0] - 1) * (p[1] - 1) * p[2];
                        break;
                    case 17:
                        (*value)[0] = b * (p[1] - 1) * (p[2]*p[2] - 1);
                        (*value)[1] = b * (p[0] + 1) * (p[2]*p[2] - 1);
                        (*value)[2] = 2 * b * (p[0] + 1) * (p[1] - 1) * p[2];
                        break;
                    case 18:
                        (*value)[0] = -b * (p[1] + 1) * (p[2]*p[2] - 1);
                        (*value)[1] = -b * (p[0] + 1) * (p[2]*p[2] - 1);
                        (*value)[0] = -2 * b * (p[0] + 1) * (p[1] + 1) * p[2];
                        break;
                    case 19:
                        (*value)[0] = b * (p[1] + 1) * (p[2]*p[2] - 1);
                        (*value)[1] = b * (p[0] - 1) * (p[2]*p[2] - 1);
                        (*value)[2] = 2 * b * (p[0] - 1) * (p[1] + 1) * p[2];
                        break;
                    default:
                        break;
                }
                break;
            }
        }
    }
}

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::getQuadratureValues(int dim,
                                          int Degree,
                                          vec2D_dbl_ptr_Type &QuadPts,
                                          vec_dbl_ptr_Type &QuadW,
                                          std::string FEType){
    double a, b, c, P1, P2;

    double b1,b2,c1,c2,d,e,f,g,h,i,j;
    if (dim==1){
        // points are for interval [0,1]
        TEUCHOS_TEST_FOR_EXCEPTION(Degree>2, std::runtime_error, "Quadrature rule in 1d only up to degree 3.");
        switch (Degree) {
            case 0:
                QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(1,0.0)));
                QuadW->resize(1);
                QuadPts->at(0).at(0) = 0.5;
                QuadW->at(0) = 1.;
                break;
            case 1:
                QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(1,0.0)));
                QuadW->resize(1);
                QuadPts->at(0).at(0) = 0.5;
                QuadW->at(0) = 1.;
                break;
            case 2:
                QuadPts.reset(new vec2D_dbl_Type(2,vec_dbl_Type(1,0.0)));
                QuadW->resize(2);
                QuadPts->at(0).at(0) = - 0.5/sqrt(3.)+0.5;
                QuadPts->at(1).at(0) = 0.5/sqrt(3.)+0.5;
                QuadW->at(0) = .5;
                QuadW->at(1) = .5;
                break;
            case 3:
                QuadPts.reset(new vec2D_dbl_Type(2,vec_dbl_Type(1,0.0)));
                QuadW->resize(2);
                QuadPts->at(0).at(0) = - 0.5/sqrt(3.)+0.5;
                QuadPts->at(1).at(0) = 0.5/sqrt(3.)+0.5;
                QuadW->at(0) = .5;
                QuadW->at(1) = .5;
                break;
            default:
                break;
        }
    }
    if (dim==2) {

        TEUCHOS_TEST_FOR_EXCEPTION(Degree>7, std::runtime_error, "Quadrature rule in 2d only up to degree 7.");
        if (Degree==3 || Degree==4)
            Degree=5;

        if (Degree==6)
            Degree=7;
        switch (Degree) {
            case 1:

                QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(2,0.0)));
                QuadW->resize(1);
                QuadPts->at(0).at(0) = 1/3.;
                QuadPts->at(0).at(1) = 1/3.;
                QuadW->at(0)	= 1/2.;
                break;

            case 2:

                QuadPts.reset(new vec2D_dbl_Type(3,vec_dbl_Type(2,0.0)));
                QuadW->resize(3);
                a = 1/6.;
                QuadPts->at(0).at(0) 	= 0.5;
                QuadPts->at(0).at(1)    = 0.5;

                QuadPts->at(1).at(0) 	= 0.;
                QuadPts->at(1).at(1) 	= 0.5;

                QuadPts->at(2).at(0) 	= 0.5;
                QuadPts->at(2).at(1) 	= 0.;

                QuadW->at(0) 		= a;
                QuadW->at(1)            = a;
                QuadW->at(2)            = a;

                break;

            case 5:
                QuadPts.reset(new vec2D_dbl_Type(7,vec_dbl_Type(2,0.0)));
                QuadW->resize(7);
                a = 0.470142064105115;
                b = 0.101286507323456;
                P1 = 0.066197076394253;
                P2 = 0.062969590272413;

                QuadPts->at(0).at(0) 	= 1/3.;
                QuadPts->at(0).at(1)    = 1/3.;

                QuadPts->at(1).at(0) 	= a;
                QuadPts->at(1).at(1) 	= a;

                QuadPts->at(2).at(0) 	= 1-2.*a;
                QuadPts->at(2).at(1) 	= a;

                QuadPts->at(3).at(0) 	= a;
                QuadPts->at(3).at(1) 	= 1-2.*a;

                QuadPts->at(4).at(0) 	= b;
                QuadPts->at(4).at(1) 	= b;

                QuadPts->at(5).at(0) 	= 1-2.*b;
                QuadPts->at(5).at(1) 	= b;

                QuadPts->at(6).at(0) 	= b;
                QuadPts->at(6).at(1) 	= 1-2.*b;

                QuadW->at(0) 			= 9/80.;
                QuadW->at(1)            = P1;
                QuadW->at(2)            = P1;
                QuadW->at(3) 			= P1;
                QuadW->at(4)            = P2;
                QuadW->at(5)            = P2;
                QuadW->at(6)            = P2;

                break;
            case 7:
                // 28 Punkte
                
                QuadPts.reset(new vec2D_dbl_Type(28,vec_dbl_Type(2,0.0)));
                QuadW.reset(new vec_dbl_Type(28,0.0));
                
                // x punkt
                QuadPts->at(0).at(0) = 0.777777777777778;
                QuadPts->at(1).at(0) = 0.111111111111111;
                QuadPts->at(2).at(0) = 0.111111111111111;
                QuadPts->at(3).at(0) = 0.666666666666667;
                QuadPts->at(4).at(0) = 0.222222222222222;
                QuadPts->at(5).at(0) = 0.111111111111111;
                QuadPts->at(6).at(0) = 0.222222222222222;
                QuadPts->at(7).at(0) = 0.111111111111111;
                QuadPts->at(8).at(0) = 0.666666666666667;
                QuadPts->at(9).at(0) = 0.555555555555556;
                QuadPts->at(10).at(0) = 0.333333333333333;
                QuadPts->at(11).at(0) = 0.111111111111111;
                QuadPts->at(12).at(0) = 0.333333333333333;
                QuadPts->at(13).at(0) = 0.111111111111111;
                QuadPts->at(14).at(0) = 0.555555555555556;
                QuadPts->at(15).at(0) = 0.555555555555556;
                QuadPts->at(16).at(0) = 0.222222222222222;
                QuadPts->at(17).at(0) = 0.222222222222222;
                QuadPts->at(18).at(0) = 0.444444444444444;
                QuadPts->at(19).at(0) = 0.444444444444444;
                QuadPts->at(20).at(0) = 0.111111111111111;
                QuadPts->at(21).at(0) = 0.444444444444444;
                QuadPts->at(22).at(0) = 0.333333333333333;
                QuadPts->at(23).at(0) = 0.222222222222222;
                QuadPts->at(24).at(0) = 0.333333333333333;
                QuadPts->at(25).at(0) = 0.222222222222222;
                QuadPts->at(26).at(0) = 0.444444444444444;
                QuadPts->at(27).at(0) = 0.333333333333333;
                
                // y punkt
                QuadPts->at(0).at(1) = 0.111111111111111;
                QuadPts->at(1).at(1) = 0.111111111111111;
                QuadPts->at(2).at(1) = 0.777777777777778;
                QuadPts->at(3).at(1) = 0.222222222222222;
                QuadPts->at(4).at(1) = 0.111111111111111;
                QuadPts->at(5).at(1) = 0.666666666666667;
                QuadPts->at(6).at(1) = 0.666666666666667;
                QuadPts->at(7).at(1) = 0.222222222222222;
                QuadPts->at(8).at(1) = 0.111111111111111;
                QuadPts->at(9).at(1) = 0.333333333333333;
                QuadPts->at(10).at(1) = 0.111111111111111;
                QuadPts->at(11).at(1) = 0.555555555555556;
                QuadPts->at(12).at(1) = 0.555555555555556;
                QuadPts->at(13).at(1) = 0.333333333333333;
                QuadPts->at(14).at(1) = 0.111111111111111;
                QuadPts->at(15).at(1) = 0.222222222222222;
                QuadPts->at(16).at(1) = 0.222222222222222;
                QuadPts->at(17).at(1) = 0.555555555555556;
                QuadPts->at(18).at(1) = 0.444444444444444;
                QuadPts->at(19).at(1) = 0.111111111111111;
                QuadPts->at(20).at(1) = 0.444444444444444;
                QuadPts->at(21).at(1) = 0.333333333333333;
                QuadPts->at(22).at(1) = 0.222222222222222;
                QuadPts->at(23).at(1) = 0.444444444444444;
                QuadPts->at(24).at(1) = 0.444444444444444;
                QuadPts->at(25).at(1) = 0.333333333333333;
                QuadPts->at(26).at(1) = 0.222222222222222;
                QuadPts->at(27).at(1) = 0.333333333333333;
                
                // Gewichte
                QuadW->at(0) 			= 0.342410714285714/2.0;
                QuadW->at(1) 			= 0.342410714285714/2.0;
                QuadW->at(2) 			= 0.342410714285714/2.0;
                QuadW->at(3) 			= -0.561160714285714/2.0;
                QuadW->at(4) 			= -0.561160714285714/2.0;
                QuadW->at(5) 			= -0.561160714285714/2.0;
                QuadW->at(6) 			= -0.561160714285714/2.0;
                QuadW->at(7) 			= -0.561160714285714/2.0;
                QuadW->at(8) 			= -0.561160714285714/2.0;
                QuadW->at(9) 			= 1.295089285714286/2.0;
                QuadW->at(10) 			= 1.295089285714286/2.0;
                QuadW->at(11) 			= 1.295089285714286/2.0;
                QuadW->at(12) 			= 1.295089285714286/2.0;
                QuadW->at(13) 			= 1.295089285714286/2.0;
                QuadW->at(14) 			= 1.295089285714286/2.0;
                QuadW->at(15) 			= 0.172767857142857/2.0;
                QuadW->at(16) 			= 0.172767857142857/2.0;
                QuadW->at(17) 			= 0.172767857142857/2.0;
                QuadW->at(18) 			= -1.354910714285714/2.0;
                QuadW->at(19) 			= -1.354910714285714/2.0;
                QuadW->at(20) 			= -1.354910714285714/2.0;
                QuadW->at(21) 			= -0.408482142857143/2.0;
                QuadW->at(22) 			= -0.408482142857143/2.0;
                QuadW->at(23) 			= -0.408482142857143/2.0;
                QuadW->at(24) 			= -0.408482142857143/2.0;
                QuadW->at(25) 			= -0.408482142857143/2.0;
                QuadW->at(26) 			= -0.408482142857143/2.0;
                QuadW->at(27) 			= 1.566517857142857/2.0;
                
                break;
                
            }
    }
    else if(dim==3){
        if (FEType.at(0)=='P') {
            if (Degree==2)
                Degree=3;
            if (Degree==4)
                Degree=5;

            TEUCHOS_TEST_FOR_EXCEPTION(Degree>6, std::runtime_error, "Tetrahedron quadrature rules only up to degree 6 available.");
            
            switch (Degree) {
                case 1:
                    QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(3,0.0)));
                    QuadW->resize(1);
                    QuadPts->at(0).at(0) 	= 0.25;
                    QuadPts->at(0).at(1) 	= 0.25;
                    QuadPts->at(0).at(2) 	= 0.25;
                    QuadW->at(0)			= 1/6.;
                    break;
                    
                case 3:
                    QuadPts.reset(new vec2D_dbl_Type(5,vec_dbl_Type(3,0.0)));
                    QuadW->resize(5);
                    a = .25;
                    b = 1./6.;
                    c = .5;
                    QuadPts->at(0).at(0) = a;
                    QuadPts->at(0).at(1) = a;
                    QuadPts->at(0).at(2) = a;
                    
                    QuadPts->at(1).at(0) = b;
                    QuadPts->at(1).at(1) = b;
                    QuadPts->at(1).at(2) = b;
                    
                    QuadPts->at(2).at(0) = b;
                    QuadPts->at(2).at(1) = b;
                    QuadPts->at(2).at(2) = c;
                    
                    QuadPts->at(3).at(0) = b;
                    QuadPts->at(3).at(1) = c;
                    QuadPts->at(3).at(2) = b;
                    
                    QuadPts->at(4).at(0) = c;
                    QuadPts->at(4).at(1) = b;
                    QuadPts->at(4).at(2) = b;
                    
                    QuadW->at(0)		 = -2./15.;
                    QuadW->at(1)		 = 3./40.;
                    QuadW->at(2)		 = 3./40.;
                    QuadW->at(3)		 = 3./40.;
                    QuadW->at(4)		 = 3./40.;
                    break;
                case 4:
                    QuadPts.reset(new vec2D_dbl_Type(11,vec_dbl_Type(3,0.0)));
                    QuadW->resize(11);
                    
                    a = .785714285714286;
                    b = .071428571428571;
                    c = .100596423833201;
                    d = .399403576166799;
                    
                    QuadPts->at(0).at(0) 	= .25;
                    QuadPts->at(0).at(1)    = .25;
                    QuadPts->at(0).at(2)    = .25;
                    
                    QuadPts->at(1).at(0) 	= a;
                    QuadPts->at(1).at(1)    = b;
                    QuadPts->at(1).at(2)    = b;
                    
                    QuadPts->at(2).at(0) 	= b;
                    QuadPts->at(2).at(1)    = b;
                    QuadPts->at(2).at(2)    = b;
                    
                    QuadPts->at(3).at(0) 	= b;
                    QuadPts->at(3).at(1)    = b;
                    QuadPts->at(3).at(2)    = a;
                    
                    QuadPts->at(4).at(0) 	= b;
                    QuadPts->at(4).at(1)    = a;
                    QuadPts->at(4).at(2)    = b;
                    
                    QuadPts->at(5).at(0) 	= c;
                    QuadPts->at(5).at(1)    = d;
                    QuadPts->at(5).at(2)    = d;
                    
                    QuadPts->at(6).at(0) 	= d;
                    QuadPts->at(6).at(1)    = c;
                    QuadPts->at(6).at(2)    = d;
                    
                    QuadPts->at(7).at(0) 	= d;
                    QuadPts->at(7).at(1)    = d;
                    QuadPts->at(7).at(2)    = c;
                    
                    QuadPts->at(8).at(0) 	= d;
                    QuadPts->at(8).at(1)    = c;
                    QuadPts->at(8).at(2)    = c;
                    
                    QuadPts->at(9).at(0) 	= c;
                    QuadPts->at(9).at(1)    = d;
                    QuadPts->at(9).at(2)    = c;
                    
                    QuadPts->at(10).at(0) 	= c;
                    QuadPts->at(10).at(1)   = c;
                    QuadPts->at(10).at(2)   = d;
                    
                    a = -.078933333333333;
                    b = .045733333333333;
                    c= .149333333333333;
                    
                    
                    QuadW->at(0) = a;
                    
                    QuadW->at(1) = b;
                    QuadW->at(2) = b;
                    QuadW->at(3) = b;
                    QuadW->at(4) = b;
                    
                    QuadW->at(5) = c;
                    QuadW->at(6) = c;
                    QuadW->at(7) = c;
                    QuadW->at(8) = c;
                    QuadW->at(9) = c;
                    QuadW->at(10) = c;
                    
                case 5:
                    QuadPts.reset(new vec2D_dbl_Type(15,vec_dbl_Type(3,0.0)));
                    QuadW->resize(15);
                    a 	= 0.25;
                    b1 	= (7.+sqrt(15.))/34.;
                    b2 	= (7.-sqrt(15.))/34.;
                    c1 	= (13.-3.*sqrt(15.))/34.;
                    c2 	= (13.+3.*sqrt(15.))/34.;
                    d 	= (5.-sqrt(15.))/20.;
                    e 	= (5.+sqrt(15.))/20.;
                    
                    QuadPts->at(0).at(0) 	= a;
                    QuadPts->at(0).at(1)    = a;
                    QuadPts->at(0).at(2)    = a;
                    
                    QuadPts->at(1).at(0) 	= b1;
                    QuadPts->at(1).at(1)    = b1;
                    QuadPts->at(1).at(2)    = b1;
                    
                    QuadPts->at(2).at(0) 	= b1;
                    QuadPts->at(2).at(1)    = b1;
                    QuadPts->at(2).at(2)    = c1;
                    
                    QuadPts->at(3).at(0) 	= b1;
                    QuadPts->at(3).at(1)    = c1;
                    QuadPts->at(3).at(2)    = b1;
                    
                    QuadPts->at(4).at(0) 	= c1;
                    QuadPts->at(4).at(1)    = b1;
                    QuadPts->at(4).at(2)    = b1;
                    
                    QuadPts->at(5).at(0) 	= b2;
                    QuadPts->at(5).at(1)    = b2;
                    QuadPts->at(5).at(2)    = b2;
                    
                    QuadPts->at(6).at(0) 	= b2;
                    QuadPts->at(6).at(1)    = b2;
                    QuadPts->at(6).at(2)    = c2;
                    
                    QuadPts->at(7).at(0) 	= b2;
                    QuadPts->at(7).at(1)    = c2;
                    QuadPts->at(7).at(2)    = b2;
                    
                    QuadPts->at(8).at(0) 	= c2;
                    QuadPts->at(8).at(1)    = b2;
                    QuadPts->at(8).at(2)    = b2;
                    
                    QuadPts->at(9).at(0) 	= d;
                    QuadPts->at(9).at(1)    = d;
                    QuadPts->at(9).at(2)    = e;
                    
                    QuadPts->at(10).at(0) 	= d;
                    QuadPts->at(10).at(1)   = e;
                    QuadPts->at(10).at(2)   = d;
                    
                    QuadPts->at(11).at(0) 	= e;
                    QuadPts->at(11).at(1)	= d;
                    QuadPts->at(11).at(2)	= d;
                    
                    QuadPts->at(12).at(0) 	= d;
                    QuadPts->at(12).at(1)	= e;
                    QuadPts->at(12).at(2)	= e;
                    
                    QuadPts->at(13).at(0) 	= e;
                    QuadPts->at(13).at(1)	= d;
                    QuadPts->at(13).at(2)	= e;
                    
                    QuadPts->at(14).at(0) 	= e;
                    QuadPts->at(14).at(1)	= e;
                    QuadPts->at(14).at(2)	= d;
                    
                    
                    P1 	= (2665.-14.*sqrt(15.))/226800.;
                    P2 	= (2665.+14.*sqrt(15.))/226800.;
                    b	= 5./567.;
                    
                    QuadW->at(0) 			= 8./405.;
                    QuadW->at(1)            = P1;
                    QuadW->at(2)            = P1;
                    QuadW->at(3) 			= P1;
                    QuadW->at(4)            = P1;
                    
                    QuadW->at(5)            = P2;
                    QuadW->at(6)            = P2;
                    QuadW->at(7)            = P2;
                    QuadW->at(8)            = P2;
                    
                    QuadW->at(9) 			= b;
                    QuadW->at(10)           = b;
                    QuadW->at(11)           = b;
                    QuadW->at(12) 			= b;
                    QuadW->at(13)           = b;
                    QuadW->at(14)           = b;
                    
                    break;
                case 6: //Keast
                    QuadPts.reset(new vec2D_dbl_Type(24,vec_dbl_Type(3,0.0)));
                    QuadW->resize(24);
                    a = .356191386222545;
                    b = .214602871259152;
                    c = .877978124396166;
                    d = .040673958534611;
                    f = .032986329573173;
                    g = .322337890142276;
                    h = .269672331458316;
                    i = .063661001875018;
                    j = .603005664791649;
                    
                    QuadPts->at(0).at(0) 	= a;
                    QuadPts->at(0).at(1)    = b;
                    QuadPts->at(0).at(2)    = b;

                    QuadPts->at(1).at(0) 	= b;
                    QuadPts->at(1).at(1)    = b;
                    QuadPts->at(1).at(2)    = b;

                    QuadPts->at(2).at(0) 	= b;
                    QuadPts->at(2).at(1)    = b;
                    QuadPts->at(2).at(2)    = a;

                    QuadPts->at(3).at(0) 	= b;
                    QuadPts->at(3).at(1)    = a;
                    QuadPts->at(3).at(2)    = b;

                    QuadPts->at(4).at(0) 	= c;
                    QuadPts->at(4).at(1)    = d;
                    QuadPts->at(4).at(2)    = d;
                    
                    QuadPts->at(5).at(0) 	= d;
                    QuadPts->at(5).at(1)    = d;
                    QuadPts->at(5).at(2)    = d;

                    QuadPts->at(6).at(0) 	= d;
                    QuadPts->at(6).at(1)    = d;
                    QuadPts->at(6).at(2)    = c;
                    
                    QuadPts->at(7).at(0) 	= d;
                    QuadPts->at(7).at(1)    = c;
                    QuadPts->at(7).at(2)    = d;

                    QuadPts->at(8).at(0) 	= f;
                    QuadPts->at(8).at(1)    = g;
                    QuadPts->at(8).at(2)    = g;

                    QuadPts->at(9).at(0) 	= g;
                    QuadPts->at(9).at(1)    = g;
                    QuadPts->at(9).at(2)    = g;

                    QuadPts->at(10).at(0) 	= g;
                    QuadPts->at(10).at(1)   = g;
                    QuadPts->at(10).at(2)   = f;

                    QuadPts->at(11).at(0) 	= g;
                    QuadPts->at(11).at(1)   = f;
                    QuadPts->at(11).at(2)   = g;
                    
                    QuadPts->at(12).at(0) 	= h;
                    QuadPts->at(12).at(1)   = i;
                    QuadPts->at(12).at(2)   = i;
                    
                    QuadPts->at(13).at(0) 	= i;
                    QuadPts->at(13).at(1)   = h;
                    QuadPts->at(13).at(2)   = i;
                    
                    QuadPts->at(14).at(0) 	= i;
                    QuadPts->at(14).at(1)   = i;
                    QuadPts->at(14).at(2)   = h;
                    
                    QuadPts->at(15).at(0) 	= j;
                    QuadPts->at(15).at(1)   = i;
                    QuadPts->at(15).at(2)   = i;

                    QuadPts->at(16).at(0) 	= i;
                    QuadPts->at(16).at(1)   = j;
                    QuadPts->at(16).at(2)   = i;

                    QuadPts->at(17).at(0) 	= i;
                    QuadPts->at(17).at(1)   = i;
                    QuadPts->at(17).at(2)   = j;
                    
                    QuadPts->at(18).at(0) 	= i;
                    QuadPts->at(18).at(1)   = h;
                    QuadPts->at(18).at(2)   = j;

                    QuadPts->at(19).at(0) 	= h;
                    QuadPts->at(19).at(1)   = j;
                    QuadPts->at(19).at(2)   = i;
                    
                    QuadPts->at(20).at(0) 	= j;
                    QuadPts->at(20).at(1)   = i;
                    QuadPts->at(20).at(2)   = h;
                    
                    QuadPts->at(21).at(0) 	= i;
                    QuadPts->at(21).at(1)   = j;
                    QuadPts->at(21).at(2)   = h;

                    QuadPts->at(22).at(0) 	= h;
                    QuadPts->at(22).at(1)   = i;
                    QuadPts->at(22).at(2)   = j;
                    
                    QuadPts->at(23).at(0) 	= j;
                    QuadPts->at(23).at(1)   = h;
                    QuadPts->at(23).at(2)   = j;
                    
                    a = .039922750258168;
                    b = .010077211055321;
                    c = .055357181543654;
                    d = .048214285714286;
                    
                    QuadW->at(0)    = a;
                    QuadW->at(1)    = a;
                    QuadW->at(2)    = a;
                    QuadW->at(3)    = a;
                    QuadW->at(4)    = b;
                    QuadW->at(5)    = b;
                    QuadW->at(6)    = b;
                    QuadW->at(7)    = b;
                    QuadW->at(8)    = c;
                    QuadW->at(9)    = c;
                    QuadW->at(10)   = c;
                    QuadW->at(11)   = c;
                    QuadW->at(12)   = d;
                    QuadW->at(13)   = d;
                    QuadW->at(14)   = d;
                    QuadW->at(15)   = d;
                    QuadW->at(16)   = d;
                    QuadW->at(17)   = d;
                    QuadW->at(18)   = d;
                    QuadW->at(19)   = d;
                    QuadW->at(20)   = d;
                    QuadW->at(21)   = d;
                    QuadW->at(22)   = d;
                    QuadW->at(23)   = d;
            }
        }
        else if(FEType.at(0)=='Q'){
            if (Degree<=3)
                Degree=3;
            else if(Degree==4 || Degree==5)
                Degree=5;
            else if(Degree==6|| Degree==7)
                Degree=7;

            TEUCHOS_TEST_FOR_EXCEPTION(Degree>7, std::logic_error, "Quadrature rules for degree > 7 not available.");
            
            switch (Degree) {
                case 1: // 1 points in each direction; order 1
                {
                    QuadPts.reset(new vec2D_dbl_Type(1,vec_dbl_Type(3,0.0)));
                    QuadW->resize(1);
                    QuadW->at(0) = 2.;
                    break;
                }
                case 3: // 2 points in each direction; order 3
                {
                    double d = 1./sqrt(3);
                    QuadPts.reset(new vec2D_dbl_Type(8,vec_dbl_Type(3,0.0)));
                    QuadW->resize(8);
                    QuadPts->at(0).at(0) 	= -d;
                    QuadPts->at(0).at(1) 	= -d;
                    QuadPts->at(0).at(2) 	= -d;
                    QuadW->at(0)			= 1.;
                    
                    QuadPts->at(1).at(0) 	= -d;
                    QuadPts->at(1).at(1) 	= -d;
                    QuadPts->at(1).at(2) 	= d;
                    QuadW->at(1)			= 1.;
                    
                    QuadPts->at(2).at(0) 	= -d;
                    QuadPts->at(2).at(1) 	= d;
                    QuadPts->at(2).at(2) 	= -d;
                    QuadW->at(2)			= 1.;
                    
                    QuadPts->at(3).at(0) 	= -d;
                    QuadPts->at(3).at(1) 	= d;
                    QuadPts->at(3).at(2) 	= d;
                    QuadW->at(3)			= 1.;
                    
                    QuadPts->at(4).at(0) 	= d;
                    QuadPts->at(4).at(1) 	= -d;
                    QuadPts->at(4).at(2) 	= -d;
                    QuadW->at(4)			= 1.;
                    
                    QuadPts->at(5).at(0) 	= d;
                    QuadPts->at(5).at(1) 	= -d;
                    QuadPts->at(5).at(2) 	= d;
                    QuadW->at(5)			= 1.;
                    
                    QuadPts->at(6).at(0) 	= d;
                    QuadPts->at(6).at(1) 	= d;
                    QuadPts->at(6).at(2) 	= -d;
                    QuadW->at(6)			= 1.;
                    
                    QuadPts->at(7).at(0) 	= d;
                    QuadPts->at(7).at(1) 	= d;
                    QuadPts->at(7).at(2) 	= d;
                    QuadW->at(7)			= 1.;
                    break;
                }
                case 5: // 3 points in each direction; order 5
                {
                    double a=sqrt(3./5);
                    double b=5./9;
                    double c=8./9;
                    std::vector<double> p(3);
                    p[0] = -a; p[1] = 0.; p[2] = a;
                    std::vector<double> w(3);
                    w[0] = b; w[1] = c; w[2] = b;
                    QuadPts.reset(new vec2D_dbl_Type(27,vec_dbl_Type(3,0.0)));
                    QuadW->resize(27);
                    int counter=0;
                    for (int i=0; i<3; i++) {
                        for (int j=0; j<3; j++) {
                            for (int k=0; k<3; k++) {
                                QuadPts->at(counter)[0] = p[k];
                                QuadPts->at(counter)[1] = p[j];
                                QuadPts->at(counter)[2] = p[i];
                                QuadW->at(counter)      = w[k]*w[j]*w[i];
                                counter++;
                            }
                        }
                    }
                    break;
                }
                case 7: // 4 points in each direction; order 7
                {
                    double aa = 2./7 * sqrt(6./5) ;
                    std::vector<double> p(4);
                    p[0] = - sqrt(3./7 + aa);
                    p[1] = - sqrt(3./7 - aa);
                    p[2] = -p[1];
                    p[3] = -p[0];
                    
                    double bb = sqrt(30.);
                    std::vector<double> w(4);
                    w[0] = ( 18. - bb ) / 36;
                    w[1] = ( 18. + bb ) / 36;
                    w[2] = w[1];
                    w[3] = w[0];
                    
                    QuadPts.reset(new vec2D_dbl_Type(64,vec_dbl_Type(3,0.0)));
                    QuadW->resize(64);

                    int counter=0;
                    for (int i=0; i<4; i++) {
                        for (int j=0; j<4; j++) {
                            for (int k=0; k<4; k++) {
                                QuadPts->at(counter)[0] = p[k];
                                QuadPts->at(counter)[1] = p[j];
                                QuadPts->at(counter)[2] = p[i];
                                QuadW->at(counter)      = w[k]*w[j]*w[i];
                                counter++;
                            }
                        }
                    }
                    break;
                }
            }
        }
    }
    
}

template <class SC, class LO, class GO, class NO>
int FE<SC,LO,GO,NO>::getPhi(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree,
                            std::string FETypeQuadPoints){

    int 			nmbLocElPts;
    int 			intFE;
    double  		value;
    vec2D_dbl_ptr_Type	QuadPts;
    if (dim==1) {
        this->getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 2;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 3;
            intFE = 2;
        }
        Phi.reset( new vec2D_dbl_Type( weightsPhi->size(), vec_dbl_Type( nmbLocElPts, 0.0 ) ) );
        for (int k=0; k<Phi->size(); k++ ){
            for (int i=0; i<Phi->at(0).size(); i++) {
                this->phi(dim,intFE,i,QuadPts->at(k),&value);
                Phi->at(k).at(i) = value;
            }
        }

    }
    else if (dim==2) {
        this->getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 3;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 6;
            intFE = 2;
        }

        Phi.reset(new vec2D_dbl_Type(weightsPhi->size(),vec_dbl_Type(nmbLocElPts,0.0)));

        for (int k=0; k<Phi->size(); k++ ){
            for (int i=0; i<Phi->at(0).size(); i++) {
                this->phi(dim,intFE,i,QuadPts->at(k),&value);
                Phi->at(k).at(i) = value;
            }
        }
    }
    else if(dim==3){
        if (FETypeQuadPoints!="")
            this->getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FETypeQuadPoints);
        else
            this->getQuadratureValues(dim, Degree, QuadPts, weightsPhi, FEType);
        
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 4;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 10;
            intFE = 2;
        }
        else if (FEType == "Q1") {
            nmbLocElPts = 8;
            intFE = 3;
        }
        else if (FEType == "Q2") {
            nmbLocElPts = 27;
            intFE = 4;
        }
        else if (FEType == "Q2-20") {
            nmbLocElPts = 20;
            intFE = 5;
        }
        else if (FEType == "P1-disc") {
            nmbLocElPts = 4;
            intFE = 1;
        }
        else if (FEType == "P1-disc-global")
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "P1-disc-global not implemented yet.");
        

        Phi.reset(new vec2D_dbl_Type(weightsPhi->size(),vec_dbl_Type(nmbLocElPts,0.0)));

        for (int k=0; k<Phi->size(); k++ ){
            for (int i=0; i<Phi->at(0).size(); i++) {
                this->phi(dim,intFE,i,QuadPts->at(k),&value);
                Phi->at(k).at(i) = value;
            }
        }
    }
    return intFE;
}
template <class SC, class LO, class GO, class NO>
int FE<SC,LO,GO,NO>::getPhiGlobal(vec2D_dbl_ptr_Type &Phi,
                            vec_dbl_ptr_Type &weightsPhi,
                            int dim,
                            std::string FEType,
                            int Degree){
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "getPhiGlobal not implemented yet.");
}
template <class SC, class LO, class GO, class NO>
int FE<SC,LO,GO,NO>::getDPhi(vec3D_dbl_ptr_Type &DPhi,
                             vec_dbl_ptr_Type &weightsDPhi,
                             int dim,
                 std::string FEType,
                 int Degree){

    int 			nmbLocElPts;
    int 			intFE;
    vec_dbl_ptr_Type 	value(new vec_dbl_Type(dim,0.0));
    vec2D_dbl_ptr_Type	QuadPts;

    if (dim==2) {
        this->getQuadratureValues(dim, Degree, QuadPts, weightsDPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 3;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 6;
            intFE = 2;
        }

        DPhi.reset(new vec3D_dbl_Type(weightsDPhi->size(),vec2D_dbl_Type(nmbLocElPts,vec_dbl_Type(2,0.0))));

        for (int k=0; k<DPhi->size(); k++ ){
            for (int i=0; i<DPhi->at(0).size(); i++) {
                this->gradPhi(dim,intFE,i,QuadPts->at(k),value);
                for (int j=0; j<2; j++) {
                    DPhi->at(k).at(i).at(j) = value->at(j);
                }
            }
        }
    }

    else if(dim==3){
    	this->getQuadratureValues(dim, Degree, QuadPts, weightsDPhi, FEType);
        if (FEType == "P0") {
            nmbLocElPts = 1;
            intFE = 0;
        }
        else if (FEType == "P1") {
            nmbLocElPts = 4;
            intFE = 1;
        }
        else if (FEType == "P2") {
            nmbLocElPts = 10;
            intFE = 2;
        }
        else if (FEType == "Q1") {
            nmbLocElPts = 8;
            intFE = 3;
        }
        else if (FEType == "Q2") {
            nmbLocElPts = 27;
            intFE = 4;
        }
        else if (FEType == "Q2-20") {
            nmbLocElPts = 20;
            intFE = 5;
        }
        else if (FEType == "P1-disc") {
            nmbLocElPts = 4;
            intFE = 6;
        }
        else if (FEType == "P1-disc-global")
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error   ,"grad of P1-disc-global not implemented yet.");

        DPhi.reset( new vec3D_dbl_Type( weightsDPhi->size(), vec2D_dbl_Type( nmbLocElPts, vec_dbl_Type(3,0.0) ) ) );
        for (int k=0; k<DPhi->size(); k++ ){
            for (int i=0; i<DPhi->at(0).size(); i++) {
                this->gradPhi(dim,intFE,i,QuadPts->at(k),value);
                for (int j=0; j<3; j++) {
                    DPhi->at(k).at(i).at(j) = value->at(j);
                }
            }
        }
    }

    return intFE;
}
*/
template <class SC, class LO, class GO, class NO>
int FE<SC,LO,GO,NO>::checkFE(int dim,
                 std::string FEType){

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


/*************************************************************
 * AceGen    6.921 MacOSX (29 Jan 19)                         *
 *           Co. J. Korelc  2013           12 Feb 19 12:07:04 *
 **************************************************************
 User     : Full professional version
 Notebook : nh3d_C
 Evaluation time                 : 6 s     Mode  : Optimal
 Number of formulae              : 181     Method: Automatic
 Subroutine                      : nh3d size: 4928
 Total size of Mathematica  code : 4928 subexpressions
 Total size of C code            : 10178 bytes */
/******************* S U B R O U T I N E *********************/
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::nh3d(double* v, double (*E), double (*Nu), double** F , double** Pmat, double**** Amat)
{
    v[356]=2e0*F[0][2];
    v[354]=2e0*F[0][1];
    v[323]=(*E)/(1e0+(*Nu));
    v[3]=((*Nu)*v[323])/(1e0-2e0*(*Nu));
    v[5]=v[323]/2e0;
    v[36]=v[5]/2e0;
    v[65]=2e0*F[0][1];
    v[86]=2e0*F[0][2];
    v[57]=2e0*F[1][0];
    v[66]=2e0*F[1][1];
    v[87]=2e0*F[1][2];
    v[58]=2e0*F[2][0];
    v[67]=2e0*F[2][1];
    v[18]=F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1];
    v[335]=(v[18]*v[18]);
    v[88]=2e0*F[2][2];
    v[24]=F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2];
    v[334]=(v[24]*v[24]);
    v[325]=v[18]*v[24];
    v[22]=F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2];
    v[15]=Power(F[0][0],2)+Power(F[1][0],2)+Power(F[2][0],2);
    v[228]=-(F[2][1]*v[18]);
    v[225]=F[2][2]*v[18];
    v[217]=-(F[1][1]*v[18]);
    v[214]=F[1][2]*v[18];
    v[194]=-(F[2][0]*v[18]);
    v[185]=-(F[1][0]*v[18]);
    v[268]=F[2][1]*v[22];
    v[264]=-(F[2][2]*v[22]);
    v[255]=F[1][1]*v[22];
    v[251]=-(F[1][2]*v[22]);
    v[190]=-(F[2][0]*v[22]);
    v[181]=-(F[1][0]*v[22]);
    v[172]=-(F[0][0]*v[22]);
    v[20]=Power(F[0][1],2)+Power(F[1][1],2)+Power(F[2][1],2);
    v[324]=-(v[20]*v[22]);
    v[327]=2e0*(v[324]+v[325]);
    v[94]=v[20]*v[88];
    v[92]=v[20]*v[87];
    v[90]=v[20]*v[86];
    v[138]=v[15]*v[20]-v[335];
    v[270]=F[2][0]*v[24];
    v[260]=-(F[2][2]*v[24]);
    v[257]=F[1][0]*v[24];
    v[247]=-(F[1][2]*v[24]);
    v[244]=F[0][0]*v[24];
    v[232]=-(F[0][2]*v[24]);
    v[222]=-(F[2][1]*v[24]);
    v[211]=-(F[1][1]*v[24]);
    v[198]=-(F[0][1]*v[24]);
    v[168]=v[18]*v[22]-v[15]*v[24];
    v[331]=2e0*v[168];
    v[329]=2e0*v[168];
    v[326]=2e0*v[168];
    v[38]=-(v[22]*v[22]);
    v[26]=Power(F[0][2],2)+Power(F[1][2],2)+Power(F[2][2],2);
    v[333]=v[20]*v[26]-v[334];
    v[351]=2e0*F[0][0]*v[333];
    v[236]=v[22]*v[24]-v[18]*v[26];
    v[332]=2e0*v[236];
    v[330]=2e0*v[236];
    v[328]=2e0*v[236];
    v[99]=v[26]*v[58];
    v[97]=v[26]*v[57];
    v[93]=v[26]*v[67];
    v[91]=v[26]*v[66];
    v[89]=v[26]*v[65];
    v[148]=v[15]*v[26]+v[38];
    v[29]=v[148]*v[20]+2e0*v[22]*v[325]-v[15]*v[334]-v[26]*v[335];
    v[336]=1e0/Power(v[29],2);
    v[32]=-v[5]+v[3]*log(sqrt(v[29]));
    v[337]=(v[3]/4e0-v[32]/2e0)*v[336];
    v[137]=v[337]*(F[2][1]*v[326]+F[2][0]*v[327]-v[335]*v[88]+v[15]*v[94]);
    v[147]=v[137]*v[138];
    v[136]=v[337]*(F[2][2]*v[329]+F[2][0]*v[330]+v[38]*v[67]+v[15]*v[93]);
    v[156]=v[136]*v[148];
    v[135]=v[337]*(F[2][2]*v[327]+F[2][1]*v[328]-v[334]*v[58]+v[20]*v[99]);
    v[165]=v[135]*v[333];
    v[134]=v[337]*(F[1][1]*v[326]+F[1][0]*v[327]-v[335]*v[87]+v[15]*v[92]);
    v[144]=v[134]*v[138];
    v[133]=v[337]*(F[1][2]*v[331]+F[1][0]*v[332]+v[38]*v[66]+v[15]*v[91]);
    v[153]=v[133]*v[148];
    v[132]=v[337]*(F[1][2]*v[327]+F[1][1]*v[328]-v[334]*v[57]+v[20]*v[97]);
    v[162]=v[132]*v[333];
    v[131]=v[337]*(F[0][0]*v[327]+F[0][1]*v[329]-v[335]*v[86]+v[15]*v[90]);
    v[130]=v[337]*(F[0][2]*v[331]+F[0][0]*v[332]+v[38]*v[65]+v[15]*v[89]);
    v[128]=v[337]*(F[0][2]*v[327]+F[0][1]*v[330]+v[351]);
    v[37]=v[32]/(2e0*v[29]);
    v[355]=v[37]*(2e0*v[172]+v[15]*v[86]);
    v[353]=v[37]*(2e0*v[232]+v[89]);
    v[352]=v[37]*(2e0*v[198]+v[90]);
    v[349]=-2e0*(F[1][0]*v[20]+v[217])*v[37];
    v[348]=-(v[37]*(2e0*v[185]+v[15]*v[66]));
    v[347]=-2e0*(F[2][0]*v[20]+v[228])*v[37];
    v[346]=-(v[37]*(2e0*v[194]+v[15]*v[67]));
    v[345]=-(v[37]*(2e0*v[251]+v[97]));
    v[344]=-(v[37]*(2e0*v[181]+v[15]*v[87]));
    v[343]=-(v[37]*(2e0*v[264]+v[99]));
    v[342]=-(v[37]*(2e0*v[190]+v[15]*v[88]));
    v[341]=-(v[37]*(2e0*v[247]+v[91]));
    v[340]=-(v[37]*(2e0*v[211]+v[92]));
    v[339]=-(v[37]*(2e0*v[260]+v[93]));
    v[338]=-(v[37]*(2e0*v[222]+v[94]));
    v[272]=v[137]*v[328]+v[37]*(2e0*v[268]+2e0*v[270]-2e0*v[18]*v[88]);
    v[267]=v[136]*v[328]+v[343];
    v[263]=v[135]*v[328]+v[339];
    v[259]=v[134]*v[328]+v[37]*(2e0*v[255]+2e0*v[257]-2e0*v[18]*v[87]);
    v[254]=v[133]*v[328]+v[345];
    v[250]=v[132]*v[328]+v[341];
    v[246]=v[131]*v[328]+v[37]*(2e0*F[0][1]*v[22]+2e0*v[244]-2e0*v[18]*v[86]);
    v[241]=v[130]*v[328]+2e0*(F[0][2]*v[22]-F[0][0]*v[26])*v[37];
    v[231]=v[137]*v[327]+v[347];
    v[227]=v[136]*v[327]+v[37]*(2e0*v[225]+2e0*v[270]-2e0*v[22]*v[67]);
    v[224]=v[135]*v[327]+v[338];
    v[301]=2e0*F[1][0]*v[165]+F[1][2]*v[224]+F[1][1]*v[263];
    v[279]=2e0*F[0][0]*v[165]+F[0][2]*v[224]+F[0][1]*v[263];
    v[220]=v[134]*v[327]+v[349];
    v[216]=v[133]*v[327]+v[37]*(2e0*v[214]+2e0*v[257]-2e0*v[22]*v[66]);
    v[213]=v[132]*v[327]+v[340];
    v[276]=2e0*F[0][0]*v[162]+F[0][2]*v[213]+F[0][1]*v[250];
    v[209]=v[131]*v[327]+2e0*(F[0][1]*v[18]-F[0][0]*v[20])*v[37];
    v[196]=v[137]*v[326]+v[346];
    v[314]=2e0*F[1][2]*v[147]+F[1][1]*v[196]+F[1][0]*v[231];
    v[296]=2e0*F[0][2]*v[147]+F[0][1]*v[196]+F[0][0]*v[231];
    v[192]=v[136]*v[326]+v[342];
    v[308]=2e0*F[1][1]*v[156]+F[1][2]*v[192]+F[1][0]*v[267];
    v[288]=2e0*F[0][1]*v[156]+F[0][2]*v[192]+F[0][0]*v[267];
    v[188]=v[135]*v[326]+v[37]*(2e0*v[225]+2e0*v[268]-2e0*v[24]*v[58]);
    v[187]=v[134]*v[326]+v[348];
    v[293]=2e0*F[0][2]*v[144]+F[0][1]*v[187]+F[0][0]*v[220];
    v[183]=v[133]*v[326]+v[344];
    v[285]=2e0*F[0][1]*v[153]+F[0][2]*v[183]+F[0][0]*v[254];
    v[179]=v[132]*v[326]+v[37]*(2e0*v[214]+2e0*v[255]-2e0*v[24]*v[57]);
    v[178]=v[131]*v[326]+2e0*(-(F[0][1]*v[15])+F[0][0]*v[18])*v[37];
    v[167]=v[137]*v[333]-v[338];
    v[303]=2e0*F[1][0]*v[167]+F[1][2]*v[231]+F[1][1]*v[272];
    v[281]=2e0*F[0][0]*v[167]+F[0][2]*v[231]+F[0][1]*v[272];
    v[166]=v[136]*v[333]-v[339];
    v[302]=2e0*F[1][0]*v[166]+F[1][2]*v[227]+F[1][1]*v[267];
    v[280]=2e0*F[0][0]*v[166]+F[0][2]*v[227]+F[0][1]*v[267];
    v[164]=v[134]*v[333]-v[340];
    v[278]=2e0*F[0][0]*v[164]+F[0][2]*v[220]+F[0][1]*v[259];
    v[163]=v[133]*v[333]-v[341];
    v[277]=2e0*F[0][0]*v[163]+F[0][2]*v[216]+F[0][1]*v[254];
    v[157]=v[137]*v[148]-v[342];
    v[309]=2e0*F[1][1]*v[157]+F[1][2]*v[196]+F[1][0]*v[272];
    v[289]=2e0*F[0][1]*v[157]+F[0][2]*v[196]+F[0][0]*v[272];
    v[155]=v[135]*v[148]-v[343];
    v[307]=2e0*F[1][1]*v[155]+F[1][2]*v[188]+F[1][0]*v[263];
    v[287]=2e0*F[0][1]*v[155]+F[0][2]*v[188]+F[0][0]*v[263];
    v[154]=v[134]*v[148]-v[344];
    v[286]=2e0*F[0][1]*v[154]+F[0][2]*v[187]+F[0][0]*v[259];
    v[284]=F[0][2]*v[179]+F[0][0]*v[250]+2e0*F[0][1]*(v[132]*v[148]-v[345]);
    v[146]=v[136]*v[138]-v[346];
    v[313]=2e0*F[1][2]*v[146]+F[1][1]*v[192]+F[1][0]*v[227];
    v[295]=2e0*F[0][2]*v[146]+F[0][1]*v[192]+F[0][0]*v[227];
    v[145]=v[135]*v[138]-v[347];
    v[312]=2e0*F[1][2]*v[145]+F[1][1]*v[188]+F[1][0]*v[224];
    v[294]=2e0*F[0][2]*v[145]+F[0][1]*v[188]+F[0][0]*v[224];
    v[292]=F[0][1]*v[183]+F[0][0]*v[216]+2e0*F[0][2]*(v[133]*v[138]-v[348]);
    v[291]=F[0][1]*v[179]+F[0][0]*v[213]+(v[132]*v[138]-v[349])*v[356];
    v[35]=v[36]+v[138]*v[37];
    v[310]=2e0*v[35];
    v[40]=v[36]+v[148]*v[37];
    v[304]=2e0*v[40];
    v[43]=v[36]+v[333]*v[37];
    v[297]=2e0*v[43];
    v[44]=v[326]*v[37];
    v[319]=2e0*F[2][1]*v[157]+F[2][2]*v[196]+F[2][0]*v[272]+v[44];
    v[306]=2e0*F[1][1]*v[154]+F[1][2]*v[187]+F[1][0]*v[259]+v[44];
    v[283]=F[0][2]*v[178]+F[0][0]*v[246]+v[354]*(v[131]*v[148]+v[355])+v[44];
    v[45]=v[327]*v[37];
    v[317]=2e0*F[2][0]*v[167]+F[2][2]*v[231]+F[2][1]*v[272]+v[45];
    v[300]=2e0*F[1][0]*v[164]+F[1][2]*v[220]+F[1][1]*v[259]+v[45];
    v[275]=F[0][2]*v[209]+F[0][1]*v[246]+2e0*F[0][0]*(v[131]*v[333]+v[352])+v[45];
    v[46]=v[328]*v[37];
    v[316]=2e0*F[2][0]*v[166]+F[2][2]*v[227]+F[2][1]*v[267]+v[46];
    v[299]=2e0*F[1][0]*v[163]+F[1][2]*v[216]+F[1][1]*v[254]+v[46];
    v[274]=F[0][1]*v[241]+2e0*F[0][0]*(v[130]*v[333]+v[353])+v[46]+F[0][2]*(v[130]*v[327]+v[37]*
                                                                            (2e0*F[0][2]*v[18]+2e0*v[244]-2e0*v[22]*v[65]));
    Pmat[0][0]=F[0][0]*v[297]+F[0][2]*v[45]+F[0][1]*v[46];
    Pmat[0][1]=F[0][1]*v[304]+F[0][2]*v[44]+F[0][0]*v[46];
    Pmat[0][2]=F[0][2]*v[310]+F[0][1]*v[44]+F[0][0]*v[45];
    Pmat[1][0]=2e0*F[1][0]*v[43]+F[1][2]*v[45]+F[1][1]*v[46];
    Pmat[1][1]=2e0*F[1][1]*v[40]+F[1][2]*v[44]+F[1][0]*v[46];
    Pmat[1][2]=2e0*F[1][2]*v[35]+F[1][1]*v[44]+F[1][0]*v[45];
    Pmat[2][0]=F[2][0]*v[297]+F[2][2]*v[45]+F[2][1]*v[46];
    Pmat[2][1]=F[2][1]*v[304]+F[2][2]*v[44]+F[2][0]*v[46];
    Pmat[2][2]=F[2][2]*v[310]+F[2][1]*v[44]+F[2][0]*v[45];
    Amat[0][0][0][0]=v[297]+v[128]*v[351]+F[0][2]*(v[128]*v[327]-v[352])+F[0][1]*(v[128]*v[328]-v[353]
                                                                                  );
    Amat[0][0][0][1]=v[274];
    Amat[0][0][0][2]=v[275];
    Amat[0][0][1][0]=v[276];
    Amat[0][0][1][1]=v[277];
    Amat[0][0][1][2]=v[278];
    Amat[0][0][2][0]=v[279];
    Amat[0][0][2][1]=v[280];
    Amat[0][0][2][2]=v[281];
    Amat[0][1][0][0]=v[274];
    Amat[0][1][0][1]=F[0][0]*v[241]+v[304]+v[130]*v[148]*v[354]+F[0][2]*(v[130]*v[326]-v[355]);
    Amat[0][1][0][2]=v[283];
    Amat[0][1][1][0]=v[284];
    Amat[0][1][1][1]=v[285];
    Amat[0][1][1][2]=v[286];
    Amat[0][1][2][0]=v[287];
    Amat[0][1][2][1]=v[288];
    Amat[0][1][2][2]=v[289];
    Amat[0][2][0][0]=v[275];
    Amat[0][2][0][1]=v[283];
    Amat[0][2][0][2]=F[0][1]*v[178]+F[0][0]*v[209]+v[310]+v[131]*v[138]*v[356];
    Amat[0][2][1][0]=v[291];
    Amat[0][2][1][1]=v[292];
    Amat[0][2][1][2]=v[293];
    Amat[0][2][2][0]=v[294];
    Amat[0][2][2][1]=v[295];
    Amat[0][2][2][2]=v[296];
    Amat[1][0][0][0]=v[276];
    Amat[1][0][0][1]=v[284];
    Amat[1][0][0][2]=v[291];
    Amat[1][0][1][0]=2e0*F[1][0]*v[162]+F[1][2]*v[213]+F[1][1]*v[250]+v[297];
    Amat[1][0][1][1]=v[299];
    Amat[1][0][1][2]=v[300];
    Amat[1][0][2][0]=v[301];
    Amat[1][0][2][1]=v[302];
    Amat[1][0][2][2]=v[303];
    Amat[1][1][0][0]=v[277];
    Amat[1][1][0][1]=v[285];
    Amat[1][1][0][2]=v[292];
    Amat[1][1][1][0]=v[299];
    Amat[1][1][1][1]=2e0*F[1][1]*v[153]+F[1][2]*v[183]+F[1][0]*v[254]+v[304];
    Amat[1][1][1][2]=v[306];
    Amat[1][1][2][0]=v[307];
    Amat[1][1][2][1]=v[308];
    Amat[1][1][2][2]=v[309];
    Amat[1][2][0][0]=v[278];
    Amat[1][2][0][1]=v[286];
    Amat[1][2][0][2]=v[293];
    Amat[1][2][1][0]=v[300];
    Amat[1][2][1][1]=v[306];
    Amat[1][2][1][2]=2e0*F[1][2]*v[144]+F[1][1]*v[187]+F[1][0]*v[220]+v[310];
    Amat[1][2][2][0]=v[312];
    Amat[1][2][2][1]=v[313];
    Amat[1][2][2][2]=v[314];
    Amat[2][0][0][0]=v[279];
    Amat[2][0][0][1]=v[287];
    Amat[2][0][0][2]=v[294];
    Amat[2][0][1][0]=v[301];
    Amat[2][0][1][1]=v[307];
    Amat[2][0][1][2]=v[312];
    Amat[2][0][2][0]=2e0*F[2][0]*v[165]+F[2][2]*v[224]+F[2][1]*v[263]+v[297];
    Amat[2][0][2][1]=v[316];
    Amat[2][0][2][2]=v[317];
    Amat[2][1][0][0]=v[280];
    Amat[2][1][0][1]=v[288];
    Amat[2][1][0][2]=v[295];
    Amat[2][1][1][0]=v[302];
    Amat[2][1][1][1]=v[308];
    Amat[2][1][1][2]=v[313];
    Amat[2][1][2][0]=v[316];
    Amat[2][1][2][1]=2e0*F[2][1]*v[156]+F[2][2]*v[192]+F[2][0]*v[267]+v[304];
    Amat[2][1][2][2]=v[319];
    Amat[2][2][0][0]=v[281];
    Amat[2][2][0][1]=v[289];
    Amat[2][2][0][2]=v[296];
    Amat[2][2][1][0]=v[303];
    Amat[2][2][1][1]=v[309];
    Amat[2][2][1][2]=v[314];
    Amat[2][2][2][0]=v[317];
    Amat[2][2][2][1]=v[319];
    Amat[2][2][2][2]=2e0*F[2][2]*v[147]+F[2][1]*v[196]+F[2][0]*v[231]+v[310];
}


/*************************************************************
 * AceGen    6.921 MacOSX (29 Jan 19)                         *
 *           Co. J. Korelc  2013           12 Feb 19 12:06:46 *
 **************************************************************
 User     : Full professional version
 Notebook : mr3d_C
 Evaluation time                 : 7 s     Mode  : Optimal
 Number of formulae              : 190     Method: Automatic
 Subroutine                      : mr3d size: 5215
 Total size of Mathematica  code : 5215 subexpressions
 Total size of C code            : 10798 bytes */

/******************* S U B R O U T I N E *********************/
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::mr3d(double* v,double (*E),double (*Nu),double (*C)
          ,double** F,double** Pmat,double**** Amat)
{
    v[366]=2e0*F[0][2];
    v[364]=2e0*F[0][1];
    v[4]=(*E)/(2e0+2e0*(*Nu));
    v[139]=((*C)*v[4])/2e0;
    v[5]=(*E)/(3e0-6e0*(*Nu));
    v[57]=2e0*F[0][0];
    v[150]=v[139]*v[57];
    v[66]=2e0*F[0][1];
    v[165]=v[139]*v[66];
    v[87]=2e0*F[0][2];
    v[167]=v[139]*v[87];
    v[58]=2e0*F[1][0];
    v[155]=v[139]*v[58];
    v[67]=2e0*F[1][1];
    v[170]=v[139]*v[67];
    v[88]=2e0*F[1][2];
    v[172]=v[139]*v[88];
    v[59]=2e0*F[2][0];
    v[159]=v[139]*v[59];
    v[68]=2e0*F[2][1];
    v[175]=v[139]*v[68];
    v[18]=F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1];
    v[345]=(v[18]*v[18]);
    v[89]=2e0*F[2][2];
    v[177]=v[139]*v[89];
    v[24]=F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2];
    v[344]=(v[24]*v[24]);
    v[335]=v[18]*v[24];
    v[22]=F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2];
    v[15]=Power(F[0][0],2)+Power(F[1][0],2)+Power(F[2][0],2);
    v[239]=-(F[2][1]*v[18]);
    v[236]=F[2][2]*v[18];
    v[228]=-(F[1][1]*v[18]);
    v[225]=F[1][2]*v[18];
    v[205]=-(F[2][0]*v[18]);
    v[196]=-(F[1][0]*v[18]);
    v[279]=F[2][1]*v[22];
    v[275]=-(F[2][2]*v[22]);
    v[266]=F[1][1]*v[22];
    v[262]=-(F[1][2]*v[22]);
    v[201]=-(F[2][0]*v[22]);
    v[192]=-(F[1][0]*v[22]);
    v[183]=-(F[0][0]*v[22]);
    v[20]=Power(F[0][1],2)+Power(F[1][1],2)+Power(F[2][1],2);
    v[334]=-(v[20]*v[22]);
    v[337]=2e0*(v[334]+v[335]);
    v[95]=v[20]*v[89];
    v[93]=v[20]*v[88];
    v[91]=v[20]*v[87];
    v[140]=v[15]*v[20]-v[345];
    v[281]=F[2][0]*v[24];
    v[271]=-(F[2][2]*v[24]);
    v[268]=F[1][0]*v[24];
    v[258]=-(F[1][2]*v[24]);
    v[255]=F[0][0]*v[24];
    v[243]=-(F[0][2]*v[24]);
    v[233]=-(F[2][1]*v[24]);
    v[222]=-(F[1][1]*v[24]);
    v[209]=-(F[0][1]*v[24]);
    v[179]=v[18]*v[22]-v[15]*v[24];
    v[341]=2e0*v[179];
    v[339]=2e0*v[179];
    v[336]=2e0*v[179];
    v[38]=-(v[22]*v[22]);
    v[26]=Power(F[0][2],2)+Power(F[1][2],2)+Power(F[2][2],2);
    v[343]=v[20]*v[26]-v[344];
    v[361]=v[343]*v[57];
    v[247]=v[22]*v[24]-v[18]*v[26];
    v[342]=2e0*v[247];
    v[340]=2e0*v[247];
    v[338]=2e0*v[247];
    v[100]=v[26]*v[59];
    v[98]=v[26]*v[58];
    v[94]=v[26]*v[68];
    v[92]=v[26]*v[67];
    v[90]=v[26]*v[66];
    v[151]=v[15]*v[26]+v[38];
    v[29]=v[151]*v[20]+2e0*v[22]*v[335]-v[15]*v[344]-v[26]*v[345];
    v[346]=1e0/Power(v[29],2);
    v[33]=-2e0*v[139]-v[4]+v[5]*log(sqrt(v[29]));
    v[347]=v[346]*(-v[33]/2e0+v[5]/4e0);
    v[138]=v[347]*(F[2][1]*v[336]+F[2][0]*v[337]-v[345]*v[89]+v[15]*v[95]);
    v[149]=v[138]*v[140];
    v[137]=v[347]*(F[2][2]*v[339]+F[2][0]*v[340]+v[38]*v[68]+v[15]*v[94]);
    v[161]=v[137]*v[151];
    v[136]=v[347]*(v[100]*v[20]+F[2][2]*v[337]+F[2][1]*v[338]-v[344]*v[59]);
    v[174]=v[136]*v[343];
    v[135]=v[347]*(F[1][1]*v[336]+F[1][0]*v[337]-v[345]*v[88]+v[15]*v[93]);
    v[146]=v[135]*v[140];
    v[134]=v[347]*(F[1][2]*v[341]+F[1][0]*v[342]+v[38]*v[67]+v[15]*v[92]);
    v[157]=v[134]*v[151];
    v[133]=v[347]*(F[1][2]*v[337]+F[1][1]*v[338]-v[344]*v[58]+v[20]*v[98]);
    v[169]=v[133]*v[343];
    v[132]=v[347]*(F[0][0]*v[337]+F[0][1]*v[339]-v[345]*v[87]+v[15]*v[91]);
    v[131]=v[347]*(F[0][2]*v[341]+F[0][0]*v[342]+v[38]*v[66]+v[15]*v[90]);
    v[129]=v[347]*(F[0][2]*v[337]+F[0][1]*v[340]+v[361]);
    v[37]=v[33]/(2e0*v[29]);
    v[365]=v[37]*(2e0*v[183]+v[15]*v[87]);
    v[363]=v[37]*(2e0*v[243]+v[90]);
    v[362]=v[37]*(2e0*v[209]+v[91]);
    v[359]=-2e0*(F[1][0]*v[20]+v[228])*v[37];
    v[358]=-(v[37]*(2e0*v[196]+v[15]*v[67]));
    v[357]=-2e0*(F[2][0]*v[20]+v[239])*v[37];
    v[356]=-(v[37]*(2e0*v[205]+v[15]*v[68]));
    v[355]=-(v[37]*(2e0*v[262]+v[98]));
    v[354]=-(v[37]*(2e0*v[192]+v[15]*v[88]));
    v[353]=-((v[100]+2e0*v[275])*v[37]);
    v[352]=-(v[37]*(2e0*v[201]+v[15]*v[89]));
    v[351]=-(v[37]*(2e0*v[258]+v[92]));
    v[350]=-(v[37]*(2e0*v[222]+v[93]));
    v[349]=-(v[37]*(2e0*v[271]+v[94]));
    v[348]=-(v[37]*(2e0*v[233]+v[95]));
    v[283]=v[138]*v[338]+v[37]*(2e0*v[279]+2e0*v[281]-2e0*v[18]*v[89]);
    v[278]=-v[159]+v[137]*v[338]+v[353];
    v[274]=-v[175]+v[136]*v[338]+v[349];
    v[270]=v[135]*v[338]+v[37]*(2e0*v[266]+2e0*v[268]-2e0*v[18]*v[88]);
    v[265]=-v[155]+v[134]*v[338]+v[355];
    v[261]=-v[170]+v[133]*v[338]+v[351];
    v[257]=v[132]*v[338]+v[37]*(2e0*F[0][1]*v[22]+2e0*v[255]-2e0*v[18]*v[87]);
    v[252]=-v[150]+v[131]*v[338]+2e0*(F[0][2]*v[22]-F[0][0]*v[26])*v[37];
    v[242]=-v[159]+v[138]*v[337]+v[357];
    v[238]=v[137]*v[337]+v[37]*(2e0*v[236]+2e0*v[281]-2e0*v[22]*v[68]);
    v[235]=-v[177]+v[136]*v[337]+v[348];
    v[312]=2e0*F[1][0]*v[174]+F[1][2]*v[235]+F[1][1]*v[274];
    v[290]=2e0*F[0][0]*v[174]+F[0][2]*v[235]+F[0][1]*v[274];
    v[231]=-v[155]+v[135]*v[337]+v[359];
    v[227]=v[134]*v[337]+v[37]*(2e0*v[225]+2e0*v[268]-2e0*v[22]*v[67]);
    v[224]=-v[172]+v[133]*v[337]+v[350];
    v[287]=2e0*F[0][0]*v[169]+F[0][2]*v[224]+F[0][1]*v[261];
    v[220]=-v[150]+v[132]*v[337]+2e0*(F[0][1]*v[18]-F[0][0]*v[20])*v[37];
    v[207]=-v[175]+v[138]*v[336]+v[356];
    v[325]=2e0*F[1][2]*v[149]+F[1][1]*v[207]+F[1][0]*v[242];
    v[307]=2e0*F[0][2]*v[149]+F[0][1]*v[207]+F[0][0]*v[242];
    v[203]=-v[177]+v[137]*v[336]+v[352];
    v[319]=2e0*F[1][1]*v[161]+F[1][2]*v[203]+F[1][0]*v[278];
    v[299]=2e0*F[0][1]*v[161]+F[0][2]*v[203]+F[0][0]*v[278];
    v[199]=v[136]*v[336]+v[37]*(2e0*v[236]+2e0*v[279]-2e0*v[24]*v[59]);
    v[198]=-v[170]+v[135]*v[336]+v[358];
    v[304]=2e0*F[0][2]*v[146]+F[0][1]*v[198]+F[0][0]*v[231];
    v[194]=-v[172]+v[134]*v[336]+v[354];
    v[296]=2e0*F[0][1]*v[157]+F[0][2]*v[194]+F[0][0]*v[265];
    v[190]=v[133]*v[336]+v[37]*(2e0*v[225]+2e0*v[266]-2e0*v[24]*v[58]);
    v[189]=-v[165]+v[132]*v[336]+2e0*(-(F[0][1]*v[15])+F[0][0]*v[18])*v[37];
    v[178]=v[177]+v[138]*v[343]-v[348];
    v[314]=2e0*F[1][0]*v[178]+F[1][2]*v[242]+F[1][1]*v[283];
    v[292]=2e0*F[0][0]*v[178]+F[0][2]*v[242]+F[0][1]*v[283];
    v[176]=v[175]+v[137]*v[343]-v[349];
    v[313]=2e0*F[1][0]*v[176]+F[1][2]*v[238]+F[1][1]*v[278];
    v[291]=2e0*F[0][0]*v[176]+F[0][2]*v[238]+F[0][1]*v[278];
    v[173]=v[172]+v[135]*v[343]-v[350];
    v[289]=2e0*F[0][0]*v[173]+F[0][2]*v[231]+F[0][1]*v[270];
    v[171]=v[170]+v[134]*v[343]-v[351];
    v[288]=2e0*F[0][0]*v[171]+F[0][2]*v[227]+F[0][1]*v[265];
    v[162]=v[138]*v[151]+v[177]-v[352];
    v[320]=2e0*F[1][1]*v[162]+F[1][2]*v[207]+F[1][0]*v[283];
    v[300]=2e0*F[0][1]*v[162]+F[0][2]*v[207]+F[0][0]*v[283];
    v[160]=v[136]*v[151]+v[159]-v[353];
    v[318]=2e0*F[1][1]*v[160]+F[1][2]*v[199]+F[1][0]*v[274];
    v[298]=2e0*F[0][1]*v[160]+F[0][2]*v[199]+F[0][0]*v[274];
    v[158]=v[135]*v[151]+v[172]-v[354];
    v[297]=2e0*F[0][1]*v[158]+F[0][2]*v[198]+F[0][0]*v[270];
    v[295]=F[0][2]*v[190]+F[0][0]*v[261]+2e0*F[0][1]*(v[133]*v[151]+v[155]-v[355]);
    v[148]=v[137]*v[140]+v[175]-v[356];
    v[324]=2e0*F[1][2]*v[148]+F[1][1]*v[203]+F[1][0]*v[238];
    v[306]=2e0*F[0][2]*v[148]+F[0][1]*v[203]+F[0][0]*v[238];
    v[147]=v[136]*v[140]+v[159]-v[357];
    v[323]=2e0*F[1][2]*v[147]+F[1][1]*v[199]+F[1][0]*v[235];
    v[305]=2e0*F[0][2]*v[147]+F[0][1]*v[199]+F[0][0]*v[235];
    v[303]=F[0][1]*v[194]+F[0][0]*v[227]+2e0*F[0][2]*(v[134]*v[140]+v[170]-v[358]);
    v[302]=F[0][1]*v[190]+F[0][0]*v[224]+(v[133]*v[140]+v[155]-v[359])*v[366];
    v[36]=v[140]*v[37]+((1e0+(*C)*(-1e0+v[15]+v[20]))*v[4])/2e0;
    v[321]=2e0*v[36];
    v[40]=v[151]*v[37]+((1e0+(*C)*(-1e0+v[15]+v[26]))*v[4])/2e0;
    v[315]=2e0*v[40];
    v[43]=v[343]*v[37]+((1e0+(*C)*(-1e0+v[20]+v[26]))*v[4])/2e0;
    v[308]=2e0*v[43];
    v[45]=-2e0*v[139]*v[24]+v[336]*v[37];
    v[330]=2e0*F[2][1]*v[162]+F[2][2]*v[207]+F[2][0]*v[283]+v[45];
    v[317]=2e0*F[1][1]*v[158]+F[1][2]*v[198]+F[1][0]*v[270]+v[45];
    v[294]=F[0][2]*v[189]+F[0][0]*v[257]+v[364]*(v[132]*v[151]+v[167]+v[365])+v[45];
    v[46]=-2e0*v[139]*v[22]+v[337]*v[37];
    v[328]=2e0*F[2][0]*v[178]+F[2][2]*v[242]+F[2][1]*v[283]+v[46];
    v[311]=2e0*F[1][0]*v[173]+F[1][2]*v[231]+F[1][1]*v[270]+v[46];
    v[286]=F[0][2]*v[220]+F[0][1]*v[257]+2e0*F[0][0]*(v[167]+v[132]*v[343]+v[362])+v[46];
    v[47]=-2e0*v[139]*v[18]+v[338]*v[37];
    v[327]=2e0*F[2][0]*v[176]+F[2][2]*v[238]+F[2][1]*v[278]+v[47];
    v[310]=2e0*F[1][0]*v[171]+F[1][2]*v[227]+F[1][1]*v[265]+v[47];
    v[285]=F[0][1]*v[252]+2e0*F[0][0]*(v[165]+v[131]*v[343]+v[363])+v[47]+F[0][2]*(v[131]*v[337]+v[37]*
                                                                                   (2e0*F[0][2]*v[18]+2e0*v[255]-2e0*v[22]*v[66]));
    Pmat[0][0]=F[0][0]*v[308]+F[0][2]*v[46]+F[0][1]*v[47];
    Pmat[0][1]=F[0][1]*v[315]+F[0][2]*v[45]+F[0][0]*v[47];
    Pmat[0][2]=F[0][2]*v[321]+F[0][1]*v[45]+F[0][0]*v[46];
    Pmat[1][0]=2e0*F[1][0]*v[43]+F[1][2]*v[46]+F[1][1]*v[47];
    Pmat[1][1]=2e0*F[1][1]*v[40]+F[1][2]*v[45]+F[1][0]*v[47];
    Pmat[1][2]=2e0*F[1][2]*v[36]+F[1][1]*v[45]+F[1][0]*v[46];
    Pmat[2][0]=F[2][0]*v[308]+F[2][2]*v[46]+F[2][1]*v[47];
    Pmat[2][1]=F[2][1]*v[315]+F[2][2]*v[45]+F[2][0]*v[47];
    Pmat[2][2]=F[2][2]*v[321]+F[2][1]*v[45]+F[2][0]*v[46];
    Amat[0][0][0][0]=v[308]+v[129]*v[361]+F[0][2]*(-v[167]+v[129]*v[337]-v[362])+F[0][1]*(-v[165]
                                                                                          +v[129]*v[338]-v[363]);
    Amat[0][0][0][1]=v[285];
    Amat[0][0][0][2]=v[286];
    Amat[0][0][1][0]=v[287];
    Amat[0][0][1][1]=v[288];
    Amat[0][0][1][2]=v[289];
    Amat[0][0][2][0]=v[290];
    Amat[0][0][2][1]=v[291];
    Amat[0][0][2][2]=v[292];
    Amat[0][1][0][0]=v[285];
    Amat[0][1][0][1]=F[0][0]*v[252]+v[315]+v[131]*v[151]*v[364]+F[0][2]*(-v[167]+v[131]*v[336]-v[365]);
    Amat[0][1][0][2]=v[294];
    Amat[0][1][1][0]=v[295];
    Amat[0][1][1][1]=v[296];
    Amat[0][1][1][2]=v[297];
    Amat[0][1][2][0]=v[298];
    Amat[0][1][2][1]=v[299];
    Amat[0][1][2][2]=v[300];
    Amat[0][2][0][0]=v[286];
    Amat[0][2][0][1]=v[294];
    Amat[0][2][0][2]=F[0][1]*v[189]+F[0][0]*v[220]+v[321]+v[132]*v[140]*v[366];
    Amat[0][2][1][0]=v[302];
    Amat[0][2][1][1]=v[303];
    Amat[0][2][1][2]=v[304];
    Amat[0][2][2][0]=v[305];
    Amat[0][2][2][1]=v[306];
    Amat[0][2][2][2]=v[307];
    Amat[1][0][0][0]=v[287];
    Amat[1][0][0][1]=v[295];
    Amat[1][0][0][2]=v[302];
    Amat[1][0][1][0]=2e0*F[1][0]*v[169]+F[1][2]*v[224]+F[1][1]*v[261]+v[308];
    Amat[1][0][1][1]=v[310];
    Amat[1][0][1][2]=v[311];
    Amat[1][0][2][0]=v[312];
    Amat[1][0][2][1]=v[313];
    Amat[1][0][2][2]=v[314];
    Amat[1][1][0][0]=v[288];
    Amat[1][1][0][1]=v[296];
    Amat[1][1][0][2]=v[303];
    Amat[1][1][1][0]=v[310];
    Amat[1][1][1][1]=2e0*F[1][1]*v[157]+F[1][2]*v[194]+F[1][0]*v[265]+v[315];
    Amat[1][1][1][2]=v[317];
    Amat[1][1][2][0]=v[318];
    Amat[1][1][2][1]=v[319];
    Amat[1][1][2][2]=v[320];
    Amat[1][2][0][0]=v[289];
    Amat[1][2][0][1]=v[297];
    Amat[1][2][0][2]=v[304];
    Amat[1][2][1][0]=v[311];
    Amat[1][2][1][1]=v[317];
    Amat[1][2][1][2]=2e0*F[1][2]*v[146]+F[1][1]*v[198]+F[1][0]*v[231]+v[321];
    Amat[1][2][2][0]=v[323];
    Amat[1][2][2][1]=v[324];
    Amat[1][2][2][2]=v[325];
    Amat[2][0][0][0]=v[290];
    Amat[2][0][0][1]=v[298];
    Amat[2][0][0][2]=v[305];
    Amat[2][0][1][0]=v[312];
    Amat[2][0][1][1]=v[318];
    Amat[2][0][1][2]=v[323];
    Amat[2][0][2][0]=2e0*F[2][0]*v[174]+F[2][2]*v[235]+F[2][1]*v[274]+v[308];
    Amat[2][0][2][1]=v[327];
    Amat[2][0][2][2]=v[328];
    Amat[2][1][0][0]=v[291];
    Amat[2][1][0][1]=v[299];
    Amat[2][1][0][2]=v[306];
    Amat[2][1][1][0]=v[313];
    Amat[2][1][1][1]=v[319];
    Amat[2][1][1][2]=v[324];
    Amat[2][1][2][0]=v[327];
    Amat[2][1][2][1]=2e0*F[2][1]*v[161]+F[2][2]*v[203]+F[2][0]*v[278]+v[315];
    Amat[2][1][2][2]=v[330];
    Amat[2][2][0][0]=v[292];
    Amat[2][2][0][1]=v[300];
    Amat[2][2][0][2]=v[307];
    Amat[2][2][1][0]=v[314];
    Amat[2][2][1][1]=v[320];
    Amat[2][2][1][2]=v[325];
    Amat[2][2][2][0]=v[328];
    Amat[2][2][2][1]=v[330];
    Amat[2][2][2][2]=2e0*F[2][2]*v[149]+F[2][1]*v[207]+F[2][0]*v[242]+v[321];
};

    
/*************************************************************
 * AceGen    6.921 MacOSX (29 Jan 19)                         *
 *           Co. J. Korelc  2013           17 Jul 19 15:09:55 *
 **************************************************************
 User     : Full professional version
 Notebook : st_venant_kirchhoff_3d
 Evaluation time                 : 3 s     Mode  : Optimal
 Number of formulae              : 91      Method: Automatic
 Subroutine                      : stvk3d size: 2846
 Total size of Mathematica  code : 2846 subexpressions
 Total size of C code            : 5830 bytes */

/******************* S U B R O U T I N E *********************/
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::stvk3d(double* v,double (*lam),double (*mue),double** F
            ,double** Pmat,double**** Amat)
{
    v[169]=Power(F[0][0],2);
    v[168]=2e0*(*mue);
    v[167]=Power(F[0][2],2);
    v[166]=F[2][2]*(*mue);
    v[165]=F[2][1]*(*mue);
    v[164]=F[2][0]*(*mue);
    v[163]=F[1][2]*(*mue);
    v[162]=F[1][1]*(*mue);
    v[161]=F[1][0]*(*mue);
    v[88]=F[0][0]*(*mue);
    v[116]=F[0][0]*v[88];
    v[70]=F[0][1]*(*lam);
    v[93]=F[0][1]*(*mue);
    v[117]=F[0][1]*v[93];
    v[71]=F[0][2]*(*lam);
    v[105]=(*mue)*v[167];
    v[72]=F[1][0]*(*lam);
    v[85]=2e0*v[161]+v[72];
    v[142]=F[1][0]*v[161];
    v[121]=F[0][0]*v[161];
    v[73]=F[1][1]*(*lam);
    v[100]=F[0][1]*v[161]+F[0][0]*v[73];
    v[82]=2e0*v[162]+v[73];
    v[143]=F[1][1]*v[162];
    v[122]=F[0][1]*v[162];
    v[108]=F[0][0]*v[162]+F[0][1]*v[72];
    v[74]=F[1][2]*(*lam);
    v[111]=F[0][2]*v[162]+F[0][1]*v[74];
    v[101]=F[0][2]*v[161]+F[0][0]*v[74];
    v[79]=2e0*v[163]+v[74];
    v[123]=v[121]+v[122]+F[0][2]*v[79];
    v[135]=F[1][2]*v[163];
    v[120]=F[0][1]*v[163]+F[0][2]*v[73];
    v[119]=F[0][0]*v[163]+F[0][2]*v[72];
    v[109]=F[0][2]*v[163];
    v[110]=v[109]+v[121]+F[0][1]*v[82];
    v[99]=v[109]+v[122]+F[0][0]*v[85];
    v[75]=F[2][0]*(*lam);
    v[86]=2e0*v[164]+v[75];
    v[156]=F[2][0]*v[164];
    v[147]=F[1][0]*v[164];
    v[126]=F[0][0]*v[164];
    v[76]=F[2][1]*(*lam);
    v[133]=F[1][1]*v[164]+F[1][0]*v[76];
    v[103]=F[0][1]*v[164]+F[0][0]*v[76];
    v[83]=2e0*v[165]+v[76];
    v[157]=F[2][1]*v[165];
    v[148]=F[1][1]*v[165];
    v[138]=F[1][0]*v[165]+F[1][1]*v[75];
    v[127]=F[0][1]*v[165];
    v[112]=F[0][0]*v[165]+F[0][1]*v[75];
    v[77]=F[2][2]*(*lam);
    v[141]=F[1][2]*v[165]+F[1][1]*v[77];
    v[134]=F[1][2]*v[164]+F[1][0]*v[77];
    v[115]=F[0][2]*v[165]+F[0][1]*v[77];
    v[104]=F[0][2]*v[164]+F[0][0]*v[77];
    v[80]=2e0*v[166]+v[77];
    v[149]=v[147]+v[148]+F[1][2]*v[80];
    v[128]=v[126]+v[127]+F[0][2]*v[80];
    v[153]=F[2][2]*v[166];
    v[146]=F[1][1]*v[166]+F[1][2]*v[76];
    v[145]=F[1][0]*v[166]+F[1][2]*v[75];
    v[139]=F[1][2]*v[166];
    v[140]=v[139]+v[147]+F[1][1]*v[83];
    v[132]=v[139]+v[148]+F[1][0]*v[86];
    v[125]=F[0][1]*v[166]+F[0][2]*v[76];
    v[124]=F[0][0]*v[166]+F[0][2]*v[75];
    v[113]=F[0][2]*v[166];
    v[114]=v[113]+v[126]+F[0][1]*v[83];
    v[102]=v[113]+v[127]+F[0][0]*v[86];
    v[24]=(-1e0+Power(F[1][0],2)+Power(F[2][0],2)+v[169])/2e0;
    v[28]=(-1e0+Power(F[0][1],2)+Power(F[1][1],2)+Power(F[2][1],2))/2e0;
    v[32]=(-1e0+Power(F[1][2],2)+Power(F[2][2],2)+v[167])/2e0;
    v[36]=(*lam)*(v[24]+v[28]+v[32]);
    v[35]=2e0*(*mue)*v[32]+v[36];
    v[37]=2e0*(*mue)*v[28]+v[36];
    v[38]=2e0*(*mue)*v[24]+v[36];
    v[39]=(F[0][0]*F[0][2]+F[1][0]*F[1][2]+F[2][0]*F[2][2])*(*mue);
    v[152]=F[2][2]*v[164]+v[39]+F[2][0]*v[77];
    v[131]=F[1][2]*v[161]+v[39]+F[1][0]*v[74];
    v[98]=v[39]+F[0][0]*v[71]+F[0][2]*v[88];
    v[40]=(F[0][1]*F[0][2]+F[1][1]*F[1][2]+F[2][1]*F[2][2])*(*mue);
    v[155]=F[2][2]*v[165]+v[40]+F[2][1]*v[77];
    v[137]=F[1][2]*v[162]+v[40]+F[1][1]*v[74];
    v[107]=v[40]+F[0][1]*v[71]+F[0][2]*v[93];
    v[41]=(F[0][0]*F[0][1]+F[1][0]*F[1][1]+F[2][0]*F[2][1])*(*mue);
    v[151]=F[2][1]*v[164]+v[41]+F[2][0]*v[76];
    v[130]=F[1][1]*v[161]+v[41]+F[1][0]*v[73];
    v[97]=v[41]+F[0][0]*v[70]+F[0][1]*v[88];
    Pmat[0][0]=F[0][0]*v[38]+F[0][2]*v[39]+F[0][1]*v[41];
    Pmat[0][1]=F[0][1]*v[37]+F[0][2]*v[40]+F[0][0]*v[41];
    Pmat[0][2]=F[0][2]*v[35]+F[0][0]*v[39]+F[0][1]*v[40];
    Pmat[1][0]=F[1][0]*v[38]+F[1][2]*v[39]+F[1][1]*v[41];
    Pmat[1][1]=F[1][1]*v[37]+F[1][2]*v[40]+F[1][0]*v[41];
    Pmat[1][2]=F[1][2]*v[35]+F[1][0]*v[39]+F[1][1]*v[40];
    Pmat[2][0]=F[2][0]*v[38]+F[2][2]*v[39]+F[2][1]*v[41];
    Pmat[2][1]=F[2][1]*v[37]+F[2][2]*v[40]+F[2][0]*v[41];
    Pmat[2][2]=F[2][2]*v[35]+F[2][0]*v[39]+F[2][1]*v[40];
    Amat[0][0][0][0]=v[105]+v[117]+((*lam)+v[168])*v[169]+v[38];
    Amat[0][0][0][1]=v[97];
    Amat[0][0][0][2]=v[98];
    Amat[0][0][1][0]=v[99];
    Amat[0][0][1][1]=v[100];
    Amat[0][0][1][2]=v[101];
    Amat[0][0][2][0]=v[102];
    Amat[0][0][2][1]=v[103];
    Amat[0][0][2][2]=v[104];
    Amat[0][1][0][0]=v[97];
    Amat[0][1][0][1]=v[105]+v[116]+v[37]+F[0][1]*(v[70]+2e0*v[93]);
    Amat[0][1][0][2]=v[107];
    Amat[0][1][1][0]=v[108];
    Amat[0][1][1][1]=v[110];
    Amat[0][1][1][2]=v[111];
    Amat[0][1][2][0]=v[112];
    Amat[0][1][2][1]=v[114];
    Amat[0][1][2][2]=v[115];
    Amat[0][2][0][0]=v[98];
    Amat[0][2][0][1]=v[107];
    Amat[0][2][0][2]=v[116]+v[117]+v[35]+F[0][2]*(F[0][2]*v[168]+v[71]);
    Amat[0][2][1][0]=v[119];
    Amat[0][2][1][1]=v[120];
    Amat[0][2][1][2]=v[123];
    Amat[0][2][2][0]=v[124];
    Amat[0][2][2][1]=v[125];
    Amat[0][2][2][2]=v[128];
    Amat[1][0][0][0]=v[99];
    Amat[1][0][0][1]=v[108];
    Amat[1][0][0][2]=v[119];
    Amat[1][0][1][0]=v[135]+v[143]+v[38]+F[1][0]*v[85];
    Amat[1][0][1][1]=v[130];
    Amat[1][0][1][2]=v[131];
    Amat[1][0][2][0]=v[132];
    Amat[1][0][2][1]=v[133];
    Amat[1][0][2][2]=v[134];
    Amat[1][1][0][0]=v[100];
    Amat[1][1][0][1]=v[110];
    Amat[1][1][0][2]=v[120];
    Amat[1][1][1][0]=v[130];
    Amat[1][1][1][1]=v[135]+v[142]+v[37]+F[1][1]*v[82];
    Amat[1][1][1][2]=v[137];
    Amat[1][1][2][0]=v[138];
    Amat[1][1][2][1]=v[140];
    Amat[1][1][2][2]=v[141];
    Amat[1][2][0][0]=v[101];
    Amat[1][2][0][1]=v[111];
    Amat[1][2][0][2]=v[123];
    Amat[1][2][1][0]=v[131];
    Amat[1][2][1][1]=v[137];
    Amat[1][2][1][2]=v[142]+v[143]+v[35]+F[1][2]*v[79];
    Amat[1][2][2][0]=v[145];
    Amat[1][2][2][1]=v[146];
    Amat[1][2][2][2]=v[149];
    Amat[2][0][0][0]=v[102];
    Amat[2][0][0][1]=v[112];
    Amat[2][0][0][2]=v[124];
    Amat[2][0][1][0]=v[132];
    Amat[2][0][1][1]=v[138];
    Amat[2][0][1][2]=v[145];
    Amat[2][0][2][0]=v[153]+v[157]+v[38]+F[2][0]*v[86];
    Amat[2][0][2][1]=v[151];
    Amat[2][0][2][2]=v[152];
    Amat[2][1][0][0]=v[103];
    Amat[2][1][0][1]=v[114];
    Amat[2][1][0][2]=v[125];
    Amat[2][1][1][0]=v[133];
    Amat[2][1][1][1]=v[140];
    Amat[2][1][1][2]=v[146];
    Amat[2][1][2][0]=v[151];
    Amat[2][1][2][1]=v[153]+v[156]+v[37]+F[2][1]*v[83];
    Amat[2][1][2][2]=v[155];
    Amat[2][2][0][0]=v[104];
    Amat[2][2][0][1]=v[115];
    Amat[2][2][0][2]=v[128];
    Amat[2][2][1][0]=v[134];
    Amat[2][2][1][1]=v[141];
    Amat[2][2][1][2]=v[149];
    Amat[2][2][2][0]=v[152];
    Amat[2][2][2][1]=v[155];
    Amat[2][2][2][2]=v[156]+v[157]+v[35]+F[2][2]*v[80];
};
/*************************************************************
 * AceGen    6.921 MacOSX (29 Jan 19)                         *
 *           Co. J. Korelc  2013           17 Jul 19 16:01:42 *
 **************************************************************
 User     : Full professional version
 Notebook : st_venant_kirchhoff_2d
 Evaluation time                 : 1 s     Mode  : Optimal
 Number of formulae              : 25      Method: Automatic
 Subroutine                      : stvk2d size: 772
 Total size of Mathematica  code : 772 subexpressions
 Total size of C code            : 1672 bytes */


/******************* S U B R O U T I N E *********************/
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::stvk2d(double* v, double (*lam),double (*mue),double** F
            ,double** Pmat,double**** Amat)
{
    v[43]=F[0][0]*F[1][0];
    v[42]=F[0][1]*F[1][1];
    v[37]=Power(F[0][0],2);
    v[12]=F[0][0]/2e0;
    v[36]=Power(F[0][1],2);
    v[34]=F[0][0]*F[0][1];
    v[11]=F[0][1]/2e0;
    v[27]=Power(F[1][0],2);
    v[31]=-1e0+v[27]+v[37];
    v[14]=F[1][0]/2e0;
    v[25]=Power(F[1][1],2);
    v[26]=-1e0+v[25]+v[36];
    v[23]=F[1][0]*F[1][1];
    v[22]=v[23]+v[34];
    v[35]=(*lam)*v[34]+(*mue)*(v[22]+v[34]);
    v[24]=(*lam)*v[23]+(*mue)*(v[22]+v[23]);
    v[21]=(*lam)*v[42]+2e0*(*mue)*(2e0*v[12]*v[14]+v[42]);
    v[20]=F[0][0]*F[1][1]*(*lam)+4e0*(*mue)*v[11]*v[14];
    v[13]=F[1][1]/2e0;
    v[30]=F[0][1]*F[1][0]*(*lam)+4e0*(*mue)*v[12]*v[13];
    v[29]=(*lam)*v[43]+2e0*(*mue)*(2e0*v[11]*v[13]+v[43]);
    v[44]=2e0*v[22];
    v[32]=((*lam)*(v[26]+v[31]))/2e0;
    Pmat[0][0]=F[0][0]*v[32]+(*mue)*(F[0][0]*v[31]+v[11]*v[44]);
    Pmat[0][1]=F[0][1]*v[32]+(*mue)*(F[0][1]*v[26]+v[12]*v[44]);
    Pmat[1][0]=F[1][0]*v[32]+(*mue)*(F[1][0]*v[31]+v[13]*v[44]);
    Pmat[1][1]=F[1][1]*v[32]+(*mue)*(F[1][1]*v[26]+v[14]*v[44]);
    Amat[0][0][0][0]=v[32]+(*lam)*v[37]+(*mue)*(v[31]+v[36]+2e0*v[37]);
    Amat[0][0][0][1]=v[35];
    Amat[0][0][1][0]=v[29];
    Amat[0][0][1][1]=v[20];
    Amat[0][1][0][0]=v[35];
    Amat[0][1][0][1]=v[32]+(*lam)*v[36]+(*mue)*(v[26]+2e0*v[36]+v[37]);
    Amat[0][1][1][0]=v[30];
    Amat[0][1][1][1]=v[21];
    Amat[1][0][0][0]=v[29];
    Amat[1][0][0][1]=v[30];
    Amat[1][0][1][0]=(*lam)*v[27]+(*mue)*(v[25]+2e0*v[27]+v[31])+v[32];
    Amat[1][0][1][1]=v[24];
    Amat[1][1][0][0]=v[20];
    Amat[1][1][0][1]=v[21];
    Amat[1][1][1][0]=v[24];
    Amat[1][1][1][1]=(*lam)*v[25]+(*mue)*(2e0*v[25]+v[26]+v[27])+v[32];
};

    
/*************************************************************
 * AceGen    6.818 Linux (13 Sep 17)                          *
 *           Co. J. Korelc  2013           21 May 19 12:35:39 *
 **************************************************************
 User     : Full professional version
 Notebook : T2_TPM_up_LE_Iso_Gal
 Evaluation time                 : 6 s     Mode  : Optimal
 Number of formulae              : 324     Method: Automatic
 Subroutine                      : SKR size: 3659
 Subroutine                      : SPP size: 2226
 Total size of Mathematica  code : 5885 subexpressions
 Total size of C code            : 16126 bytes */
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::SMTSetElSpecBiot(ElementSpec *es,int *idata,int ic,int ng, vec_dbl_Type& paraVec)
{
    int dim = domainVec_[0]->getDimension();
    int intc,nd,i;FILE *SMSFile;
    static int pn[24]={1, 4, 2, 5, 3, 6, 0, 1, 4, 6, -1, 4, 2, 5, -1, 6, 4, 5, -1, 6, 5, 3, -1, 0};
    static int dof[9]={2, 2, 2, 2, 2, 2, 1, 1, 1};
    static int nsto[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    static int ndat[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};

    static char *nid[]={"D","D","D","D","D","D",
        "p","p","p"};

    
    //Name der Daten in es->Data
    static char *gdcs[]={   "$[Gamma]$NM -Newmark Parameter",
                            "$[Beta]$NM -Newmark Parameter",
                            "ns0s -initial volume fraction solid",
                            "kL -Darcy parameter",
                            "E -Youngs modulus",
                            "$[Nu]$ -Poissons Ratio"};
        
    //soll in es->Data
//    static double defd[]={  gamma,
//                            beta,
//                            ns_init_solid,
//                            darcyPara,
//                            youngsModulus,
//                            poissonRatio,
//                            0e0 /*this is not used*/};
    static char *gpcs[]={""};
    static char *npcs[]={"p","u1","u2","u3","$[Sigma]$11","$[Sigma]$22",
        "$[Sigma]$33","$[Sigma]$12","$[Sigma]$23","$[Sigma]$31","seepage1","seepage2"};
    static char *sname[]={""};
    static char *idname[]={""};
    static int idindex[1];
    static char *rdname[]={""};
    static char *cswitch[]={""};
    static int iswitch[1]={0};
    static double dswitch[1]={0e0};
    static int rdindex[1];
    static int nspecs[9];
    static double version[3]={6.818,6.818,11.1};
    static double pnweights[9]={1e0,1e0,1e0,1e0,1e0,1e0,
        0e0,0e0,0e0};
    //bestimmt die Koordinaten fuer das Refernezelement
    static double rnodes[27]={  1e0,0e0,0e0,
                                0e0,1e0,0e0,
                                0e0,0e0,0e0,
                                0.5e0,0e0,0e0,
                                0.5e0,0.5e0,0e0,
                                0e0,0.5e0,0e0,
                                1e0,0e0,0e0,
                                0e0,1e0,0e0,
                                0e0,0e0,0e0 };
    
    es->ReferenceNodes=rnodes;
    if(ng==-1) es->Data= &paraVec[0];
    es->id.NoDomainData=6;
    es->Code="T2_TPM_up_LE_Iso_Gal";
    es->Version=version;
    es->MainTitle="";
    es->SubTitle="";
    es->SubSubTitle="";
    es->Bibliography="";
    //WorkingVectoSize fehlt hier. Wird fuer v[WorkingVectorSize] benoetigt.
    es->id.NoDimensions = dim ;
    es->id.NoDOFGlobal=15; //anpassen
    es->id.NoDOFCondense=0;
    es->id.NoNodes=9;//anpassen
    es->id.NoSegmentPoints=23;//was tut das?
    es->Segments=pn;
    es->PostNodeWeights=pnweights;
    es->id.NoIntSwitch=0;//was tut das?
    es->IntSwitch=iswitch;//was tut das?
    es->id.DemoLimitation=0;//was tut das?
    es->id.NoDoubleSwitch=0;//was tut das?
    es->DoubleSwitch=dswitch;//was tut das?
    es->id.NoCharSwitch=0;//was tut das?
    es->CharSwitch=cswitch;//was tut das?
    es->DOFGlobal=dof;
    es->NodeID=nid;
    es->id.NoGPostData=0;
    es->id.NoNPostData=12;//was tut das?
    es->id.SymmetricTangent=0;
    es->id.CreateDummyNodes=0;
    es->id.PostIterationCall=0;//was tut das?
    es->id.DOFScaling=0;//was tut das?
    es->Topology="XX";
    es->DomainDataNames=gdcs;//warum? dies sind die Parameter
    es->GPostNames=gpcs;
    es->NPostNames=npcs;
    es->AdditionalNodes="{}&";
    es->AdditionalGraphics="{}&";
    es->MMAInitialisation=MMAInitialisationCode;
    es->MMANextStep="";
    es->MMAStepBack="";
    es->MMAPreIteration="";
    es->IDataNames=idname;
    es->IDataIndex=idindex;
    es->RDataNames=rdname;
    es->RDataIndex=rdindex;
    es->id.NoIData=0;
    es->id.NoRData=0;
    es->id.ShapeSensitivity=0;
    es->id.EBCSensitivity=0;
    es->id.SensitivityOrder=0;
    es->id.NoSensNames=0;
    es->SensitivityNames=sname;
    es->NodeSpecs=nspecs;
    //es->user.SPP=SPP; //not used atm
//    es->user.SKR=SKR; //do we need this?
    
    es->id.DefaultIntegrationCode=35;
    if(ic==-1){intc=35;} else {intc=ic;};
    es->id.IntCode=intc;
    // we should know the suitable length of idata
    es->IntPoints = SMTMultiIntPoints(&intc,idata,&es->id.NoIntPoints,
                      &es->id.NoIntPointsA,&es->id.NoIntPointsB,&es->id.NoIntPointsC,1);
    es->id.NoAdditionalData=(int)(0);
    //Laenge von NoTimeStorage muss zu hp bzw. ht passen, siehe ed->hp und ed->ht in SKR
    es->id.NoTimeStorage=(int)(24);
    es->id.NoElementData=(int)(0);
    
    es->NoNodeStorage=nsto;es->NoNodeData=ndat;
    
};

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::SMTSetElSpecBiotStVK(ElementSpec *es,int *idata,int ic,int ng, vec_dbl_Type& paraVec)
{
    
  int intc,nd,i;FILE *SMSFile;

  static int pn[9]={1, 2, 3, 0, 1, 2, 3, -1, 0};
  static int dof[9]={2, 2, 2, 2, 2, 2, 1, 1, 1};
  static int nsto[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};

  static int ndat[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};

  static char *nid[]={"D","D","D","D","D","D",
                       "p","p","p"};
  static char *gdcs[]={ "E -Youngs modulus",
                        "$[Nu]$ -Poissons Ratio",
                        "ns0s -initial volume fraction solid",
                        "kL -Darcy parameter",
                        "$[Gamma]$NM -Newmark Parameter",
                        "$[Beta]$NM -Newmark Parameter"};
  
  static char *gpcs[]={"p","u1","u2","u3","$[Sigma]$11","$[Sigma]$22",
                       "$[Sigma]$33","$[Sigma]$12","$[Sigma]$23","$[Sigma]$31","seepage1","seepage2"};
  static char *npcs[]={""};
  static char *sname[]={""};
  static char *idname[]={""};
  static int idindex[1];
  static char *rdname[]={""};
  static char *cswitch[]={""};
  static int iswitch[1]={0};
  static double dswitch[1]={0e0};
  static int rdindex[1];
  static int nspecs[9];
  static double version[3]={6.818,6.818,11.1};
  static double pnweights[9]={1e0,1e0,1e0,0e0,0e0,0e0,
  0e0,0e0,0e0};
  static double rnodes[27]={1e0,0e0,0e0,0e0,1e0,0e0,
  0e0,0e0,0e0,0.5e0,0e0,0e0,
  0.5e0,0.5e0,0e0,0e0,0.5e0,0e0,
  1e0,0e0,0e0,0e0,1e0,0e0,
  0e0,0e0,0e0};
  es->ReferenceNodes=rnodes;
  if(ng==-1) es->Data= &paraVec[0];
  es->id.NoDomainData=6;
  es->Code="T2T1_up_nichtlinear_iso_gal_2020_02_05";es->Version=version;
  es->MainTitle="";
  es->SubTitle="";
  es->SubSubTitle="";
  es->Bibliography="";
  es->id.NoDimensions=2;es->id.NoDOFGlobal=15;es->id.NoDOFCondense=0;es->id.NoNodes=9;
  es->id.NoSegmentPoints=8;es->Segments=pn;es->PostNodeWeights=pnweights;
  es->id.NoIntSwitch=0;es->IntSwitch=iswitch;es->id.DemoLimitation=0;
  es->id.NoDoubleSwitch=0;es->DoubleSwitch=dswitch;
  es->id.NoCharSwitch=0;es->CharSwitch=cswitch;
  es->DOFGlobal=dof;es->NodeID=nid;es->id.NoGPostData=12;es->id.NoNPostData=0;
  es->id.SymmetricTangent=0;es->id.CreateDummyNodes=0;es->id.PostIterationCall=0;es->id.DOFScaling=0;
  es->Topology="XX";es->DomainDataNames=gdcs;es->GPostNames=gpcs;es->NPostNames=npcs;
  es->AdditionalNodes="{}&";
  es->AdditionalGraphics="{}&";
  es->MMAInitialisation=MMAInitialisationCode;
  es->MMANextStep="";
  es->MMAStepBack="";
  es->MMAPreIteration="";
  es->IDataNames=idname;es->IDataIndex=idindex;es->RDataNames=rdname;es->RDataIndex=rdindex;
  es->id.NoIData=0;es->id.NoRData=0;
  es->id.ShapeSensitivity=0; es->id.EBCSensitivity=0;es->id.SensitivityOrder=0;
  es->id.NoSensNames=0;es->SensitivityNames=sname;es->NodeSpecs=nspecs;
//  es->user.SPP=SPP;
//    es->user.SKR=SKR;

  es->id.DefaultIntegrationCode=35;
  if(ic==-1){intc=35;} else {intc=ic;};
  es->id.IntCode=intc;
  es->IntPoints = SMTMultiIntPoints(&intc,idata,&es->id.NoIntPoints,
                      &es->id.NoIntPointsA,&es->id.NoIntPointsB,&es->id.NoIntPointsC,1);

  es->id.NoAdditionalData=(int)(0);
  es->id.NoTimeStorage=(int)(24);
  es->id.NoElementData=(int)(0);

  es->NoNodeStorage=nsto;es->NoNodeData=ndat;

};

/*************************************************************
* AceGen    7.114 Linux (9 Jul 20)                           *
*           Co. J. Korelc  2020           4 Sep 20 16:26:52  *
**************************************************************
User     : Full professional version
Notebook : O2O1_up_lin_el_iso_gal_Newmark_mesh_koeln_1D_Ehlers_2020_09_04
Evaluation time                 : 23 s    Mode  : Optimal
Number of formulae              : 521     Method: Automatic
Subroutine                      : SKR size: 7982
Subroutine                      : SPP size: 2615
Total size of Mathematica  code : 10597 subexpressions
Total size of C code            : 31924 bytes */
//Added 07.09.2020 CH
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::SMTSetElSpecBiot3D(ElementSpec *es,int *idata,int ic,int ng, vec_dbl_Type& paraVec)
{
    int intc,nd,i;FILE *SMSFile;
    static int pn[33]={1, 2, 4, 0, 1, 4, 3, 0, 1, 2, 3, 0, 2, 3, 4, 0, 1, 2, 4, -1, 1, 4, 3, -2, 1, 2, 3, -3, 2, 3, 4, -4, 0};
    static int dof[14]={3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1};
    static int nsto[14]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    static int ndat[14]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    //static int esdat[1]={0}; only used in new release
    
    static char *nid[]={"D","D","D","D","D","D",
        "D","D","D","D","p","p",
        "p","p"};
    static char *gdcs[]={"E -Youngs modulus","$[Nu]$ -Poissons Ratio","bx -body force","by -body force","bz -body force","ns0s -initial volume fraction solid",
        "kL -Darcy parameter","$[Rho]$SR -effective density solid","$[Rho]$FR -effective density solid","g -gravity","$[Gamma]$NM -Newmark Parameter","$[Beta]$NM -Newmark Parameter"};
    static double defd[]={200e0,0e0,0e0,0e0,0e0,0.67e0,
        0.1e-3,2000e0,1000e0,0.9810000000000001e1,0.5e0,0.25e0,
        0e0};
    static char *gpcs[]={"p","u1","u2","u3","$[Sigma]$11","$[Sigma]$22",
        "$[Sigma]$33","$[Sigma]$12","$[Sigma]$23","$[Sigma]$31","seepage1","seepage2",
        "seepage3"};
    static char *npcs[]={""};
    static char *sname[]={""};
    static char *idname[]={""};
    static int idindex[1];
    static char *rdname[]={""};
    static char *cswitch[]={""};
    static int iswitch[1]={0};
    static double dswitch[1]={0e0};
    static int rdindex[1];
    static int nspecs[14];
    static double version[3]={7.114,7.114,11.1};
    static double pnweights[14]={1e0,1e0,1e0,1e0,0e0,0e0,
        0e0,0e0,0e0,0e0,0e0,0e0,
        0e0,0e0};
    static double rnodes[42]={1e0,0e0,0e0,0e0,1e0,0e0,
        0e0,0e0,1e0,0e0,0e0,0e0,
        0.5e0,0.5e0,0e0,0e0,0.5e0,0.5e0,
        0.5e0,0e0,0.5e0,0.5e0,0e0,0e0,
        0e0,0.5e0,0e0,0e0,0e0,0.5e0,
        1e0,0e0,0e0,0e0,1e0,0e0,
        0e0,0e0,1e0,0e0,0e0,0e0};
    es->ReferenceNodes=rnodes;
    
    if(ng==-1)
        es->Data=&paraVec[0];
    //es->Data=defd;
    es->id.NoDomainData=12;
    es->Code="O2O1_up_lin_el_iso_gal_Newmark_mesh_koeln_1D_Ehlers_2020_09_04";es->Version=version;
    es->MainTitle="";
    es->SubTitle="";
    es->SubSubTitle="";
    es->Bibliography="";
    es->id.NoDimensions=3;es->id.NoDOFGlobal=34;es->id.NoDOFCondense=0;es->id.NoNodes=14;
    es->id.NoSegmentPoints=32;es->Segments=pn;es->PostNodeWeights=pnweights;
    es->id.NoIntSwitch=0;es->IntSwitch=iswitch;
    //es->id.LocalReKe=0; LocalReKe is new and had a different name before (see sms.h(pp)). Was changed with new AceGen/AceFEM release.
    es->id.NoDoubleSwitch=0;es->DoubleSwitch=dswitch;
    es->id.NoCharSwitch=0;es->id.WorkingVectorSize=1817;es->CharSwitch=cswitch;
    es->DOFGlobal=dof;es->NodeID=nid;es->id.NoGPostData=13;es->id.NoNPostData=0;
    es->id.SymmetricTangent=0;es->id.PostIterationCall=0;es->id.DOFScaling=0;
    es->Topology="XX";es->DomainDataNames=gdcs;es->GPostNames=gpcs;es->NPostNames=npcs;
    es->AdditionalNodes="{}&";
    es->AdditionalGraphics="{Null,Null,Null}";
    es->MMAInitialisation=MMAInitialisationCode;
    es->MMANextStep="";
    es->MMAStepBack="";
    es->MMAPreIteration="";
    es->IDataNames=idname;es->IDataIndex=idindex;es->RDataNames=rdname;es->RDataIndex=rdindex;
    es->id.NoIData=0;es->id.NoRData=0;
    es->id.ShapeSensitivity=0; es->id.EBCSensitivity=0;es->id.SensitivityOrder=0;
    es->id.NoSensNames=0;es->SensitivityNames=sname;es->NodeSpecs=nspecs;
//    es->user.SPP=SPP;
//    es->user.SKR=SKR;
    
    es->id.DefaultIntegrationCode=40;
    if(ic==-1){intc=40;} else {intc=ic;};
    es->id.IntCode=intc;
    es->IntPoints = SMTMultiIntPoints(&intc,idata,&es->id.NoIntPoints,
                      &es->id.NoIntPointsA,&es->id.NoIntPointsB,&es->id.NoIntPointsC,1);
    es->id.NoAdditionalData=(int)(0);
    es->id.NoTimeStorage=(int)(60);
    es->id.NoElementData=(int)(0);
    
    
    es->NoNodeStorage=nsto;es->NoNodeData=ndat;
    
    //es->ExtraSensitivityData=esdat;
    //ExtraSensitivityData is new. Was introduced with new AceGen/AceFEM release.
};



template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::SKR_Biot(double* v,ElementSpec *es,ElementData *ed,NodeSpec **ns
                               ,NodeData **nd,double *rdata,int *idata,double *p,double **s)
{
    int i107,i109,i215,i230;
    v[5162]=0e0;
    v[5163]=0e0;
    v[5164]=0e0;
    v[5165]=0e0;
    v[5166]=0e0;
    v[5167]=0e0;
    v[5168]=0e0;
    v[5169]=1e0;
    v[5170]=0e0;
    v[5171]=-1e0;
    v[5172]=0e0;
    v[5173]=0e0;
    v[5174]=0e0;
    v[5175]=0e0;
    v[5176]=0e0;
    v[5132]=0e0;
    v[5133]=0e0;
    v[5134]=0e0;
    v[5135]=0e0;
    v[5136]=0e0;
    v[5137]=0e0;
    v[5138]=0e0;
    v[5139]=1e0;
    v[5140]=0e0;
    v[5141]=0e0;
    v[5142]=0e0;
    v[5143]=-1e0;
    v[5144]=0e0;
    v[5145]=0e0;
    v[5146]=0e0;
    v[5102]=0e0;
    v[5103]=0e0;
    v[5104]=0e0;
    v[5105]=0e0;
    v[5106]=0e0;
    v[5107]=0e0;
    v[5108]=1e0;
    v[5109]=0e0;
    v[5110]=-1e0;
    v[5111]=0e0;
    v[5112]=0e0;
    v[5113]=0e0;
    v[5114]=0e0;
    v[5115]=0e0;
    v[5116]=0e0;
    v[5072]=0e0;
    v[5073]=0e0;
    v[5074]=0e0;
    v[5075]=0e0;
    v[5076]=0e0;
    v[5077]=0e0;
    v[5078]=1e0;
    v[5079]=0e0;
    v[5080]=0e0;
    v[5081]=0e0;
    v[5082]=-1e0;
    v[5083]=0e0;
    v[5084]=0e0;
    v[5085]=0e0;
    v[5086]=0e0;
    v[5057]=0e0;
    v[5058]=0e0;
    v[5059]=0e0;
    v[5060]=0e0;
    v[5061]=0e0;
    v[5062]=0e0;
    v[5063]=0e0;
    v[5064]=0e0;
    v[5065]=0e0;
    v[5066]=0e0;
    v[5067]=0e0;
    v[5068]=0e0;
    v[5069]=1e0;
    v[5070]=0e0;
    v[5071]=-1e0;
    v[5042]=0e0;
    v[5043]=0e0;
    v[5044]=0e0;
    v[5045]=0e0;
    v[5046]=0e0;
    v[5047]=0e0;
    v[5048]=0e0;
    v[5049]=0e0;
    v[5050]=0e0;
    v[5051]=0e0;
    v[5052]=0e0;
    v[5053]=0e0;
    v[5054]=0e0;
    v[5055]=1e0;
    v[5056]=-1e0;
    
    //es->Data sind Materialparameter
    v[191]=es->Data[5];
    v[194]=es->Data[4]/(2e0*(1e0+v[191]));
    v[192]=(2e0*v[191]*v[194])/(1e0-2e0*v[191]);
    v[186]=es->Data[3];
    v[183]=es->Data[2];
    // nd sind Nodalwerte, Anzahl an structs in nd sollte den Knoten entsprechen, bei P2-P1 in 2D also 9
    // In X stehen die Koordinaten, X[0] ist x-Koordinate, X[1] ist y-Koordinate, etc.
    v[1]=nd[0]->X[0];
    v[2]=nd[0]->X[1];
    v[3]=nd[1]->X[0];
    v[4]=nd[1]->X[1];
    v[5]=nd[2]->X[0];
    v[6]=nd[2]->X[1];
    v[7]=nd[3]->X[0];
    v[8]=nd[3]->X[1];
    v[9]=nd[4]->X[0];
    v[10]=nd[4]->X[1];
    v[11]=nd[5]->X[0];
    v[12]=nd[5]->X[1];
    // at ist die Loesung im letzten Newtonschritt.
    v[19]=nd[0]->at[0];
    v[20]=nd[0]->at[1];
    v[21]=nd[1]->at[0];
    v[22]=nd[1]->at[1];
    v[23]=nd[2]->at[0];
    v[24]=nd[2]->at[1];
    v[25]=nd[3]->at[0];
    v[26]=nd[3]->at[1];
    v[27]=nd[4]->at[0];
    v[28]=nd[4]->at[1];
    v[29]=nd[5]->at[0];
    v[30]=nd[5]->at[1];
    v[31]=nd[6]->at[0];
    v[32]=nd[7]->at[0];
    v[33]=nd[8]->at[0];
    v[172]=v[32]-v[33];
    v[171]=v[31]-v[33];
    // ap ist die Loesung im letzten Zeitschritt.
    v[490]=-nd[0]->ap[0]+v[19];
    v[497]=-nd[0]->ap[1]+v[20];
    v[492]=-nd[1]->ap[0]+v[21];
    v[498]=-nd[1]->ap[1]+v[22];
    v[493]=-nd[2]->ap[0]+v[23];
    v[499]=-nd[2]->ap[1]+v[24];
    v[494]=-nd[3]->ap[0]+v[25];
    v[500]=-nd[3]->ap[1]+v[26];
    v[495]=-nd[4]->ap[0]+v[27];
    v[501]=-nd[4]->ap[1]+v[28];
    v[496]=-nd[5]->ap[0]+v[29];
    v[488]=-nd[5]->ap[1]+v[30];
    //rdata ist die Zeitschrittweite, RD_TimeIncrement wird in sms.h definiert, entsprechend wird auch die Laenge von rdata dort definiert. Standard 400, aber auch nicht gesetzt. Wert muss selber initialisiert werden; eventuell kuerzer moeglich.
    v[49]=rdata[RD_TimeIncrement];
    //hp:history previous (timestep); previous solution
    v[50]=ed->hp[0];
    v[51]=ed->hp[1];
    v[52]=ed->hp[2];
    v[53]=ed->hp[3];
    v[54]=ed->hp[4];
    v[55]=ed->hp[5];
    v[56]=ed->hp[6];
    v[57]=ed->hp[7];
    v[58]=ed->hp[8];
    v[59]=ed->hp[9];
    v[60]=ed->hp[10];
    v[61]=ed->hp[11];
    v[62]=ed->hp[12];
    v[63]=ed->hp[13];
    v[64]=ed->hp[14];
    v[65]=ed->hp[15];
    v[66]=ed->hp[16];
    v[67]=ed->hp[17];
    v[68]=ed->hp[18];
    v[69]=ed->hp[19];
    v[70]=ed->hp[20];
    v[71]=ed->hp[21];
    v[72]=ed->hp[22];
    v[73]=ed->hp[23];
    v[75]=es->Data[1];
    v[489]=1e0/(Power(v[49],2)*v[75]);
    v[491]=-((v[49]*v[49])*(0.5e0-v[75]));
    v[79]=es->Data[0]/(v[49]*v[75]);
    v[486]=-(v[49]*v[79]);
    v[487]=(1e0+v[486]/2e0)*v[49];
    v[77]=1e0+v[486];
    v[76]=v[487]*v[62]+v[50]*v[77]+v[490]*v[79];
    v[80]=v[487]*v[68]+v[56]*v[77]+v[497]*v[79];
    v[81]=v[487]*v[63]+v[51]*v[77]+v[492]*v[79];
    v[82]=v[487]*v[69]+v[57]*v[77]+v[498]*v[79];
    v[83]=v[487]*v[64]+v[52]*v[77]+v[493]*v[79];
    v[84]=v[487]*v[70]+v[58]*v[77]+v[499]*v[79];
    v[85]=v[487]*v[65]+v[53]*v[77]+v[494]*v[79];
    v[86]=v[487]*v[71]+v[59]*v[77]+v[500]*v[79];
    v[87]=v[487]*v[66]+v[54]*v[77]+v[495]*v[79];
    v[88]=v[487]*v[72]+v[60]*v[77]+v[501]*v[79];
    v[89]=v[487]*v[67]+v[55]*v[77]+v[496]*v[79];
    v[90]=v[487]*v[73]+v[61]*v[77]+v[488]*v[79];
    
    ed->ht[0]=v[76];
    ed->ht[1]=v[81];
    ed->ht[2]=v[83];
    ed->ht[3]=v[85];
    ed->ht[4]=v[87];
    ed->ht[5]=v[89];
    ed->ht[6]=v[80];
    ed->ht[7]=v[82];
    ed->ht[8]=v[84];
    ed->ht[9]=v[86];
    ed->ht[10]=v[88];
    ed->ht[11]=v[90];
    ed->ht[12]=v[489]*(v[490]-v[49]*v[50]+v[491]*v[62]);
    ed->ht[13]=v[489]*(v[492]-v[49]*v[51]+v[491]*v[63]);
    ed->ht[14]=v[489]*(v[493]-v[49]*v[52]+v[491]*v[64]);
    ed->ht[15]=v[489]*(v[494]-v[49]*v[53]+v[491]*v[65]);
    ed->ht[16]=v[489]*(v[495]-v[49]*v[54]+v[491]*v[66]);
    ed->ht[17]=v[489]*(v[496]-v[49]*v[55]+v[491]*v[67]);
    ed->ht[18]=v[489]*(v[497]-v[49]*v[56]+v[491]*v[68]);
    ed->ht[19]=v[489]*(v[498]-v[49]*v[57]+v[491]*v[69]);
    ed->ht[20]=v[489]*(v[499]-v[49]*v[58]+v[491]*v[70]);
    ed->ht[21]=v[489]*(v[500]-v[49]*v[59]+v[491]*v[71]);
    ed->ht[22]=v[489]*(v[501]-v[49]*v[60]+v[491]*v[72]);
    ed->ht[23]=v[489]*(v[488]-v[49]*v[61]+v[491]*v[73]);
    //es->id.NoIntPoints wird in SMTSetElSpec gesetzt
    for(i107=1;i107<=es->id.NoIntPoints;i107++){
        i109=4*(-1+i107);
        v[108]=es->IntPoints[i109];
        v[137]=4e0*v[108];
        v[132]=-1e0+v[137];
        v[110]=es->IntPoints[1+i109];
        v[136]=4e0*v[110];
        v[133]=-1e0+v[136];
        v[119]=-1e0+v[108]+v[110];
        v[5192]=0e0;
        v[5193]=0e0;
        v[5194]=0e0;
        v[5195]=0e0;
        v[5196]=0e0;
        v[5197]=0e0;
        v[5198]=0e0;
        v[5199]=0e0;
        v[5200]=0e0;
        v[5201]=0e0;
        v[5202]=0e0;
        v[5203]=0e0;
        v[5204]=-v[108];
        v[5205]=-v[110];
        v[5206]=v[119];
        v[139]=-4e0*v[119];
        v[140]=-v[137]+v[139];
        v[138]=-v[136]+v[139];
        v[135]=-1e0+2e0*v[108]+2e0*v[110]+2e0*v[119];
        v[5177]=0e0;
        v[5178]=v[132];
        v[5179]=0e0;
        v[5180]=0e0;
        v[5181]=0e0;
        v[5182]=v[135];
        v[5183]=0e0;
        v[5184]=0e0;
        v[5185]=0e0;
        v[5186]=0e0;
        v[5187]=0e0;
        v[5188]=v[140];
        v[5189]=0e0;
        v[5190]=0e0;
        v[5191]=0e0;
        v[5147]=0e0;
        v[5148]=0e0;
        v[5149]=0e0;
        v[5150]=v[133];
        v[5151]=0e0;
        v[5152]=v[135];
        v[5153]=0e0;
        v[5154]=0e0;
        v[5155]=0e0;
        v[5156]=v[138];
        v[5157]=0e0;
        v[5158]=0e0;
        v[5159]=0e0;
        v[5160]=0e0;
        v[5161]=0e0;
        v[5117]=v[132];
        v[5118]=0e0;
        v[5119]=0e0;
        v[5120]=0e0;
        v[5121]=v[135];
        v[5122]=0e0;
        v[5123]=0e0;
        v[5124]=0e0;
        v[5125]=0e0;
        v[5126]=0e0;
        v[5127]=v[140];
        v[5128]=0e0;
        v[5129]=0e0;
        v[5130]=0e0;
        v[5131]=0e0;
        v[5087]=0e0;
        v[5088]=0e0;
        v[5089]=v[133];
        v[5090]=0e0;
        v[5091]=v[135];
        v[5092]=0e0;
        v[5093]=0e0;
        v[5094]=0e0;
        v[5095]=v[138];
        v[5096]=0e0;
        v[5097]=0e0;
        v[5098]=0e0;
        v[5099]=0e0;
        v[5100]=0e0;
        v[5101]=0e0;
        v[164]=v[135]*v[84];
        v[161]=v[135]*v[83];
        v[153]=v[135]*v[24];
        v[154]=v[153]+v[133]*v[22]+v[138]*v[28]+v[137]*(v[26]-v[30]);
        v[152]=v[153]+v[132]*v[20]+v[136]*(v[26]-v[28])+v[140]*v[30];
        v[150]=v[135]*v[23];
        v[151]=v[150]+v[133]*v[21]+v[138]*v[27]+v[137]*(v[25]-v[29]);
        v[149]=v[150]+v[132]*v[19]+v[136]*(v[25]-v[27])+v[140]*v[29];
        v[145]=v[135]*v[6];
        v[146]=v[10]*v[138]+v[145]+v[133]*v[4]+v[137]*(-v[12]+v[8]);
        v[144]=v[12]*v[140]+v[145]+v[132]*v[2]+v[136]*(-v[10]+v[8]);
        v[142]=v[135]*v[5];
        v[143]=v[142]+v[133]*v[3]+v[137]*(-v[11]+v[7])+v[138]*v[9];
        v[141]=v[1]*v[132]+v[11]*v[140]+v[142]+v[136]*(v[7]-v[9]);
        v[155]=-(v[143]*v[144])+v[141]*v[146];
        v[525]=-(v[143]/v[155]);
        v[529]=v[140]*v[525];
        v[524]=v[141]/v[155];
        v[527]=v[138]*v[524];
        v[523]=v[146]/v[155];
        v[528]=v[140]*v[523];
        v[522]=-(v[144]/v[155]);
        v[526]=v[138]*v[522];
        v[519]=1e0/Power(v[155],2);
        v[503]=v[137]/v[155];
        v[502]=v[136]/v[155];
        v[267]=v[143]*v[502];
        v[266]=v[141]*v[503];
        v[263]=v[146]*v[502];
        v[262]=v[144]*v[503];
        v[148]=es->IntPoints[3+i109]*fabs(v[155]);
        v[176]=(v[146]*v[149]-v[144]*v[151])/v[155];
        v[179]=(-(v[143]*v[152])+v[141]*v[154])/v[155];
        v[195]=v[176]+v[179];
        v[505]=v[192]*v[195]-v[108]*v[31]-v[110]*v[32]+v[119]*v[33];
        v[185]=1e0-v[183]*(1e0+v[195]);
        v[504]=v[186]/(v[155]*v[185]);
        v[510]=v[185]*v[504];
        v[198]=-(((v[143]*v[149]-v[141]*v[151]-v[146]*v[152]+v[144]*v[154])*v[194])/v[155]);
        v[204]=2e0*v[176]*v[194]+v[505];
        v[205]=2e0*v[179]*v[194]+v[505];
        v[208]=(-(v[144]*(v[161]+v[133]*v[81]+v[137]*v[85]+v[138]*v[87]-v[137]*v[89]))+v[146]*(v[161]
                                                                                               +v[132]*v[76]+v[136]*v[85]-v[136]*v[87]+v[140]*v[89])+v[141]*(v[164]+v[133]*v[82]+v[137]*v[86]
                                                                                                                                                             +v[138]*v[88]-v[137]*v[90])-v[143]*(v[164]+v[132]*v[80]+v[136]*v[86]-v[136]*v[88]+v[140]*v[90]))
        /v[155];
        v[210]=(-(v[146]*v[171])+v[144]*v[172])*v[185]*v[504];
        v[211]=(v[143]*v[171]-v[141]*v[172])*v[185]*v[504];
        v[217]=(v[146]*v[198]-v[143]*v[205])/v[155];
        v[225]=v[136]*v[217];
        v[218]=(-(v[144]*v[198])+v[141]*v[205])/v[155];
        v[226]=v[137]*v[218];
        v[219]=(-(v[143]*v[198])+v[146]*v[204])/v[155];
        v[223]=v[136]*v[219];
        v[220]=(v[141]*v[198]-v[144]*v[204])/v[155];
        v[224]=v[137]*v[220];
        v[221]=(v[146]*v[210]-v[143]*v[211])/v[155];
        v[222]=(-(v[144]*v[210])+v[141]*v[211])/v[155];
        v[5023]=v[132]*v[219];
        v[5024]=v[132]*v[217];
        v[5025]=v[133]*v[220];
        v[5026]=v[133]*v[218];
        v[5027]=v[135]*(v[219]+v[220]);
        v[5028]=v[135]*(v[217]+v[218]);
        v[5029]=v[223]+v[224];
        v[5030]=v[225]+v[226];
        v[5031]=v[138]*v[220]-v[223];
        v[5032]=v[138]*v[218]-v[225];
        v[5033]=v[140]*v[219]-v[224];
        v[5034]=v[140]*v[217]-v[226];
        v[5035]=-(v[108]*v[208])+v[221];
        v[5036]=-(v[110]*v[208])+v[222];
        v[5037]=v[119]*v[208]-v[221]-v[222];
//        for(i215=1;i215<=15;i215++){
//        
//            for(i230=1;i230<=15;i230++){
//                std::cout << " pre s[i215-1][i230-1]:"<< s[i215-1][i230-1] << std::endl;
//            }
//            
//        }
        
        for(i215=1;i215<=15;i215++){
            v[233]=v[5041+i215];
            v[234]=v[5056+i215];
            v[509]=v[5191+i215]*v[79];
            v[261]=v[135]*v[509];
            v[260]=v[509]/v[155];
            v[246]=((v[141]*v[233]-v[143]*v[234])*v[510])/v[155];
            v[247]=((-(v[144]*v[233])+v[146]*v[234])*v[510])/v[155];
            v[238]=v[137]*v[5071+i215]+v[5086+i215];
            v[239]=v[136]*v[5101+i215]+v[5116+i215];
            v[240]=(-(v[144]*v[238])+v[146]*v[239])/v[155];
            v[241]=v[137]*v[5131+i215]+v[5146+i215];
            v[242]=v[136]*v[5161+i215]+v[5176+i215];
            v[243]=(v[141]*v[241]-v[143]*v[242])/v[155];
            
            v[256]=v[194]*(v[141]*v[238]-v[143]*v[239]-v[144]*v[241]+v[146]*v[242])*v[519];
            v[245]=v[143]*v[246]-v[146]*v[247];
            v[248]=-(v[141]*v[246])+v[144]*v[247];
            v[249]=-v[240]-v[243];
            v[252]=-(v[192]*v[249]);
            v[520]=(2e0*v[194]*v[243]+v[252])/v[155];
            v[521]=(2e0*v[194]*v[240]+v[252])/v[155];
            v[255]=v[146]*v[256]-v[143]*v[520];
            v[268]=v[136]*v[255];
            v[257]=-(v[144]*v[256])+v[141]*v[520];
            v[269]=v[137]*v[257];
            v[258]=-(v[143]*v[256])+v[146]*v[521];
            v[264]=v[136]*v[258];
            v[259]=v[141]*v[256]-v[144]*v[521];
            v[265]=v[137]*v[259];
            
            v[5207]=v[132]*(v[258]+v[146]*v[260]);
            v[5208]=v[132]*(v[255]-v[143]*v[260]);

            v[5209]=v[133]*(v[259]-v[144]*v[260]);
            v[5210]=v[133]*(v[257]+v[141]*v[260]);
            v[5211]=v[135]*(v[258]+v[259])+v[261]*(v[522]+v[523]);
            v[5212]=v[135]*(v[255]+v[257])+v[261]*(v[524]+v[525]);
            v[5213]=v[264]+v[265]+(-v[262]+v[263])*v[509];
            v[5214]=v[268]+v[269]+(v[266]-v[267])*v[509];
            v[5215]=v[138]*v[259]-v[264]+v[509]*(-v[263]+v[526]);
            v[5216]=v[138]*v[257]-v[268]+v[509]*(v[267]+v[527]);
            v[5217]=v[140]*v[258]-v[265]+v[509]*(v[262]+v[528]);
            v[5218]=v[140]*v[255]-v[269]+v[509]*(-v[266]+v[529]);
            v[5219]=v[245]+v[108]*v[249];
            v[5220]=v[248]+v[110]*v[249];
            v[5221]=-v[245]-v[248]-v[119]*v[249];
            p[i215-1]+=v[148]*v[5022+i215];
            for(i230=1;i230<=15;i230++){
                s[i215-1][i230-1]+=v[148]*v[5206+i230];

            };/* end for */
        };/* end for */
    };/* end for */
};

template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::SKR_Biot_StVK(double* v,ElementSpec *es,ElementData *ed,NodeSpec **ns
                                    ,NodeData **nd,double *rdata,int *idata,double *p,double **s)
{
    int i114,i116,i233,i250,b297;
    FILE *SMSFile;
    v[5132]=0e0;
    v[5133]=0e0;
    v[5134]=0e0;
    v[5135]=0e0;
    v[5136]=0e0;
    v[5137]=0e0;
    v[5138]=0e0;
    v[5139]=1e0;
    v[5140]=0e0;
    v[5141]=-1e0;
    v[5142]=0e0;
    v[5143]=0e0;
    v[5144]=0e0;
    v[5145]=0e0;
    v[5146]=0e0;
    v[5102]=0e0;
    v[5103]=0e0;
    v[5104]=0e0;
    v[5105]=0e0;
    v[5106]=0e0;
    v[5107]=0e0;
    v[5108]=0e0;
    v[5109]=1e0;
    v[5110]=0e0;
    v[5111]=0e0;
    v[5112]=0e0;
    v[5113]=-1e0;
    v[5114]=0e0;
    v[5115]=0e0;
    v[5116]=0e0;
    v[5072]=0e0;
    v[5073]=0e0;
    v[5074]=0e0;
    v[5075]=0e0;
    v[5076]=0e0;
    v[5077]=0e0;
    v[5078]=1e0;
    v[5079]=0e0;
    v[5080]=-1e0;
    v[5081]=0e0;
    v[5082]=0e0;
    v[5083]=0e0;
    v[5084]=0e0;
    v[5085]=0e0;
    v[5086]=0e0;
    v[5042]=0e0;
    v[5043]=0e0;
    v[5044]=0e0;
    v[5045]=0e0;
    v[5046]=0e0;
    v[5047]=0e0;
    v[5048]=1e0;
    v[5049]=0e0;
    v[5050]=0e0;
    v[5051]=0e0;
    v[5052]=-1e0;
    v[5053]=0e0;
    v[5054]=0e0;
    v[5055]=0e0;
    v[5056]=0e0;
    v[5177]=0e0;
    v[5178]=0e0;
    v[5179]=0e0;
    v[5180]=0e0;
    v[5181]=0e0;
    v[5182]=0e0;
    v[5183]=0e0;
    v[5184]=0e0;
    v[5185]=0e0;
    v[5186]=0e0;
    v[5187]=0e0;
    v[5188]=0e0;
    v[5189]=1e0;
    v[5190]=0e0;
    v[5191]=-1e0;
    v[5162]=0e0;
    v[5163]=0e0;
    v[5164]=0e0;
    v[5165]=0e0;
    v[5166]=0e0;
    v[5167]=0e0;
    v[5168]=0e0;
    v[5169]=0e0;
    v[5170]=0e0;
    v[5171]=0e0;
    v[5172]=0e0;
    v[5173]=0e0;
    v[5174]=0e0;
    v[5175]=1e0;
    v[5176]=-1e0;
    v[1]=nd[0]->X[0];
    v[2]=nd[0]->X[1];
    v[3]=nd[1]->X[0];
    v[4]=nd[1]->X[1];
    v[5]=nd[2]->X[0];
    v[6]=nd[2]->X[1];
    v[7]=nd[3]->X[0];
    v[8]=nd[3]->X[1];
    v[9]=nd[4]->X[0];
    v[10]=nd[4]->X[1];
    v[11]=nd[5]->X[0];
    v[12]=nd[5]->X[1];
    v[19]=nd[0]->at[0];
    v[20]=nd[0]->at[1];
    v[21]=nd[1]->at[0];
    v[22]=nd[1]->at[1];
    v[23]=nd[2]->at[0];
    v[24]=nd[2]->at[1];
    v[25]=nd[3]->at[0];
    v[26]=nd[3]->at[1];
    v[27]=nd[4]->at[0];
    v[28]=nd[4]->at[1];
    v[29]=nd[5]->at[0];
    v[30]=nd[5]->at[1];
    v[31]=nd[6]->at[0];
    v[32]=nd[7]->at[0];
    v[33]=nd[8]->at[0];
    v[179]=v[32]-v[33];
    v[178]=v[31]-v[33];
    v[34]=nd[0]->ap[0];
    v[574]=v[19]-v[34];
    v[35]=nd[0]->ap[1];
    v[576]=v[20]-v[35];
    v[36]=nd[1]->ap[0];
    v[578]=v[21]-v[36];
    v[37]=nd[1]->ap[1];
    v[579]=v[22]-v[37];
    v[38]=nd[2]->ap[0];
    v[580]=v[23]-v[38];
    v[39]=nd[2]->ap[1];
    v[581]=v[24]-v[39];
    v[40]=nd[3]->ap[0];
    v[582]=v[25]-v[40];
    v[41]=nd[3]->ap[1];
    v[583]=v[26]-v[41];
    v[42]=nd[4]->ap[0];
    v[584]=v[27]-v[42];
    v[43]=nd[4]->ap[1];
    v[585]=v[28]-v[43];
    v[44]=nd[5]->ap[0];
    v[586]=v[29]-v[44];
    v[45]=nd[5]->ap[1];
    v[587]=v[30]-v[45];
    v[49]=es->Data[0];
    v[50]=es->Data[1];
    v[571]=v[49]/(1e0+v[50]);
    v[51]=(v[50]*v[571])/(1e0-2e0*v[50]);
    v[53]=v[571]/2e0;
    v[54]=es->Data[2];
    v[55]=es->Data[3];
    v[56]=rdata[RD_TimeIncrement];
    v[57]=ed->hp[0];
    v[58]=ed->hp[1];
    v[59]=ed->hp[2];
    v[60]=ed->hp[3];
    v[61]=ed->hp[4];
    v[62]=ed->hp[5];
    v[63]=ed->hp[6];
    v[64]=ed->hp[7];
    v[65]=ed->hp[8];
    v[66]=ed->hp[9];
    v[67]=ed->hp[10];
    v[68]=ed->hp[11];
    v[69]=ed->hp[12];
    v[70]=ed->hp[13];
    v[71]=ed->hp[14];
    v[72]=ed->hp[15];
    v[73]=ed->hp[16];
    v[74]=ed->hp[17];
    v[75]=ed->hp[18];
    v[76]=ed->hp[19];
    v[77]=ed->hp[20];
    v[78]=ed->hp[21];
    v[79]=ed->hp[22];
    v[80]=ed->hp[23];
    v[82]=es->Data[5];
    v[575]=1e0/(Power(v[56],2)*v[82]);
    v[577]=-((v[56]*v[56])*(0.5e0-v[82]));
    v[86]=es->Data[4]/(v[56]*v[82]);
    v[572]=-(v[56]*v[86]);
    v[573]=v[56]*(1e0+v[572]/2e0);
    v[84]=1e0+v[572];
    v[83]=v[573]*v[69]+v[57]*v[84]+v[574]*v[86];
    v[87]=v[573]*v[75]+v[63]*v[84]+v[576]*v[86];
    v[88]=v[573]*v[70]+v[58]*v[84]+v[578]*v[86];
    v[89]=v[573]*v[76]+v[64]*v[84]+v[579]*v[86];
    v[90]=v[573]*v[71]+v[59]*v[84]+v[580]*v[86];
    v[91]=v[573]*v[77]+v[65]*v[84]+v[581]*v[86];
    v[92]=v[573]*v[72]+v[60]*v[84]+v[582]*v[86];
    v[93]=v[573]*v[78]+v[66]*v[84]+v[583]*v[86];
    v[94]=v[573]*v[73]+v[61]*v[84]+v[584]*v[86];
    v[95]=v[573]*v[79]+v[67]*v[84]+v[585]*v[86];
    v[96]=v[573]*v[74]+v[62]*v[84]+v[586]*v[86];
    v[97]=v[573]*v[80]+v[68]*v[84]+v[587]*v[86];
    v[98]=v[575]*(-(v[56]*v[57])+v[574]+v[577]*v[69]);
    v[102]=v[575]*(v[576]-v[56]*v[63]+v[577]*v[75]);
    v[103]=v[575]*(v[578]-v[56]*v[58]+v[577]*v[70]);
    v[104]=v[575]*(v[579]-v[56]*v[64]+v[577]*v[76]);
    v[105]=v[575]*(v[580]-v[56]*v[59]+v[577]*v[71]);
    v[106]=v[575]*(v[581]-v[56]*v[65]+v[577]*v[77]);
    v[107]=v[575]*(v[582]-v[56]*v[60]+v[577]*v[72]);
    v[108]=v[575]*(v[583]-v[56]*v[66]+v[577]*v[78]);
    v[109]=v[575]*(v[584]-v[56]*v[61]+v[577]*v[73]);
    v[110]=v[575]*(v[585]-v[56]*v[67]+v[577]*v[79]);
    v[111]=v[575]*(v[586]-v[56]*v[62]+v[577]*v[74]);
    v[112]=v[575]*(v[587]-v[56]*v[68]+v[577]*v[80]);
    ed->ht[0]=v[83];
    ed->ht[1]=v[88];
    ed->ht[2]=v[90];
    ed->ht[3]=v[92];
    ed->ht[4]=v[94];
    ed->ht[5]=v[96];
    ed->ht[6]=v[87];
    ed->ht[7]=v[89];
    ed->ht[8]=v[91];
    ed->ht[9]=v[93];
    ed->ht[10]=v[95];
    ed->ht[11]=v[97];
    ed->ht[12]=v[98];
    ed->ht[13]=v[103];
    ed->ht[14]=v[105];
    ed->ht[15]=v[107];
    ed->ht[16]=v[109];
    ed->ht[17]=v[111];
    ed->ht[18]=v[102];
    ed->ht[19]=v[104];
    ed->ht[20]=v[106];
    ed->ht[21]=v[108];
    ed->ht[22]=v[110];
    ed->ht[23]=v[112];
    for(i114=1;i114<=es->id.NoIntPoints;i114++){
     i116=4*(-1+i114);
     v[115]=es->IntPoints[i116];
     v[144]=4e0*v[115];
     v[139]=-1e0+v[144];
     v[117]=es->IntPoints[1+i116];
     v[143]=4e0*v[117];
     v[140]=-1e0+v[143];
     v[126]=-1e0+v[115]+v[117];
     v[5192]=0e0;
     v[5193]=0e0;
     v[5194]=0e0;
     v[5195]=0e0;
     v[5196]=0e0;
     v[5197]=0e0;
     v[5198]=0e0;
     v[5199]=0e0;
     v[5200]=0e0;
     v[5201]=0e0;
     v[5202]=0e0;
     v[5203]=0e0;
     v[5204]=-v[115];
     v[5205]=-v[117];
     v[5206]=v[126];
     v[146]=-4e0*v[126];
     v[147]=-v[144]+v[146];
     v[145]=-v[143]+v[146];
     v[142]=-1e0+2e0*v[115]+2e0*v[117]+2e0*v[126];
     v[5147]=0e0;
     v[5148]=v[139];
     v[5149]=0e0;
     v[5150]=0e0;
     v[5151]=0e0;
     v[5152]=v[142];
     v[5153]=0e0;
     v[5154]=0e0;
     v[5155]=0e0;
     v[5156]=0e0;
     v[5157]=0e0;
     v[5158]=v[147];
     v[5159]=0e0;
     v[5160]=0e0;
     v[5161]=0e0;
     v[5117]=0e0;
     v[5118]=0e0;
     v[5119]=0e0;
     v[5120]=v[140];
     v[5121]=0e0;
     v[5122]=v[142];
     v[5123]=0e0;
     v[5124]=0e0;
     v[5125]=0e0;
     v[5126]=v[145];
     v[5127]=0e0;
     v[5128]=0e0;
     v[5129]=0e0;
     v[5130]=0e0;
     v[5131]=0e0;
     v[5087]=v[139];
     v[5088]=0e0;
     v[5089]=0e0;
     v[5090]=0e0;
     v[5091]=v[142];
     v[5092]=0e0;
     v[5093]=0e0;
     v[5094]=0e0;
     v[5095]=0e0;
     v[5096]=0e0;
     v[5097]=v[147];
     v[5098]=0e0;
     v[5099]=0e0;
     v[5100]=0e0;
     v[5101]=0e0;
     v[5057]=0e0;
     v[5058]=0e0;
     v[5059]=v[140];
     v[5060]=0e0;
     v[5061]=v[142];
     v[5062]=0e0;
     v[5063]=0e0;
     v[5064]=0e0;
     v[5065]=v[145];
     v[5066]=0e0;
     v[5067]=0e0;
     v[5068]=0e0;
     v[5069]=0e0;
     v[5070]=0e0;
     v[5071]=0e0;
     v[171]=v[142]*v[91];
     v[168]=v[142]*v[90];
     v[160]=v[142]*v[24];
     v[161]=v[160]+v[140]*v[22]+v[145]*v[28]+v[144]*(v[26]-v[30]);
     v[159]=v[160]+v[139]*v[20]+v[143]*(v[26]-v[28])+v[147]*v[30];
     v[157]=v[142]*v[23];
     v[158]=v[157]+v[140]*v[21]+v[145]*v[27]+v[144]*(v[25]-v[29]);
     v[156]=v[157]+v[139]*v[19]+v[143]*(v[25]-v[27])+v[147]*v[29];
     v[152]=v[142]*v[6];
     v[153]=v[10]*v[145]+v[152]+v[140]*v[4]+v[144]*(-v[12]+v[8]);
     v[151]=v[12]*v[147]+v[152]+v[139]*v[2]+v[143]*(-v[10]+v[8]);
     v[149]=v[142]*v[5];
     v[150]=v[149]+v[140]*v[3]+v[144]*(-v[11]+v[7])+v[145]*v[9];
     v[148]=v[1]*v[139]+v[11]*v[147]+v[149]+v[143]*(v[7]-v[9]);
     v[162]=-(v[150]*v[151])+v[148]*v[153];
     v[610]=-(v[150]/v[162]);
     v[614]=v[147]*v[610];
     v[609]=v[148]/v[162];
     v[612]=v[145]*v[609];
     v[608]=v[153]/v[162];
     v[613]=v[147]*v[608];
     v[607]=-(v[151]/v[162]);
     v[611]=v[145]*v[607];
     v[589]=v[144]/v[162];
     v[588]=v[143]/v[162];
     v[292]=v[150]*v[588];
     v[291]=v[148]*v[589];
     v[288]=v[153]*v[588];
     v[287]=v[151]*v[589];
     v[155]=es->IntPoints[3+i116]*fabs(v[162]);
     v[163]=(v[153]*v[156]-v[151]*v[158])/v[162];
     v[164]=(-(v[150]*v[156])+v[148]*v[158])/v[162];
     v[165]=(v[153]*v[159]-v[151]*v[161])/v[162];
     v[166]=(-(v[150]*v[159])+v[148]*v[161])/v[162];
     v[190]=1e0-v[163];
     v[191]=1e0-v[166];
     v[195]=0.5e0*(-1e0+(v[165]*v[165])+(v[190]*v[190]));
     v[200]=0.5e0*(-1e0+(v[164]*v[164])+(v[191]*v[191]));
     v[591]=-(v[115]*v[31])-v[117]*v[32]+v[126]*v[33]+(v[195]+v[200])*v[51];
     v[206]=1e0-(1e0+v[163]+v[166])*v[54];
     v[590]=v[55]/(v[162]*v[206]);
     v[604]=v[206]*v[590];
     v[222]=2e0*v[195]*v[53]+v[591];
     v[223]=2e0*v[200]*v[53]+v[591];
     v[226]=(-(v[151]*(v[168]+v[140]*v[88]+v[144]*v[92]+v[145]*v[94]-v[144]*v[96]))+v[153]*(v[168]
      +v[139]*v[83]+v[143]*v[92]-v[143]*v[94]+v[147]*v[96])+v[148]*(v[171]+v[140]*v[89]+v[144]*v[93]
      +v[145]*v[95]-v[144]*v[97])-v[150]*(v[171]+v[139]*v[87]+v[143]*v[93]-v[143]*v[95]+v[147]*v[97]))
      /v[162];
     v[228]=(-(v[153]*v[178])+v[151]*v[179])*v[206]*v[590];
     v[229]=(v[150]*v[178]-v[148]*v[179])*v[206]*v[590];
     v[236]=-1e0*(v[164]*v[190]+v[165]*v[191])*v[53];
     v[237]=(v[153]*v[228]-v[150]*v[229])/v[162];
     v[238]=(-(v[151]*v[228])+v[148]*v[229])/v[162];
     v[239]=(-(v[150]*v[223])+v[153]*v[236])/v[162];
     v[245]=v[143]*v[239];
     v[240]=(v[148]*v[223]-v[151]*v[236])/v[162];
     v[246]=v[144]*v[240];
     v[241]=(v[153]*v[222]-v[150]*v[236])/v[162];
     v[243]=v[143]*v[241];
     v[242]=(-(v[151]*v[222])+v[148]*v[236])/v[162];
     v[244]=v[144]*v[242];
     v[5023]=v[139]*v[241];
     v[5024]=v[139]*v[239];
     v[5025]=v[140]*v[242];
     v[5026]=v[140]*v[240];
     v[5027]=v[142]*(v[241]+v[242]);
     v[5028]=v[142]*(v[239]+v[240]);
     v[5029]=v[243]+v[244];
     v[5030]=v[245]+v[246];
     v[5031]=v[145]*v[242]-v[243];
     v[5032]=v[145]*v[240]-v[245];
     v[5033]=v[147]*v[241]-v[244];
     v[5034]=v[147]*v[239]-v[246];
     v[5035]=-(v[115]*v[226])+v[237];
     v[5036]=-(v[117]*v[226])+v[238];
     v[5037]=v[126]*v[226]-v[237]-v[238];
     for(i233=1;i233<=15;i233++){
      v[253]=v[5161+i233];
      v[254]=v[5176+i233];
      v[595]=v[5191+i233]*v[86];
      v[286]=v[142]*v[595];
      v[285]=v[595]/v[162];
      v[256]=v[144]*v[5041+i233]+v[5056+i233];
      v[257]=v[143]*v[5071+i233]+v[5086+i233];
      v[258]=(-(v[151]*v[256])+v[153]*v[257])/v[162];
      v[259]=v[144]*v[5101+i233]+v[5116+i233];
      v[260]=v[143]*v[5131+i233]+v[5146+i233];
      v[261]=(v[148]*v[259]-v[150]*v[260])/v[162];
      v[269]=((v[148]*v[253]-v[150]*v[254])*v[604])/v[162];
      v[270]=((-(v[151]*v[253])+v[153]*v[254])*v[604])/v[162];
      v[264]=(2e0*(v[148]*v[256]-v[150]*v[257]-v[151]*v[259]+v[153]*v[260])*v[53])/v[162];
      v[606]=0.5e0*v[264];
      v[605]=0.5e0*v[264];
      v[268]=v[150]*v[269]-v[153]*v[270];
      v[271]=-(v[148]*v[269])+v[151]*v[270];
      v[272]=v[258]+v[261];
      v[274]=v[272]*v[51];
      v[273]=v[274]+2e0*v[261]*v[53];
      v[275]=v[274]+2e0*v[258]*v[53];
      v[276]=1e0*v[165]*v[275]-v[191]*v[605];
      v[277]=1e0*v[164]*v[273]-v[190]*v[606];
      v[278]=-1e0*v[191]*v[273]+v[165]*v[605];
      v[280]=-1e0*v[190]*v[275]+v[164]*v[606];
      v[281]=(v[153]*v[276]-v[150]*v[278])/v[162];
      v[293]=v[143]*v[281];
      v[282]=(-(v[151]*v[276])+v[148]*v[278])/v[162];
      v[294]=v[144]*v[282];
      v[283]=(-(v[150]*v[277])+v[153]*v[280])/v[162];
      v[289]=v[143]*v[283];
      v[284]=(v[148]*v[277]-v[151]*v[280])/v[162];
      v[290]=v[144]*v[284];
      v[5207]=v[139]*(v[283]+v[153]*v[285]);
      v[5208]=v[139]*(v[281]-v[150]*v[285]);
      v[5209]=v[140]*(v[284]-v[151]*v[285]);
      v[5210]=v[140]*(v[282]+v[148]*v[285]);
      v[5211]=v[142]*(v[283]+v[284])+v[286]*(v[607]+v[608]);
      v[5212]=v[142]*(v[281]+v[282])+v[286]*(v[609]+v[610]);
      v[5213]=v[289]+v[290]+(-v[287]+v[288])*v[595];
      v[5214]=v[293]+v[294]+(v[291]-v[292])*v[595];
      v[5215]=v[145]*v[284]-v[289]+v[595]*(-v[288]+v[611]);
      v[5216]=v[145]*v[282]-v[293]+v[595]*(v[292]+v[612]);
      v[5217]=v[147]*v[283]-v[290]+v[595]*(v[287]+v[613]);
      v[5218]=v[147]*v[281]-v[294]+v[595]*(-v[291]+v[614]);
      v[5219]=v[268]-v[115]*v[272];
      v[5220]=v[271]-v[117]*v[272];
      v[5221]=-v[268]-v[271]+v[126]*v[272];
      p[i233-1]+=v[155]*v[5022+i233];
      for(i250=1;i250<=15;i250++){
       s[i233-1][i250-1]+=v[155]*v[5206+i250];
      };/* end for */
     };/* end for */
    };/* end for */
    if(idata[ID_CurrentElement]==1e0){
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","ElementNumber: ",(double)idata[ID_CurrentElement]
      ,"   Time: ",(double)rdata[RD_Time],"   TimeSteps: ",(double)idata[ID_Step],"   IterationSteps: ",
      (double)idata[ID_Iteration]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g ","Emodul: ",(double)v[49],"   Poissonsratio: ",(double)v[50]
      ,"   InitialSolidVolumeFraction: ",(double)v[54]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g ","TimeIncrement :",(double)rdata[RD_TimeIncrement]
      ,"   LoadIncrement: ",(double)rdata[RD_MultiplierIncrement]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[0],"  CoorX: ",(double
      )nd[0]->X[0],"  CoorY: ",(double)nd[0]->X[1],"  NumberDOF: ",(double)nd[0]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[0]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[0]->DOF[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[1],"  CoorX: ",(double
      )nd[1]->X[0],"  CoorY: ",(double)nd[1]->X[1],"  NumberDOF: ",(double)nd[1]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[1]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[1]->DOF[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[2],"  CoorX: ",(double
      )nd[2]->X[0],"  CoorY: ",(double)nd[2]->X[1],"  NumberDOF: ",(double)nd[2]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[2]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[2]->DOF[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[3],"  CoorX: ",(double
      )nd[3]->X[0],"  CoorY: ",(double)nd[3]->X[1],"  NumberDOF: ",(double)nd[3]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[3]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[3]->DOF[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[4],"  CoorX: ",(double
      )nd[4]->X[0],"  CoorY: ",(double)nd[4]->X[1],"  NumberDOF: ",(double)nd[4]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[4]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[4]->DOF[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[5],"  CoorX: ",(double
      )nd[5]->X[0],"  CoorY: ",(double)nd[5]->X[1],"  NumberDOF: ",(double)nd[5]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[5]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[5]->DOF[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[6],"  CoorX: ",(double
      )nd[6]->X[0],"  CoorY: ",(double)nd[6]->X[1],"  NumberDOF: ",(double)nd[6]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[6]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[7],"  CoorX: ",(double
      )nd[7]->X[0],"  CoorY: ",(double)nd[7]->X[1],"  NumberDOF: ",(double)nd[7]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[7]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[8],"  CoorX: ",(double
      )nd[8]->X[0],"  CoorY: ",(double)nd[8]->X[1],"  NumberDOF: ",(double)nd[8]->id.NoDOF);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g ","DOFposition: ",(double)nd[8]->DOF[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g "
      ,"hp or h1 field: ",(double)v[57],(double)v[58],(double)v[59],(double)v[60],(double)v[61],(double
      )v[62],(double)v[63],(double)v[64],(double)v[65],(double)v[66],(double)v[67],(double)v[68],(double
      )v[69],(double)v[70],(double)v[71],(double)v[72],(double)v[73],(double)v[74],(double)v[75],(double
      )v[76],(double)v[77],(double)v[78],(double)v[79],(double)v[80]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g "
      ,"ht or h2 field: ",(double)v[83],(double)v[88],(double)v[90],(double)v[92],(double)v[94],(double
      )v[96],(double)v[87],(double)v[89],(double)v[91],(double)v[93],(double)v[95],(double)v[97],(double
      )v[98],(double)v[103],(double)v[105],(double)v[107],(double)v[109],(double)v[111],(double)v[102],
      (double)v[104],(double)v[106],(double)v[108],(double)v[110],(double)v[112]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ","ap or up field: ",(double
      )v[34],(double)v[35],(double)v[36],(double)v[37],(double)v[38],(double)v[39],(double)v[40],(double
      )v[41],(double)v[42],(double)v[43],(double)v[44],(double)v[45],(double)nd[6]->ap[0],(double)nd[7]
      ->ap[0],(double)nd[8]->ap[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g ","at or ul field: ",(double
      )v[19],(double)v[20],(double)v[21],(double)v[22],(double)v[23],(double)v[24],(double)v[25],(double
      )v[26],(double)v[27],(double)v[28],(double)v[29],(double)v[30],(double)v[31],(double)v[32],(double
      )v[33]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[0],"  Ux: ",(double)nd[0]
      ->at[0],"  Uy: ",(double)nd[0]->at[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[1],"  Ux: ",(double)nd[1]
      ->at[0],"  Uy: ",(double)nd[1]->at[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[2],"  Ux: ",(double)nd[2]
      ->at[0],"  Uy: ",(double)nd[2]->at[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[3],"  Ux: ",(double)nd[3]
      ->at[0],"  Uy: ",(double)nd[3]->at[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[4],"  Ux: ",(double)nd[4]
      ->at[0],"  Uy: ",(double)nd[4]->at[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g %s %g ","Nodes: ",(double)ed->Nodes[5],"  Ux: ",(double)nd[5]
      ->at[0],"  Uy: ",(double)nd[5]->at[1]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g ","Nodes: ",(double)ed->Nodes[6],"  p: ",(double)nd[6]->at[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g ","Nodes: ",(double)ed->Nodes[7],"  p: ",(double)nd[7]->at[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s %g %s %g ","Nodes: ",(double)ed->Nodes[8],"  p: ",(double)nd[8]->at[0]);
     fclose(SMSFile);};};
     ++idata[ID_NoMessages];
     if(1){SMSFile=fopen("myoutput.dat","a");if(SMSFile!=NULL){
     fprintf(SMSFile,"\n%s ","%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
     fclose(SMSFile);};};
    } else {
    };
};


// Added 07.09.2020 CH
template <class SC, class LO, class GO, class NO>
void FE<SC,LO,GO,NO>::SKR_Biot3D(double* v, ElementSpec *es, ElementData *ed, NodeSpec **ns, NodeData **nd, double *rdata, int *idata, double *p, double **s)
{
    int i254,i256,i438,i455;
    v[1514]=0e0;
    v[1515]=0e0;
    v[1516]=0e0;
    v[1517]=0e0;
    v[1518]=0e0;
    v[1519]=0e0;
    v[1520]=0e0;
    v[1521]=0e0;
    v[1522]=0e0;
    v[1523]=0e0;
    v[1524]=0e0;
    v[1525]=0e0;
    v[1526]=0e0;
    v[1527]=0e0;
    v[1528]=0e0;
    v[1529]=0e0;
    v[1530]=0e0;
    v[1531]=0e0;
    v[1532]=0e0;
    v[1533]=0e0;
    v[1534]=0e0;
    v[1535]=0e0;
    v[1536]=0e0;
    v[1537]=0e0;
    v[1538]=0e0;
    v[1539]=0e0;
    v[1540]=0e0;
    v[1541]=0e0;
    v[1542]=0e0;
    v[1543]=0e0;
    v[1544]=1e0;
    v[1545]=0e0;
    v[1546]=0e0;
    v[1547]=-1e0;
    v[1480]=0e0;
    v[1481]=0e0;
    v[1482]=0e0;
    v[1483]=0e0;
    v[1484]=0e0;
    v[1485]=0e0;
    v[1486]=0e0;
    v[1487]=0e0;
    v[1488]=0e0;
    v[1489]=0e0;
    v[1490]=0e0;
    v[1491]=0e0;
    v[1492]=0e0;
    v[1493]=0e0;
    v[1494]=0e0;
    v[1495]=0e0;
    v[1496]=0e0;
    v[1497]=0e0;
    v[1498]=0e0;
    v[1499]=0e0;
    v[1500]=0e0;
    v[1501]=0e0;
    v[1502]=0e0;
    v[1503]=0e0;
    v[1504]=0e0;
    v[1505]=0e0;
    v[1506]=0e0;
    v[1507]=0e0;
    v[1508]=0e0;
    v[1509]=0e0;
    v[1510]=0e0;
    v[1511]=1e0;
    v[1512]=0e0;
    v[1513]=-1e0;
    v[1446]=0e0;
    v[1447]=0e0;
    v[1448]=0e0;
    v[1449]=0e0;
    v[1450]=0e0;
    v[1451]=0e0;
    v[1452]=0e0;
    v[1453]=0e0;
    v[1454]=0e0;
    v[1455]=0e0;
    v[1456]=0e0;
    v[1457]=0e0;
    v[1458]=0e0;
    v[1459]=0e0;
    v[1460]=0e0;
    v[1461]=0e0;
    v[1462]=0e0;
    v[1463]=0e0;
    v[1464]=0e0;
    v[1465]=0e0;
    v[1466]=0e0;
    v[1467]=0e0;
    v[1468]=0e0;
    v[1469]=0e0;
    v[1470]=0e0;
    v[1471]=0e0;
    v[1472]=0e0;
    v[1473]=0e0;
    v[1474]=0e0;
    v[1475]=0e0;
    v[1476]=0e0;
    v[1477]=0e0;
    v[1478]=1e0;
    v[1479]=-1e0;
    v[1]=nd[0]->X[0];
    v[2]=nd[0]->X[1];
    v[3]=nd[0]->X[2];
    v[4]=nd[1]->X[0];
    v[5]=nd[1]->X[1];
    v[6]=nd[1]->X[2];
    v[7]=nd[2]->X[0];
    v[8]=nd[2]->X[1];
    v[9]=nd[2]->X[2];
    v[10]=nd[3]->X[0];
    v[11]=nd[3]->X[1];
    v[12]=nd[3]->X[2];
    v[13]=nd[4]->X[0];
    v[14]=nd[4]->X[1];
    v[15]=nd[4]->X[2];
    v[16]=nd[5]->X[0];
    v[17]=nd[5]->X[1];
    v[18]=nd[5]->X[2];
    v[19]=nd[6]->X[0];
    v[20]=nd[6]->X[1];
    v[21]=nd[6]->X[2];
    v[22]=nd[7]->X[0];
    v[23]=nd[7]->X[1];
    v[24]=nd[7]->X[2];
    v[25]=nd[8]->X[0];
    v[26]=nd[8]->X[1];
    v[27]=nd[8]->X[2];
    v[28]=nd[9]->X[0];
    v[29]=nd[9]->X[1];
    v[30]=nd[9]->X[2];
    v[43]=nd[0]->at[0];
    v[44]=nd[0]->at[1];
    v[45]=nd[0]->at[2];
    v[46]=nd[1]->at[0];
    v[47]=nd[1]->at[1];
    v[48]=nd[1]->at[2];
    v[49]=nd[2]->at[0];
    v[50]=nd[2]->at[1];
    v[51]=nd[2]->at[2];
    v[52]=nd[3]->at[0];
    v[53]=nd[3]->at[1];
    v[54]=nd[3]->at[2];
    v[55]=nd[4]->at[0];
    v[56]=nd[4]->at[1];
    v[57]=nd[4]->at[2];
    v[58]=nd[5]->at[0];
    v[59]=nd[5]->at[1];
    v[60]=nd[5]->at[2];
    v[61]=nd[6]->at[0];
    v[62]=nd[6]->at[1];
    v[63]=nd[6]->at[2];
    v[64]=nd[7]->at[0];
    v[65]=nd[7]->at[1];
    v[66]=nd[7]->at[2];
    v[67]=nd[8]->at[0];
    v[68]=nd[8]->at[1];
    v[69]=nd[8]->at[2];
    v[70]=nd[9]->at[0];
    v[71]=nd[9]->at[1];
    v[72]=nd[9]->at[2];
    v[73]=nd[10]->at[0];
    v[74]=nd[11]->at[0];
    v[75]=nd[12]->at[0];
    v[76]=nd[13]->at[0];
    v[388]=v[75]-v[76];
    v[387]=v[74]-v[76];
    v[386]=v[73]-v[76];
    v[112]=es->Data[1];
    v[962]=es->Data[0]/(1e0+v[112]);
    v[113]=(v[112]*v[962])/(1e0-2e0*v[112]);
    v[115]=v[962]/2e0;
    v[116]=es->Data[2];
    v[117]=es->Data[3];
    v[118]=es->Data[4];
    v[119]=es->Data[5];
    v[121]=es->Data[7];
    v[122]=es->Data[8];
    v[1012]=v[121]-v[122];
    v[971]=es->Data[6]/(es->Data[9]*v[122]);
    v[433]=-(v[122]*v[971]);
    v[995]=v[118]*v[433];
    v[994]=v[117]*v[433];
    v[993]=v[116]*v[433];
    v[963]=-(es->Data[10]/es->Data[11]);
    v[190]=1e0+v[963];
    v[127]=rdata[RD_TimeIncrement];
    v[964]=v[127]*(1e0+v[963]/2e0);
    v[189]=-(v[963]/v[127]);
    v[188]=ed->hp[0]*v[190]+v[189]*(-nd[0]->ap[0]+v[43])+ed->hp[30]*v[964];
    v[192]=ed->hp[10]*v[190]+v[189]*(-nd[0]->ap[1]+v[44])+ed->hp[40]*v[964];
    v[193]=ed->hp[20]*v[190]+v[189]*(-nd[0]->ap[2]+v[45])+ed->hp[50]*v[964];
    v[194]=ed->hp[1]*v[190]+v[189]*(-nd[1]->ap[0]+v[46])+ed->hp[31]*v[964];
    v[195]=ed->hp[11]*v[190]+v[189]*(-nd[1]->ap[1]+v[47])+ed->hp[41]*v[964];
    v[196]=ed->hp[21]*v[190]+v[189]*(-nd[1]->ap[2]+v[48])+ed->hp[51]*v[964];
    v[197]=ed->hp[2]*v[190]+v[189]*(-nd[2]->ap[0]+v[49])+ed->hp[32]*v[964];
    v[198]=ed->hp[12]*v[190]+v[189]*(-nd[2]->ap[1]+v[50])+ed->hp[42]*v[964];
    v[199]=ed->hp[22]*v[190]+v[189]*(-nd[2]->ap[2]+v[51])+ed->hp[52]*v[964];
    v[200]=ed->hp[3]*v[190]+v[189]*(-nd[3]->ap[0]+v[52])+ed->hp[33]*v[964];
    v[201]=ed->hp[13]*v[190]+v[189]*(-nd[3]->ap[1]+v[53])+ed->hp[43]*v[964];
    v[202]=ed->hp[23]*v[190]+v[189]*(-nd[3]->ap[2]+v[54])+ed->hp[53]*v[964];
    v[203]=ed->hp[4]*v[190]+v[189]*(-nd[4]->ap[0]+v[55])+ed->hp[34]*v[964];
    v[204]=ed->hp[14]*v[190]+v[189]*(-nd[4]->ap[1]+v[56])+ed->hp[44]*v[964];
    v[205]=ed->hp[24]*v[190]+v[189]*(-nd[4]->ap[2]+v[57])+ed->hp[54]*v[964];
    v[206]=ed->hp[5]*v[190]+v[189]*(-nd[5]->ap[0]+v[58])+ed->hp[35]*v[964];
    v[207]=ed->hp[15]*v[190]+v[189]*(-nd[5]->ap[1]+v[59])+ed->hp[45]*v[964];
    v[208]=ed->hp[25]*v[190]+v[189]*(-nd[5]->ap[2]+v[60])+ed->hp[55]*v[964];
    v[209]=ed->hp[6]*v[190]+v[189]*(-nd[6]->ap[0]+v[61])+ed->hp[36]*v[964];
    v[210]=ed->hp[16]*v[190]+v[189]*(-nd[6]->ap[1]+v[62])+ed->hp[46]*v[964];
    v[211]=ed->hp[26]*v[190]+v[189]*(-nd[6]->ap[2]+v[63])+ed->hp[56]*v[964];
    v[212]=ed->hp[7]*v[190]+v[189]*(-nd[7]->ap[0]+v[64])+ed->hp[37]*v[964];
    v[213]=ed->hp[17]*v[190]+v[189]*(-nd[7]->ap[1]+v[65])+ed->hp[47]*v[964];
    v[214]=ed->hp[27]*v[190]+v[189]*(-nd[7]->ap[2]+v[66])+ed->hp[57]*v[964];
    v[215]=ed->hp[8]*v[190]+v[189]*(-nd[8]->ap[0]+v[67])+ed->hp[38]*v[964];
    v[216]=ed->hp[18]*v[190]+v[189]*(-nd[8]->ap[1]+v[68])+ed->hp[48]*v[964];
    v[217]=ed->hp[28]*v[190]+v[189]*(-nd[8]->ap[2]+v[69])+ed->hp[58]*v[964];
    v[218]=ed->hp[9]*v[190]+v[189]*(-nd[9]->ap[0]+v[70])+ed->hp[39]*v[964];
    v[219]=ed->hp[19]*v[190]+v[189]*(-nd[9]->ap[1]+v[71])+ed->hp[49]*v[964];
    v[220]=ed->hp[29]*v[190]+v[189]*(-nd[9]->ap[2]+v[72])+ed->hp[59]*v[964];
    for(i254=1;i254<=es->id.NoIntPoints;i254++){
     i256=4*(-1+i254);
     v[255]=es->IntPoints[i256];
     v[286]=4e0*v[255];
     v[374]=-(v[214]*v[286]);
     v[367]=-(v[213]*v[286]);
     v[360]=-(v[212]*v[286]);
     v[281]=-1e0+v[286];
     v[257]=es->IntPoints[1+i256];
     v[285]=4e0*v[257];
     v[333]=-(v[285]*v[69]);
     v[326]=-(v[285]*v[68]);
     v[319]=-(v[285]*v[67]);
     v[310]=-(v[27]*v[285]);
     v[303]=-(v[26]*v[285]);
     v[296]=-(v[25]*v[285]);
     v[282]=-1e0+v[285];
     v[258]=es->IntPoints[2+i256];
     v[287]=4e0*v[258];
     v[371]=-(v[220]*v[287]);
     v[364]=-(v[219]*v[287]);
     v[357]=-(v[218]*v[287]);
     v[331]=-(v[287]*v[72]);
     v[324]=-(v[287]*v[71]);
     v[317]=-(v[287]*v[70]);
     v[308]=-(v[287]*v[30]);
     v[301]=-(v[287]*v[29]);
     v[294]=-(v[28]*v[287]);
     v[283]=-1e0+v[287];
     v[260]=v[255]*(-1e0+2e0*v[255]);
     v[261]=v[257]*(-1e0+2e0*v[257]);
     v[262]=v[258]*(-1e0+2e0*v[258]);
     v[263]=1e0-v[255]-v[257]-v[258];
     v[1548]=0e0;
     v[1549]=0e0;
     v[1550]=0e0;
     v[1551]=0e0;
     v[1552]=0e0;
     v[1553]=0e0;
     v[1554]=0e0;
     v[1555]=0e0;
     v[1556]=0e0;
     v[1557]=0e0;
     v[1558]=0e0;
     v[1559]=0e0;
     v[1560]=0e0;
     v[1561]=0e0;
     v[1562]=0e0;
     v[1563]=0e0;
     v[1564]=0e0;
     v[1565]=0e0;
     v[1566]=0e0;
     v[1567]=0e0;
     v[1568]=0e0;
     v[1569]=0e0;
     v[1570]=0e0;
     v[1571]=0e0;
     v[1572]=0e0;
     v[1573]=0e0;
     v[1574]=0e0;
     v[1575]=0e0;
     v[1576]=0e0;
     v[1577]=0e0;
     v[1578]=v[255];
     v[1579]=v[257];
     v[1580]=v[258];
     v[1581]=v[263];
     v[288]=-4e0*v[263];
     v[291]=-v[287]-v[288];
     v[290]=-v[285]-v[288];
     v[289]=-v[286]-v[288];
     v[284]=1e0+v[288];
     v[1412]=0e0;
     v[1413]=0e0;
     v[1414]=v[281];
     v[1415]=0e0;
     v[1416]=0e0;
     v[1417]=0e0;
     v[1418]=0e0;
     v[1419]=0e0;
     v[1420]=0e0;
     v[1421]=0e0;
     v[1422]=0e0;
     v[1423]=v[284];
     v[1424]=0e0;
     v[1425]=0e0;
     v[1426]=v[285];
     v[1427]=0e0;
     v[1428]=0e0;
     v[1429]=0e0;
     v[1430]=0e0;
     v[1431]=0e0;
     v[1432]=v[287];
     v[1433]=0e0;
     v[1434]=0e0;
     v[1435]=v[289];
     v[1436]=0e0;
     v[1437]=0e0;
     v[1438]=-v[285];
     v[1439]=0e0;
     v[1440]=0e0;
     v[1441]=-v[287];
     v[1442]=0e0;
     v[1443]=0e0;
     v[1444]=0e0;
     v[1445]=0e0;
     v[1378]=0e0;
     v[1379]=0e0;
     v[1380]=0e0;
     v[1381]=0e0;
     v[1382]=0e0;
     v[1383]=v[282];
     v[1384]=0e0;
     v[1385]=0e0;
     v[1386]=0e0;
     v[1387]=0e0;
     v[1388]=0e0;
     v[1389]=v[284];
     v[1390]=0e0;
     v[1391]=0e0;
     v[1392]=v[286];
     v[1393]=0e0;
     v[1394]=0e0;
     v[1395]=v[287];
     v[1396]=0e0;
     v[1397]=0e0;
     v[1398]=0e0;
     v[1399]=0e0;
     v[1400]=0e0;
     v[1401]=-v[286];
     v[1402]=0e0;
     v[1403]=0e0;
     v[1404]=v[290];
     v[1405]=0e0;
     v[1406]=0e0;
     v[1407]=-v[287];
     v[1408]=0e0;
     v[1409]=0e0;
     v[1410]=0e0;
     v[1411]=0e0;
     v[1344]=0e0;
     v[1345]=0e0;
     v[1346]=0e0;
     v[1347]=0e0;
     v[1348]=0e0;
     v[1349]=0e0;
     v[1350]=0e0;
     v[1351]=0e0;
     v[1352]=v[283];
     v[1353]=0e0;
     v[1354]=0e0;
     v[1355]=v[284];
     v[1356]=0e0;
     v[1357]=0e0;
     v[1358]=0e0;
     v[1359]=0e0;
     v[1360]=0e0;
     v[1361]=v[285];
     v[1362]=0e0;
     v[1363]=0e0;
     v[1364]=v[286];
     v[1365]=0e0;
     v[1366]=0e0;
     v[1367]=-v[286];
     v[1368]=0e0;
     v[1369]=0e0;
     v[1370]=-v[285];
     v[1371]=0e0;
     v[1372]=0e0;
     v[1373]=v[291];
     v[1374]=0e0;
     v[1375]=0e0;
     v[1376]=0e0;
     v[1377]=0e0;
     v[1310]=0e0;
     v[1311]=v[281];
     v[1312]=0e0;
     v[1313]=0e0;
     v[1314]=0e0;
     v[1315]=0e0;
     v[1316]=0e0;
     v[1317]=0e0;
     v[1318]=0e0;
     v[1319]=0e0;
     v[1320]=v[284];
     v[1321]=0e0;
     v[1322]=0e0;
     v[1323]=v[285];
     v[1324]=0e0;
     v[1325]=0e0;
     v[1326]=0e0;
     v[1327]=0e0;
     v[1328]=0e0;
     v[1329]=v[287];
     v[1330]=0e0;
     v[1331]=0e0;
     v[1332]=v[289];
     v[1333]=0e0;
     v[1334]=0e0;
     v[1335]=-v[285];
     v[1336]=0e0;
     v[1337]=0e0;
     v[1338]=-v[287];
     v[1339]=0e0;
     v[1340]=0e0;
     v[1341]=0e0;
     v[1342]=0e0;
     v[1343]=0e0;
     v[1276]=0e0;
     v[1277]=0e0;
     v[1278]=0e0;
     v[1279]=0e0;
     v[1280]=v[282];
     v[1281]=0e0;
     v[1282]=0e0;
     v[1283]=0e0;
     v[1284]=0e0;
     v[1285]=0e0;
     v[1286]=v[284];
     v[1287]=0e0;
     v[1288]=0e0;
     v[1289]=v[286];
     v[1290]=0e0;
     v[1291]=0e0;
     v[1292]=v[287];
     v[1293]=0e0;
     v[1294]=0e0;
     v[1295]=0e0;
     v[1296]=0e0;
     v[1297]=0e0;
     v[1298]=-v[286];
     v[1299]=0e0;
     v[1300]=0e0;
     v[1301]=v[290];
     v[1302]=0e0;
     v[1303]=0e0;
     v[1304]=-v[287];
     v[1305]=0e0;
     v[1306]=0e0;
     v[1307]=0e0;
     v[1308]=0e0;
     v[1309]=0e0;
     v[1242]=0e0;
     v[1243]=0e0;
     v[1244]=0e0;
     v[1245]=0e0;
     v[1246]=0e0;
     v[1247]=0e0;
     v[1248]=0e0;
     v[1249]=v[283];
     v[1250]=0e0;
     v[1251]=0e0;
     v[1252]=v[284];
     v[1253]=0e0;
     v[1254]=0e0;
     v[1255]=0e0;
     v[1256]=0e0;
     v[1257]=0e0;
     v[1258]=v[285];
     v[1259]=0e0;
     v[1260]=0e0;
     v[1261]=v[286];
     v[1262]=0e0;
     v[1263]=0e0;
     v[1264]=-v[286];
     v[1265]=0e0;
     v[1266]=0e0;
     v[1267]=-v[285];
     v[1268]=0e0;
     v[1269]=0e0;
     v[1270]=v[291];
     v[1271]=0e0;
     v[1272]=0e0;
     v[1273]=0e0;
     v[1274]=0e0;
     v[1275]=0e0;
     v[1208]=v[281];
     v[1209]=0e0;
     v[1210]=0e0;
     v[1211]=0e0;
     v[1212]=0e0;
     v[1213]=0e0;
     v[1214]=0e0;
     v[1215]=0e0;
     v[1216]=0e0;
     v[1217]=v[284];
     v[1218]=0e0;
     v[1219]=0e0;
     v[1220]=v[285];
     v[1221]=0e0;
     v[1222]=0e0;
     v[1223]=0e0;
     v[1224]=0e0;
     v[1225]=0e0;
     v[1226]=v[287];
     v[1227]=0e0;
     v[1228]=0e0;
     v[1229]=v[289];
     v[1230]=0e0;
     v[1231]=0e0;
     v[1232]=-v[285];
     v[1233]=0e0;
     v[1234]=0e0;
     v[1235]=-v[287];
     v[1236]=0e0;
     v[1237]=0e0;
     v[1238]=0e0;
     v[1239]=0e0;
     v[1240]=0e0;
     v[1241]=0e0;
     v[1174]=0e0;
     v[1175]=0e0;
     v[1176]=0e0;
     v[1177]=v[282];
     v[1178]=0e0;
     v[1179]=0e0;
     v[1180]=0e0;
     v[1181]=0e0;
     v[1182]=0e0;
     v[1183]=v[284];
     v[1184]=0e0;
     v[1185]=0e0;
     v[1186]=v[286];
     v[1187]=0e0;
     v[1188]=0e0;
     v[1189]=v[287];
     v[1190]=0e0;
     v[1191]=0e0;
     v[1192]=0e0;
     v[1193]=0e0;
     v[1194]=0e0;
     v[1195]=-v[286];
     v[1196]=0e0;
     v[1197]=0e0;
     v[1198]=v[290];
     v[1199]=0e0;
     v[1200]=0e0;
     v[1201]=-v[287];
     v[1202]=0e0;
     v[1203]=0e0;
     v[1204]=0e0;
     v[1205]=0e0;
     v[1206]=0e0;
     v[1207]=0e0;
     v[1140]=0e0;
     v[1141]=0e0;
     v[1142]=0e0;
     v[1143]=0e0;
     v[1144]=0e0;
     v[1145]=0e0;
     v[1146]=v[283];
     v[1147]=0e0;
     v[1148]=0e0;
     v[1149]=v[284];
     v[1150]=0e0;
     v[1151]=0e0;
     v[1152]=0e0;
     v[1153]=0e0;
     v[1154]=0e0;
     v[1155]=v[285];
     v[1156]=0e0;
     v[1157]=0e0;
     v[1158]=v[286];
     v[1159]=0e0;
     v[1160]=0e0;
     v[1161]=-v[286];
     v[1162]=0e0;
     v[1163]=0e0;
     v[1164]=-v[285];
     v[1165]=0e0;
     v[1166]=0e0;
     v[1167]=v[291];
     v[1168]=0e0;
     v[1169]=0e0;
     v[1170]=0e0;
     v[1171]=0e0;
     v[1172]=0e0;
     v[1173]=0e0;
     v[370]=v[202]*v[284];
     v[977]=-(v[217]*v[285])+v[370];
     v[363]=v[201]*v[284];
     v[976]=-(v[216]*v[285])+v[363];
     v[356]=v[200]*v[284];
     v[975]=-(v[215]*v[285])+v[356];
     v[330]=v[284]*v[54];
     v[965]=v[330]-v[286]*v[66];
     v[335]=v[333]+v[283]*v[51]+v[285]*v[60]+v[286]*v[63]+v[291]*v[72]+v[965];
     v[332]=v[331]+v[282]*v[48]+v[286]*v[57]+v[287]*v[60]+v[290]*v[69]+v[965];
     v[329]=v[330]+v[331]+v[333]+v[281]*v[45]+v[285]*v[57]+v[287]*v[63]+v[289]*v[66];
     v[323]=v[284]*v[53];
     v[966]=v[323]-v[286]*v[65];
     v[328]=v[326]+v[283]*v[50]+v[285]*v[59]+v[286]*v[62]+v[291]*v[71]+v[966];
     v[325]=v[324]+v[282]*v[47]+v[286]*v[56]+v[287]*v[59]+v[290]*v[68]+v[966];
     v[322]=v[323]+v[324]+v[326]+v[281]*v[44]+v[285]*v[56]+v[287]*v[62]+v[289]*v[65];
     v[316]=v[284]*v[52];
     v[967]=v[316]-v[286]*v[64];
     v[321]=v[319]+v[283]*v[49]+v[285]*v[58]+v[286]*v[61]+v[291]*v[70]+v[967];
     v[318]=v[317]+v[282]*v[46]+v[286]*v[55]+v[287]*v[58]+v[290]*v[67]+v[967];
     v[315]=v[316]+v[317]+v[319]+v[281]*v[43]+v[285]*v[55]+v[287]*v[61]+v[289]*v[64];
     v[307]=v[12]*v[284];
     v[968]=-(v[24]*v[286])+v[307];
     v[312]=v[18]*v[285]+v[21]*v[286]+v[291]*v[30]+v[310]+v[283]*v[9]+v[968];
     v[309]=v[15]*v[286]+v[18]*v[287]+v[27]*v[290]+v[308]+v[282]*v[6]+v[968];
     v[306]=v[15]*v[285]+v[21]*v[287]+v[24]*v[289]+v[281]*v[3]+v[307]+v[308]+v[310];
     v[300]=v[11]*v[284];
     v[969]=-(v[23]*v[286])+v[300];
     v[305]=v[17]*v[285]+v[20]*v[286]+v[29]*v[291]+v[303]+v[283]*v[8]+v[969];
     v[302]=v[14]*v[286]+v[17]*v[287]+v[26]*v[290]+v[301]+v[282]*v[5]+v[969];
     v[340]=-(v[305]*v[309])+v[302]*v[312];
     v[299]=v[2]*v[281]+v[14]*v[285]+v[20]*v[287]+v[23]*v[289]+v[300]+v[301]+v[303];
     v[342]=-(v[302]*v[306])+v[299]*v[309];
     v[341]=v[305]*v[306]-v[299]*v[312];
     v[293]=v[10]*v[284];
     v[970]=-(v[22]*v[286])+v[293];
     v[298]=v[16]*v[285]+v[19]*v[286]+v[28]*v[291]+v[296]+v[283]*v[7]+v[970];
     v[295]=v[13]*v[286]+v[16]*v[287]+v[25]*v[290]+v[294]+v[282]*v[4]+v[970];
     v[348]=-(v[298]*v[302])+v[295]*v[305];
     v[344]=v[298]*v[309]-v[295]*v[312];
     v[292]=v[1]*v[281]+v[13]*v[285]+v[19]*v[287]+v[22]*v[289]+v[293]+v[294]+v[296];
     v[350]=-(v[295]*v[299])+v[292]*v[302];
     v[349]=v[298]*v[299]-v[292]*v[305];
     v[346]=v[295]*v[306]-v[292]*v[309];
     v[345]=-(v[298]*v[306])+v[292]*v[312];
     v[336]=v[292]*v[340]+v[295]*v[341]+v[298]*v[342];
     v[981]=v[115]/(v[336]*v[336]);
     v[972]=v[115]/v[336];
     v[480]=v[971]/v[336];
     v[521]=v[350]/v[336];
     v[520]=v[349]/v[336];
     v[1027]=v[520]+v[521];
     v[1021]=v[287]*v[520]+v[285]*v[521];
     v[519]=v[348]/v[336];
     v[1030]=v[519]+v[521];
     v[1024]=v[287]*v[519]+v[286]*v[521];
     v[1018]=v[285]*v[519]+v[286]*v[520];
     v[985]=v[519]+v[520];
     v[1015]=v[521]+v[985];
     v[518]=v[346]/v[336];
     v[517]=v[345]/v[336];
     v[1026]=v[517]+v[518];
     v[1020]=v[287]*v[517]+v[285]*v[518];
     v[516]=v[344]/v[336];
     v[1029]=v[516]+v[518];
     v[1023]=v[287]*v[516]+v[286]*v[518];
     v[1017]=v[285]*v[516]+v[286]*v[517];
     v[984]=v[516]+v[517];
     v[1014]=v[518]+v[984];
     v[515]=v[342]/v[336];
     v[514]=v[341]/v[336];
     v[1025]=v[514]+v[515];
     v[1019]=v[287]*v[514]+v[285]*v[515];
     v[513]=v[340]/v[336];
     v[1028]=v[513]+v[515];
     v[1022]=v[287]*v[513]+v[286]*v[515];
     v[1016]=v[285]*v[513]+v[286]*v[514];
     v[983]=v[513]+v[514];
     v[1013]=v[515]+v[983];
     v[264]=v[263]*(-1e0+2e0*v[263]);
     v[265]=4e0*v[255]*v[257];
     v[266]=4e0*v[257]*v[258];
     v[267]=4e0*v[255]*v[258];
     v[268]=4e0*v[255]*v[263];
     v[269]=4e0*v[257]*v[263];
     v[270]=4e0*v[258]*v[263];
     v[1582]=0e0;
     v[1583]=0e0;
     v[1584]=v[260];
     v[1585]=0e0;
     v[1586]=0e0;
     v[1587]=v[261];
     v[1588]=0e0;
     v[1589]=0e0;
     v[1590]=v[262];
     v[1591]=0e0;
     v[1592]=0e0;
     v[1593]=v[264];
     v[1594]=0e0;
     v[1595]=0e0;
     v[1596]=v[265];
     v[1597]=0e0;
     v[1598]=0e0;
     v[1599]=v[266];
     v[1600]=0e0;
     v[1601]=0e0;
     v[1602]=v[267];
     v[1603]=0e0;
     v[1604]=0e0;
     v[1605]=v[268];
     v[1606]=0e0;
     v[1607]=0e0;
     v[1608]=v[269];
     v[1609]=0e0;
     v[1610]=0e0;
     v[1611]=v[270];
     v[1612]=0e0;
     v[1613]=0e0;
     v[1614]=0e0;
     v[1615]=0e0;
     v[1616]=0e0;
     v[1617]=v[260];
     v[1618]=0e0;
     v[1619]=0e0;
     v[1620]=v[261];
     v[1621]=0e0;
     v[1622]=0e0;
     v[1623]=v[262];
     v[1624]=0e0;
     v[1625]=0e0;
     v[1626]=v[264];
     v[1627]=0e0;
     v[1628]=0e0;
     v[1629]=v[265];
     v[1630]=0e0;
     v[1631]=0e0;
     v[1632]=v[266];
     v[1633]=0e0;
     v[1634]=0e0;
     v[1635]=v[267];
     v[1636]=0e0;
     v[1637]=0e0;
     v[1638]=v[268];
     v[1639]=0e0;
     v[1640]=0e0;
     v[1641]=v[269];
     v[1642]=0e0;
     v[1643]=0e0;
     v[1644]=v[270];
     v[1645]=0e0;
     v[1646]=0e0;
     v[1647]=0e0;
     v[1648]=0e0;
     v[1649]=0e0;
     v[1650]=v[260];
     v[1651]=0e0;
     v[1652]=0e0;
     v[1653]=v[261];
     v[1654]=0e0;
     v[1655]=0e0;
     v[1656]=v[262];
     v[1657]=0e0;
     v[1658]=0e0;
     v[1659]=v[264];
     v[1660]=0e0;
     v[1661]=0e0;
     v[1662]=v[265];
     v[1663]=0e0;
     v[1664]=0e0;
     v[1665]=v[266];
     v[1666]=0e0;
     v[1667]=0e0;
     v[1668]=v[267];
     v[1669]=0e0;
     v[1670]=0e0;
     v[1671]=v[268];
     v[1672]=0e0;
     v[1673]=0e0;
     v[1674]=v[269];
     v[1675]=0e0;
     v[1676]=0e0;
     v[1677]=v[270];
     v[1678]=0e0;
     v[1679]=0e0;
     v[1680]=0e0;
     v[1681]=0e0;
     v[1682]=0e0;
     v[1683]=0e0;
     v[314]=es->IntPoints[3+i256]*fabs(v[336]);
     v[393]=(v[315]*v[340]+v[318]*v[341]+v[321]*v[342])/v[336];
     v[396]=(v[322]*v[344]+v[325]*v[345]+v[328]*v[346])/v[336];
     v[398]=(v[329]*v[348]+v[332]*v[349]+v[335]*v[350])/v[336];
     v[974]=v[393]+v[396]+v[398];
     v[973]=-(v[255]*v[73])-v[257]*v[74]-v[258]*v[75]-v[263]*v[76]+v[113]*v[974];
     v[402]=(v[322]*v[340]+v[325]*v[341]+v[328]*v[342]+v[315]*v[344]+v[318]*v[345]+v[321]*v[346]
      )*v[972];
     v[403]=(v[329]*v[340]+v[332]*v[341]+v[335]*v[342]+v[315]*v[348]+v[318]*v[349]+v[321]*v[350]
      )*v[972];
     v[406]=(v[329]*v[344]+v[332]*v[345]+v[335]*v[346]+v[322]*v[348]+v[325]*v[349]+v[328]*v[350]
      )*v[972];
     v[408]=2e0*v[115]*v[393]+v[973];
     v[409]=2e0*v[115]*v[396]+v[973];
     v[410]=2e0*v[115]*v[398]+v[973];
     v[417]=-(v[119]*(-1e0+v[974]));
     v[426]=v[122]*(-1e0+v[417])-v[121]*v[417];
     v[425]=v[116]*v[426];
     v[427]=v[117]*v[426];
     v[428]=v[118]*v[426];
     v[430]=(v[341]*(v[194]*v[282]+v[203]*v[286]+v[206]*v[287]+v[215]*v[290]+v[356]+v[357]+v[360])
      +v[345]*(v[195]*v[282]+v[204]*v[286]+v[207]*v[287]+v[216]*v[290]+v[363]+v[364]+v[367])+v[349]*
      (v[196]*v[282]+v[205]*v[286]+v[208]*v[287]+v[217]*v[290]+v[370]+v[371]+v[374])+v[340]*
      (v[188]*v[281]+v[203]*v[285]+v[209]*v[287]+v[212]*v[289]+v[357]+v[975])+v[342]*(v[197]*v[283]
      +v[206]*v[285]+v[209]*v[286]+v[218]*v[291]+v[360]+v[975])+v[344]*(v[192]*v[281]+v[204]*v[285]
      +v[210]*v[287]+v[213]*v[289]+v[364]+v[976])+v[346]*(v[198]*v[283]+v[207]*v[285]+v[210]*v[286]
      +v[219]*v[291]+v[367]+v[976])+v[348]*(v[193]*v[281]+v[205]*v[285]+v[211]*v[287]+v[214]*v[289]
      +v[371]+v[977])+v[350]*(v[199]*v[283]+v[208]*v[285]+v[211]*v[286]+v[220]*v[291]+v[374]+v[977]))
      /v[336];
     v[432]=(v[340]*v[386]+v[341]*v[387]+v[342]*v[388])*v[480]+v[993];
     v[434]=(v[344]*v[386]+v[345]*v[387]+v[346]*v[388])*v[480]+v[994];
     v[435]=(v[348]*v[386]+v[349]*v[387]+v[350]*v[388])*v[480]+v[995];
     v[440]=(v[340]*v[432]+v[344]*v[434]+v[348]*v[435])/v[336];
     v[441]=(v[341]*v[432]+v[345]*v[434]+v[349]*v[435])/v[336];
     v[442]=(v[342]*v[432]+v[346]*v[434]+v[350]*v[435])/v[336];
     v[443]=(v[340]*v[403]+v[344]*v[406]+v[348]*v[410])/v[336];
     v[444]=(v[341]*v[403]+v[345]*v[406]+v[349]*v[410])/v[336];
     v[445]=(v[342]*v[403]+v[346]*v[406]+v[350]*v[410])/v[336];
     v[980]=v[444]+v[445];
     v[446]=(v[340]*v[402]+v[348]*v[406]+v[344]*v[409])/v[336];
     v[447]=(v[341]*v[402]+v[349]*v[406]+v[345]*v[409])/v[336];
     v[448]=(v[342]*v[402]+v[350]*v[406]+v[346]*v[409])/v[336];
     v[979]=v[447]+v[448];
     v[449]=(v[344]*v[402]+v[348]*v[403]+v[340]*v[408])/v[336];
     v[450]=(v[345]*v[402]+v[349]*v[403]+v[341]*v[408])/v[336];
     v[451]=(v[346]*v[402]+v[350]*v[403]+v[342]*v[408])/v[336];
     v[978]=v[450]+v[451];
     v[1102]=v[260]*v[425]+v[281]*v[449];
     v[1103]=v[260]*v[427]+v[281]*v[446];
     v[1104]=v[260]*v[428]+v[281]*v[443];
     v[1105]=v[261]*v[425]+v[282]*v[450];
     v[1106]=v[261]*v[427]+v[282]*v[447];
     v[1107]=v[261]*v[428]+v[282]*v[444];
     v[1108]=v[262]*v[425]+v[283]*v[451];
     v[1109]=v[262]*v[427]+v[283]*v[448];
     v[1110]=v[262]*v[428]+v[283]*v[445];
     v[1111]=v[264]*v[425]+v[284]*(v[449]+v[978]);
     v[1112]=v[264]*v[427]+v[284]*(v[446]+v[979]);
     v[1113]=v[264]*v[428]+v[284]*(v[443]+v[980]);
     v[1114]=v[265]*v[425]+v[285]*v[449]+v[286]*v[450];
     v[1115]=v[265]*v[427]+v[285]*v[446]+v[286]*v[447];
     v[1116]=v[265]*v[428]+v[285]*v[443]+v[286]*v[444];
     v[1117]=v[266]*v[425]+v[287]*v[450]+v[285]*v[451];
     v[1118]=v[266]*v[427]+v[287]*v[447]+v[285]*v[448];
     v[1119]=v[266]*v[428]+v[287]*v[444]+v[285]*v[445];
     v[1120]=v[267]*v[425]+v[287]*v[449]+v[286]*v[451];
     v[1121]=v[267]*v[427]+v[287]*v[446]+v[286]*v[448];
     v[1122]=v[267]*v[428]+v[287]*v[443]+v[286]*v[445];
     v[1123]=v[268]*v[425]+v[289]*v[449]-v[286]*v[978];
     v[1124]=v[268]*v[427]+v[289]*v[446]-v[286]*v[979];
     v[1125]=v[268]*v[428]+v[289]*v[443]-v[286]*v[980];
     v[1126]=v[269]*v[425]+v[290]*v[450]-v[285]*(v[449]+v[451]);
     v[1127]=v[269]*v[427]+v[290]*v[447]-v[285]*(v[446]+v[448]);
     v[1128]=v[269]*v[428]+v[290]*v[444]-v[285]*(v[443]+v[445]);
     v[1129]=v[270]*v[425]-v[287]*(v[449]+v[450])+v[291]*v[451];
     v[1130]=v[270]*v[427]-v[287]*(v[446]+v[447])+v[291]*v[448];
     v[1131]=v[270]*v[428]-v[287]*(v[443]+v[444])+v[291]*v[445];
     v[1132]=v[255]*v[430]+v[440];
     v[1133]=v[257]*v[430]+v[441];
     v[1134]=v[258]*v[430]+v[442];
     v[1135]=v[263]*v[430]-v[440]-v[441]-v[442];
     for(i438=1;i438<=34;i438++){
      v[458]=v[1139+i438];
      v[459]=v[1173+i438];
      v[460]=v[1207+i438];
      v[461]=v[1241+i438];
      v[462]=v[1275+i438];
      v[463]=v[1309+i438];
      v[464]=v[1343+i438];
      v[465]=v[1377+i438];
      v[466]=v[1411+i438];
      v[467]=v[1445+i438];
      v[468]=v[1479+i438];
      v[469]=v[1513+i438];
      v[470]=v[1547+i438];
      v[992]=-(v[287]*v[470]);
      v[991]=-(v[285]*v[470]);
      v[988]=-(v[286]*v[470]);
      v[986]=v[189]*v[470];
      v[512]=v[284]*v[986];
      v[511]=v[470]*v[521];
      v[510]=v[470]*v[518];
      v[509]=v[470]*v[515];
      v[508]=v[470]*v[520];
      v[507]=v[470]*v[517];
      v[506]=v[470]*v[514];
      v[505]=v[470]*v[519];
      v[504]=v[470]*v[516];
      v[503]=v[470]*v[513];
      v[471]=(v[342]*v[458]+v[341]*v[459]+v[340]*v[460])/v[336];
      v[472]=(v[346]*v[461]+v[345]*v[462]+v[344]*v[463])/v[336];
      v[497]=(v[346]*v[458]+v[345]*v[459]+v[344]*v[460]+v[342]*v[461]+v[341]*v[462]+v[340]*v[463]
       )*v[981];
      v[474]=(v[350]*v[464]+v[349]*v[465]+v[348]*v[466])/v[336];
      v[492]=(v[350]*v[461]+v[349]*v[462]+v[348]*v[463]+v[346]*v[464]+v[345]*v[465]+v[344]*v[466]
       )*v[981];
      v[493]=(v[350]*v[458]+v[349]*v[459]+v[348]*v[460]+v[342]*v[464]+v[341]*v[465]+v[340]*v[466]
       )*v[981];
      v[477]=(v[350]*v[467]+v[349]*v[468]+v[348]*v[469])/v[336];
      v[478]=(v[346]*v[467]+v[345]*v[468]+v[344]*v[469])/v[336];
      v[479]=(v[342]*v[467]+v[341]*v[468]+v[340]*v[469])/v[336];
      v[481]=(v[348]*v[477]+v[344]*v[478]+v[340]*v[479])*v[480];
      v[482]=(v[349]*v[477]+v[345]*v[478]+v[341]*v[479])*v[480];
      v[483]=(v[350]*v[477]+v[346]*v[478]+v[342]*v[479])*v[480];
      v[485]=v[471]+v[472]+v[474];
      v[982]=v[1012]*v[119]*(v[118]*v[1581+i438]+v[117]*v[1615+i438]+v[116]*v[1649+i438])
       +v[113]*v[485];
      v[486]=2e0*v[115]*v[474]+v[982];
      v[489]=2e0*v[115]*v[472]+v[982];
      v[490]=2e0*v[115]*v[471]+v[982];
      v[491]=v[344]*v[492]+v[340]*v[493]+v[486]*v[519];
      v[494]=v[345]*v[492]+v[341]*v[493]+v[486]*v[520];
      v[495]=v[346]*v[492]+v[342]*v[493]+v[486]*v[521];
      v[990]=v[494]+v[495];
      v[496]=v[348]*v[492]+v[340]*v[497]+v[489]*v[516];
      v[498]=v[349]*v[492]+v[341]*v[497]+v[489]*v[517];
      v[499]=v[350]*v[492]+v[342]*v[497]+v[489]*v[518];
      v[989]=v[498]+v[499];
      v[500]=v[348]*v[493]+v[344]*v[497]+v[490]*v[513];
      v[501]=v[349]*v[493]+v[345]*v[497]+v[490]*v[514];
      v[502]=v[350]*v[493]+v[346]*v[497]+v[490]*v[515];
      v[987]=v[501]+v[502];
      v[1684]=v[281]*(v[500]+v[189]*v[503]);
      v[1685]=v[281]*(v[496]+v[189]*v[504]);
      v[1686]=v[281]*(v[491]+v[189]*v[505]);
      v[1687]=v[282]*(v[501]+v[189]*v[506]);
      v[1688]=v[282]*(v[498]+v[189]*v[507]);
      v[1689]=v[282]*(v[494]+v[189]*v[508]);
      v[1690]=v[283]*(v[502]+v[189]*v[509]);
      v[1691]=v[283]*(v[499]+v[189]*v[510]);
      v[1692]=v[283]*(v[495]+v[189]*v[511]);
      v[1693]=v[1013]*v[512]+v[284]*(v[500]+v[987]);
      v[1694]=v[1014]*v[512]+v[284]*(v[496]+v[989]);
      v[1695]=v[1015]*v[512]+v[284]*(v[491]+v[990]);
      v[1696]=v[285]*v[500]+v[286]*v[501]+v[1016]*v[986];
      v[1697]=v[285]*v[496]+v[286]*v[498]+v[1017]*v[986];
      v[1698]=v[285]*v[491]+v[286]*v[494]+v[1018]*v[986];
      v[1699]=v[287]*v[501]+v[285]*v[502]+v[1019]*v[986];
      v[1700]=v[287]*v[498]+v[285]*v[499]+v[1020]*v[986];
      v[1701]=v[287]*v[494]+v[285]*v[495]+v[1021]*v[986];
      v[1702]=v[287]*v[500]+v[286]*v[502]+v[1022]*v[986];
      v[1703]=v[287]*v[496]+v[286]*v[499]+v[1023]*v[986];
      v[1704]=v[287]*v[491]+v[286]*v[495]+v[1024]*v[986];
      v[1705]=v[289]*v[500]-v[286]*v[987]+v[189]*(v[289]*v[503]+v[1025]*v[988]);
      v[1706]=v[289]*v[496]+v[189]*(v[289]*v[504]+v[1026]*v[988])-v[286]*v[989];
      v[1707]=v[289]*v[491]+v[189]*(v[289]*v[505]+v[1027]*v[988])-v[286]*v[990];
      v[1708]=v[290]*v[501]-v[285]*(v[500]+v[502])+v[189]*(v[290]*v[506]+v[1028]*v[991]);
      v[1709]=v[290]*v[498]-v[285]*(v[496]+v[499])+v[189]*(v[290]*v[507]+v[1029]*v[991]);
      v[1710]=v[290]*v[494]-v[285]*(v[491]+v[495])+v[189]*(v[290]*v[508]+v[1030]*v[991]);
      v[1711]=-(v[287]*(v[500]+v[501]))+v[291]*v[502]+v[189]*(v[291]*v[509]+v[983]*v[992]);
      v[1712]=-(v[287]*(v[496]+v[498]))+v[291]*v[499]+v[189]*(v[291]*v[510]+v[984]*v[992]);
      v[1713]=-(v[287]*(v[491]+v[494]))+v[291]*v[495]+v[189]*(v[291]*v[511]+v[985]*v[992]);
      v[1714]=v[481]-v[255]*v[485];
      v[1715]=v[482]-v[257]*v[485];
      v[1716]=v[483]-v[258]*v[485];
      v[1717]=-v[481]-v[482]-v[483]-v[263]*v[485];
      p[i438-1]+=v[1101+i438]*v[314];
      for(i455=1;i455<=34;i455++){
       s[i438-1][i455-1]+=v[1683+i455]*v[314];
      };/* end for */
     };/* end for */
    };/* end for */
    };

    /******************* S U B R O U T I N E *********************/
    void SPP(double v[1817],ElementSpec *es,ElementData *ed,NodeSpec **ns
         ,NodeData **nd,double *rdata,int *idata,double **gpost,double **npost)
    {
    int i777,i779;
    v[524]=nd[0]->X[0];
    v[525]=nd[0]->X[1];
    v[526]=nd[0]->X[2];
    v[527]=nd[1]->X[0];
    v[528]=nd[1]->X[1];
    v[529]=nd[1]->X[2];
    v[530]=nd[2]->X[0];
    v[531]=nd[2]->X[1];
    v[532]=nd[2]->X[2];
    v[533]=nd[3]->X[0];
    v[534]=nd[3]->X[1];
    v[535]=nd[3]->X[2];
    v[536]=nd[4]->X[0];
    v[537]=nd[4]->X[1];
    v[538]=nd[4]->X[2];
    v[539]=nd[5]->X[0];
    v[540]=nd[5]->X[1];
    v[541]=nd[5]->X[2];
    v[542]=nd[6]->X[0];
    v[543]=nd[6]->X[1];
    v[544]=nd[6]->X[2];
    v[545]=nd[7]->X[0];
    v[546]=nd[7]->X[1];
    v[547]=nd[7]->X[2];
    v[548]=nd[8]->X[0];
    v[549]=nd[8]->X[1];
    v[550]=nd[8]->X[2];
    v[551]=nd[9]->X[0];
    v[552]=nd[9]->X[1];
    v[553]=nd[9]->X[2];
    v[566]=nd[0]->at[0];
    v[567]=nd[0]->at[1];
    v[568]=nd[0]->at[2];
    v[569]=nd[1]->at[0];
    v[570]=nd[1]->at[1];
    v[571]=nd[1]->at[2];
    v[572]=nd[2]->at[0];
    v[573]=nd[2]->at[1];
    v[574]=nd[2]->at[2];
    v[575]=nd[3]->at[0];
    v[576]=nd[3]->at[1];
    v[577]=nd[3]->at[2];
    v[578]=nd[4]->at[0];
    v[579]=nd[4]->at[1];
    v[580]=nd[4]->at[2];
    v[581]=nd[5]->at[0];
    v[582]=nd[5]->at[1];
    v[583]=nd[5]->at[2];
    v[584]=nd[6]->at[0];
    v[585]=nd[6]->at[1];
    v[586]=nd[6]->at[2];
    v[587]=nd[7]->at[0];
    v[588]=nd[7]->at[1];
    v[589]=nd[7]->at[2];
    v[590]=nd[8]->at[0];
    v[591]=nd[8]->at[1];
    v[592]=nd[8]->at[2];
    v[593]=nd[9]->at[0];
    v[594]=nd[9]->at[1];
    v[595]=nd[9]->at[2];
    v[596]=nd[10]->at[0];
    v[597]=nd[11]->at[0];
    v[598]=nd[12]->at[0];
    v[599]=nd[13]->at[0];
    v[911]=v[598]-v[599];
    v[910]=v[597]-v[599];
    v[909]=v[596]-v[599];
    v[635]=es->Data[1];
    v[1031]=es->Data[0]/(1e0+v[635]);
    v[636]=(v[1031]*v[635])/(1e0-2e0*v[635]);
    v[638]=v[1031]/2e0;
    v[642]=es->Data[5];
    v[643]=es->Data[6];
    for(i777=1;i777<=es->id.NoIntPoints;i777++){
     i779=4*(-1+i777);
     v[778]=es->IntPoints[i779];
     v[809]=4e0*v[778];
     v[804]=-1e0+v[809];
     v[780]=es->IntPoints[1+i779];
     v[808]=4e0*v[780];
     v[856]=-(v[592]*v[808]);
     v[849]=-(v[591]*v[808]);
     v[842]=-(v[590]*v[808]);
     v[833]=-(v[550]*v[808]);
     v[826]=-(v[549]*v[808]);
     v[819]=-(v[548]*v[808]);
     v[805]=-1e0+v[808];
     v[781]=es->IntPoints[2+i779];
     v[810]=4e0*v[781];
     v[854]=-(v[595]*v[810]);
     v[847]=-(v[594]*v[810]);
     v[840]=-(v[593]*v[810]);
     v[831]=-(v[553]*v[810]);
     v[824]=-(v[552]*v[810]);
     v[817]=-(v[551]*v[810]);
     v[806]=-1e0+v[810];
     v[783]=v[778]*(-1e0+2e0*v[778]);
     v[784]=v[780]*(-1e0+2e0*v[780]);
     v[785]=v[781]*(-1e0+2e0*v[781]);
     v[786]=1e0-v[778]-v[780]-v[781];
     v[811]=-4e0*v[786];
     v[814]=-v[810]-v[811];
     v[813]=-v[808]-v[811];
     v[812]=-v[809]-v[811];
     v[807]=1e0+v[811];
     v[853]=v[577]*v[807];
     v[1032]=-(v[589]*v[809])+v[853];
     v[858]=v[1032]+v[574]*v[806]+v[583]*v[808]+v[586]*v[809]+v[595]*v[814]+v[856];
     v[855]=v[1032]+v[571]*v[805]+v[580]*v[809]+v[583]*v[810]+v[592]*v[813]+v[854];
     v[852]=v[568]*v[804]+v[580]*v[808]+v[586]*v[810]+v[589]*v[812]+v[853]+v[854]+v[856];
     v[846]=v[576]*v[807];
     v[1033]=-(v[588]*v[809])+v[846];
     v[851]=v[1033]+v[573]*v[806]+v[582]*v[808]+v[585]*v[809]+v[594]*v[814]+v[849];
     v[848]=v[1033]+v[570]*v[805]+v[579]*v[809]+v[582]*v[810]+v[591]*v[813]+v[847];
     v[845]=v[567]*v[804]+v[579]*v[808]+v[585]*v[810]+v[588]*v[812]+v[846]+v[847]+v[849];
     v[839]=v[575]*v[807];
     v[1034]=-(v[587]*v[809])+v[839];
     v[844]=v[1034]+v[572]*v[806]+v[581]*v[808]+v[584]*v[809]+v[593]*v[814]+v[842];
     v[841]=v[1034]+v[569]*v[805]+v[578]*v[809]+v[581]*v[810]+v[590]*v[813]+v[840];
     v[838]=v[566]*v[804]+v[578]*v[808]+v[584]*v[810]+v[587]*v[812]+v[839]+v[840]+v[842];
     v[830]=v[535]*v[807];
     v[1035]=-(v[547]*v[809])+v[830];
     v[835]=v[1035]+v[532]*v[806]+v[541]*v[808]+v[544]*v[809]+v[553]*v[814]+v[833];
     v[832]=v[1035]+v[529]*v[805]+v[538]*v[809]+v[541]*v[810]+v[550]*v[813]+v[831];
     v[829]=v[526]*v[804]+v[538]*v[808]+v[544]*v[810]+v[547]*v[812]+v[830]+v[831]+v[833];
     v[823]=v[534]*v[807];
     v[1036]=-(v[546]*v[809])+v[823];
     v[828]=v[1036]+v[531]*v[806]+v[540]*v[808]+v[543]*v[809]+v[552]*v[814]+v[826];
     v[825]=v[1036]+v[528]*v[805]+v[537]*v[809]+v[540]*v[810]+v[549]*v[813]+v[824];
     v[863]=-(v[828]*v[832])+v[825]*v[835];
     v[822]=v[525]*v[804]+v[537]*v[808]+v[543]*v[810]+v[546]*v[812]+v[823]+v[824]+v[826];
     v[865]=-(v[825]*v[829])+v[822]*v[832];
     v[864]=v[828]*v[829]-v[822]*v[835];
     v[816]=v[533]*v[807];
     v[1037]=-(v[545]*v[809])+v[816];
     v[821]=v[1037]+v[530]*v[806]+v[539]*v[808]+v[542]*v[809]+v[551]*v[814]+v[819];
     v[818]=v[1037]+v[527]*v[805]+v[536]*v[809]+v[539]*v[810]+v[548]*v[813]+v[817];
     v[871]=-(v[821]*v[825])+v[818]*v[828];
     v[867]=v[821]*v[832]-v[818]*v[835];
     v[815]=v[524]*v[804]+v[536]*v[808]+v[542]*v[810]+v[545]*v[812]+v[816]+v[817]+v[819];
     v[873]=-(v[818]*v[822])+v[815]*v[825];
     v[872]=v[821]*v[822]-v[815]*v[828];
     v[869]=v[818]*v[829]-v[815]*v[832];
     v[868]=-(v[821]*v[829])+v[815]*v[835];
     v[859]=v[815]*v[863]+v[818]*v[864]+v[821]*v[865];
     v[1048]=v[873]/v[859];
     v[1047]=v[872]/v[859];
     v[1046]=v[871]/v[859];
     v[1045]=v[869]/v[859];
     v[1044]=v[868]/v[859];
     v[1043]=v[867]/v[859];
     v[1042]=v[865]/v[859];
     v[1041]=v[864]/v[859];
     v[1040]=v[863]/v[859];
     v[787]=v[786]*(-1e0+2e0*v[786]);
     v[788]=4e0*v[778]*v[780];
     v[789]=4e0*v[780]*v[781];
     v[790]=4e0*v[778]*v[781];
     v[791]=4e0*v[778]*v[786];
     v[792]=4e0*v[780]*v[786];
     v[793]=4e0*v[781]*v[786];
     v[803]=v[596]*v[778]+v[597]*v[780]+v[598]*v[781]+v[599]*v[786];
     v[916]=(v[838]*v[863]+v[841]*v[864]+v[844]*v[865])/v[859];
     v[919]=(v[845]*v[867]+v[848]*v[868]+v[851]*v[869])/v[859];
     v[921]=(v[852]*v[871]+v[855]*v[872]+v[858]*v[873])/v[859];
     v[1038]=v[916]+v[919]+v[921];
     v[1039]=v[1038]*v[636]-v[803];
     v[1049]=-(v[643]/(1e0+(-1e0+v[1038])*v[642]));
     gpost[i777-1][0]=v[803];
     gpost[i777-1][1]=v[566]*v[783]+v[569]*v[784]+v[572]*v[785]+v[575]*v[787]+v[578]*v[788]
      +v[581]*v[789]+v[584]*v[790]+v[587]*v[791]+v[590]*v[792]+v[593]*v[793];
     gpost[i777-1][2]=v[567]*v[783]+v[570]*v[784]+v[573]*v[785]+v[576]*v[787]+v[579]*v[788]
      +v[582]*v[789]+v[585]*v[790]+v[588]*v[791]+v[591]*v[792]+v[594]*v[793];
     gpost[i777-1][3]=v[568]*v[783]+v[571]*v[784]+v[574]*v[785]+v[577]*v[787]+v[580]*v[788]
      +v[583]*v[789]+v[586]*v[790]+v[589]*v[791]+v[592]*v[792]+v[595]*v[793];
     gpost[i777-1][4]=v[1039]+v[1031]*v[916];
     gpost[i777-1][5]=v[1039]+2e0*v[638]*v[919];
     gpost[i777-1][6]=v[1039]+v[1031]*v[921];
     gpost[i777-1][7]=v[638]*(v[1043]*v[838]+v[1044]*v[841]+v[1045]*v[844]+v[1040]*v[845]
      +v[1041]*v[848]+v[1042]*v[851]);
     gpost[i777-1][8]=v[638]*(v[1046]*v[845]+v[1047]*v[848]+v[1048]*v[851]+v[1043]*v[852]
      +v[1044]*v[855]+v[1045]*v[858]);
     gpost[i777-1][9]=v[638]*(v[1046]*v[838]+v[1047]*v[841]+v[1048]*v[844]+v[1040]*v[852]
      +v[1041]*v[855]+v[1042]*v[858]);
     gpost[i777-1][10]=v[1049]*(v[1040]*v[909]+v[1041]*v[910]+v[1042]*v[911]);
     gpost[i777-1][11]=v[1049]*(v[1043]*v[909]+v[1044]*v[910]+v[1045]*v[911]);
     gpost[i777-1][12]=v[1049]*(v[1046]*v[909]+v[1047]*v[910]+v[1048]*v[911]);
    };/* end for */
};
}
#endif
