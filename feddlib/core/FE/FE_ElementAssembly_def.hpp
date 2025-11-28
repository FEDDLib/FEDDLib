#ifndef FE_ELEMENTASSEMBLY_DEF_hpp
#define FE_ELEMENTASSEMBLY_DEF_hpp

#include <string>
#include "feddlib/core/core_config.h"
#ifdef FEDD_HAVE_ACEGENINTERFACE
#include <aceinterface.hpp>
#endif

/*!
 Definition of FE_ElementAssembly

 @brief  FE_ElementAssembly
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */



using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::outArg;

namespace FEDD {


template <class SC, class LO, class GO, class NO>
FE_ElementAssembly<SC,LO,GO,NO>::FE_ElementAssembly(bool saveAssembly):
domainVec_(0),
setZeros_(false),
myeps_(),
saveAssembly_(saveAssembly)
{
}

template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::addFE(DomainConstPtr_Type domain){
    
    if (saveAssembly_){
        DomainPtr_Type domainNC = Teuchos::rcp_const_cast<Domain_Type>( domain );
        domainNC->initializeFEData();
    }
    domainVec_.push_back(domain);

}

template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::doSetZeros(double eps){

    setZeros_ = true;
    myeps_ = eps;

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
void FE_ElementAssembly<SC,LO,GO,NO>::globalAssembly(std::string ProblemType,
                        int dim,                    
                        int degree,                       
                        MultiVectorPtr_Type sol_rep,
                        BlockMatrixPtr_Type &A,
                        BlockMultiVectorPtr_Type &resVec,
                        ParameterListPtr_Type params,
                        std::string assembleMode,
                        bool callFillComplete,
                        int FELocExternal){

    int problemSize = domainVec_.size(); // Problem size, as in number of subproblems

    // Depending on problem size we extract necessary information 

	/// Tupel construction follows follwing pattern:
	/// std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
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
void FE_ElementAssembly<SC,LO,GO,NO>::addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement_vec element, tuple_disk_vec_ptr_Type problemDisk){
		
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
void FE_ElementAssembly<SC,LO,GO,NO>::addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec,tuple_disk_vec_ptr_Type problemDisk, FiniteElement_vec element)
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
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyLinearElasticity(int dim,
	                                    std::string FEType,
	                                    int degree,
										int dofs,
										MultiVectorPtr_Type d_rep,
	                                    BlockMatrixPtr_Type &A,
										BlockMultiVectorPtr_Type &resVec,
 										ParameterListPtr_Type params,
 										bool reAssemble,
 										std::string assembleMode,
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
	/// std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
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
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyNonLinearElasticity(int dim,
	                                    std::string FEType,
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
	/// std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numNodes=6;
	if(dim==3){
		numNodes=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type displacement ("Displacement",FEType,dofs,numNodes);
	problemDisk->push_back(displacement);

    int neoHookeNum = params->sublist("Parameter").get("Neo-Hooke Modell",1);

    std::string nonLinElasModell = "NonLinearElasticity2";
    if(neoHookeNum == 1)
        nonLinElasModell = "NonLinearElasticity";

    //std::cout << " ######## Assembly Modell: " << nonLinElasModell << " ############ " <<  std::endl;


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
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyNonLinearElasticity(int dim,
	                                    std::string FEType,
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
	/// std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
	int numNodes=6;
	if(dim==3){
		numNodes=10;
	}
	tuple_disk_vec_ptr_Type problemDisk = Teuchos::rcp(new tuple_disk_vec_Type(0));
	tuple_ssii_Type displacement ("Displacement",FEType,dofs,numNodes);
	problemDisk->push_back(displacement);

    int neoHookeNum = params->sublist("Parameter").get("Neo-Hooke Modell",1);

    std::string nonLinElasModell = "NonLinearElasticity2";
    if(neoHookeNum == 1)
        nonLinElasModell = "NonLinearElasticity";

    //std::cout << " ######## Assembly Modell: " << nonLinElasModell << " ############ " <<  std::endl;

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
void FE_ElementAssembly<SC,LO,GO,NO>::addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock, int dofs){

    Teuchos::ArrayRCP<SC>  resArray_block = res->getBlockNonConst(0)->getDataNonConst(0);

	vec_LO_Type nodeList_block = elementBlock.getVectorNodeList();

	for(int i=0; i< nodeList_block.size() ; i++){
		for(int d=0; d<dofs; d++)
			resArray_block[nodeList_block[i]*dofs+d] += (*rhsVec)[i*dofs+d];
	}
}

// Check the order of chemistry and solid in system matrix
template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyAceDeformDiffu(int dim,
                        std::string FETypeChem,
                        std::string FETypeSolid,
                        int degree,
                        int dofsChem,
                        int dofsSolid,
                        MultiVectorPtr_Type c_rep,
                        MultiVectorPtr_Type d_rep,
                        BlockMatrixPtr_Type &A,
                        BlockMultiVectorPtr_Type &resVec,
                        ParameterListPtr_Type params,
                        std::string assembleMode,
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
	/// std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
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

	std::string SCIModel = params->sublist("Parameter").get("Structure Model","SCI_simple");

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
        std::cout << " Determinante " << detB << std::endl;*/
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
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyAceDeformDiffuBlock(int dim,
                        std::string FETypeChem,
                        std::string FETypeSolid,
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
                        std::string assembleMode,
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
	/// std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
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
	
	std::string SCIModel = params->sublist("Parameter").get("Structure Model","SCI_simple");

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
void FE_ElementAssembly<SC,LO,GO,NO>::addFeBlockMatrix(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element1, FiniteElement element2, MapConstPtr_Type mapFirstRow,MapConstPtr_Type mapSecondRow, tuple_disk_vec_ptr_Type problemDisk){
		
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
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyNavierStokes(int dim,
	                                    std::string FETypeVelocity,
	                                    std::string FETypePressure,
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
 										std::string assembleMode,
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
	/// std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e. "P2"), int: Degrees of Freedom per Node, int: Number of Nodes per element)
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
        if(params->sublist("Material").get("Newtonian",true) == false)
	 	    initAssembleFEElements("GeneralizedNewtonian",problemDisk,elements, params,pointsRep,domainVec_.at(FElocVel)->getElementMap()); // In cas of non Newtonian Fluid
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
            
            //AssembleFENavierStokesPtr_Type elTmp = Teuchos::rcp_dynamic_cast<AssembleFENavierStokes_Type>(assemblyFEElements_[T] ); // Why should we need a pointer we can directly call it using assemblyFEElements_[T]? Because assemblyFixedPoint is not a function of base class

            if(params->sublist("Material").get("Newtonian",true) == false)
            {
                AssembleFEGeneralizedNewtonianPtr_Type elTmp = Teuchos::rcp_dynamic_cast<AssembleFEGeneralizedNewtonian_Type>( assemblyFEElements_[T] );
                elTmp->assembleFixedPoint();
                elementMatrix =  elTmp->getFixedPointMatrix(); 
            }
            else // Newtonian Case
            {
              AssembleFENavierStokesPtr_Type elTmp = Teuchos::rcp_dynamic_cast<AssembleFENavierStokes_Type>( assemblyFEElements_[T] );
              elTmp->assembleFixedPoint();
              elementMatrix =  elTmp->getFixedPointMatrix(); 
            }
            

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
 \brief  Method to loop over all assembleFESpecific elements and set the defined linearization 
@param[in] string linearization e.g. "Picard" or "Newton"
*/
template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::changeLinearizationFE(std::string linearization)
{
    for (UN T=0; T<assemblyFEElements_.size(); T++) // For each assembledFEElement change Linearization
    {	
        assemblyFEElements_[T]->changeLinearization(linearization);
    }
}


/*!
 \brief Postprocessing: Using a converged velocity solution -> compute averaged viscosity inside an element at center of mass
@param[in] dim Dimension
@param[in] FEType FE Discretization
@param[in] degree Degree of basis function
@param[in] repeated solution fields for u and p
@param[in] parameter lists
*/


template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::computeSteadyViscosityFE_CM(int dim,
	                                    std::string FETypeVelocity,
	                                    std::string FETypePressure,
										int dofsVelocity,
										int dofsPressure,
										MultiVectorPtr_Type u_rep,
										MultiVectorPtr_Type p_rep,
 										ParameterListPtr_Type params){
	

    UN FElocVel = checkFE(dim,FETypeVelocity); // Checks for different domains which belongs to a certain fetype
    UN FElocPres = checkFE(dim,FETypePressure); // Checks for different domains which belongs to a certain fetype

	ElementsPtr_Type elements = domainVec_.at(FElocVel)->getElementsC();
	ElementsPtr_Type elementsPres = domainVec_.at(FElocPres)->getElementsC();

	vec_dbl_Type solution(0);
	vec_dbl_Type solution_u;
	vec_dbl_Type solution_p;
    vec_dbl_Type solution_viscosity;

    // We have to compute viscosity solution in each element 
	MultiVectorPtr_Type Sol_viscosity = Teuchos::rcp( new MultiVector_Type( domainVec_.at(FElocVel)->getElementMap(), 1 ) ); //
    BlockMultiVectorPtr_Type visco_output = Teuchos::rcp( new BlockMultiVector_Type(1) );
    visco_output->addBlock(Sol_viscosity,0);
   

	for (UN T=0; T<assemblyFEElements_.size(); T++) {
       
		vec_dbl_Type solution(0);
		solution_u = getSolution(elements->getElement(T).getVectorNodeList(), u_rep,dofsVelocity); // get the solution inside an element on the nodes
		solution_p = getSolution(elementsPres->getElement(T).getVectorNodeList(), p_rep,dofsPressure);

		solution.insert( solution.end(), solution_u.begin(), solution_u.end() ); // here we insert the solution
		solution.insert( solution.end(), solution_p.begin(), solution_p.end() );

		assemblyFEElements_[T]->updateSolution(solution); // here we update the value of the solutions inside an element
 
        assemblyFEElements_[T]->computeLocalconstOutputField(); //  we compute the viscosity inside an element
        solution_viscosity = assemblyFEElements_[T]->getLocalconstOutputField();

        Teuchos::ArrayRCP<SC>  resArray_block = visco_output->getBlockNonConst(0)->getDataNonConst(0); // First 
        resArray_block[T] = solution_viscosity[0]; // although it is a vector it only has one entry because we compute the value in the center of the element
          
	} // end loop over all elements
    // We could also instead of just overwrite it add an aditional block such that we could also compute other output fields and save it in there
    this->const_output_fields= visco_output;


}



/*!

 \brief Inserting local rhsVec into global residual Mv;


@param[in] res BlockMultiVector of residual vec; Repeated distribution; 2 blocks
@param[in] rhsVec sorted the same way as residual vec
@param[in] element of block1

*/

template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::addFeBlockMv(BlockMultiVectorPtr_Type &res, vec_dbl_ptr_Type rhsVec, FiniteElement elementBlock1,FiniteElement elementBlock2, int dofs1, int dofs2 ){

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
TODO: column indices pre determine

@param[in] &A Global Block Matrix
@param[in] elementMatrix Stiffness matrix of one element
@param[in] element Corresponding finite element
@param[in] map Map that corresponds to repeated nodes of first block
@param[in] map Map that corresponds to repeated nodes of second block

*/
template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::addFeBlock(BlockMatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type mapRow, int row, int column, tuple_disk_vec_ptr_Type problemDisk){
		
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
void FE_ElementAssembly<SC,LO,GO,NO>::initAssembleFEElements(std::string elementType,tuple_disk_vec_ptr_Type problemDisk,ElementsPtr_Type elements, ParameterListPtr_Type params,vec2D_dbl_ptr_Type pointsRep, MapConstPtr_Type elementMap){
    
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
vec2D_dbl_Type FE_ElementAssembly<SC,LO,GO,NO>::getCoordinates(vec_LO_Type localIDs, vec2D_dbl_ptr_Type points){

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
vec_dbl_Type FE_ElementAssembly<SC,LO,GO,NO>::getSolution(vec_LO_Type localIDs, MultiVectorPtr_Type u_rep, int dofsVelocity){

    Teuchos::ArrayRCP<SC>  uArray = u_rep->getDataNonConst(0);
	
	vec_dbl_Type solution(0);
	for(int i=0; i < localIDs.size() ; i++){
		for(int d=0; d<dofsVelocity; d++)
			solution.push_back(uArray[localIDs[i]*dofsVelocity+d]);
	}

    return solution;
}


    
template <class SC, class LO, class GO, class NO>
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyEmptyMatrix(MatrixPtr_Type &A){
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
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyLaplaceAssFE(int dim,
                                        std::string FEType,
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
void FE_ElementAssembly<SC,LO,GO,NO>::assemblyNonlinearLaplace(
    int dim, std::string FEType, int degree, MultiVectorPtr_Type u_rep,
    BlockMatrixPtr_Type &A, BlockMultiVectorPtr_Type &resVec,
    ParameterListPtr_Type params, std::string assembleMode, bool callFillComplete,
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
    // std::string: Physical Entity (i.e. Velocity) , std::string: Discretisation (i.e.
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


/*!
 \brief Checks which domain corresponds to certain FE Type and dimension
 */
template <class SC, class LO, class GO, class NO>
int FE_ElementAssembly<SC,LO,GO,NO>::checkFE(int dim,
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

}
#endif
