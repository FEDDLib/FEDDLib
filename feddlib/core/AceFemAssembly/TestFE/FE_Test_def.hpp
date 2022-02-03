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
	                                    MatrixPtr_Type &A,
	                                    bool callFillComplete,
	                                    int FELocExternal){

	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblem.xml");
	
    UN FEloc = checkFE(dim,FEType);

	ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

	int dofsElement = elements->getElement(0).getVectorNodeList().size();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

	MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

	vec2D_dbl_Type nodes;

	//SmallMatrixPtr_Type elementMatrix =Teuchos::rcp( new SmallMatrix_Type( dofsElement));

	for (UN T=0; T<elements->numberElements(); T++) {
		
		nodes = getCoordinates(elements->getElement(T).getVectorNodeList(), pointsRep);

		AssembleFEFactory<SC,LO,GO,NO> assembleFEFactory;// = new AssembleFEFactory<SC,LO,GO,NO>();

		AssembleFEPtr_Type assemblyFE = assembleFEFactory.build("Laplace",elements->getElement(T).getFlag(),nodes, params);

		SmallMatrixPtr_Type elementMatrix = assemblyFE->assembleJacobian();

		addFeMatrix(A,elementMatrix, elements->getElement(T), map);
		
	}
	if (callFillComplete)
	    A->fillComplete();

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

	ParameterListPtr_Type params = Teuchos::getParametersFromXmlFile("parametersProblem.xml");

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

    for (UN T=0; T<elements->numberElements(); T++) {

		nodes = getCoordinates(elements->getElement(T).getVectorNodeList(), pointsRep);

		AssembleFEFactory<SC,LO,GO,NO> assembleFEFactory;// = new AssembleFEFactory<SC,LO,GO,NO>();

		AssembleFEPtr_Type assemblyFE = assembleFEFactory.build("Laplace",elements->getElement(T).getFlag(),nodes, params);

		assemblyFE->addRHSFunc(func);

		vec_dbl_Type elementVec = assemblyFE->assembleRHS();

		addFeVector(a, elementVec, elements->getElement(T)); // if they are both multivectors its actually super simple! Import entries and add
	}
}


/*!

 \brief Inserting element stiffness matrices into global stiffness matrix

@param[in] &A Global Matrix
@param[in] elementMatrix Stiffness matrix of one element
@param[in] element Corresponding finite element
@param[in] map Map that corresponds to repeated nodes

*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::addFeMatrix(MatrixPtr_Type &A, SmallMatrixPtr_Type elementMatrix, FiniteElement element, MapConstPtr_Type map){

		int numNodes = elementMatrix->size();
		Teuchos::Array<SC> value( numNodes, 0. );
        Teuchos::Array<GO> columnIndices( numNodes, 0 );

		//Teuchos::ArrayView<const LO> indices;
		//Teuchos::ArrayView<const SC> values;

		for (UN i=0; i < numNodes ; i++) {
			GO row = map->getGlobalElement( element.getNode(i) );
			//A->getGlobalRowView(element.getNode(i), indices,valuesRow);
			for(UN j=0; j < numNodes ; j++){
			    
			    columnIndices[j] = map->getGlobalElement( element.getNode(j) );
				
				value[j] = (*elementMatrix)[i][j];
			    				    
			}

		  	A->insertGlobalValues( row, columnIndices(), value() ); // Automatically adds entries if a value already exists
		}
}

/*!

 \brief Inserting element rhs vector into global rhs vector

@param[in] &a Global Matrix
@param[in] elementVector Stiffness matrix of one element
@param[in] element Corresponding finite element

*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::addFeVector(MultiVectorPtr_Type &a, vec_dbl_Type elementVector, FiniteElement element){

        Teuchos::ArrayRCP<SC>  globalVec = a->getDataNonConst(0);
        //Teuchos::ArrayRCP<SC>  elementVec = elementVector->getDataNonConst(0);
		int numNodes = elementVector.size();

		for(UN j=0; j < element.getVectorNodeList().size() ; j++){
			globalVec[element.getNode(j)] += elementVector[j];
		    				    
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

};
#endif







