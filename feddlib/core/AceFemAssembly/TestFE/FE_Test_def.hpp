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
	
	ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

	vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

	MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

	SmallMatrixPtr_Type elementMatrix(elements->getElement(0).size());

	for (UN T=0; T<elements->numberElements(); T++) {
		
		AssemblyFE_Type assemblyFe = new AssemblyFEFactory("Laplace",elements->getElement(T).getFlag(),elements->getElement(T).getVectorNodeList(), params);

		assemblyFE.assemblyLaplace(elementMatrix);

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
    ElementsPtr_Type elements = domainVec_.at(FEloc)->getElementsC();

    vec2D_dbl_ptr_Type pointsRep = domainVec_.at(FEloc)->getPointsRepeated();

    MapConstPtr_Type map = domainVec_.at(FEloc)->getMapRepeated();

    Teuchos::ArrayRCP< SC > valuesRhs = a->getDataNonConst(0);
    int parameters;
    double x;

    std::vector<double> valueFunc(dim);
    SC* paras = &(funcParameter[0]);
    
    func( &x, &valueFunc[0], paras );
    SC value;

    for (UN T=0; T<elements->numberElements(); T++) {

		AssemblyFE_Type assemblyFe = new AssemblyFEFactory("Laplace",elements->getElement(T).getFlag(),elements->getElement(T).getVectorNodeList(), params);

		assemblyFE.assembleRhs(elementVec,func, funcParameter);

		addFeVec(valuesRhs,elementVec, elements->getElement(T), map);
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
void FE_Test<SC,LO,GO,NO>::addFeMatrix(SmallMatrixPtr_Type &A, SmallMatrix<SC> elementMatrix, FiniteElement element, MapConstPtr_Type map){

		int numNodes = elementMatrix[0].size();
		Teuchos::Array<SC> value( numNodes, 0. );
        Teuchos::Array<GO> columnIndices( numNodes.size(), 0 );

		for (UN i=0; i < numNodes ; i++) {
			for(UN j=0; j < numNodes ; j++){
			    
			    columnIndices[j] = map->getGlobalElement( element.getNode(j) );
				value[j] = (*elementMatrix)[i][j]
			    				    
			}
			GO row = map->getGlobalElement( element.getNode(i) );
		    A->addGlobalValues( row, columnIndices(), value() ); // Check if this functions (addGlobalValues) exists for Matrix

		}

}

/*!

 \brief Inserting element rhs vector into global rhs vector

@param[in] &a Global Matrix
@param[in] elementVector Stiffness matrix of one element
@param[in] element Corresponding finite element

*/
template <class SC, class LO, class GO, class NO>
void FE_Test<SC,LO,GO,NO>::addFeVector(VectorPtr_Type &a, vec_dbl_Type elementVector, FiniteElement element){

		int numNodes = elementVector.size();
		
		for(UN j=0; j < numNodes ; j++){
		    
			a[element.getNode(j)] += elementVector[j];
		    				    
		}

}

};
#endif







