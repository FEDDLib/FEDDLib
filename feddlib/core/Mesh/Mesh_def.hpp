#ifndef Mesh_def_hpp
#define Mesh_def_hpp

#include "Mesh_decl.hpp"

/*!
Definition of Mesh

@brief  Mesh
@author Christian Hochmuth
@version 1.0
@copyright CH
*/
using Teuchos::reduceAll;
using Teuchos::REDUCE_SUM;
using Teuchos::outArg;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;


using namespace std;
namespace FEDD {
template <class SC, class LO, class GO, class NO>
Mesh<SC,LO,GO,NO>::Mesh():
numElementsGlob_(0),
mapUnique_(),
mapRepeated_(),
pointsRep_(),
pointsUni_(),
bcFlagRep_(),
bcFlagUni_(),
surfaceElements_(),
elementMap_(),
comm_(),
edgeElements_(),
pointsRepRef_(),
pointsUniRef_(),
mapUniqueP2Map_(),
mapRepeatedP2Map_(),
AABBTree_()
{

    surfaceElements_.reset(new Elements());
    
    elementsC_.reset(new Elements());    

    edgeElements_ = Teuchos::rcp( new EdgeElements_Type() );
    
    FEType_ = "P1"; // We generally assume the mesh to be p1. In case of P1 or Q2 the FEType is allways adjusted

    volumeID_=0;
}

template <class SC, class LO, class GO, class NO>
Mesh<SC,LO,GO,NO>::Mesh(CommConstPtrConst_Type& comm):
numElementsGlob_(0),
mapUnique_(),
mapRepeated_(),
pointsRep_(),
pointsUni_(),
bcFlagRep_(),
bcFlagUni_(),
surfaceElements_(),
elementMap_(),
edgeMap_(),
comm_(comm),
pointsRepRef_(),
pointsUniRef_(),
mapUniqueP2Map_(),
volumeID_(0),
mapRepeatedP2Map_(),
AABBTree_()
{
    AABBTree_.reset(new AABBTree_Type());
    surfaceElements_.reset(new Elements());

    elementsC_.reset(new Elements());

    edgeElements_ = Teuchos::rcp( new EdgeElements_Type() );

    FEType_ = "P1";

    volumeID_=0;

}

template <class SC, class LO, class GO, class NO>
Mesh<SC,LO,GO,NO>::~Mesh(){

}

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::setElementFlags(std::string type){

    ElementsPtr_Type elements = this->getElementsC();
//    this->elementFlag_.reset( new vec_int_Type( elements->numberElements(), 0 ) );
    if (type == "TPM_square") {
        double xRef, yRef;

        for (int i=0; i<elements->numberElements(); i++) {
            xRef = ( this->pointsRep_->at( elements->getElement(i).getNode(0) )[0] + this->pointsRep_->at( elements->getElement(i).getNode(1) )[0] + this->pointsRep_->at( elements->getElement(i).getNode(2) )[0] ) / 3.;
            yRef = ( this->pointsRep_->at( elements->getElement(i).getNode(0) )[1] + this->pointsRep_->at( elements->getElement(i).getNode(1) )[1] + this->pointsRep_->at( elements->getElement(i).getNode(2) )[1] ) / 3.;
            if ( xRef>=0.3  && xRef<=0.7) {
                if ( yRef>= 0.6) {
                    elements->getElement(i).setFlag(1);
//                    this->elementFlag_->at(i) = 1;
                }
            }
        }
    }
    else if (type == "Excavation1"){
    }
    else{
    
    }

}

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::setParameterList( ParameterListPtr_Type& pL ) {
    pList_ = pL;
}

template <class SC, class LO, class GO, class NO>
ParameterListConstPtr_Type Mesh<SC,LO,GO,NO>::getParameterList( ) const{
    return pList_;
}
    
template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type Mesh<SC,LO,GO,NO>::getElementsFlag() const{
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "we are not using the correct flags here. use the flags of elementC_." );
    vec_int_ptr_Type tmp;
    return tmp;
}

template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::MapConstPtr_Type Mesh<SC,LO,GO,NO>::getMapUnique() const{

    return mapUnique_;
}

template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::MapConstPtr_Type Mesh<SC,LO,GO,NO>::getMapRepeated() const{

    return mapRepeated_;
}

template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::MapConstPtr_Type Mesh<SC,LO,GO,NO>::getMapUniqueP2() const{

    return mapUniqueP2Map_;
}

template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::MapConstPtr_Type Mesh<SC,LO,GO,NO>::getMapRepeatedP2() const{

    return mapRepeatedP2Map_;
}

template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::MapConstPtr_Type Mesh<SC,LO,GO,NO>::getElementMap(){
    TEUCHOS_TEST_FOR_EXCEPTION( elementMap_.is_null(), std::runtime_error, "Element map of mesh does not exist." );
    return elementMap_;
}


// edgeMap
template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::MapConstPtr_Type Mesh<SC,LO,GO,NO>::getEdgeMap(){
    TEUCHOS_TEST_FOR_EXCEPTION( edgeMap_.is_null(), std::runtime_error, "Element map of mesh does not exist." );
    return edgeMap_;
}

template <class SC, class LO, class GO, class NO>
vec2D_dbl_ptr_Type Mesh<SC,LO,GO,NO>::getPointsRepeated() const{

    return pointsRep_;
}

template <class SC, class LO, class GO, class NO>
vec2D_dbl_ptr_Type Mesh<SC,LO,GO,NO>::getPointsUnique() const{

    return pointsUni_;
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type Mesh<SC,LO,GO,NO>::getBCFlagRepeated() const{

    return bcFlagRep_;
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type Mesh<SC,LO,GO,NO>::getBCFlagUnique() const{

    return bcFlagUni_;
}

template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::ElementsPtr_Type Mesh<SC,LO,GO,NO>::getElementsC(){
    return elementsC_;
}

template <class SC, class LO, class GO, class NO>
typename Mesh<SC,LO,GO,NO>::ElementsPtr_Type Mesh<SC,LO,GO,NO>::getSurfaceElements(){
    return surfaceElements_;
}

template <class SC, class LO, class GO, class NO>
int Mesh<SC,LO,GO,NO>::getDimension(){

    return dim_;
}

template <class SC, class LO, class GO, class NO>
GO Mesh<SC,LO,GO,NO>::getNumElementsGlobal(){

    return numElementsGlob_;
}

template <class SC, class LO, class GO, class NO>
LO Mesh<SC,LO,GO,NO>::getNumElements(){
    TEUCHOS_TEST_FOR_EXCEPTION( this->elementsC_.is_null(), std::runtime_error ,"Elements do not exist." );
    return this->elementsC_->numberElements();
}

template <class SC, class LO, class GO, class NO>
LO Mesh<SC,LO,GO,NO>::getNumPoints(std::string type){
    if (!type.compare("Unique"))
    return pointsUni_->size();
    else if(!type.compare("Repeated"))
    return pointsRep_->size();

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Select valid map type: unique or repeated.");
    return 0;
}

template <class SC, class LO, class GO, class NO>
int Mesh<SC,LO,GO,NO>::getOrderElement(){

    switch (dim_) {
        case 2:
            if ( !FEType_.compare("P1") )
                return 3;
            else if ( !FEType_.compare("P1-disc") || !FEType_.compare("P1-disc-global") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "P1-disc only available in 3D.");
            }
            else if( !FEType_.compare("P2") )
                return 6;
            else if( !FEType_.compare("P2-CR") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "P2-CR only available in 3D.");
            }
            else if( !FEType_.compare("Q2-20") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Q2-20 only available in 3D.");
            }
            else if( !FEType_.compare("Q2") ){
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Q2 only available in 3D.");
            }
            break;
        case 3:
            if ( !FEType_.compare("P1") )
                return 4;
            if ( !FEType_.compare("P1-disc") || !FEType_.compare("P1-disc-global") )
                return 4;
            else if( !FEType_.compare("P2") )
                return 10;
            else if( !FEType_.compare("P2-CR") )
                return 15;
            else if( !FEType_.compare("Q2-20") )
                return 20;
            else if( !FEType_.compare("Q2") )
                return 27;
            break;
        default:
            return -1;
            break;
    }
    return -1;
}

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::setReferenceConfiguration()
{
    // Bemerkung: Repeated und Unique sind unterschiedlich lang!!! => zwei Schleifen

    // Setze zunaechst alles auf Null, andernfalls kann man nicht drauf zugreifen
    //    vec2D_dbl_ptr_Type zeroRep(new vec2D_dbl_Type(pointsRep_->size(),vec_dbl_Type(pointsRep_->at(0).size(),0.0)));
    //    vec2D_dbl_ptr_Type zeroUni(new vec2D_dbl_Type(pointsUni_->size(),vec_dbl_Type(pointsUni_->at(0).size(),0.0)));


    pointsRepRef_.reset( new vec2D_dbl_Type() );
    pointsUniRef_.reset( new vec2D_dbl_Type() );
    // zeroRep und zeroUni leben nur hier drinnen, weswegen wir Pointer gleich Pointer setzen koennen.
    // TODO: *PointsRepRef_ = *zeroRep funktioniert nicht.
    *pointsRepRef_ = *pointsRep_;
    *pointsUniRef_ = *pointsUni_;

    //    // Repeated
    //    for(int i = 0; i < pointsRep_->size(); i++)
    //    {
    //        for(int j = 0; j < pointsRep_->at(0).size(); j++)
    //        {
    //            pointsRepRef_->at(i).at(j) = PointsRep_->at(i).at(j);
    //        }
    //    }
    //
    //    // Unique
    //    for(int i = 0; i < PointsUni_->size(); i++)
    //    {
    //        for(int j = 0; j < PointsUni_->at(0).size(); j++)
    //        {
    //            PointsUniRef_->at(i).at(j) = PointsUni_->at(i).at(j);
    //        }
    //    }
}

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::moveMesh( MultiVectorPtr_Type displacementUnique, MultiVectorPtr_Type displacementRepeated )
{
    // Bemerkung: Repeated und Unique sind unterschiedlich lang!!! => zwei Schleifen
    TEUCHOS_TEST_FOR_EXCEPTION (displacementRepeated.is_null(), std::runtime_error," displacementRepeated in moveMesh is null.")
    TEUCHOS_TEST_FOR_EXCEPTION (displacementUnique.is_null(), std::runtime_error," displacementRepeated in moveMesh is null.")
    // Repeated
    Teuchos::ArrayRCP<const SC> values = displacementRepeated->getData(0); //only 1 MV
    for(int i = 0; i < pointsRepRef_->size(); i++)
    {
        for(int j = 0; j < pointsRepRef_->at(0).size(); j++)
        {
            // Sortierung von DisplacementRepeated ist x-y-x-y-x-y-x-y bzw. x-y-z-x-y-z-x-y-z
            // Achtung: DisplacementRepeated ist ein Pointer der mit (*) dereferenziert werden muss.
            // Operator[] kann nicht auf einen Pointer angewendet werden!!!
            // Es sei denn es ist ein Array.
            pointsRep_->at(i).at(j) = pointsRepRef_->at(i).at(j) + values[dim_*i+j];
        }
    }

    // Unique
    values = displacementUnique->getData(0); //only 1 MV
    for(int i = 0; i < pointsUniRef_->size(); i++)
    {
        for(int j = 0; j < pointsUniRef_->at(0).size(); j++)
        {
            // Sortierung von DisplacementRepeated ist x-y-x-y-x-y-x-y bzw. x-y-z-x-y-z-x-y-z
            // Erklaerung: DisplacementUnique ist ein Vector-Pointer, wo in jedem Eintrag ein MultiVector-Pointer drin steht (std::vector<MultiVector_ptr_Type>)
            // Greife mit ->at auf den Eintrag des Vektors zu (hier nur ein Eintrag vorhanden), dereferenziere den damit erhaltenen MultiVector-Pointer (als Referenz) um einen
            // MultiVector zu erhalten.
            // Greife dann mit [] auf das entsprechende Array (double *&) im MultiVector zu (hier gibt es nur einen)
            // und anschliessend mit [] auf den Wert des Arrays.
            // Beachte falls x ein Array ist (also z.B. double *), dann ist x[i] := *(x+i)!!!
            // Liefert also direkt den Wert und keinen Pointer auf einen double.
            // Achtung: MultiVector[] liefert double* wohingegen MultiVector() Epetra_Vector* zurueck liefert
            pointsUni_->at(i).at(j) = pointsUniRef_->at(i).at(j) + values[dim_*i+j];
        }
    }
}

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::create_AABBTree(){
    if (AABBTree_.is_null()){
        AABBTree_.reset(new AABBTree_Type());
    }
    AABBTree_->createTreeFromElements(
        getElementsC(),
        getPointsRepeated()
    );
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type Mesh<SC,LO,GO,NO>::findElemsForPoints(
    vec2D_dbl_ptr_Type queryPoints
){
    int numPoints = queryPoints->size();

    // Return vector. -1 means that point is in no elem, otherwise entry is the
    // elem the point is in
    vec_int_ptr_Type pointToElem(
        new vec_int_Type(
            numPoints,
            -1
        )
    );

    // Create tree if it is empty
    if (AABBTree_->isEmpty()){
        AABBTree_->createTreeFromElements(
            getElementsC(),
            getPointsRepeated(),
            false
        );
    }

    // Query the AABBTree
    map<int, list<int> > treeToItem;
    map<int, list<int> > itemToTree;
    tie(treeToItem, itemToTree) = AABBTree_->scanTree(queryPoints, false);

    // FIXME: put this in a function of AABBTree?
    // unnest the returned answer for each query_point
    int point = -1;
    bool found = false;
    list<int> rectangles;
    list<int> elements;
    for (auto keyValue: itemToTree){
        // FIXME: put this in a function of AABBTree?
        // rectangles is a list<int> of all rectangles point is in
        // find the element(s) that is/are in all of said rectangles.
        // If there is only one element that is the element the point is in,
        // if not we have to query all remaining elements
        point = keyValue.first;
        rectangles = keyValue.second;


        // query all remaining elements
        for(auto rectangle: rectangles){
            elements = AABBTree_->getElements(rectangle);
            for(auto element: elements){
                found = isPointInElem(queryPoints->at(point), element);
                if( found ){
                    pointToElem->at(point) = element;
                    break;
                }
            }
            if( found ){
                // we already found the element, no need to check additional rectangles
                break;
            }
        }

    }
    return pointToElem;
}

template <class SC, class LO, class GO, class NO>
vec_dbl_Type Mesh<SC,LO,GO,NO>::getBaryCoords(vec_dbl_Type point, int element){
    vec_int_Type localNodes = elementsC_->getElement(element).getVectorNodeList();
    vec2D_dbl_Type coords(
        localNodes.size(),
        vec_dbl_Type(2, 0.0) // FIXME: this depends on the dimension
    );
    int node = 0;
    for (int localNode=0; localNode<localNodes.size(); localNode++){
        node = localNodes.at(localNode); // get global node_id
        coords.at(localNode) = pointsRep_->at(node); //get global coordinates
    }

    double px, py, x1, x2, x3, y1, y2, y3;
    px = point.at(0); py = point.at(1);
    x1 = coords.at(0).at(0); y1 = coords.at(0).at(1);
    x2 = coords.at(1).at(0); y2 = coords.at(1).at(1);
    x3 = coords.at(2).at(0); y3 = coords.at(2).at(1);

    // baryzentric coordinates
    double det_T = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);

    vec_dbl_Type baryCoords(
        3,
        0.0
    );
    baryCoords[0] =(y2 - y3) * (px - x3) + (x3 - x2) * (py - y3);
    baryCoords[0] = baryCoords[0] / det_T;

    baryCoords[1] = (y3 -y1) * (px - x3) + (x1 - x3) * (py - y3);
    baryCoords[1] = baryCoords[1] / det_T;

    baryCoords[2] = 1 - baryCoords[1] - baryCoords[0];
    return baryCoords;
}


template <class SC, class LO, class GO, class NO>
bool Mesh<SC,LO,GO,NO>:: isPointInElem(vec_dbl_Type point, int element){
    // FIXME: This is only valid for a triangle
    vec_dbl_Type baryCoords;
    baryCoords = getBaryCoords(point, element);

    if(baryCoords[0] >= 0 && baryCoords[1] >= 0 && baryCoords[2] >= 0){
        return true;
    }
    return false;
}

template <class SC, class LO, class GO, class NO>
vec2D_int_ptr_Type Mesh<SC,LO,GO,NO>::getElements(){
    this->elementsVec_ = Teuchos::rcp( new vec2D_int_Type( this->elementsC_->numberElements() ) );
    for (int i=0; i<this->elementsVec_->size(); i++)
        this->elementsVec_->at(i) = this->elementsC_->getElement(i).getVectorNodeList();

    return this->elementsVec_;
}

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::correctNormalDirections(){

    int outwardNormals = 0;
    int inwardNormals = 0;
    for (UN T=0; T<elementsC_->numberElements(); T++) {
        FiniteElement fe = elementsC_->getElement( T );
        ElementsPtr_Type subEl = fe.getSubElements(); // might be null
        for (int surface=0; surface<fe.numSubElements(); surface++) {
            FiniteElement feSub = subEl->getElement( surface  );
            vec_int_Type nodeListElement = fe.getVectorNodeList();
            if(subEl->getDimension() == dim_-1 ){
                vec_int_Type nodeList = feSub.getVectorNodeListNonConst();
                int numNodes_T = nodeList.size();
                
                vec_dbl_Type v_E(dim_,1.);
                double norm_v_E=1.;
                LO id0 = nodeList[0];

                Helper::computeSurfaceNormal(dim_, pointsRep_,nodeList,v_E,norm_v_E);

                std::sort(nodeList.begin(), nodeList.end());
                std::sort(nodeListElement.begin(), nodeListElement.end());

                std::vector<int> v_symDifference;

                std::set_symmetric_difference(
                    nodeList.begin(), nodeList.end(),
                    nodeListElement.begin(), nodeListElement.end(),
                    std::back_inserter(v_symDifference));

                LO id1 = v_symDifference[0]; // This is a node that is not part of the surface, i.e. 4th element node in 3D

                vec_dbl_Type p0(dim_,0.);
                for(int i=0; i< dim_; i++)
                    p0[i] = pointsRep_->at(id1)[i] - pointsRep_->at(id0)[i];

                double sum = 0.;
                for(int i=0; i< dim_; i++)
                    sum += p0[i] * v_E[i];
                
                if(sum<=0){
                    outwardNormals++;
                }
                if(sum>0){
                    inwardNormals++;
                }
                if(sum>0) // if the sum is greater than 0, the normal is in inward direction and thus is flipped.
                    flipSurface(subEl,surface);
                    

            }
        }
    }
    reduceAll<int, int> (*this->getComm(), REDUCE_SUM, inwardNormals, outArg (inwardNormals));
    reduceAll<int, int> (*this->getComm(), REDUCE_SUM, outwardNormals, outArg (outwardNormals));

    if(this->getComm()->getRank() == 0){
        cout << " ############################################ " << endl;
        cout << " Mesh Orientation Statistic " << endl;
        cout << " Number of outward normals " << outwardNormals << endl;
        cout << " Number of inward normals " << inwardNormals << endl;
        cout << " ############################################ " << endl;
    }

}
template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::correctElementOrientation(){

    int posDet=0;
    int negDet=0;

    SC detB;
    SmallMatrix<SC> B(dim_);
    SmallMatrix<SC> Binv(dim_);

    for (UN T=0; T<elementsC_->numberElements(); T++) {
        Helper::buildTransformation(elementsC_->getElement(T).getVectorNodeList(), pointsRep_, B);
        detB = B.computeInverse(Binv);

        if(detB<0){
            negDet++;
        }
        if(detB>0){
            posDet++;
        }
        if(detB<0) // If the determinant is smaller than zero we flip the element
            flipElement(elementsC_,T); 

    }
    cout << " Finished " << endl;
    reduceAll<int, int> (*this->getComm(), REDUCE_SUM, negDet, outArg (negDet));
    reduceAll<int, int> (*this->getComm(), REDUCE_SUM, posDet, outArg (posDet));

    if(this->getComm()->getRank() == 0){
        cout << " ############################################ " << endl;
        cout << " Mesh Orientation Statistic " << endl;
        cout << " Number of positive dets " << posDet << endl;
        cout << " Number of negative dets " << negDet << endl;
        cout << " ############################################ " << endl;
    }

}

// We allways want a outward normal direction
// Assumptions: We are flipping a surface which is a subelement. Subelement are generally element that are on the boundary layers of the domain. Thus, they are unique. (There are no two identical triangles in two elements as subelements)
// Question: Easiest way to flip the surface without redoing whole dim-element (surface being dim-1-element)
template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::flipSurface(ElementsPtr_Type subEl, int surfaceNumber){

    vec_LO_Type surfaceElements_vec = subEl->getElement(surfaceNumber).getVectorNodeList();

    if(dim_ == 2){
        if(FEType_ == "P1"){
            LO id1,id2;
            id1= surfaceElements_vec[0];
            id2= surfaceElements_vec[1];
           
            surfaceElements_vec[0] = id2;
            surfaceElements_vec[1] = id1;
        }
        else if(FEType_ == "P2"){
            LO id1,id2,id3;
            id1= surfaceElements_vec[0];
            id2= surfaceElements_vec[1];
            id3= surfaceElements_vec[2];

            surfaceElements_vec[0] = id2;
            surfaceElements_vec[1] = id1;
            surfaceElements_vec[2] = id3;
        }
    }
    else if(dim_ == 3){

        if(FEType_ == "P1"){
            LO id1,id2,id3,id4,id5,id6;
            id1= surfaceElements_vec[0];
            id2= surfaceElements_vec[1];
            id3= surfaceElements_vec[2];
           
            surfaceElements_vec[0] = id1;
            surfaceElements_vec[1] = id3;
            surfaceElements_vec[2] = id2;           
        }
        else if(FEType_ == "P2"){
            LO id1,id2,id3,id4,id5,id6;
            id1= surfaceElements_vec[0];
            id2= surfaceElements_vec[1];
            id3= surfaceElements_vec[2];
            id4= surfaceElements_vec[3];
            id5= surfaceElements_vec[4];
            id6= surfaceElements_vec[5];

            surfaceElements_vec[0] = id1;
            surfaceElements_vec[1] = id3;
            surfaceElements_vec[2] = id2;
            surfaceElements_vec[3] = id6;
            surfaceElements_vec[4] = id5;
            surfaceElements_vec[5] = id4;
        }
        else    
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "We can only flip normals for P1 or P2 elements. Invalid " << FEType_ << " " );
      
    }  
    FiniteElement feFlipped(surfaceElements_vec,subEl->getElement(surfaceNumber).getFlag());
    subEl->switchElement(surfaceNumber,feFlipped); // We can switch the current element with the newly defined element which has just a different node ordering. 
 

}

// We allways want a positive determinant
template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::flipElement(ElementsPtr_Type elements, int elementNumber){

    vec_LO_Type surfaceElements_vec = elements->getElement(elementNumber).getVectorNodeList();
    if(dim_ == 2){
        if(FEType_ == "P1"){
            LO id1,id2,id3;
            id1= surfaceElements_vec[0];
            id2= surfaceElements_vec[1];
            id3= surfaceElements_vec[2];

            surfaceElements_vec[0] = id1;
            surfaceElements_vec[1] = id3;
            surfaceElements_vec[2] = id2;

        }
        else if(FEType_ == "P2"){
            LO id1,id2,id3,id4,id5,id6;
            id1= surfaceElements_vec[0];
            id2= surfaceElements_vec[1];
            id3= surfaceElements_vec[2];
            id4= surfaceElements_vec[3];
            id5= surfaceElements_vec[4];
            id6= surfaceElements_vec[5];

            surfaceElements_vec[0] = id1;
            surfaceElements_vec[1] = id3;
            surfaceElements_vec[2] = id2;
            surfaceElements_vec[3] = id6;
            surfaceElements_vec[4] = id5;
            surfaceElements_vec[5] = id4;
        }
    }
    else if(dim_ == 3){

        if(FEType_ == "P1"){
            LO id1,id2,id3,id4;
            id1= surfaceElements_vec[0];
            id2= surfaceElements_vec[1];
            id3= surfaceElements_vec[2];
            id4= surfaceElements_vec[3];

            surfaceElements_vec[0] = id1;
            surfaceElements_vec[1] = id3;
            surfaceElements_vec[2] = id2;      
            surfaceElements_vec[3] = id4;           
     
        }
        else if(FEType_ == "P2"){
            LO id1,id2,id3,id4,id5,id6, id7,id8,id9,id0;
            id0= surfaceElements_vec[0];
            id1= surfaceElements_vec[1];
            id2= surfaceElements_vec[2];
            id3= surfaceElements_vec[3];
            id4= surfaceElements_vec[4];
            id5= surfaceElements_vec[5];
            id6= surfaceElements_vec[6];
            id7= surfaceElements_vec[7];
            id8= surfaceElements_vec[8];
            id9= surfaceElements_vec[9];

            surfaceElements_vec[0] = id0;
            surfaceElements_vec[1] = id2;
            surfaceElements_vec[2] = id1;
            surfaceElements_vec[3] = id3;

            surfaceElements_vec[4] = id6;
            surfaceElements_vec[5] = id5;
            surfaceElements_vec[6] = id4;
            surfaceElements_vec[7] = id7;
            surfaceElements_vec[8] = id9;
            surfaceElements_vec[9] = id8;
        }
        else    
            TEUCHOS_TEST_FOR_EXCEPTION( true, std::runtime_error, "We can only flip normals for P1 or P2 elements. Invalid " << FEType_ << " " );
      
    }  
    FiniteElement feFlipped(surfaceElements_vec,elements->getElement(elementNumber).getFlag()); 
    ElementsPtr_Type subElements = elements->getElement(elementNumber).getSubElements();
    for(int T = 0; T<elements->getElement(elementNumber).numSubElements(); T++){
        if(!feFlipped.subElementsInitialized())
            feFlipped.initializeSubElements( this->FEType_, this->dim_ -1) ;
        feFlipped.addSubElement(subElements->getElement(T));
    }   
    elements->switchElement(elementNumber,feFlipped); // We can switch the current element with the newly defined element which has just a different node ordering. 

}

/*!

 \brief Building edgeMap after refinement.

@param[in] mapGlobalProc Map of global processor numbers
@param[in] mapProc Map of local processor number
*/

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::buildEdgeMap(){

    // -----------------
    int maxRank = std::get<1>(this->rankRange_);
    const int myRank = this->comm_->getRank();
    // We need this to transfer some information
    vec_GO_Type globalProcs(0);
    for (int i=0; i<= maxRank; i++)
        globalProcs.push_back(i);

    Teuchos::ArrayView<GO> globalProcArray = Teuchos::arrayViewFromVector( globalProcs);

    vec_GO_Type localProc(0);
    localProc.push_back(this->comm_->getRank());
    Teuchos::ArrayView<GO> localProcArray = Teuchos::arrayViewFromVector( localProc);

    MapPtr_Type mapGlobalProc =
        Teuchos::RCP( new Map_Type( this->getElementMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), globalProcArray, 0, this->comm_) );

    MapPtr_Type mapProc =
        Teuchos::RCP( new Map_Type( this->getElementMap()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), localProcArray, 0, this->comm_) );

    // ------------------

    vec2D_int_Type interfaceEdgesLocalId(1,vec_int_Type(1));

    MultiVectorLOPtr_Type exportLocalEntry = Teuchos::rcp( new MultiVectorLO_Type( mapProc, 1 ) );

    // (A) First we determine a Map only for the interface Nodes
    // This will reduce the size of the Matrix we build later significantly if only look at the interface edges
    int numEdges= this->edgeElements_->numberElements();
    vec2D_GO_Type inzidenzIndices(0,vec_GO_Type(2)); // Vector that stores global IDs of each edge (in Repeated Sense)
    vec_LO_Type localEdgeIndex(0); // stores the local ID of edges in question 
    vec_GO_Type id(2);
    int edgesUnique=0;
    EdgeElementsPtr_Type edgeElements = this->edgeElements_; // Edges

    vec2D_dbl_ptr_Type points = this->pointsRep_;

    int interfaceNum=0;
    for(int i=0; i<numEdges; i++ ){
        //if(edgeElements->getElement(i).isInterfaceElement()){
            id[0] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(0)); 
            id[1] = this->mapRepeated_->getGlobalElement(edgeElements->getElement(i).getNode(1));
            
            sort(id.begin(),id.end());
            inzidenzIndices.push_back(id);

            localEdgeIndex.push_back(i);
            interfaceNum++;
        //}

        //else{
        //	edgesUnique++;
        //}


    }


    // This Matrix is row based, where the row is based on mapInterfaceNodesUnqiue
    // We then add a '1' Entry when two global Node IDs form an edge
    MatrixPtr_Type inzidenzMatrix = Teuchos::RCP( new Matrix_Type(this->mapUnique_, 40 ) );
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
    exportLocalEntry->putScalar(   edgesUnique );

    MultiVectorLOPtr_Type newEdgesUniqueGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
    newEdgesUniqueGlobal->putScalar(   0 ); 
    newEdgesUniqueGlobal->importFromVector( exportLocalEntry, true, "Insert");
    // offset EdgesUnique for proc and globally
    Teuchos::ArrayRCP< const SC > newEdgesList = newEdgesUniqueGlobal->getData(0);

    GO procOffsetEdges=0;
    for(int i=0; i< myRank; i++)
        procOffsetEdges= procOffsetEdges + newEdgesList[i];

    // global IDs for map
    vec_GO_Type vecGlobalIDsEdges(this->edgeElements_->numberElements()); 

    // Step 1: adding unique global edge IDs
    int count=0;
    for(int i=0; i< this->edgeElements_->numberElements(); i++){
        //if(!this->edgeElements_->getElement(i).isInterfaceElement()){
            vecGlobalIDsEdges.at(i) = procOffsetEdges+count;
            count++;
        //}
    }	
    
    // Now we add the repeated ids, by first turning interfaceEdgesTag into a map
    // Offset for interface IDS:

    GO offsetInterface=0;
    for(int i=0; i< maxRank+1; i++)
            offsetInterface=  offsetInterface + newEdgesList[i];
    
    //Now we count the row entries on each processor an set global IDs
    Teuchos::ArrayView<const LO> indices;
    Teuchos::ArrayView<const SC> values;
    vec2D_GO_Type inzidenzIndicesUnique(0,vec_GO_Type(2)); // Vector that stores only both global IDs if the first is part of my unique Interface Nodes
    MapConstPtr_Type colMap = inzidenzMatrix->getMap("col");
    MapConstPtr_Type rowMap = inzidenzMatrix->getMap("row");
    int numRows = rowMap->getNodeNumElements();
    int uniqueEdges =0;
    for(int i=0; i<numRows; i++ ){
        inzidenzMatrix->getLocalRowView(i, indices,values); 
        uniqueEdges = uniqueEdges+indices.size();
        vec_GO_Type edgeTmp(2);
        for(int j=0; j<indices.size(); j++){
            edgeTmp[0] = rowMap->getGlobalElement(i);
            edgeTmp[1] = colMap->getGlobalElement(indices[j]);
            inzidenzIndicesUnique.push_back(edgeTmp);
        }
    }

    exportLocalEntry->putScalar( uniqueEdges );
    MultiVectorLOPtr_Type newEdgesInterfaceGlobal= Teuchos::rcp( new MultiVectorLO_Type( mapGlobalProc, 1 ) );
    newEdgesInterfaceGlobal->putScalar(   0 ); 
    newEdgesInterfaceGlobal->importFromVector( exportLocalEntry, true, "Insert");

    // offset EdgesUnique for proc and globally
    Teuchos::ArrayRCP< const SC > numUniqueInterface = newEdgesInterfaceGlobal->getData(0);

    procOffsetEdges=0;
    for(int i=0; i< myRank; i++)
        procOffsetEdges= procOffsetEdges + numUniqueInterface[i];

    int numInterfaceEdges=0;
    
    vec_GO_Type uniqueInterfaceIDsList_(inzidenzIndicesUnique.size());
    for(int i=0; i< uniqueInterfaceIDsList_.size(); i++)
        uniqueInterfaceIDsList_[i] = procOffsetEdges + i;

    MatrixPtr_Type indMatrix = Teuchos::rcp( new Matrix_Type(this->mapUnique_, 40 ) );

    for(int i=0; i<inzidenzIndicesUnique.size(); i++ ){
        index[0] = inzidenzIndicesUnique[i][0];
        col[0] = inzidenzIndicesUnique[i][1];
        Teuchos::Array<SC> value2(1,uniqueInterfaceIDsList_[i]);
        indMatrix->insertGlobalValues(index[0], col(), value2());
    }
    indMatrix->fillComplete(); 
    MatrixPtr_Type importMatrix = Teuchos::rcp( new Matrix_Type(this->mapRepeated_, 40 ) );
    
    importMatrix->importFromVector(indMatrix,false,"Insert");
    importMatrix->fillComplete(); 		
    
    // Determine global indices
    GO edgeID=0;
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
                edgeID = indicesTmp[entry][1];
                vecGlobalIDsEdges.at(localEdgeIndex[i]) = offsetInterface + edgeID;
                found =true;
            }
        }
        if(found == false)
            cout << " Asking for row " << rowMap->getLocalElement(inzidenzIndices[i][0]) << " for Edge [" << inzidenzIndices[i][0] << ",  " << inzidenzIndices[i][1] << "], on Proc " << myRank << " but no Value found " <<endl;
        }


    Teuchos::RCP<std::vector<GO>> edgesGlobMapping = Teuchos::rcp( new vector<GO>( vecGlobalIDsEdges ) );
    Teuchos::ArrayView<GO> edgesGlobMappingArray = Teuchos::arrayViewFromVector( *edgesGlobMapping);

    // Based on repeated edge map we should be able to identify the interface edges! 

    this->edgeMap_.reset(new Map<LO,GO,NO>(this->getMapRepeated()->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesGlobMappingArray, 0, this->comm_) );
    //this->edgeMap_->print();
    MapConstPtr_Type edgeMapUnique = this->edgeMap_->buildUniqueMap( this->rankRange_ );

    MultiVectorPtr_Type repeatedEdges = Teuchos::rcp( new MultiVector_Type( this->edgeMap_, 1 ) );
    repeatedEdges->putScalar(  1);

    MultiVectorPtr_Type uniqueEdgesVec = Teuchos::rcp( new MultiVector_Type( edgeMapUnique, 1 ) );
    uniqueEdgesVec->putScalar(   0 ); 

    // Adding all edges up to determine multiplicity
    uniqueEdgesVec->exportFromVector( repeatedEdges, false, "Add");

    // Distributing it back to repeated edges
    repeatedEdges->putScalar(0); 
    repeatedEdges->importFromVector(uniqueEdgesVec,false,"Insert");
    
	Teuchos::ArrayRCP< SC > repeatedEdgesArray  = repeatedEdges->getDataNonConst(0);
    for(int E = 0; E < edgeElements_->numberElements() ; E++)
        if(repeatedEdgesArray[E] > 1)
            edgeElements->getElement(E).setInterfaceElement(true);

}

// We will build th eedge list based on the elements. We will ignore the P2 Points here.
template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::buildEdges(ElementsPtr_Type elements){

    TEUCHOS_TEST_FOR_EXCEPTION( elements.is_null() , std::runtime_error ,"Elements null.");
    //CommConstPtr_Type comm = this->getComm();
    //bool verbose(comm->getRank()==0);
    FEDD_TIMER_START(edgesTotal,"Mesh : Building edges from scratch.");

    FEDD_TIMER_START(edgeList,"Mesh : 0. Building edge list");

    MapConstPtr_Type repeatedMap =  this->getMapRepeated();
    // Building local edges with repeated node list
    vec2D_int_Type localEdgeIndices(0);
    setLocalEdgeIndices( localEdgeIndices );
    EdgeElementsPtr_Type edges = Teuchos::rcp( new EdgeElements_Type() );
    for (int i=0; i<elements->numberElements(); i++) {
        for (int j=0; j<localEdgeIndices.size(); j++) {
            int id1 = elements->getElement( i ).getNode( localEdgeIndices[j][0] );
            int id2 = elements->getElement( i ).getNode( localEdgeIndices[j][1] );
            vec_int_Type edgeVec(2);
            if (id1<id2){
                edgeVec[0] = id1;
                edgeVec[1] = id2;
            }
            else{
                edgeVec[0] = id2;
                edgeVec[1] = id1;
            }
            int flag=0.;
            if(min(bcFlagRep_->at(id1),bcFlagRep_->at(id2)) ==0 || max(bcFlagRep_->at(id1),bcFlagRep_->at(id2)) ==10 || max(bcFlagRep_->at(id1),bcFlagRep_->at(id2)) ==15)
                flag = min(bcFlagRep_->at(id1),bcFlagRep_->at(id2));
            else
                flag=max(bcFlagRep_->at(id1),bcFlagRep_->at(id2));

            FiniteElement edge( edgeVec,flag);
            edges->addEdge( edge, i ); // Adding Edge with local Element ID
            if ( !elements->getElement( i ).subElementsInitialized() )
                elements->getElement( i ).initializeSubElements("P1",1);
            elements->getElement( i ).addSubElement(edge);
        }
    }
    FEDD_TIMER_STOP(edgeList);
    vec2D_GO_Type combinedEdgeElements;
    FEDD_TIMER_START(edgeListUniqueTimer,"Mesh : 1. Make Edge List Unique");
    edges->sortUniqueAndSetGlobalIDsParallel(this->elementMap_,combinedEdgeElements);    
    edgeElements_ = edges;
    FEDD_TIMER_STOP(edgeListUniqueTimer);

	// ------------------------------------------------------------------------------------------------------
    // Part VIII: Updating the EdgeMap 
    // The general idea is to indentifiy an edge by their global node indices. 
    // This way two edges on different processors can be indentified and given the same global ID
    // First we determine the global IDs of non interface edges 
    // Then we determine those of interface edges with the procedure discribed above
    // ------------------------------------------------------------------------------------------------------
    FEDD_TIMER_START(edgeMapTimer,"Mesh: 2. Creating EdgeMap");		
	this->buildEdgeMap(); // This function is new to the Mesh_def Class
	FEDD_TIMER_STOP(edgeMapTimer);	


	// Part IX: Updating elementsOfEdgeGlobal and elementsOfEdgeLocal
    // this step is only necessary if we have more than 1 processor, as the edgeElements function work serially
    // we started the setup before (sortUniqueAndSetGlobalIDsParallel) an now finalize it with the information of other processors
    // the edges on the interface need the global element number of the neighbouring processor

    FEDD_TIMER_START(elementsOfEdgeTimer," Mesh: 3. Updating ElementsOfEdgeLocal/Global");		
    this->edgeElements_->setElementsEdges( combinedEdgeElements );

    this->edgeElements_->setUpElementsOfEdge( this->elementMap_, this->edgeMap_);

    int maxRank = std::get<1>(this->rankRange_);
    this->updateElementsOfEdgesLocalAndGlobal(maxRank);

    FEDD_TIMER_STOP(elementsOfEdgeTimer);

    FEDD_TIMER_STOP(edgesTotal);

    //for(int E = 0; E< this->edgeElements_->numberElements(); E++)
    //    if(this->edgeElements_->getElementsOfEdgeGlobal().at(E).size()>1)
    //        cout << " Edge " << this->getMapRepeated()->getGlobalElement(this->edgeElements_->getElement(E).getNode(0)) << " " << this->getMapRepeated()->getGlobalElement(this->edgeElements_->getElement(E).getNode(1)) << " has global element neighbour " << this->edgeElements_->getElementsOfEdgeGlobal().at(E).at(0) << " and " << this->edgeElements_->getElementsOfEdgeGlobal().at(E).at(1) << endl;
    //edges->setElementsEdges( combinedEdgeElements );
    Teuchos::TimeMonitor::report(cout);//,"Mesh Refinement");

}

/*!

 \brief Updating ElementsOfEdgesLocal and ElementsOfEdgesGlobal.

@param[in] maxRank The maximal processor rank.
@param[in] edgeMap Map of global edge ids.

*/

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::updateElementsOfEdgesLocalAndGlobal(int maxRank){

	if(maxRank >0 && this->dim_ == 2){
		vec_GO_Type edgesInterfaceGlobalID(0);
		LO id=0;
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(this->edgeElements_->getElement(i).isInterfaceElement() ){		
				this->edgeElements_->setElementsOfEdgeLocalEntry(i,-1);
				edgesInterfaceGlobalID.push_back(this->edgeMap_->getGlobalElement(i)); // extracting the global IDs of the new interfaceEdges
			}
			
		}

		// communticating elements across interface
		Teuchos::ArrayView<GO> edgesInterfaceGlobalID_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalID);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalID_, 0, this->comm_) );
		//mapGlobalInterface->print();

		// Global IDs of Procs
		// Setting newPoints as to be communicated Values
		MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP < SC > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			interfaceElementsEntries[i] = this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).at(0);
		}

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );

		MultiVectorLOPtr_Type isInterfaceElement_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_imp->putScalar(   0 ); 
		isInterfaceElement_imp->importFromVector( interfaceElements, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar(   0 ); 
		isInterfaceElement_exp->exportFromVector( interfaceElements, false, "Insert");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar(   0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_imp, false, "Insert");

		isInterfaceElement2_imp->exportFromVector(isInterfaceElement_exp, false, "Insert");

		interfaceElementsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			this->edgeElements_->setElementsOfEdgeGlobalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),interfaceElementsEntries[i]);
			}

	}

	// Contrary to the 2D case in 3D an edge can be part of more than two elements
	// We need to determine how many elements are connected to an edge and import the global IDs from different Processors
	if(maxRank >0 && this->dim_ == 3){

		// First we determine the interface edges, which entries in elementsOfEdgesGlobal/Local have to be completed
		vec_GO_Type edgesInterfaceGlobalID(0);
		vec_int_Type numberElements(0); // represents number of elements the edge is connected to on my Processor
		
		for(int i=0; i< this->edgeElements_->numberElements(); i++){
			if(this->edgeElements_->getElement(i).isInterfaceElement() ){
				edgesInterfaceGlobalID.push_back(this->edgeMap_->getGlobalElement(i)); // extracting the global IDs of the new interfaceEdges
			}	
		}
		sort(edgesInterfaceGlobalID.begin(), edgesInterfaceGlobalID.end());

		for(int i=0; i< edgesInterfaceGlobalID.size(); i++){
			numberElements.push_back(this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).size());
		}
			
		
		// from this we build a map
		Teuchos::ArrayView<GO> edgesInterfaceGlobalID_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalID);

		MapPtr_Type mapGlobalInterface =
			Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalID_, 0, this->comm_) );

		// As edges can be part of multiple elements on different processors we collect the number of elements connected to the edge in total
		MultiVectorLOPtr_Type numberInterfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP < SC > numberInterfaceElementsEntries  = numberInterfaceElements->getDataNonConst(0);

		for(int i=0; i< numberInterfaceElementsEntries.size(); i++)
			numberInterfaceElementsEntries[i] = numberElements[i];

		MapConstPtr_Type mapGlobalInterfaceUnique = mapGlobalInterface->buildUniqueMap( this->rankRange_ );
	
		// Element are unique to each processor. This means that the number we collect is the number of elements that are connected to my edge on other processors.
		// With the following communication we add up all the entries for a certain global Edge ID
		// Potential causes of error:
		//  - if an edge is not identified as an interface edge, of course it will not import nor export its interface Information, making itself and others incomplete

		MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
		isInterfaceElement_exp->putScalar(   0 ); 
		isInterfaceElement_exp->exportFromVector( numberInterfaceElements, false, "Add");

		MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		isInterfaceElement2_imp->putScalar(   0 ); 
		isInterfaceElement2_imp->importFromVector(isInterfaceElement_exp, true, "Insert");

		Teuchos::ArrayRCP < SC > numberInterfaceElementsImportEntries  = isInterfaceElement2_imp->getDataNonConst(0);

		vec_int_Type missingEntries(numberInterfaceElementsEntries.size());
		// With this number we can complete the elementsOfEdgeLocal List with -1 for the elements not on our processor
		for(int i=0; i<numberInterfaceElementsEntries.size() ; i++){
			for(int j=0; j< numberInterfaceElementsImportEntries[i] - numberInterfaceElementsEntries[i];j++){
				this->edgeElements_->setElementsOfEdgeLocalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),-1);
				missingEntries[i] = numberInterfaceElementsImportEntries[i] -numberInterfaceElementsEntries[i];
			}
		}

		// Next we need to identify the global Element IDs of those missing entries and communicate them
		// Hey i got the global Elements ... of edge ... -> exchange
		// Elements are uniquely distributed -> you cannot import an element you already have
		// I need x number of entries -> import all i need, export all i have 
		// Global IDs of Procs
		// Setting newPoints as to be communicated Values

		// Communicating max number of necessary values:
		vec_int_Type::iterator it;
		it = max_element(numberElements.begin(), numberElements.end());
		int myNumberElementsMax = numberElements.at(distance(numberElements.begin(), it)); // accumulate(errorElement.begin(), errorElement.end(),0.0);

		reduceAll<int, int> (*this->comm_, REDUCE_MAX,  myNumberElementsMax , outArg ( myNumberElementsMax));

		MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
		Teuchos::ArrayRCP < SC > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

		vec2D_int_Type importElements(this->edgeElements_->getElementsOfEdgeGlobal().size(),vec_int_Type( 0));

		// We extended this function to also consider the ranks. Before we would only exchange the information, but in case one processor received information from more than one processor at 
		// the same time some of the information would get lost. Now we only send the information one processor holds to the other at the same time and move through the processor only destributing their
		// edgeOfElementGlobal Information in a circle like order
		for(int k=0; k< maxRank+1 ; k++){
			
			vec_GO_Type edgesInterfaceGlobalIDProc;				
			if(this->comm_->getRank() == k ){
				edgesInterfaceGlobalIDProc = edgesInterfaceGlobalID; // extracting the global IDs of the new interfaceEdges
			}
					

			// from this we build a map
			Teuchos::ArrayView<GO> edgesInterfaceGlobalIDProc_ = Teuchos::arrayViewFromVector( edgesInterfaceGlobalIDProc);

			MapPtr_Type mapGlobalInterfaceProcs =
				Teuchos::rcp( new Map_Type( this->edgeMap_->getUnderlyingLib(), Teuchos::OrdinalTraits<GO>::invalid(), edgesInterfaceGlobalIDProc_, 0, this->comm_) );

			for(int j=0; j< myNumberElementsMax; j++){
				MultiVectorLOPtr_Type interfaceElements = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceProcs, 1 ) );
				Teuchos::ArrayRCP < SC > interfaceElementsEntries  = interfaceElements->getDataNonConst(0);

				for(int i=0; i< interfaceElementsEntries.size() ; i++){		
					if(numberElements[i] > j && this->comm_->getRank() == k )
						interfaceElementsEntries[i] = this->edgeElements_->getElementsOfEdgeGlobal(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i])).at(j);
					else
						interfaceElementsEntries[i] = -1; 
				}

				MultiVectorLOPtr_Type isInterfaceElement_exp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterfaceUnique, 1 ) );
				isInterfaceElement_exp->putScalar(   -1 ); 
				isInterfaceElement_exp->exportFromVector( interfaceElements, false, "Insert");

				if(this->comm_->getRank() == k && mapGlobalInterfaceUnique->getNodeNumElements() > 0){
					Teuchos::ArrayRCP < SC > interfaceElementsEntries_exp  = isInterfaceElement_exp->getDataNonConst(0);
					for(int i=0; i<  interfaceElementsEntries_exp.size() ; i++){
						LO id = mapGlobalInterface->getLocalElement(mapGlobalInterfaceUnique->getGlobalElement(i));
						interfaceElementsEntries_exp[i] = interfaceElementsEntries[id];
					}
									
				}

			
				MultiVectorLOPtr_Type isInterfaceElement2_imp = Teuchos::rcp( new MultiVectorLO_Type( mapGlobalInterface, 1 ) );
				isInterfaceElement2_imp->putScalar(   0 ); 
				isInterfaceElement2_imp->importFromVector(isInterfaceElement_exp, false, "Insert");

				interfaceElementsEntries  = isInterfaceElement2_imp->getDataNonConst(0);

				for(int i=0; i< interfaceElementsEntries.size() ; i++){
					if(this->comm_->getRank() != k && interfaceElementsEntries[i] != -1)
						importElements[i].push_back( interfaceElementsEntries[i]);
				}
				
			}

		}
		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			sort(importElements[i].begin(),importElements[i].end());
			importElements[i].erase( unique(importElements[i].begin(), importElements[i].end() ), importElements[i].end() );
			if(importElements[i].size() != missingEntries[i])
			   cout << " On Processor " << this->comm_->getRank() << " uneven number for edge imported: " << importElements[i].size() << " missing " << missingEntries[i] << " " << edgesInterfaceGlobalID[i] << endl; // " something went wrong while updating elementsOfEdgesGlobal as the imported entries do not match the supposed number of imported entries. Please check." << endl; 	
		}

		for(int i=0; i< interfaceElementsEntries.size() ; i++){
			for(int j=0; j < importElements[i].size(); j++){
				if(importElements[i][j] != -1)
					this->edgeElements_->setElementsOfEdgeGlobalEntry(this->edgeMap_->getLocalElement(edgesInterfaceGlobalID[i]),importElements[i][j]);
			}
		}

	}

}

template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::setLocalEdgeIndices(vec2D_int_Type &localEdgeIndices ){
    if ( dim_ == 2 ) {
        localEdgeIndices.resize(3,vec_int_Type(2,-1));
        localEdgeIndices.at(0).at(0) = 0;
        localEdgeIndices.at(0).at(1) = 1;
        localEdgeIndices.at(1).at(0) = 0;
        localEdgeIndices.at(1).at(1) = 2;
        localEdgeIndices.at(2).at(0) = 1;
        localEdgeIndices.at(2).at(1) = 2;
    }
    else if( dim_ == 3) {
        localEdgeIndices.resize(6,vec_int_Type(2,-1));
        localEdgeIndices.at(0).at(0) = 0;
        localEdgeIndices.at(0).at(1) = 1;
        localEdgeIndices.at(1).at(0) = 0;
        localEdgeIndices.at(1).at(1) = 2;
        localEdgeIndices.at(2).at(0) = 1;
        localEdgeIndices.at(2).at(1) = 2;
        localEdgeIndices.at(3).at(0) = 0;
        localEdgeIndices.at(3).at(1) = 3;
        localEdgeIndices.at(4).at(0) = 1;
        localEdgeIndices.at(4).at(1) = 3;
        localEdgeIndices.at(5).at(0) = 2;
        localEdgeIndices.at(5).at(1) = 3;
        
    }
}

/*!
 \brief function, that determines volume of tetrahedra.

@param[in] elements Elements.
@param[in] edgeElements Edges.
@param[in] points Points.

@param[out] volumeTetrahedra Volume of tetrahedras
*/
template <class SC, class LO, class GO, class NO>
vec_dbl_Type Mesh<SC,LO,GO,NO>::determineVolTet(ElementsPtr_Type elements,vec2D_dbl_ptr_Type points){
	
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

		
		volumeTetrahedra[k] = fabs(p3[0] * v_K[0] + p3[1] * v_K[1] +p3[2] * v_K[2]) / 6. ;
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
void Mesh<SC,LO,GO,NO>::calcDiamTetraeder(){
	
    vec_dbl_Type volTet = determineVolTet(this->elementsC_, this->pointsRep_);

	// vec_dbl_Type diamElements(elements->numberElements());
    double diamElement=0;
	double a,b,c,A,B,C;

	vec2D_dbl_Type p(4,vec_dbl_Type(this->dim_));
	for(int k=0; k< this->elementsC_->numberElements() ; k++){	
		
		// Calculating edges of Tetraeder
		vec_LO_Type nodeList = 	this->elementsC_->getElement(k).getVectorNodeListNonConst();
		
		for(int i=0; i<4; i++){
			p[i] = this->pointsRep_->at(nodeList[i]);
		}

		a = sqrt(pow(p[0][0] - p[1][0],2)+ pow(p[0][1] - p[1][1],2) +pow(p[0][2] - p[1][2],2));
		b = sqrt(pow(p[0][0] - p[2][0],2)+ pow(p[0][1] - p[2][1],2) +pow(p[0][2] - p[2][2],2));
		c = sqrt(pow(p[0][0] - p[3][0],2)+ pow(p[0][1] - p[3][1],2) +pow(p[0][2] - p[3][2],2));

		A = sqrt(pow(p[3][0] - p[2][0],2)+ pow(p[3][1] - p[2][1],2) +pow(p[3][2] - p[2][2],2));
		B = sqrt(pow(p[1][0] - p[3][0],2)+ pow(p[1][1] - p[3][1],2) +pow(p[1][2] - p[3][2],2));
		C = sqrt(pow(p[1][0] - p[2][0],2)+ pow(p[1][1] - p[2][1],2) +pow(p[1][2] - p[2][2],2));		


		diamElement = sqrt(((a*A+b*B+c*C)*(-a*A+b*B+c*C)*(a*A-b*B+c*C)*(a*A+b*B-c*C))) / (12*volTet[k]) ;	

        this->elementsC_->getElement(k).setDiamElement(diamElement);
	}
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
void Mesh<SC,LO,GO,NO>::calcRhoTetraeder(){
	
    vec_dbl_Type volTet = determineVolTet(this->elementsC_, this->pointsRep_);

    double rhoElement =0.;
	
	vec_LO_Type surfaceOfEl(4);
	
    vec2D_dbl_ptr_Type points = this->pointsRep_;
	vec_dbl_Type vecTmp(3),vecTmp1(3),vecTmp2(3);
    double lengthA, lengthB,lengthC,s1,A1,A2,A3,A4;
	for(int k=0; k< this->elementsC_->numberElements() ; k++){	
		// Calculating edges of Tetraeder
        // ----------
        // areas
        // ----------
        FiniteElement elementTmp = this->elementsC_->getElement(k);
        vecTmp[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(1)).at(0);
		vecTmp[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(1)).at(1);
		vecTmp[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(1)).at(2);
		lengthA = sqrt(pow(vecTmp[0],2)+pow(vecTmp[1],2)+pow(vecTmp[2],2));
	
		vecTmp1[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp1[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		vecTmp1[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(2)).at(2);
		lengthB = sqrt(pow(vecTmp1[0],2)+pow(vecTmp1[1],2)+pow(vecTmp1[2],2));

		vecTmp2[0] = points->at(elementTmp.getNode(1)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp2[1] = points->at(elementTmp.getNode(1)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		vecTmp2[2] = points->at(elementTmp.getNode(1)).at(2) - points->at(elementTmp.getNode(2)).at(2);
		lengthC = sqrt(pow(vecTmp2[0],2)+pow(vecTmp2[1],2)+pow(vecTmp2[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		A1 = sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));

        // ----------
        vecTmp[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(1)).at(0);
		vecTmp[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(1)).at(1);
		vecTmp[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(1)).at(2);
		lengthA = sqrt(pow(vecTmp[0],2)+pow(vecTmp[1],2)+pow(vecTmp[2],2));
	
		vecTmp1[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(3)).at(0);
		vecTmp1[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(3)).at(1);
		vecTmp1[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(3)).at(2);
		lengthB = sqrt(pow(vecTmp1[0],2)+pow(vecTmp1[1],2)+pow(vecTmp1[2],2));

		vecTmp2[0] = points->at(elementTmp.getNode(1)).at(0) - points->at(elementTmp.getNode(3)).at(0);
		vecTmp2[1] = points->at(elementTmp.getNode(1)).at(1) - points->at(elementTmp.getNode(3)).at(1);
		vecTmp2[2] = points->at(elementTmp.getNode(1)).at(2) - points->at(elementTmp.getNode(3)).at(2);
		lengthC = sqrt(pow(vecTmp2[0],2)+pow(vecTmp2[1],2)+pow(vecTmp2[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		A2 = sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));

        // ----------

        vecTmp[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		vecTmp[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(2)).at(2);
		lengthA = sqrt(pow(vecTmp[0],2)+pow(vecTmp[1],2)+pow(vecTmp[2],2));
	
		vecTmp1[0] = points->at(elementTmp.getNode(0)).at(0) - points->at(elementTmp.getNode(3)).at(0);
		vecTmp1[1] = points->at(elementTmp.getNode(0)).at(1) - points->at(elementTmp.getNode(3)).at(1);
		vecTmp1[2] = points->at(elementTmp.getNode(0)).at(2) - points->at(elementTmp.getNode(3)).at(2);
		lengthB = sqrt(pow(vecTmp1[0],2)+pow(vecTmp1[1],2)+pow(vecTmp1[2],2));

		vecTmp2[0] = points->at(elementTmp.getNode(2)).at(0) - points->at(elementTmp.getNode(3)).at(0);
		vecTmp2[1] = points->at(elementTmp.getNode(2)).at(1) - points->at(elementTmp.getNode(3)).at(1);
		vecTmp2[2] = points->at(elementTmp.getNode(2)).at(2) - points->at(elementTmp.getNode(3)).at(2);
		lengthC = sqrt(pow(vecTmp2[0],2)+pow(vecTmp2[1],2)+pow(vecTmp2[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		A3 = sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));

        // ----------
        vecTmp[0] = points->at(elementTmp.getNode(1)).at(0) - points->at(elementTmp.getNode(2)).at(0);
		vecTmp[1] = points->at(elementTmp.getNode(1)).at(1) - points->at(elementTmp.getNode(2)).at(1);
		vecTmp[2] = points->at(elementTmp.getNode(1)).at(2) - points->at(elementTmp.getNode(2)).at(2);
		lengthA = sqrt(pow(vecTmp[0],2)+pow(vecTmp[1],2)+pow(vecTmp[2],2));
	
		vecTmp1[0] = points->at(elementTmp.getNode(1)).at(0) - points->at(elementTmp.getNode(3)).at(0);
		vecTmp1[1] = points->at(elementTmp.getNode(1)).at(1) - points->at(elementTmp.getNode(3)).at(1);
		vecTmp1[2] = points->at(elementTmp.getNode(1)).at(2) - points->at(elementTmp.getNode(3)).at(2);
		lengthB = sqrt(pow(vecTmp1[0],2)+pow(vecTmp1[1],2)+pow(vecTmp1[2],2));

		vecTmp2[0] = points->at(elementTmp.getNode(2)).at(0) - points->at(elementTmp.getNode(3)).at(0);
		vecTmp2[1] = points->at(elementTmp.getNode(2)).at(1) - points->at(elementTmp.getNode(3)).at(1);
		vecTmp2[2] = points->at(elementTmp.getNode(2)).at(2) - points->at(elementTmp.getNode(3)).at(2);
		lengthC = sqrt(pow(vecTmp2[0],2)+pow(vecTmp2[1],2)+pow(vecTmp2[2],2));

		s1 = (lengthA+lengthB+lengthC)/2.;
		A4 = sqrt(s1*(s1-lengthA)*(s1-lengthB)*(s1-lengthC));


		rhoElement = (6*volTet[k]) / (A1 + A2 +A3 +A4); // DURCHMESSER

        this->elementsC_->getElement(k).setRhoElement(rhoElement);

	}

}


template <class SC, class LO, class GO, class NO>
void Mesh<SC,LO,GO,NO>::determineLongestEdge( ){
	// We have to determine which edge is longer, as we use the opposite node of the longer edge for the element construction
    vec2D_dbl_ptr_Type points = this->pointsRep_;

    vec_dbl_Type P1(this->dim_),P2(this->dim_);
	double maxLength=0.0;
	int maxEntry=0;
	LO p1ID,p2ID;
    vec2D_LO_Type edgeVec = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
	for(int k=0; k< this->elementsC_->numberElements() ; k++){	

        FiniteElement elementTmp = this->elementsC_->getElement(k);
        
        double max_length=0.;
        for(int i=0; i<edgeVec.size();i++){
            p1ID =this->elementsC_->getElement(k).getNode(edgeVec[i][0]);
            p2ID =this->elementsC_->getElement(k).getNode(edgeVec[i][1]);
            P1 = points->at(p1ID);
            P2 = points->at(p2ID);
            double sum=0;
            for(int j=0; j< P1.size();j++)
                sum += pow(P1[j]-P2[j],2);
            sum =  sqrt(sum);
            if( sum > max_length)
                max_length=  sum;
        }    
	
        this->elementsC_->getElement(k).setLongestEdgeLength(max_length);

		
	}
		
	
}

}
#endif
