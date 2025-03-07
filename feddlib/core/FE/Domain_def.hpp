#ifndef Domain_def_hpp
#define Domain_def_hpp
#include "Domain_decl.hpp"

/*!
 Definition of Domain

 @brief  Domain
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */

namespace FEDD {
template <class SC, class LO, class GO, class NO>
Domain<SC,LO,GO,NO>::Domain():
comm_(),
mesh_(),
mapVecFieldUnique_(),
mapVecFieldRepeated_(),
geometries2DVec_(),
geometries3DVec_(),
distancesToInterface_(),
interfaceMapUnique_(),
interfaceMapVecFieldUnique_(),
globalInterfaceMapUnique_(),
globalInterfaceMapVecFieldUnique_(),
partialGlobalInterfaceVecFieldMap_(),
otherGlobalInterfaceMapUnique_(),
otherGlobalInterfaceMapVecFieldUnique_(),
otherPartialGlobalInterfaceVecFieldMap_()
{

}

template <class SC, class LO, class GO, class NO>
Domain<SC,LO,GO,NO>::Domain(CommConstPtr_Type comm):
comm_(comm),
mesh_(),
mapVecFieldUnique_(),
mapVecFieldRepeated_(),
geometries2DVec_(),
geometries3DVec_(),
distancesToInterface_(),
interfaceMapUnique_(),
interfaceMapVecFieldUnique_(),
globalInterfaceMapUnique_(),
globalInterfaceMapVecFieldUnique_(),
partialGlobalInterfaceVecFieldMap_(),
otherGlobalInterfaceMapUnique_(),
otherGlobalInterfaceMapVecFieldUnique_(),
otherPartialGlobalInterfaceVecFieldMap_()
{

}

template <class SC, class LO, class GO, class NO>
Domain<SC,LO,GO,NO>::Domain(CommConstPtr_Type comm, int dimension):
comm_(comm),
mesh_(),
dim_(dimension),
mapVecFieldUnique_(),
mapVecFieldRepeated_(),
geometries2DVec_(),
geometries3DVec_(),
distancesToInterface_(),
interfaceMapUnique_(),
interfaceMapVecFieldUnique_(),
globalInterfaceMapUnique_(),
globalInterfaceMapVecFieldUnique_(),
partialGlobalInterfaceVecFieldMap_(),
otherGlobalInterfaceMapUnique_(),
otherGlobalInterfaceMapVecFieldUnique_(),
otherPartialGlobalInterfaceVecFieldMap_()
{

}

// Constructor for structured 2D Meshes.  
template <class SC, class LO, class GO, class NO>
Domain<SC,LO,GO,NO>::Domain(vec_dbl_Type coor, double l, double h, CommConstPtr_Type comm):
comm_(comm),
mesh_(),
mapVecFieldUnique_(),
mapVecFieldRepeated_(),
geometries2DVec_(),
geometries3DVec_(),
distancesToInterface_(),
interfaceMapUnique_(),
interfaceMapVecFieldUnique_(),
globalInterfaceMapUnique_(),
globalInterfaceMapVecFieldUnique_(),
partialGlobalInterfaceVecFieldMap_(),
    otherGlobalInterfaceMapUnique_(),
    otherGlobalInterfaceMapVecFieldUnique_(),
    otherPartialGlobalInterfaceVecFieldMap_()
{
    coorRec	= coor;
    length 	= l;
	height 	= h;
    width = -1;
    // Available 2D geometries 
    geometries2DVec_.reset(new string_vec_Type(0));
    geometries2DVec_->push_back("Square");
    geometries2DVec_->push_back("BFS");
    geometries2DVec_->push_back("SquareTPM");
    geometries2DVec_->push_back("structuredMiniTest");
//    geometries2DVec->push_back("REC");
}

// Constructor for 3D structured meshes
template <class SC, class LO, class GO, class NO>
Domain<SC,LO,GO,NO>::Domain(vec_dbl_Type coor, double l, double w, double h, CommConstPtr_Type comm):
comm_(comm),
mesh_(),
mapVecFieldUnique_(),
mapVecFieldRepeated_(),
geometries2DVec_(),
geometries3DVec_(),
distancesToInterface_(),
interfaceMapUnique_(),
interfaceMapVecFieldUnique_(),
globalInterfaceMapUnique_(),
globalInterfaceMapVecFieldUnique_(),
partialGlobalInterfaceVecFieldMap_()
{
    coorRec	= coor;
    length 	= l;
    width 	= w;
    height	= h;
    // Different geometries available as geometries. 
    geometries3DVec_.reset(new string_vec_Type(0));
    geometries3DVec_->push_back("Square"); // for 3D this is synonymous to Cube with 6-Element subcube structure
    geometries3DVec_->push_back("BFS"); // Backward-facing step geometry
    geometries3DVec_->push_back("Square5Element"); // this is a cube with different 5-Element per subcube structure


}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::info() const{

    LO minNumberNodes;
    LO maxNumberNodes;
    LO numberNodes = this->getMapUnique()->getNodeNumElements();
    Teuchos::reduceAll( *this->comm_, Teuchos::REDUCE_MIN, numberNodes, Teuchos::ptrFromRef(minNumberNodes) );
    Teuchos::reduceAll( *this->comm_, Teuchos::REDUCE_MAX, numberNodes, Teuchos::ptrFromRef(maxNumberNodes) );
    
    bool verbose(comm_->getRank()==0);
    if (verbose) {
        std::cout << "\t### Domain ###" << std::endl;
        std::cout << "\t### Dimension: "<< dim_ << std::endl;
        std::cout << "\t### Mesh type: "<< meshType_ << std::endl;
        std::cout << "\t### Mesh flags: "<< flagsOption_ << std::endl;
        std::cout << "\t### Subdomains: "<< n_ << std::endl;
        std::cout << "\t### H/h: "<< m_ << std::endl;
        std::cout << "\t### FE type: "<< FEType_ << std::endl;
        std::cout << "\t### Number Nodes: "<< mesh_->getMapUnique()->getMaxAllGlobalIndex()+1 << std::endl;
        std::cout << "\t### Minimum number of nodes: "<< minNumberNodes << "  Maximum number of nodes: " << maxNumberNodes << std::endl;
        std::cout << "\t### Empty ranks for coarse solves: "<< numProcsCoarseSolve_ << std::endl;
    }
}
    
template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::initializeFEData(){

    this->getElementsC()->initializeFEData( this->getPointsRepeated() );
}
    
template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type Domain<SC,LO,GO,NO>::getElementsFlag() const{
    return mesh_->getElementsFlag();
}


template <class SC, class LO, class GO, class NO>
LO Domain<SC,LO,GO,NO>::getApproxEntriesPerRow() const{
    if (this->dim_ == 2) {
        if ( this->FEType_ == "P1" ) {
            return 44;
        }
        else if ( this->FEType_ == "P2" ) {
            return 60;
        }
        else {
            return 60;
        }
    } else {
        if ( this->FEType_ == "P1" ) {
            return 400;
        }
        else if ( this->FEType_ == "P2" ) {
            return 460;
        }
        else {
            return 400;
        }
    }
}



template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::buildMesh(int flagsOption , std::string meshType, int dim, std::string FEType, int N, int M, int numProcsCoarseSolve){

    int geoNumber = checkGeomentry(meshType, dim);

#ifdef ASSERTS_WARNINGS
    MYASSERT(geoNumber!=-1, "Geometry not known for this Dimension.")
#endif

    MeshStrPtr_Type meshStructured = Teuchos::rcp(new MeshStr_Type(comm_));
    n_ = N;
    m_ = M;
    dim_ = dim;
    FEType_ = FEType;
    numProcsCoarseSolve_ = numProcsCoarseSolve;
    meshType_ = meshType;
    flagsOption_ = flagsOption;

    switch (dim) {
        case 2:
            switch (geoNumber) {
                case 0:
                    meshStructured->setGeometry2DRectangle(coorRec, length, height);
                    meshStructured->buildMesh2D(FEType, n_, m_, numProcsCoarseSolve);
                    break;
                case 1:
                    meshStructured->setGeometry2DRectangle(coorRec, length, height);
                    meshStructured->buildMesh2DBFS(FEType, n_, m_, numProcsCoarseSolve);
                    break;
                case 2:
                    meshStructured->setGeometry2DRectangle(coorRec, length, height);
                    meshStructured->buildMesh2DTPM(FEType, n_, m_, numProcsCoarseSolve);
                    break;
                case 3:
                    meshStructured->setGeometry2DRectangle(coorRec, length, height);
                    meshStructured->buildMesh2DMiniTPM(FEType, n_, m_, numProcsCoarseSolve);
                    break;
                default:
                    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Select valid mesh. Structured types are 'structured' and 'structured_bfs' in 2D. TPM test meshes also available.");
                    break;
            }

            break;
        case 3:
            switch (geoNumber) {
                case 0:
                    meshStructured->setGeometry3DBox(coorRec, length, width, height);
                    meshStructured->buildMesh3D( FEType, n_, m_, numProcsCoarseSolve);
                    break;
                case 1:
                    meshStructured->setGeometry3DBox(coorRec, length, width, height);
                    meshStructured->buildMesh3DBFS(	FEType, n_, m_, numProcsCoarseSolve);
                    break;
                case 2:
                    meshStructured->setGeometry3DBox(coorRec, length, width, height);
                    meshStructured->buildMesh3D5Elements(	FEType, n_, m_, numProcsCoarseSolve);
                break;
                default:
                    TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"Select valid mesh. Structured types are 'structured' and 'structured_bfs' in 3D." );
                    break;
            }
            break;
        default:
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Select valid mesh dimension. 2 or 3 dimensional meshes can be constructed.");
            break;
    }
    meshStructured->buildElementMap();
    meshStructured->setStructuredMeshFlags(flagsOption,FEType);
    meshStructured->buildSurfaces(flagsOption,FEType);
    
    mesh_ = meshStructured;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::initializeUnstructuredMesh(int dimension, string feType, int volumeID){
    
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp(new MeshUnstr_Type(comm_, volumeID));
    mesh_ = meshUnstructured;
    mesh_->dim_ = dimension;
    FEType_ = feType;
    meshType_ = "unstructured";
    numProcsCoarseSolve_ = 0;
    n_ = comm_->getSize() - numProcsCoarseSolve_;
    m_ = -1;
    flagsOption_ = -1;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::readMeshSize(string filename, string delimiter){
    
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( mesh_ );
    TEUCHOS_TEST_FOR_EXCEPTION( meshUnstructured.is_null(), std::runtime_error, "Unstructured Mesh is null." );
    
    meshUnstructured->setMeshFileName( filename, delimiter );
    meshUnstructured->readMeshSize( );
    
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::partitionMesh( bool partitionDistance ){

    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( mesh_ );
    
    if (partitionDistance)
        partitionDistanceToInterface();

    mesh_ = meshUnstructured;

}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::partitionDistanceToInterface( ){

    vec_dbl_Type tmp = *distancesToInterface_;

    distancesToInterface_.reset( new vec_dbl_Type ( mesh_->getMapRepeated()->getNodeNumElements() ) );

    for (UN i=0; i<distancesToInterface_->size(); i++) {
        GO index = mesh_->getMapRepeated()->getGlobalElement(i);
        distancesToInterface_->at(i) = tmp[index];
    }
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::readAndPartitionMesh( std::string filename, std::string delimiter, int dim, std::string FEType, int volumeID ){

    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp(new MeshUnstr_Type(comm_, volumeID));

    dim_ = dim;
    FEType_ = FEType;

    mesh_ = meshUnstructured;
    mesh_->dim_ = dim_;
    numProcsCoarseSolve_ = 0;
    n_ = comm_->getSize() - numProcsCoarseSolve_;
    m_ = -1;
    flagsOption_ = -1;
    meshType_ = "unstructured";
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::buildP2ofP1Domain( DomainPtr_Type domainP1 ){ //P1 mesh must be parallel

    MeshUnstrPtr_Type meshUnstructuredP1 = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP1->mesh_ );
    
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp( new MeshUnstr_Type( comm_, meshUnstructuredP1->volumeID_) );

    n_ = domainP1->n_;
    m_ = domainP1->m_;
    dim_ = domainP1->dim_;
    FEType_ = "P2";
    numProcsCoarseSolve_ = 0;
    flagsOption_ = -1;
    meshType_ = "unstructured";

    meshUnstructured->buildP2ofP1MeshEdge( meshUnstructuredP1 );
    mesh_ = meshUnstructured;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::initWithDomain(DomainPtr_Type domainP1){ 

	n_ = domainP1->n_;
	m_ = domainP1->m_;
    dim_ = domainP1->dim_;
	FEType_ = domainP1->FEType_;

    numProcsCoarseSolve_ = 0;
    flagsOption_ = -1;

    meshType_ = domainP1->meshType_;
	
    mesh_ = domainP1->mesh_;
}

/*template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::initMeshRef( DomainPtr_Type domainP1 ){ 
	// Initialize MeshRefinementType as through other function like meshPartitioner and buildP2OfP1 Mesh meshUnstr Type is required

	MeshUnstrPtr_Type meshUnstr = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainP1->mesh_ , true);
	MeshUnstrRefPtr_Type meshUnstrRefTmp = Teuchos::rcp( new MeshUnstrRef_Type( comm_, meshUnstr->volumeID_, meshUnstr ) );
	mesh_ = meshUnstrRefTmp;

}*/

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::setMesh(MeshUnstrPtr_Type meshUnstr){ 

    mesh_ = meshUnstr;
}


template <class SC, class LO, class GO, class NO>
void Domain<SC, LO, GO, NO>::initDummyMesh(MapPtr_Type map)
{
    MeshUnstrPtr_Type outputMesh = Teuchos::rcp( new MeshUnstr_Type( comm_) );

    outputMesh->dim_ = this->dim_ ;
	outputMesh->FEType_ = this->FEType_ ;

	outputMesh->mapUnique_ = map;
	outputMesh->mapRepeated_ = map;
	
    mesh_ = outputMesh;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::exportMesh(bool exportEdges, bool exportSurfaces, string exportMesh){ 

    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( mesh_ );

    meshUnstructured->exportMesh(this->getMapUnique() , this->getMapRepeated(), exportEdges, exportSurfaces, exportMesh);
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::preProcessMesh(bool correctSurfaceNormals, bool correctElementOrientation){ 

    if(correctSurfaceNormals)
        mesh_->correctNormalDirections();

    if(correctElementOrientation)
        mesh_->correctElementOrientation();


}


template <class SC, class LO, class GO, class NO>
UN Domain<SC,LO,GO,NO>::getDimension() const{

    return dim_;
}

template <class SC, class LO, class GO, class NO>
std::string Domain<SC,LO,GO,NO>::getFEType() const{

    return FEType_;
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::CommConstPtr_Type Domain<SC,LO,GO,NO>::getComm() const{
    return comm_;
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getMapUnique() const{

    return mesh_->getMapUnique();
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getMapRepeated() const{

    return mesh_->getMapRepeated();
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getMapUniqueP2() const{

    return mesh_->getMapUniqueP2();
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getMapRepeatedP2() const{

    return mesh_->getMapRepeatedP2();
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getElementMap() const{

    return mesh_->getElementMap();
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getEdgeMap() const{

    return mesh_->getEdgeMap();
}

template <class SC, class LO, class GO, class NO>
vec2D_dbl_ptr_Type Domain<SC,LO,GO,NO>::getPointsRepeated() const{

    return mesh_->getPointsRepeated();
}

template <class SC, class LO, class GO, class NO>
vec2D_dbl_ptr_Type Domain<SC,LO,GO,NO>::getPointsUnique() const{

    return mesh_->getPointsUnique();
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type Domain<SC,LO,GO,NO>::getBCFlagRepeated() const{

    return mesh_->getBCFlagRepeated();
}

template <class SC, class LO, class GO, class NO>
vec_int_ptr_Type Domain<SC,LO,GO,NO>::getBCFlagUnique() const{

    return mesh_->getBCFlagUnique();
}

template <class SC, class LO, class GO, class NO>
vec2D_int_ptr_Type Domain<SC,LO,GO,NO>::getElements() const{
    TEUCHOS_TEST_FOR_EXCEPTION(mesh_->getElementsC().is_null(), std::runtime_error, "Elements is null for this mesh.");
    return mesh_->getElements();
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::ElementsPtr_Type Domain<SC,LO,GO,NO>::getElementsC() const{
    TEUCHOS_TEST_FOR_EXCEPTION(mesh_->getElementsC().is_null(), std::runtime_error, "Elements is null for this mesh.");
    return mesh_->getElementsC();
}


template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getMapVecFieldUnique() const{
    if ( mapVecFieldUnique_.is_null() ) {
        MapConstPtr_Type mapTmp = this->getMapUnique();
        mapVecFieldUnique_ = mapTmp->buildVecFieldMap(dim_);
    }
    return mapVecFieldUnique_;
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getMapVecFieldRepeated() const{
    if ( mapVecFieldRepeated_.is_null() ) {
        MapConstPtr_Type mapTmp = this->getMapRepeated();
        mapVecFieldRepeated_ = mapTmp->buildVecFieldMap(dim_);
    }
    return mapVecFieldRepeated_;
}

template <class SC, class LO, class GO, class NO>
GO Domain<SC,LO,GO,NO>::getNumElementsGlobal() const{

    return mesh_->getNumElementsGlobal();
}

template <class SC, class LO, class GO, class NO>
LO Domain<SC,LO,GO,NO>::getNumElements() const{

    return mesh_->getNumElements();
}

template <class SC, class LO, class GO, class NO>
LO Domain<SC,LO,GO,NO>::getNumPoints(std::string type) const{

    return mesh_->getNumPoints(type);
}

template <class SC, class LO, class GO, class NO>
int Domain<SC,LO,GO,NO>::checkGeomentry(std::string meshType, int dim) const{

    int notfoundLabel;
    switch (dim) {
        case 2:
            for (int i = 0; i<geometries2DVec_->size(); i++) {
                notfoundLabel = meshType.compare(geometries2DVec_->at(i));
                if (notfoundLabel==0) {
                    return i;
                }
            }
            break;
        case 3:
            for (int i = 0; i<geometries3DVec_->size(); i++) {
                notfoundLabel = meshType.compare(geometries3DVec_->at(i));
                if (notfoundLabel==0) {
                    return i;
                }
            }
            break;
        default:
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Select valid geometry.");
            break;
    }
    return -1;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::identifyInterfaceParallelAndDistance( DomainPtr_Type domainOther, vec_int_Type interfaceID_vec ){
    
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( mesh_ );
    MeshUnstrPtr_Type meshUnstructuredOther = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( domainOther->mesh_ );
    meshUnstructured->buildMeshInterfaceParallelAndDistance( meshUnstructuredOther, interfaceID_vec, distancesToInterface_ );
    
}
    
// Abstand zum Interface berechnen
template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::calculateDistancesToInterface()
{
    // IdentifyInterface() ist nur fuer unstrukturierte Gitter programmiert
    // Mesh_->MeshInterface_ gibt es somit nicht.
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( mesh_ );

    // Fuer jede Interface-Flag gibt es "aussen" einen Eintrag im Vektor.
    // Wenn das Interface also fuer 2 Flags bestimmt wird, dann ist IndicesGlobalMatched_.size() = 2.
    // Danach gibt es zwei Zeilenvektoren mit der globalen ID fuer this (thisMesh) und otherMesh,
    // also z.B. fuer Fluid (this) und Struktur (other).
    // Das heisst wir haben IndicesGlobalMatched_.at(0).size() = 2.
    // Letztlich haben wir NumberOfInterfaceNodesWithFlag = IndicesGlobalMatched_.at(0).at(0).size().
    // In den Eintragen IndicesGlobalMatched_.at(0).at(0).at(i) steht die globale ID von Interfaceknoten i.
    vec3D_GO_ptr_Type indicesGlobalMatched = meshUnstructured->meshInterface_->getIndicesGlobalMatched();

    // SourceNodes (z.B. Fluidgebiet) aus dem Gitter ziehen
    // Gefaehrlich, da Veraenderungen an SourceNodesRep auch PointsRep_ veraendern;
    // Sonst auf 0.0 setzen und dann Werte reinschreiben.
    // Am besten waere es Getter zu haben mit return vec2D_dbl_ptr_Type bzw. const vec2vec2D_dbl_ptr_Type, damit nichts passieren kann.
    vec2D_dbl_ptr_Type sourceNodesRep = meshUnstructured->getPointsRepeated();

    // Zaehle mit Hilfe von indicesGlobalMatched wie viele Interfaceknoten es insgesamt gibt.
    int numberInterfaceNodes = 0;
    for(int i = 0; i < indicesGlobalMatched->size(); i++) // Schleife ueber jede Flag
    {
        for(int j = 0; j < indicesGlobalMatched->at(i).at(0).size(); j++) // Schleife ueber jeden Interfaceknoten mit der Flag
        {
            numberInterfaceNodes = numberInterfaceNodes + 1;
        }
    }

    // EndNodes (Interfaceknoten) aus dem Gitter ziehen; mit Hilfe von IndicesGlobalMatched
    vec2D_dbl_ptr_Type endNodesRep(new vec2D_dbl_Type( numberInterfaceNodes, vec_dbl_Type( dim_, 0.0 ) ) );
    int counter = 0; // Die Zeile in die bei EndNodesRep hineingeschrieben wird
    int globalIDOfInterfaceNode;
    for(int i = 0; i < indicesGlobalMatched->size(); i++)
    {
        for(int j = 0; j < indicesGlobalMatched->at(i).at(0).size(); j++)
        {
            // ->at(i).at(0) ist this(z.B. Fluid) und ->at(i).at(1) waere other (z.B. Struktur).
            globalIDOfInterfaceNode = indicesGlobalMatched->at(i).at(0).at(j);

            for(int k = 0; k < dim_; k++) // Schleife ueber x- und y-Koordinate und ggf. z
            {
                endNodesRep->at(counter).at(k) = meshUnstructured->getPointsRepeated()->at( globalIDOfInterfaceNode ).at( k );
            }

            counter = counter + 1;
        }
    }

    // Distanz zum Interface auf eine unrealistisch grosse Zahl setzen.
    // In DistancesToInterface_->at(i) steht die Distanz der globalen Knoten-ID i zum Interface
    vec_dbl_ptr_Type tempVec( new vec_dbl_Type( meshUnstructured->getPointsRepeated()->size(), 1000.0 ) );
    distancesToInterface_ = tempVec;

    // Berechne nun den Abstand von den SourceNodesRep zu den EndNodesRep
    double distance = 0.0;
    for(int i = 0; i < sourceNodesRep->size(); i++)
    {
        for(int j = 0; j < endNodesRep->size(); j++)
        {
            for(int k = 0; k < dim_; k++)
            {
                distance = distance + pow( sourceNodesRep->at(i).at(k) - endNodesRep->at(j).at(k), 2.0 );
            }

            // Noch die Wurzel ziehen
            distance = sqrt(distance);

            if(distancesToInterface_->at(i) > distance)
            {
                distancesToInterface_->at(i) = distance;
            }

            // Reseten
            distance = 0.0;

        }
    }
}
    
template <class SC, class LO, class GO, class NO>
vec_dbl_ptr_Type Domain<SC,LO,GO,NO>::getDistancesToInterface() const{
    return distancesToInterface_;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::setReferenceConfiguration()
{
    mesh_->setReferenceConfiguration();
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MeshPtr_Type Domain<SC,LO,GO,NO>::getMesh(){
    return mesh_;
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MeshConstPtr_Type Domain<SC,LO,GO,NO>::getMesh() const{
    MeshConstPtr_Type meshConst = mesh_;
    return meshConst;
}
    
template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::moveMesh(MultiVectorPtr_Type displacementUnique, MultiVectorPtr_Type displacementRepeated)
{
    mesh_->moveMesh(displacementUnique, displacementRepeated);
}

//only need in geometry problem.
template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::buildUniqueInterfaceMaps()
{
    // Mesh_ umcasten  in unstructured und dann meshUnstructured nutzen,
    // da nur unstructured Attribut MeshInterface_ besitzt.
    // Jeder Prozessor kennt also das komplette matched Interface
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( this->getMesh() );

    vec3D_GO_ptr_Type indicesGlobalMatchedOrigin = meshUnstructured->getMeshInterface()->getIndicesGlobalMatchedOrigin();

    // lokale Interface ID, die ich dem Knoten zuschreibe; Wir zaehlen von 0 bis #\Gamma-1
    GO localInterfaceID = 0; // long long

    // In diesen beiden Vektoren stehen die Interface IDs in der Interface-Nummerierung,
    // die der Prozessor haelt. Unique!!!
    vec_GO_Type vecInterfaceMap; // vec_long ist long long wg. 64

    // indicesGlobalMatchedOrigin sind die IndicesGlobalMatched von this-Sicht aus.
    // D.h. in indicesGlobalMatchedOrigin.at(0).at(0) ist Fluid GID und
    // indicesGlobalMatchedOrigin.at(0).at(1) ist Struktur GID aus Fluid-Sicht.

    // Fluid und Struktur (bzw. this und other) haben gleich viele Interfaceknoten, deswegen
    // muessen wir hier nicht zwischen differenzieren
    // ACHTUNG: indicesGlobalMatchedOrigin ist hier nicht partitioniert!!!!
    for(int i = 0; i < indicesGlobalMatchedOrigin->size(); i++) // Schleife ueber jede flag
    {
        for(int j = 0; j < indicesGlobalMatchedOrigin->at(i).at(0).size(); j++) // GIDs innerhalb der flag (vom Fluid)
        {
            // Wir muessen long long anstatt int nutzen, da MyGID auf 64 gestellt ist/ genutzt wird
            GO globalIDOfInterfaceNode = indicesGlobalMatchedOrigin->at(i).at(0).at(j);
            // liefert true, falls Proz. GID besitzt
            if( this->getMapUnique()->getLocalElement(globalIDOfInterfaceNode) != Teuchos::OrdinalTraits<LO>::invalid())
            {
                // Schreibe die lokale Interface ID hinein
                vecInterfaceMap.push_back(localInterfaceID);
            }

            localInterfaceID = localInterfaceID + 1;

        }
    }

    // Am Ende steht in localInterfaceID wie viele Interface-Knoten es insgesamt gibt
    GO numberInterfaceNodes = localInterfaceID; // long long wg. 64
    
    // Baue nun die InterfaceMap (node)
    std::string ulib = this->getMapUnique()->getUnderlyingLib();
    Teuchos::ArrayView<GO> vecInterfaceMapArray =  Teuchos::arrayViewFromVector( vecInterfaceMap );
    interfaceMapUnique_ = Teuchos::rcp(new Map_Type( ulib, numberInterfaceNodes, vecInterfaceMapArray, 0, comm_ ) ); //maybe numberInterfaceNodes instead of -1
    // dof-Map bauen
    interfaceMapVecFieldUnique_ = interfaceMapUnique_->buildVecFieldMap(dim_/*dofs*/);
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getInterfaceMapUnique() const
{
    return interfaceMapUnique_;
}


template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getInterfaceMapVecFieldUnique() const
{
    return interfaceMapVecFieldUnique_;
}


template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getGlobalInterfaceMapUnique() const
{
    return globalInterfaceMapUnique_;
}


template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getGlobalInterfaceMapVecFieldUnique() const
{
    return globalInterfaceMapVecFieldUnique_;
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MapConstPtr_Type Domain<SC,LO,GO,NO>::getOtherGlobalInterfaceMapVecFieldUnique() const
{
    return otherGlobalInterfaceMapVecFieldUnique_;
}
    
template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::setPartialCoupling(int flag, std::string type){

    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( this->getMesh() );
    meshUnstructured->getMeshInterface()->setPartialCoupling(flag, type);
    
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::buildInterfaceMaps()
{
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp_dynamic_cast<MeshUnstr_Type>( this->getMesh() );
    MeshInterfacePtr_Type interface = meshUnstructured->getMeshInterface();
    vec3D_GO_ptr_Type indicesMatchedUni = interface->getIndicesGlobalMatchedUnique();
    
    vec3D_GO_ptr_Type indicesMatchedGlobalSerial = interface->getIndicesGlobalMatchedOrigin();
    
    MapConstPtr_Type mapUni = this->getMapUnique();
    vec_int_ptr_Type flagPointsUni = this->getBCFlagUnique();
    vec_GO_Type vecGlobalInterfaceID(0);
    vec_GO_Type vecOtherGlobalInterfaceID(0);
    vec_GO_Type vecInterfaceID(0);
    vec_int_Type vecInterfaceFlag(0);
    LO localID = 0;

    for(int i = 0; i < indicesMatchedGlobalSerial->size(); i++) // Schleife ueber jede flag
    {
        for(int j = 0; j < indicesMatchedGlobalSerial->at(i).at(0).size(); j++) {// GIDs innerhalb der flag (vom Fluid)
            LO index = mapUni->getLocalElement( indicesMatchedGlobalSerial->at(i).at(0).at(j) );
            if ( index != Teuchos::OrdinalTraits<LO>::invalid() ){
                vecGlobalInterfaceID.push_back( indicesMatchedGlobalSerial->at(i).at(0).at(j) );
                vecOtherGlobalInterfaceID.push_back( indicesMatchedGlobalSerial->at(i).at(1).at(j) );
                vecInterfaceID.push_back( (GO) localID );
                vecInterfaceFlag.push_back( (*flagPointsUni)[index] );
            }
            localID++;
        }
    }
//    std::sort( vecGlobalInterfaceID.begin(), vecGlobalInterfaceID.end() );
    // Baue nun die InterfaceMap fuer Fluid oder Struktur
    std::string ulib = this->getMapUnique()->getUnderlyingLib();
    Teuchos::ArrayView<GO> vecInterfaceGlobalMapArray =  Teuchos::arrayViewFromVector( vecGlobalInterfaceID );
    Teuchos::ArrayView<GO> vecInterfaceMapArray =  Teuchos::arrayViewFromVector( vecInterfaceID );
    Teuchos::ArrayView<GO> vecOtherInterfaceGlobalMapArray =  Teuchos::arrayViewFromVector( vecOtherGlobalInterfaceID );
    
    globalInterfaceMapUnique_ = Teuchos::rcp(new Map_Type( ulib, -1, vecInterfaceGlobalMapArray, 0, comm_ ) );
    interfaceMapUnique_ = Teuchos::rcp(new Map_Type( ulib, -1, vecInterfaceMapArray, 0, comm_ ) );

    otherGlobalInterfaceMapUnique_ = Teuchos::rcp(new Map_Type( ulib, -1, vecOtherInterfaceGlobalMapArray, 0, comm_ ) );

    
    if ( interface->sizePartialCoupling() == 0 ) {
        globalInterfaceMapVecFieldUnique_ = globalInterfaceMapUnique_->buildVecFieldMap(dim_);
        interfaceMapVecFieldUnique_ = interfaceMapUnique_->buildVecFieldMap(dim_);
        otherGlobalInterfaceMapVecFieldUnique_ = otherGlobalInterfaceMapUnique_->buildVecFieldMap(dim_);
    }
    // we add the dofs which are not part of the interface as dummy interface dofs, which do nothing (unit matrix in the diagonal coupling block). We need to do this for the monolithic use of FROSch, since we need a constant number of dofs for a node for each block
    else{
        if (this->comm_->getRank() == 0)
            std::cout << "-- ### Warning! Unique and unique-vector-field map might not be compatiable due to partial coupling conditions. ### --"  << std::endl;
        int numDofs = dim_;
        // we use nodewise ordering
        Teuchos::ArrayView<const GO> elementListGlobal = globalInterfaceMapUnique_->getNodeElementList();
        Teuchos::ArrayView<const GO> otherElementListGlobal = otherGlobalInterfaceMapUnique_->getNodeElementList();
        
        Teuchos::Array<GO> elementListFieldGlobal( 0 );
        Teuchos::Array<GO> elListFieldPartial( 0 );
    
        Teuchos::Array<GO> otherElementListFieldGlobal( 0 );
        Teuchos::Array<GO> otherElListFieldPartial( 0 );
        
        for (UN i=0; i<elementListGlobal.size(); i++) {
            int loc = interface->isPartialCouplingFlag( vecInterfaceFlag[i] );
            if ( loc > -1 ) {
                std::string partialType = interface->getPartialCouplingType( loc );
                if (partialType == "X_Y" && dim_ == 3) {
                    elementListFieldGlobal.push_back (numDofs * elementListGlobal[i] + 0 );
                    elementListFieldGlobal.push_back (numDofs * elementListGlobal[i] + 1 );
                    elementListFieldGlobal.push_back (numDofs * elementListGlobal[i] + 2 );
                    elListFieldPartial.push_back (numDofs * elementListGlobal[i] + 0 );
                    elListFieldPartial.push_back (numDofs * elementListGlobal[i] + 1 );
                    
                    otherElementListFieldGlobal.push_back (numDofs * otherElementListGlobal[i] + 0 );
                    otherElementListFieldGlobal.push_back (numDofs * otherElementListGlobal[i] + 1 );
                    otherElementListFieldGlobal.push_back (numDofs * otherElementListGlobal[i] + 2 );
                    otherElListFieldPartial.push_back (numDofs * otherElementListGlobal[i] + 0 );
                    otherElListFieldPartial.push_back (numDofs * otherElementListGlobal[i] + 1 );

                }
                else
                    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Implement other partial coupling types.");
            }
            else{
                for (UN dof=0; dof<numDofs; dof++){
                    elementListFieldGlobal.push_back(  numDofs * elementListGlobal[i] + dof );
                    elListFieldPartial.push_back ( numDofs * elementListGlobal[i] + dof );
                    otherElementListFieldGlobal.push_back(  numDofs * otherElementListGlobal[i] + dof );
                    otherElListFieldPartial.push_back ( numDofs * otherElementListGlobal[i] + dof );
                }
            }
        }
        globalInterfaceMapVecFieldUnique_ = Teuchos::rcp(new Map_Type( globalInterfaceMapUnique_->getUnderlyingLib(), -1, elementListFieldGlobal(), 0/*index base*/, this->getComm() ) );
        otherGlobalInterfaceMapVecFieldUnique_ = Teuchos::rcp(new Map_Type( globalInterfaceMapUnique_->getUnderlyingLib(), -1, otherElementListFieldGlobal(), 0/*index base*/, this->getComm() ) );
        // This is only a temporary map. We need to make sure that the dummy values have a higher GID, as the real coupling values. Otherwise we might get a problem during the assembly of the coupling matrices
        interfaceMapVecFieldUnique_ = Teuchos::rcp(new Map_Type( globalInterfaceMapUnique_->getUnderlyingLib(), -1, elementListFieldGlobal.size(), 0/*index base*/, this->getComm() ) );
        
        partialGlobalInterfaceVecFieldMap_ = Teuchos::rcp(new Map_Type( interfaceMapUnique_->getUnderlyingLib(), -1, elListFieldPartial(), 0/*index base*/, this->getComm() ) );
        otherPartialGlobalInterfaceVecFieldMap_ = Teuchos::rcp(new Map_Type( interfaceMapUnique_->getUnderlyingLib(), -1, otherElListFieldPartial(), 0/*index base*/, this->getComm() ) );
        
    }

}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::toNodeID(UN dim, GO dofID, GO& nodeID, LO& localDofNumber )
{
    nodeID = (GO) (dofID/dim);
    localDofNumber = (LO) (dofID%dim);
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::toDofID(UN dim, GO nodeID, LO localDofNumber, GO& dofID )
{
    dofID = (GO) ( dim * nodeID + localDofNumber);
}

template <class SC, class LO, class GO, class NO>
vec_long_Type Domain<SC,LO,GO,NO>::getLocalInterfaceIDInGlobal() const
{
    return vecLocalInterfaceIDinGlobal_;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC,LO,GO,NO>::setDummyInterfaceDomain(DomainPtr_Type domain)
{
    MeshUnstrPtr_Type meshUnstructured = Teuchos::rcp(new MeshUnstr_Type(comm_, 10/*default volume flag*/));
    mesh_ = meshUnstructured;
    
    this->dim_ = domain->getDimension();
    this->mesh_->mapUnique_ = Teuchos::rcp_const_cast<Map_Type>( domain->getInterfaceMapUnique() );
    this->mesh_->mapRepeated_ = Teuchos::rcp_const_cast<Map_Type>( domain->getInterfaceMapUnique() );
    this->mapVecFieldRepeated_ = Teuchos::rcp_const_cast<Map_Type>( domain->getInterfaceMapVecFieldUnique() );
    this->mapVecFieldUnique_ = Teuchos::rcp_const_cast<Map_Type>( domain->getInterfaceMapVecFieldUnique() );
    
    MapConstPtr_Type mapGI = domain->getGlobalInterfaceMapUnique();
    MapConstPtr_Type mapD = domain->getMapUnique();
    this->mesh_->pointsRep_.reset(new std::vector<std::vector<SC> >(mapGI->getNodeNumElements(), std::vector<SC>(this->dim_, 0.0 ) ) );
    this->mesh_->pointsUni_.reset(new std::vector<std::vector<SC> >(mapGI->getNodeNumElements(), std::vector<SC>(this->dim_, 0.0 ) ) );
    vec2D_dbl_ptr_Type pointsUni = domain->getPointsUnique();
    for (int i=0; i<mapGI->getNodeNumElements(); i++) {
        LO index = mapD->getLocalElement( mapGI->getGlobalElement( i ) );
        for (int j=0; j<this->dim_; j++){
            (*this->mesh_->pointsRep_)[i][j] = (*pointsUni)[index][j];
            (*this->mesh_->pointsUni_)[i][j] = (*pointsUni)[index][j];
        }
    }
    
}
//template <class SC, class LO, class GO, class NO>
//void Domain<SC,LO,GO,NO>::setDummyInterfaceMesh( MeshPtr_Type mesh )
//{
//    this->mesh_ = mesh;
//}
template <class SC, class LO, class GO, class NO>
int Domain<SC,LO,GO,NO>::findInPointsUnique(const vec_dbl_Type& x) const{

    double eps = 0.0000001;
    
    vec2D_dbl_ptr_Type points = this->getPointsUnique();
    
    if (this->getDimension()==2) {
        auto iterator = std::find_if( points->begin(), points->end(),
                                     [&] (const vector<double>& a){
                                         if (a[0] >= x[0]-eps && a[0] <= x[0]+eps
                                             && a[1] >= x[1]-eps && a[1] <= x[1]+eps)
                                             return true;
                                         else
                                             return false;
                                     }
                                     );
        if ( iterator != points->end() )
            return std::distance(points->begin(),iterator);
        else
             return -1;
    }
    else if(this->getDimension()==3) {
        auto iterator = std::find_if(points->begin(),points->end(),
                                     [&] (const vector<double>& a){
                                         if (a[0] >= x[0]-eps && a[0] <= x[0]+eps
                                             && a[1] >= x[1]-eps && a[1] <= x[1]+eps
                                             && a[2] >= x[2]-eps && a[2] <= x[2]+eps)
                                             return true;
                                         else
                                             return false;
                                     }
                                     );
        if ( iterator != points->end() )
            return std::distance(points->begin(),iterator);
        else
            return -1;
    }
    return -1;
}

template <class SC, class LO, class GO, class NO>
typename Domain<SC,LO,GO,NO>::MultiVectorPtr_Type Domain<SC,LO,GO,NO>::getNodeListMV() const{
    
    MapConstPtr_Type map = this->getMapRepeated();
    MultiVectorPtr_Type nodeList = Teuchos::rcp( new MultiVector_Type ( map, dim_ ) );
    vec2D_dbl_ptr_Type pointsRepeated = this->getPointsRepeated();
    for (int i=0; i<nodeList->getLocalLength(); i++) {
        for (int j=0; j<dim_; j++) {
            Teuchos::ArrayRCP< SC > values = nodeList->getDataNonConst( j );
            values[i] = (*pointsRepeated)[i][j];
        }
    }
    return nodeList;
}

template <class SC, class LO, class GO, class NO>
void Domain<SC, LO, GO, NO>::exportNodeFlags(string name)
{
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(this->getMapUnique()));
        vec_int_ptr_Type BCFlags = this->getBCFlagUnique(); // Unique flags at points

        Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
        for(int i=0; i< entries.size(); i++){
            entries[i] = BCFlags->at(i);
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

        exPara->setup("Mesh_Node_Flags_"+name,this->getMesh(), this->FEType_);

        exPara->addVariable(exportSolutionConst, "Flags", "Scalar", 1,this->getMapUnique()); 
        exPara->save(0.0);

        exPara->closeExporter();
} 

template <class SC, class LO, class GO, class NO>
void Domain<SC, LO, GO, NO>::exportSurfaceNormals(string name)
{
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(this->getMapVecFieldUnique()));
        exportSolution->putScalar(0.);
        ElementsPtr_Type elementsC = this->getElementsC(); // Unique flags at points

        Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);

        // We iterate over all elements and compute the surface normals to export them
        for (UN T=0; T<elementsC->numberElements(); T++) {
            FiniteElement fe = elementsC->getElement( T );
            ElementsPtr_Type subEl = fe.getSubElements(); // might be null
            for (int surface=0; surface<fe.numSubElements(); surface++) {
                FiniteElement feSub = subEl->getElement( surface  );
                vec_int_Type nodeListElement = fe.getVectorNodeList();
                if(subEl->getDimension() == dim_-1 ){
                    vec_int_Type nodeList = feSub.getVectorNodeListNonConst();
                    
                    vec_dbl_Type v_E(dim_,1.);
                    double norm_v_E=1.;

                    Helper::computeSurfaceNormal(this->dim_,this->getPointsRepeated(),nodeList,v_E,norm_v_E);

                    for(int j=0; j< nodeList.size(); j++){
                        for(int i=0; i< this->dim_; i++){
                            if(this->getMapUnique()->getLocalElement(this->getMapRepeated()->getGlobalElement(nodeList[j])) != -1 ){  // We only write when we have the node in unique distribution, because we are lazy. Otherwise we would do an import etc. But this will not give us further information.
                                LO id = this->getMapUnique()->getLocalElement(this->getMapRepeated()->getGlobalElement(nodeList[j]));
                                entries[id*this->dim_+i] = 1./norm_v_E * v_E[i];
                            }
                        }
                    }
                }
            }
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

        exPara->setup("Mesh_Surface_Directions_"+name,this->getMesh(), this->FEType_);

        exPara->addVariable(exportSolutionConst, "SurfaceNormals", "Vector", this->dim_,this->getMapUnique()); 
        exPara->save(0.0);

        exPara->closeExporter();
} 

template <class SC, class LO, class GO, class NO>
void Domain<SC, LO, GO, NO>::exportElementOrientation(string name)
{
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(this->getElementMap()));
        exportSolution->putScalar(0.);
        ElementsPtr_Type elementsC = this->getElementsC(); // Unique flags at points

        Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);

        SC detB;
        SmallMatrix<SC> B(this->dim_);
        SmallMatrix<SC> Binv(this->dim_);
        // We iterate over all elements and compute the surface normals to export them
        for (UN T=0; T<elementsC->numberElements(); T++) {
           Helper::buildTransformation(elementsC->getElement(T).getVectorNodeList(), this->getPointsRepeated(), B);
           detB = B.computeInverse(Binv);
           entries[T] = detB;
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

        exPara->setup("Mesh_Element_Orientation_"+name,this->getMesh(), "P0");

        exPara->addVariable(exportSolutionConst, "VolElement", "Scalar", 1,this->getElementMap()); 
        exPara->save(0.0);

        exPara->closeExporter();
} 


template <class SC, class LO, class GO, class NO>
void Domain<SC, LO, GO, NO>::exportElementFlags(string name)
{
        Teuchos::RCP<ExporterParaView<SC,LO,GO,NO> > exPara(new ExporterParaView<SC,LO,GO,NO>());

        Teuchos::RCP<MultiVector<SC,LO,GO,NO> > exportSolution(new MultiVector<SC,LO,GO,NO>(this->getElementMap()));

        Teuchos::ArrayRCP< SC > entries  = exportSolution->getDataNonConst(0);
        ElementsPtr_Type elements = this->getElementsC(); // element list

        for(int i=0; i< entries.size(); i++){
            entries[i] = elements->getElement(i).getFlag(); // element flags
        }

        Teuchos::RCP<const MultiVector<SC,LO,GO,NO> > exportSolutionConst = exportSolution;

        exPara->setup("Mesh_Element_Flags_"+name,this->getMesh(), "P0");

        exPara->addVariable(exportSolutionConst, "Flags", "Scalar", 1,this->getElementMap(), this->getElementMap());

        exPara->save(0.0);

        exPara->closeExporter();

}   
    
}
#endif