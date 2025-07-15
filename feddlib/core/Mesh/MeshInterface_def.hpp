#ifndef MESHINTERFACE_def_hpp
#define MESHINTERFACE_def_hpp
#include "MeshInterface_decl.hpp"

/*!
 Definition of MeshInterface

 @brief  MeshInterface
 @author Christian Hochmuth
 @version 1.0
 @copyright CH
 */


//search all code for these functions and move them to Tools file.
template <typename T>
std::vector<T> sort_from_ref(
                        std::vector<T> const& in,
                        std::vector<int> const& reference
                        ) {
    std::vector<T> ret(in.size());

    int const size = in.size();
    for (int i = 0; i < size; ++i)
        ret[i] = in[reference[i]];

    return ret;
};

template <typename T>
std::vector<T> sort_from_ref(
                        std::vector<T> const& in,
                        std::vector<long long> const& reference
                        ) {
    std::vector<T> ret(in.size());

    int const size = in.size();
    for (long long i = 0; i < size; ++i)
        ret[i] = in[reference[i]];

    return ret;
};

namespace FEDD {

template <class SC, class LO, class GO, class NO>
MeshInterface<SC,LO,GO,NO>::MeshInterface():
indicesGlobalMatched_(),
indicesGlobalMatchedOrigin_(),
indicesGlobalMatchedUnique_(),
isPartitioned_(false),
partialCouplingFlag_(0),
partialCouplingType_(0),
comm_()
{

}

template <class SC, class LO, class GO, class NO>
MeshInterface<SC,LO,GO,NO>::MeshInterface(CommConstPtr_Type comm):
indicesGlobalMatched_(),
indicesGlobalMatchedOrigin_(),
indicesGlobalMatchedUnique_(),
isPartitioned_(false),
partialCouplingFlag_(0),
partialCouplingType_(0),
comm_(comm)
{
    
}
    
template <class SC, class LO, class GO, class NO>
void MeshInterface<SC,LO,GO,NO>::partitionMeshInterface(MapPtr_Type mapRepeated, MapPtr_Type mapUnique){

    // TODO: irgenwann mal in indicesGlobalMatchedOrigin_ abaendern
    vec3D_GO_ptr_Type indicesSerial = indicesGlobalMatchedOrigin_;

    // Das ist der repated Vektor
    indicesGlobalMatched_.reset(new vec3D_GO_Type ( indicesSerial->size() , vec2D_GO_Type( 2 , vec_GO_Type(0) ) ) );

    // Das ist der unique Vektor
    indicesGlobalMatchedUnique_.reset(new vec3D_GO_Type ( indicesSerial->size() , vec2D_GO_Type( 2 , vec_GO_Type(0) ) ) );

    for (int flag=0; flag<indicesSerial->size(); flag++) {
        for (int i=0; i<indicesSerial->at(flag).at(0).size(); i++) {
            // Fuer den repated Vektor
            GO indexThis = mapRepeated->getLocalElement( indicesSerial->at(flag).at(0).at(i) );
            if ( indexThis != Teuchos::OrdinalTraits<LO>::invalid() ) {
                indicesGlobalMatched_->at(flag).at(0).push_back( indicesSerial->at(flag).at(0).at(i) );
                indicesGlobalMatched_->at(flag).at(1).push_back( indicesSerial->at(flag).at(1).at(i) );
            }

            GO indexThisUni = mapUnique->getLocalElement( indicesSerial->at(flag).at(0).at(i) );
            // Fuer den unique Vektor
            if ( indexThisUni != Teuchos::OrdinalTraits<LO>::invalid() ) {
                indicesGlobalMatchedUnique_->at(flag).at(0).push_back( indicesSerial->at(flag).at(0).at(i) );
                indicesGlobalMatchedUnique_->at(flag).at(1).push_back( indicesSerial->at(flag).at(1).at(i) );
            }

        }
    }

    isPartitioned_ = true;
}


template <class SC, class LO, class GO, class NO>
void MeshInterface<SC,LO,GO,NO>::determineInterface( vec2D_dbl_ptr_Type pointsRepThis,  vec2D_dbl_ptr_Type pointsRepOther, vec_int_ptr_Type flagThis, vec_int_ptr_Type flagOther, vec_int_Type relevant_flag_vec ){

    indicesGlobalMatched_.reset(new vec3D_GO_Type( relevant_flag_vec.size(), vec2D_GO_Type(2) ) );

    // Damit wir auf den nicht-partionierten Vektor zugreifen koennen.
    indicesGlobalMatchedOrigin_.reset(new vec3D_GO_Type( relevant_flag_vec.size(), vec2D_GO_Type(2) ) );

    // Points N x dim vectors
    for (int flag=0; flag<relevant_flag_vec.size(); flag++) {
        vec2D_dbl_Type pointThis(0);
        vec2D_dbl_Type pointOther(0);
        vec_GO_Type indexGlobalThis(0);
        vec_GO_Type indexGlobalOther(0);
        vec_int_Type indexThis(0);
        vec_int_Type indexOther(0);
        int counter = 0;
        for (int i=0; i<pointsRepThis->size(); i++) {
            if (flagThis->at(i) == relevant_flag_vec.at(flag)) {
                pointThis.push_back( pointsRepThis->at(i) );
                indexGlobalThis.push_back(i);
                indexThis.push_back(counter);
                counter++;
            }
        }
        counter = 0;
        for (int i=0; i<pointsRepOther->size(); i++) {
            if (flagOther->at(i) == relevant_flag_vec.at(flag)) {
                pointOther.push_back( pointsRepOther->at(i) );
                indexGlobalOther.push_back(i);
                indexOther.push_back(counter);
                counter++;
            }
        }

        // std::cout << "serial determineInterface - number flag nodes this:" << indexThis.size()<< " other:" << indexOther.size() << std::endl;
        
        sort(indexThis.begin(), indexThis.end(),
             [&](const int& a, const int& b) {
                 return  pointThis[a] < pointThis[b];
             }
             );

        sort(indexOther.begin(), indexOther.end(),
             [&](const int& a, const int& b) {
                 return  pointOther[a] < pointOther[b];
             }
             );

        indexGlobalThis         = sort_from_ref(indexGlobalThis, indexThis);
        indexGlobalOther   		= sort_from_ref(indexGlobalOther, indexOther);

        TEUCHOS_TEST_FOR_EXCEPTION( indexGlobalThis.size()!=indexGlobalOther.size(), std::logic_error, "Interfaces do not have the same length!");


        indicesGlobalMatched_->at(flag).at(0) = indexGlobalThis;
        indicesGlobalMatched_->at(flag).at(1) = indexGlobalOther;

        indicesGlobalMatchedOrigin_->at(flag).at(0) = indexGlobalThis;
        indicesGlobalMatchedOrigin_->at(flag).at(1) = indexGlobalOther;

    }

}
    
template <class SC, class LO, class GO, class NO>
int MeshInterface<SC,LO,GO,NO>::isPartialCouplingFlag(int flag){

    auto it = std::find( partialCouplingFlag_.begin(), partialCouplingFlag_.end(), flag );
    if ( it!=partialCouplingFlag_.end() )
        return std::distance( partialCouplingFlag_.begin(), it );
    else
        return -1;
}
    
template <class SC, class LO, class GO, class NO>
void MeshInterface<SC,LO,GO,NO>::setPartialCoupling(int flag, std::string type){
    partialCouplingFlag_.push_back(flag);
    partialCouplingType_.push_back(type);
    
}

template <class SC, class LO, class GO, class NO>
int MeshInterface<SC,LO,GO,NO>::getPartialCouplingFlag(int i){
    TEUCHOS_TEST_FOR_EXCEPTION(sizePartialCoupling()-1 < i, std::runtime_error, "There is no partial coupling for this index");
    return partialCouplingFlag_[i];
}

template <class SC, class LO, class GO, class NO>
std::string MeshInterface<SC,LO,GO,NO>::getPartialCouplingType(int i){
    TEUCHOS_TEST_FOR_EXCEPTION(sizePartialCoupling()-1 < i, std::runtime_error, "There is no partial coupling for this index");
    return partialCouplingType_[i];
}
    
template <class SC, class LO, class GO, class NO>
int MeshInterface<SC,LO,GO,NO>::sizePartialCoupling(){
    return partialCouplingFlag_.size();
}
    
template <class SC, class LO, class GO, class NO>
void MeshInterface<SC,LO,GO,NO>::determineInterfaceParallelAndDistance( vec2D_dbl_ptr_Type pointsUniThis,  vec2D_dbl_ptr_Type pointsUniOther, vec_int_ptr_Type flagUniThis, vec_int_ptr_Type flagUniOther, vec_int_Type relevant_flag_vec, MapConstPtr_Type mapUniThis, MapConstPtr_Type mapUniOther, vec_dbl_ptr_Type &distancesToInterface, vec2D_dbl_ptr_Type pointsRepThis, int dim ) {
    
    SC eps100 = 100. * Teuchos::ScalarTraits<SC>::eps();
    GO invalid = Teuchos::OrdinalTraits<GO>::invalid();

    indicesGlobalMatched_.reset(new vec3D_GO_Type( relevant_flag_vec.size(), vec2D_GO_Type(2) ) );
    
    // Damit wir auf den nicht-partionierten Vektor zugreifen koennen.
    // Build in parallel first and then build this global index vector
    // Adjust size
    indicesGlobalMatchedOrigin_.reset(new vec3D_GO_Type( relevant_flag_vec.size(), vec2D_GO_Type(2) ) );
    
    // send data for the relevant flag
    for (int flagIndex=0; flagIndex<relevant_flag_vec.size(); flagIndex++) {
        int flag = relevant_flag_vec[flagIndex];
        //Warning !!! The following information is local!
        vec_GO_Type indexGlobalCommThis(0);
        vec_GO_Type indexGlobalCommOther(0);
        
        for (int i=0; i<pointsUniThis->size(); i++) {
            if ( flagUniThis->at(i) == flag ) {
                indexGlobalCommThis.push_back( mapUniThis->getGlobalElement( i ) );
            }
        }

        for (int i=0; i<pointsUniOther->size(); i++) {
            if ( flagUniOther->at(i) == flag ) {
                indexGlobalCommOther.push_back( mapUniOther->getGlobalElement( i ) );
            }
        }
        
        //Communicate everything with MultiVectors.
        //First determine the length of the new global vector
        GO numInterfaceGlobalThis = 0;
        GO numInterfaceGlobalOther = 0;

        Teuchos::reduceAll( *this->comm_, Teuchos::REDUCE_SUM, (GO) indexGlobalCommThis.size(), Teuchos::ptrFromRef( numInterfaceGlobalThis ) );

        Teuchos::reduceAll( *this->comm_, Teuchos::REDUCE_SUM, (GO) indexGlobalCommOther.size(), Teuchos::ptrFromRef( numInterfaceGlobalOther ) );
        
        MapPtr_Type mapThis = Teuchos::rcp( new Map_Type( -1, Teuchos::arrayViewFromVector( indexGlobalCommThis ), 0, this->comm_ ) );

        MapPtr_Type mapOther = Teuchos::rcp( new Map_Type(  -1, Teuchos::arrayViewFromVector( indexGlobalCommOther ), 0, this->comm_ ) );
        
        // std::cout << "numInterfaceGlobalThis:" << numInterfaceGlobalThis << std::endl;
        // std::cout << "numInterfaceGlobalOther:" << numInterfaceGlobalOther << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION( numInterfaceGlobalThis != numInterfaceGlobalOther, std::runtime_error, "DetermineInterfaceInParallel failed. ThisMesh and OtherMesh seem to have different numbers of interface nodes." );
        
        std::vector<GO> gatherAllIndices(numInterfaceGlobalThis);
        std::iota ( std::begin( gatherAllIndices ), std::end( gatherAllIndices ), 0 );

        MapPtr_Type linearMapThis = Teuchos::rcp( new Map_Type(  numInterfaceGlobalThis, indexGlobalCommThis.size(), 0, this->comm_ ) );
        MapPtr_Type linearMapOther = Teuchos::rcp( new Map_Type(numInterfaceGlobalThis, indexGlobalCommOther.size(), 0, this->comm_ ) );
        MapPtr_Type gatherAllMap = Teuchos::rcp( new Map_Type(  invalid, Teuchos::arrayViewFromVector( gatherAllIndices ), 0, this->comm_ ) );

        // We would like to use the Teuchos version of MPI_Allgatherv, which does not exist. Therefore we gatherv on a root and broadcast afterwards
        // Gather local lengths first
        vec_int_Type localSizeThis( this->comm_->getSize(),0 );
        vec_int_Type localSizeOther( this->comm_->getSize(),0 );
        int root = 0;

        int sizeThis = indexGlobalCommThis.size();
        int sizeOther = indexGlobalCommOther.size();
        Teuchos::gather<int, int>( &sizeThis, 1, &localSizeThis[0], 1, root, *this->comm_ );
        Teuchos::gather<int, int>( &sizeOther, 1, &localSizeOther[0], 1, root, *this->comm_ );

        GO* sendThis = NULL;
        if (indexGlobalCommThis.size()>0)
            sendThis = &indexGlobalCommThis[0];
        GO* sendOther = NULL;
        if (indexGlobalCommOther.size()>0)
            sendOther = &indexGlobalCommOther[0];
        
        vec_GO_Type gatheredThis( numInterfaceGlobalThis , -1);
        vec_GO_Type gatheredOther( numInterfaceGlobalThis, -1);
        
        vec_int_Type displacementsThis( ((int) this->comm_->getRank()==root ) * this->comm_->getSize(),0 );
        vec_int_Type displacementsOther( ((int) this->comm_->getRank()==root ) * this->comm_->getSize(),0 );
        int* displThis = NULL;
        if (displacementsThis.size()>0)
            displThis = &displacementsThis[0];
        int* displOther = NULL;
        if (displacementsOther.size()>0)
            displOther = &displacementsOther[0];
        
        for (int i=1; i<displacementsThis.size(); i++)
            displacementsThis[i] = displacementsThis[i-1] + localSizeThis[i-1];
        for (int i=1; i<displacementsOther.size(); i++)
            displacementsOther[i] = displacementsOther[i-1] + localSizeOther[i-1];
        
        Teuchos::gatherv<int,GO>( sendThis, indexGlobalCommThis.size(), &gatheredThis[0], &localSizeThis[0], displThis, root, *this->comm_ );
        Teuchos::gatherv<int,GO>( sendOther, indexGlobalCommOther.size(), &gatheredOther[0], &localSizeOther[0], displOther, root, *this->comm_ );
        
        //Now we broadcast from root
        Teuchos::broadcast<int,GO>( *this->comm_, root, Teuchos::arrayViewFromVector( gatheredThis ) );
        Teuchos::broadcast<int,GO>( *this->comm_, root, Teuchos::arrayViewFromVector( gatheredOther ) );

        MapPtr_Type mapAllThis = Teuchos::rcp( new Map_Type(invalid, Teuchos::arrayViewFromVector( gatheredThis ), 0, this->comm_ ) );
        MapPtr_Type mapAllOther = Teuchos::rcp( new Map_Type( invalid, Teuchos::arrayViewFromVector( gatheredOther ), 0, this->comm_ ) );

        bool meshOnRank = false;
        if (pointsUniThis->size() > 0)
            meshOnRank = true;
        
        MultiVectorPtr_Type mvLocalThis;
        MultiVectorPtr_Type mvLocalOther;

        mvLocalThis = Teuchos::rcp( new MultiVector_Type( mapThis, dim + 1 ) );
        mvLocalOther = Teuchos::rcp( new MultiVector_Type( mapOther, dim + 1 ) );

        // Fill local vectors with node coordinates and global index; last column is the global index
        // This
        for (int j=0; j<dim; j++) {
            Teuchos::ArrayRCP< SC > data =  mvLocalThis->getDataNonConst(j);
            for (int i=0; i<data.size(); i++) {
                LO index = mapUniThis->getLocalElement( indexGlobalCommThis[i] );
                data[i] = pointsUniThis->at( index )[j];
            }
        }

        Teuchos::ArrayRCP< SC > dataThis =  mvLocalThis->getDataNonConst( dim );
        for (int i=0; i<dataThis.size(); i++)
            dataThis[i] = indexGlobalCommThis[i];

        
        // Other
        for (int j=0; j<dim; j++) {
            Teuchos::ArrayRCP< SC > data =  mvLocalOther->getDataNonConst(j);
            for (int i=0; i<data.size(); i++) {
                LO index = mapUniOther->getLocalElement( indexGlobalCommOther[i] );
                data[i] = pointsUniOther->at( index )[j];
            }
        }

        Teuchos::ArrayRCP< SC > dataOther =  mvLocalOther->getDataNonConst( dim );
        for (int i=0; i<dataOther.size(); i++)
            dataOther[i] = indexGlobalCommOther[i];


        MultiVectorPtr_Type mvGlobalThis;
        MultiVectorPtr_Type mvGlobalOther;

        mvGlobalThis = Teuchos::rcp( new MultiVector_Type( mapAllThis, dim + 1 ) );
        mvGlobalOther = Teuchos::rcp( new MultiVector_Type( mapAllOther, dim + 1 ) );
        // Communicate, we might want to use gatherv and broadcast here aswell
        mvGlobalThis->exportFromVector( mvLocalThis, true, "Insert", "Reverse" );
        mvGlobalOther->exportFromVector( mvLocalOther, true, "Insert", "Reverse" );

        // now we can go through the same data on every rank like in the serial case
        // copy to std::vectors first; we should try to avoid this in a optimized implementation
        // Points N x dim vectors
        
        vec_int_Type indexThis( numInterfaceGlobalThis );
        vec_int_Type indexOther( numInterfaceGlobalThis );
        std::iota ( std::begin( indexThis ), std::end( indexThis ), 0 );
        std::iota ( std::begin( indexOther ), std::end( indexOther ), 0 );

        vec2D_dbl_Type pointThis( numInterfaceGlobalThis, vec_dbl_Type(dim) );
        vec2D_dbl_Type pointOther( numInterfaceGlobalThis, vec_dbl_Type(dim) );

        vec_GO_Type indexGlobalThis( numInterfaceGlobalThis );
        vec_GO_Type indexGlobalOther( numInterfaceGlobalThis );
        

        int counter = 0;
        {
            Teuchos::ArrayRCP< const SC > dataX = mvGlobalThis->getData(0);
            Teuchos::ArrayRCP< const SC > dataY = mvGlobalThis->getData(1);
            Teuchos::ArrayRCP< const SC > dataZ;
            Teuchos::ArrayRCP< const SC > dataGlob;
            if (dim==3) {
                dataZ = mvGlobalThis->getData(2);
                dataGlob = mvGlobalThis->getData(dim);
            }
            else if (dim==2)
                dataGlob = mvGlobalThis->getData(dim);
        
            for (int i=0; i<pointThis.size(); i++) {
                pointThis[i][0] = dataX[i];
                pointThis[i][1] = dataY[i];
                if (dim==3)
                    pointThis[i][2] = dataZ[i];
                indexGlobalThis[i] = (GO) dataGlob[i] + eps100;
            }
        }
        {
            Teuchos::ArrayRCP< const SC > dataX = mvGlobalOther->getData(0);
            Teuchos::ArrayRCP< const SC > dataY = mvGlobalOther->getData(1);
            Teuchos::ArrayRCP< const SC > dataZ;
            Teuchos::ArrayRCP< const SC > dataGlob;
            if (dim==3) {
                dataZ = mvGlobalOther->getData(2);
                dataGlob = mvGlobalOther->getData(dim);
            }
            else if (dim==2)
                dataGlob = mvGlobalOther->getData(dim);
            
            for (int i=0; i<pointOther.size(); i++) {
                pointOther[i][0] = dataX[i];
                pointOther[i][1] = dataY[i];
                if (dim==3)
                    pointOther[i][2] = dataZ[i];
                indexGlobalOther[i] = (GO) dataGlob[i] + eps100;
         
            }
        }
        
        sort(indexThis.begin(), indexThis.end(),
             [&](const int& a, const int& b) {
                 return  pointThis[a] < pointThis[b];
             }
             );
        
        sort(indexOther.begin(), indexOther.end(),
             [&](const int& a, const int& b) {
                 return  pointOther[a] < pointOther[b];
             }
             );
        
        indexGlobalThis = sort_from_ref(indexGlobalThis, indexThis);
        indexGlobalOther = sort_from_ref(indexGlobalOther, indexOther);
        
        
        indicesGlobalMatched_->at(flagIndex).at(0) = indexGlobalThis;
        indicesGlobalMatched_->at(flagIndex).at(1) = indexGlobalOther;
//        std::cout << "indexGlobalThis->size():" << indexGlobalThis.size() << std::endl;
//        std::cout << "indexGlobalOther->size():" << indexGlobalOther.size() << std::endl;
        indicesGlobalMatchedOrigin_->at(flagIndex).at(0) = indexGlobalThis;
        indicesGlobalMatchedOrigin_->at(flagIndex).at(1) = indexGlobalOther;
     
        //determine distance for this flag.
        calculateDistancesToInterfaceParallel( distancesToInterface, pointThis, pointsRepThis );
    }
}
 
    
template <class SC, class LO, class GO, class NO>
void MeshInterface<SC,LO,GO,NO>::calculateDistancesToInterfaceParallel( vec_dbl_ptr_Type &distancesToInterface, vec2D_dbl_Type &pointThis/*global interface, every proc has same information*/, vec2D_dbl_ptr_Type sourceNodesRep /*partitioned points*/)
{
    
    int dim = -1;
    if (pointThis.size() > 0)
        dim = pointThis[0].size();

    // See comments for calculateDistancesToInterface() in class Domain.
    // This is the parallel version
        
    // Not partitioned interfaces
    vec3D_GO_ptr_Type indicesGlobalMatched = this->indicesGlobalMatched_;
    
    // SourceNodes (z.B. Fluidgebiet) aus dem Gitter ziehen
    // Gefaehrlich, da Veraenderungen an SourceNodesRep auch PointsRep_ veraendern;
    // Sonst auf 0.0 setzen und dann Werte reinschreiben.
    // Am besten waere es Getter zu haben mit return vec2D_dbl_ptr_Type bzw. const vec2vec2D_dbl_ptr_Type, damit nichts passieren kann.
    
    // Zaehle mit Hilfe von indicesGlobalMatched wie viele Interfaceknoten es insgesamt gibt.
    // Here we should loop over all flags in the case of more than one flag.
    int numberInterfaceNodes = indicesGlobalMatched->at(0).at(0).size();
    
    // EndNodes (Interfaceknoten) aus dem Gitter ziehen; mit Hilfe von IndicesGlobalMatched
    vec2D_dbl_ptr_Type endNodesRep = Teuchos::rcpFromRef(pointThis);
    
    // Distanz zum Interface auf eine unrealistisch grosse Zahl setzen.
    // In DistancesToInterface_->at(i) steht die Distanz der globalen Knoten-ID i zum Interface
    
    if ( distancesToInterface.is_null() )
        distancesToInterface = Teuchos::rcp( new vec_dbl_Type( sourceNodesRep->size(), 1000.0 ) );
    
    // Berechne nun den Abstand von den sourceNodesRep zu den EndNodesRep
    double distance = 0.0;
    for(int i = 0; i < sourceNodesRep->size(); i++)
    {
        for(int j = 0; j < endNodesRep->size(); j++)
        {
            for(int k = 0; k < dim; k++)
            {
                distance = distance + std::pow( sourceNodesRep->at(i).at(k) - endNodesRep->at(j).at(k), 2.0 );
            }
            
            // Noch die Wurzel ziehen
            distance = std::sqrt(distance);
            
            if(distancesToInterface->at(i) > distance)
                distancesToInterface->at(i) = distance;
            
            // Reseten
            distance = 0.0;
        }
    }
}

    
template <class SC, class LO, class GO, class NO>
void MeshInterface<SC,LO,GO,NO>::buildFromOtherInterface( Teuchos::RCP<MeshInterface> otherMeshInterface){

    TEUCHOS_TEST_FOR_EXCEPTION( otherMeshInterface->indicesGlobalMatched_->size()<=0, std::logic_error, "MeshInterface given is empty.");
    TEUCHOS_TEST_FOR_EXCEPTION( otherMeshInterface->indicesGlobalMatched_->at(0).size()<=0, std::logic_error, "MeshInterface given is empty.");
    TEUCHOS_TEST_FOR_EXCEPTION( otherMeshInterface->isPartitioned_, std::logic_error, "Other MeshInterface already partitioned.");
    TEUCHOS_TEST_FOR_EXCEPTION( isPartitioned_, std::logic_error, "This MeshInterface already partitioned.");

    indicesGlobalMatched_.reset(new vec3D_GO_Type( otherMeshInterface->indicesGlobalMatched_->size(), vec2D_GO_Type ( 2 , vec_GO_Type( 0 ) ) ) );

    indicesGlobalMatchedOrigin_.reset(new vec3D_GO_Type( otherMeshInterface->indicesGlobalMatchedOrigin_->size(), vec2D_GO_Type ( 2 , vec_GO_Type( 0 ) ) ) );
    
    for (int i=0; i<otherMeshInterface->indicesGlobalMatched_->size(); i++) {
        indicesGlobalMatched_->at(i).at(0).resize( otherMeshInterface->indicesGlobalMatched_->at(i).at(0).size() );
        indicesGlobalMatched_->at(i).at(1).resize( otherMeshInterface->indicesGlobalMatched_->at(i).at(1).size() );
        indicesGlobalMatchedOrigin_->at(i).at(0).resize( otherMeshInterface->indicesGlobalMatchedOrigin_->at(i).at(0).size() );
        indicesGlobalMatchedOrigin_->at(i).at(1).resize( otherMeshInterface->indicesGlobalMatchedOrigin_->at(i).at(1).size() );
        for (int j=0; j<otherMeshInterface->indicesGlobalMatched_->at(i).at(0).size(); j++) {
            indicesGlobalMatched_->at(i).at(0).at(j) = otherMeshInterface->indicesGlobalMatched_->at(i).at(1).at(j);
            indicesGlobalMatched_->at(i).at(1).at(j) = otherMeshInterface->indicesGlobalMatched_->at(i).at(0).at(j);

            indicesGlobalMatchedOrigin_->at(i).at(0).at(j) = otherMeshInterface->indicesGlobalMatchedOrigin_->at(i).at(1).at(j);
            indicesGlobalMatchedOrigin_->at(i).at(1).at(j) = otherMeshInterface->indicesGlobalMatchedOrigin_->at(i).at(0).at(j);

        }
    }
}
    
template <class SC, class LO, class GO, class NO>
void MeshInterface<SC,LO,GO,NO>::print(CommConstPtr_Type comm){

    for (int i=0; i<indicesGlobalMatched_->size(); i++) {
        for (int j=0; j<indicesGlobalMatched_->at(i).at(0).size(); j++) {
            std::cout <<  comm->getRank()<<" Matched IDs for flag " << i << " :" << indicesGlobalMatched_->at(i).at(0).at(j) << " - " << indicesGlobalMatched_->at(i).at(0).at(j) << std::endl;
        }
    }
}

template <class SC, class LO, class GO, class NO>
typename MeshInterface<SC,LO,GO,NO>::vec3D_GO_ptr_Type MeshInterface<SC,LO,GO,NO>::getIndicesGlobalMatched(){
    return indicesGlobalMatched_;
}

template <class SC, class LO, class GO, class NO>
typename MeshInterface<SC,LO,GO,NO>::vec3D_GO_ptr_Type MeshInterface<SC,LO,GO,NO>::getIndicesGlobalMatchedOrigin()
{
    return indicesGlobalMatchedOrigin_;
}

template <class SC, class LO, class GO, class NO>
typename MeshInterface<SC,LO,GO,NO>::vec3D_GO_ptr_Type MeshInterface<SC,LO,GO,NO>::getIndicesGlobalMatchedUnique()
{
    return indicesGlobalMatchedUnique_;
}


}
#endif
